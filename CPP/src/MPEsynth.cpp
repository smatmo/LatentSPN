/*
MPEsynth.cpp

Generate SPNs with Poon & Domingos structure and perform MPE inference by
- max-backtracking in the original SPN (corresponds to MPE inference in the augmented SPN with deterministic twin-weights)
- max-backtracking in the original SPN with correction weigths corresponding to uniform twin-weights in the augmented SPN
- exhaustive enumeration

Robert Peharz
October 2016
*/


#include "Network.hpp"
#include "DataSource.hpp"
#include "PDStruct.hpp"
#include "Utils.hpp"
#include "FileIO.hpp"
#include "EM.hpp"
#include "Utils.hpp"
#include <iostream>
#include <string>
#include <random>
#include <iomanip>
#include <sys/types.h>
#include <sys/stat.h>


using namespace SPN;
using namespace std;

//
// draw sample from a Dirichlet
//
DoubleVec Dirichlet(double alpha, unsigned K, std::default_random_engine &generator) {
	std::gamma_distribution<double> distribution(alpha, 1.0);

	DoubleVec ret(K, 0.0);
	double sum = 0.0;
	for (unsigned k = 0; k < K; k++) {
		double val = distribution(generator);
		ret[k] = val;
		sum += val;
	}
	for (unsigned k = 0; k < K; k++)
		ret[k] /= sum;
	return ret;
}


//
// increase state counter for exhaustive enumeration
//
bool incStateCounter(UnsignedVec &state, const UnsignedVec &numStates) {
	if (state.size() != numStates.size())
		throw std::runtime_error("incStateCounter: (state.size() != maxState.size()");

	unsigned k = 0;
	for (k = 0; k < state.size(); k++) {
		state[k]++;
		if (state[k] == numStates[k])
			state[k] = 0;
		else
			break;
	}
	return k == state.size();
}


DoubleVec Unsigned2DoubleVec(const UnsignedVec &uv) {
	DoubleVec ret;
	ret.resize(uv.size());
	for (size_t k = 0; k < uv.size(); k++)
		ret[k] = uv[k];
	return ret;
}

//
// generate a PD structure
//
Network* makeFullPDStructure(NetworkFactory &factory,
	unsigned width, unsigned height, 
	UnsignedVec coarseRes, UnsignedVec coarseResX, UnsignedVec coarseResY,
	size_t maxNumSumsPerRectangle, 
	size_t numDistPerPixel,	NodeType inputType) {

	std::vector<PDRectangle* > rectangles;
	std::vector<PDPartition* > partitions;
	std::map<PDRectangle, PDRectangle* > rectangleMap;
	std::map<unsigned, std::vector<PDRectangle*> > rectanglesOfLevel;

	///////////////////////////////////////////////////////////////////
	// Determine Structure                                           //
	// we start with a single rectangle and decompose it further     //
	///////////////////////////////////////////////////////////////////
	PDRectangle* root = new PDRectangle;
	root->x_ = 0;
	root->y_ = 0;
	root->width_ = width;
	root->height_ = height;

	// rectangles will contain a list of pointers to all generated rectangles
	rectangles.push_back(root);
	rectangleMap[*root] = root;

	for (unsigned k = 0; k < rectangles.size(); k++) {
		// if rectangle was already split, or if it is of unit size, don't split
		if ((rectangles[k]->childPartitions_.empty() == false)
			|| ((rectangles[k]->width_ == 1) && (rectangles[k]->height_ == 1)))
			continue;

		PDRectangle *curRectangle = rectangles[k];

		//
		// deltaX, deltaY denote the used resolutions 
		//
		unsigned deltaX = 1;
		unsigned deltaY = 1;

		if ((coarseResX.empty()) || (coarseResY.empty())) {
			for (unsigned l = 0; l < coarseRes.size(); l++) {
				size_t idx = coarseRes.size() - 1 - l;
				if (
					(coarseRes[idx] <= curRectangle->width_) &&
					(coarseRes[idx] <= curRectangle->height_) &&
					((coarseRes[idx] != curRectangle->height_) ||
					(coarseRes[idx] != curRectangle->width_))) {
					deltaX = coarseRes[idx];
					deltaY = coarseRes[idx];
					break;
				}
			}
		}
		else {
			for (unsigned l = 0; l < coarseResX.size(); l++) {
				size_t idx = coarseResX.size() - 1 - l;
				if (coarseResX[idx] < curRectangle->width_) {
					deltaX = coarseResX[idx];
					break;
				}
			}
			for (unsigned l = 0; l < coarseResY.size(); l++) {
				size_t idx = coarseResY.size() - 1 - l;
				if (coarseResY[idx] < curRectangle->height_) {
					deltaY = coarseResY[idx];
					break;
				}
			}
		}

		// split into smaller rectangles -- x dimension
		PDRectangle r1, r2;
		PDRectangle *r1ptr;
		PDRectangle *r2ptr;

		r1.y_ = curRectangle->y_;
		r1.height_ = curRectangle->height_;
		r1.x_ = curRectangle->x_;
		r2.y_ = curRectangle->y_;
		r2.height_ = curRectangle->height_;

		for (unsigned x = deltaX; x < curRectangle->width_; x = x + deltaX) {
			r1.width_ = x;
			r2.x_ = curRectangle->x_ + x;
			r2.width_ = curRectangle->width_ - x;

			if (rectangleMap.count(r1) > 0)
				r1ptr = rectangleMap[r1];
			else {
				r1ptr = new PDRectangle;
				*r1ptr = r1;
				rectangleMap[r1] = r1ptr;
				rectangles.push_back(r1ptr);
			}
			if (rectangleMap.count(r2) > 0)
				r2ptr = rectangleMap[r2];
			else {
				r2ptr = new PDRectangle;
				*r2ptr = r2;
				rectangleMap[r2] = r2ptr;
				rectangles.push_back(r2ptr);
			}

			PDPartition *p = new PDPartition;
			partitions.push_back(p);
			p->parentRectangle_ = rectangles[k];
			p->childRectangles_.push_back(r1ptr);
			p->childRectangles_.push_back(r2ptr);
			curRectangle->childPartitions_.push_back(p);
			r1ptr->parentPartitions_.push_back(p);
			r2ptr->parentPartitions_.push_back(p);
		}

		// split into smaller rectangles -- y dimension
		r1.y_ = curRectangle->y_;
		r1.width_ = curRectangle->width_;
		r1.x_ = curRectangle->x_;
		r2.x_ = curRectangle->x_;
		r2.width_ = curRectangle->width_;

		for (unsigned y = deltaY; y < curRectangle->height_; y = y + deltaY) {
			r1.height_ = y;
			r2.y_ = curRectangle->y_ + y;
			r2.height_ = curRectangle->height_ - y;

			if (rectangleMap.count(r1) > 0)
				r1ptr = rectangleMap[r1];
			else {
				r1ptr = new PDRectangle;
				*r1ptr = r1;
				rectangleMap[r1] = r1ptr;
				rectangles.push_back(r1ptr);
			}
			if (rectangleMap.count(r2) > 0)
				r2ptr = rectangleMap[r2];
			else {
				r2ptr = new PDRectangle;
				*r2ptr = r2;
				rectangleMap[r2] = r2ptr;
				rectangles.push_back(r2ptr);
			}

			PDPartition *p = new PDPartition;
			partitions.push_back(p);
			p->parentRectangle_ = rectangles[k];
			p->childRectangles_.push_back(r1ptr);
			p->childRectangles_.push_back(r2ptr);
			curRectangle->childPartitions_.push_back(p);
			r1ptr->parentPartitions_.push_back(p);
			r2ptr->parentPartitions_.push_back(p);
		}
	}

	//
	// set levels
	//
	root->setLevelRecursive(0);
	int numLevels = -1;
	for (unsigned k = 0; k < rectangles.size(); k++)
		if ((rectangles[k]->level_ + 1) > numLevels)
			numLevels = rectangles[k]->level_ + 1;

	std::cout << "#rectangle levels " << numLevels << std::endl;

	//
	// set all rectangles of size 1x1 (pixels) to the highest level
	//
	for (unsigned k = 0; k < rectangles.size(); k++)
		if ((rectangles[k]->width_ == 1) && (rectangles[k]->height_ == 1))
			rectangles[k]->level_ = numLevels - 1;

	//
	// organize level map of rectangles
	//
	for (unsigned k = 0; k < rectangles.size(); k++)
		rectanglesOfLevel[rectangles[k]->level_].push_back(rectangles[k]);

	std::cout << "#rectangles " << rectangles.size() << std::endl;

	// 
	Network *spn = factory.generateNetwork();

	// generate input layer
	ValLayer *inputLayer = (ValLayer*)factory.generateLayer(Layer::VAL_LAYER, width * height);

	// generate input layer
	Layer *distLayer;
	if (inputType == Node::GAUSSIAN_NODE)
		distLayer = factory.generateLayer(Layer::GAUSSIAN_LAYER, 0);
	else if (inputType == Node::INDICATOR_NODE)
		distLayer = factory.generateLayer(Layer::INDICATOR_LAYER, 0);
	else throw std::runtime_error("makeFullStructure(NetworkFactory &factory, NodeType inputType): wrong input type");

	for (unsigned k = 0; k < rectanglesOfLevel[numLevels - 1].size(); k++) {
		unsigned pixX = rectanglesOfLevel[numLevels - 1][k]->x_;
		unsigned pixY = rectanglesOfLevel[numLevels - 1][k]->y_;
		unsigned inputIdx = (unsigned)(pixY * width + pixX);

		for (size_t g = 0; g < numDistPerPixel; g++) {
			Node *distNode = factory.generateNode(inputType);
			distNode->connectChild((*inputLayer)(inputIdx));
			rectanglesOfLevel[numLevels - 1][k]->nodePtrs_[g] = distNode;
			distLayer->insertNode(distNode);
		}
	}

	spn->addLayer(inputLayer);
	spn->addLayer(distLayer);

	for (unsigned k = 1; k < (unsigned)numLevels; k++) {
		unsigned curLevelIdx = (unsigned)numLevels - 1 - k;

		SumLayer* sumLayer = (SumLayer*)factory.generateLayer(Layer::SUM_LAYER, 0);
		ProdLayer* prodLayer = (ProdLayer*)factory.generateLayer(Layer::PROD_LAYER, 0);

		for (unsigned l = 0; l < rectanglesOfLevel[curLevelIdx].size(); l++) {
			PDRectangle *curRectangle = rectanglesOfLevel[curLevelIdx][l];

			if ((curRectangle->width_ == 1) && (curRectangle->height_ == 1))
				throw std::runtime_error("makeFullStructure(NetworkFactory &factory, NodeType inputType): pixel in higher layer");;

			for (unsigned p = 0; p < curRectangle->childPartitions_.size(); p++) {
				PDPartition *curPartition = curRectangle->childPartitions_[p];
				PDRectangle *r1 = curPartition->childRectangles_[0];
				PDRectangle *r2 = curPartition->childRectangles_[1];

				for (auto s1 = r1->nodePtrs_.begin(); s1 != r1->nodePtrs_.end(); s1++) {
					for (auto s2 = r2->nodePtrs_.begin(); s2 != r2->nodePtrs_.end(); s2++) {
						Node *curProd = factory.generateNode(Node::PROD_NODE);
						curProd->connectChild(s1->second);
						curProd->connectChild(s2->second);
						curRectangle->prodNodePtrs_.push_back(curProd);
						prodLayer->insertNode(curProd);
					}
				}
			}

			unsigned numSums = 1;
			if ((k + 1) < (unsigned)numLevels)
				numSums = maxNumSumsPerRectangle;

			for (unsigned s = 0; s < numSums; s++) {
				Node *curSum = factory.generateNode(Node::SUM_NODE);
				for (unsigned p = 0; p < curRectangle->prodNodePtrs_.size(); p++)
					curSum->connectChild(curRectangle->prodNodePtrs_[p]);
				curRectangle->nodePtrs_[s] = curSum;
				sumLayer->insertNode(curSum);
			}
		}

		spn->addLayer(prodLayer);
		spn->addLayer(sumLayer);
	}

	return spn;
}


int main(int argc, char** argv) {

	//
	// get directory where to write results to
	//
	std::string resultPath;

	if (argc > 1) {
		resultPath = argv[1];
		struct stat info;

		if (stat(resultPath.c_str(), &info) != 0) {
			cout << resultPath << " not accessible" << endl;
			return -1;
		}
		else if ((info.st_mode & S_IFDIR) == 0) {
			cout << resultPath << " is not a directory" << endl;
			return -1;
		}

		resultPath = resultPath + "/";
	}
	else
		resultPath = "";
	

	unsigned numRuns = 100;
	std::default_random_engine generator;
	NetworkFactory factory;


	//
	// spn is the Poon & Domingos structure
	// copyspn is simulating the augmented SPN using correction weights
	//
	Network *spn;
	Network *copyspn;
	
	
	//
	// rectangle lengths to iterate over
	//
	UnsignedVec DRange;
	DRange.push_back(2);
	DRange.push_back(3);
	DRange.push_back(4);
	//DRange.push_back(5);

	
	//
	// Dirichlet prior params to iterate over
	//
	DoubleVec alphaRange;
	alphaRange.push_back(0.5);
	alphaRange.push_back(1);
	alphaRange.push_back(2.0);


	//
	// seed for the random generator
	//
	UnsignedVec seedVec(3,0);



	//
	// main loop: iterate over lengths/prior params
	//
	for (size_t DRangeCount = 0; DRangeCount < DRange.size(); DRangeCount++) {
		for (size_t aRangeCount = 0; aRangeCount < alphaRange.size(); aRangeCount++) {

			unsigned D = DRange[DRangeCount];
			double alpha = alphaRange[aRangeCount];
			
			UnsignedVec coarseRes = UnsignedVec(1, 1);
			UnsignedVec coarseResX, coarseResY;
			unsigned width = D;
			unsigned height = D;
			unsigned numDist = 2;
			unsigned numSums = 2;

			//
			// generate Poon & Domingos structure
			//
			spn = makeFullPDStructure(factory, width, height, coarseRes, coarseResX, coarseResY, numSums, numDist, Node::INDICATOR_NODE);
			
			//
			// set indicator values
			//
			if (spn->getLayer(0)->getType() != Layer::VAL_LAYER)
				throw std::runtime_error("first layer ought be a val layer");

			for (size_t k = 0; k < spn->getLayer(0)->getNumNodes(); k++) {
				ValNode* val = (ValNode*)spn->getLayer(0)->getNode(k);
				if (val->getNumParents() != 2)
					throw std::runtime_error("number of parents should be 2");
				for (size_t p = 0; p < 2; p++) {
					if (val->getParent(p)->getType() != Node::INDICATOR_NODE)
						throw std::runtime_error("parent should be indicator node");
					((IndicatorNode*)val->getParent(p))->setIndValue(val, (double)p);
				}
			}

			//
			// check structure
			//
			spn->performStructureCheck();
			if (!spn->isFeedforward())
				throw std::runtime_error("not feedforward");
			if (!spn->isComplete())
				throw std::runtime_error("not complete");
			if (!spn->isDecomposable())
				throw std::runtime_error("not decomposable");
						
			//
			// copy (and check structure of the copy)
			//
			copyspn = factory.copyNetwork(spn);			
			copyspn->performStructureCheck();
			if (!copyspn->isFeedforward())
				throw std::runtime_error("not feedforward");
			if (!copyspn->isComplete())
				throw std::runtime_error("not complete");
			if (!copyspn->isDecomposable())
				throw std::runtime_error("not decomposable");

			//
			// open file
			//
			string fileName = resultPath + "MPEresult_D" +
					unsigned2str((unsigned)DRangeCount) + "_alpha" +
					unsigned2str((unsigned)aRangeCount) + ".txt";
			
			std::ofstream outputFile;
			outputFile.open(fileName.c_str());
			if (!outputFile.is_open())
				throw std::runtime_error("unable to open file.");

			outputFile << std::fixed << std::setprecision(16);


			//
			// loop over runs 
			//
			for (unsigned i = 0; i < numRuns; i++) {

				std::cout << "D " << D << "   alpha " << alpha << "   run " << i << endl;

				//
				// random seed
				//
				seedVec[0] = DRangeCount;
				seedVec[1] = aRangeCount;
				seedVec[2] = i;
				std::seed_seq seed(seedVec.begin(), seedVec.end());
				generator.seed(seed);

				//
				// init sum-weights
				//
				for (size_t k = 0; k < spn->getNumLayers(); k++) {
					Layer* curLayer = spn->getLayer(k);
					Layer* copyCurLayer = copyspn->getLayer(k);
					if (curLayer->getType() == Layer::SUM_LAYER) {
						if ((copyCurLayer->getType() != Layer::SUM_LAYER) || (copyCurLayer->getNumNodes() != curLayer->getNumNodes()))
							throw std::runtime_error("copy inconsistency");

						for (size_t n = 0; n < curLayer->getNumNodes(); n++) {
							SumNode* sum = (SumNode*)curLayer->getNode(n);
							SumNode* copysum = (SumNode*)copyCurLayer->getNode(n);

							DoubleVec v = Dirichlet(alpha, (unsigned)sum->getNumChildren(), generator);
							for (size_t wc = 0; wc < v.size(); wc++)
								v[wc] = log(v[wc]);

							sum->setLogWeights(v);
							copysum->setLogWeights(v);
						}
					}
				}

				// set correction weights stemming from augmentation
				copyspn->setMPEcorrectionWeights();
				
				// "input" only consisting of nans, i.e. missing data
				DoubleVec input(D*D, SPN::nan);


				///////////////////////////////////////////
				//
				// MPE by max back-tracking
				//
				///////////////////////////////////////////

				//
				// Max backtracking in the original SPN, corresponds to MPE in 
				// augmented SPN with deterministic twin-weights
				//
				spn->setMaxMode();
				spn->setVals(input);
				spn->eval();
				spn->maxBacktracking();
				DoubleVec MPE_aug_det = spn->getVals();


				//
				// Max backtracking in SPN with correction weights, corresponds to MPE in 
				// augmented SPN with uniform twin-weights
				//
				copyspn->setMaxMode();
				copyspn->setVals(input);
				copyspn->eval();
				copyspn->maxBacktracking();
				DoubleVec MPE_aug_uni = copyspn->getVals();
				

				///////////////////////////////////////////
				//
				// MPE inference by exhaustive search
				//
				///////////////////////////////////////////

				UnsignedVec numStates(D*D, 2);
				UnsignedVec states(D*D, 0);


				//
				// exhaustive MPE in original SPN, i.e. MAP in augmented SPN (LVs marginalized)
				//
				
				// init state counter
				for (size_t k = 0; k < states.size(); k++)
					states[k] = 0;

				DoubleVec MPE_exhaustive_orig;
				double maxLL = -inf;
				bool overflow = false;
				unsigned long count = 0;

				// run over all state combinations, find maximum
				spn->setSumMode();
				while (!overflow) {
					spn->setVals(Unsigned2DoubleVec(states));
					spn->eval();
					double tmpLL = spn->getTopLayer()->getNode(0)->getLogVal();
					if (tmpLL > maxLL){
						maxLL = tmpLL;
						MPE_exhaustive_orig = Unsigned2DoubleVec(states);
					}
					overflow = incStateCounter(states, numStates);
					count++;
				}
				

				//
				// exhaustive MPE in augmented SPN, uniform twin-weights
				//

				// init state counter
				for (size_t k = 0; k < states.size(); k++)
					states[k] = 0;

				DoubleVec MPE_exhaustive_aug_uni;
				maxLL = -inf;
				overflow = false;
				count = 0;
				
				// run over all state combinations, find maximum
				copyspn->setMaxMode();
				while (!overflow) {
					copyspn->setVals(Unsigned2DoubleVec(states));
					copyspn->eval();
					double tmpLL = copyspn->getTopLayer()->getNode(0)->getLogVal();
					if (tmpLL > maxLL){
						maxLL = tmpLL;
						MPE_exhaustive_aug_uni = Unsigned2DoubleVec(states);
					}
					overflow = incStateCounter(states, numStates);
					count++;
				}


				//
				// exhaustive MPE in augmented SPN, deterministic twin-weights
				//

				// init state counter
				for (size_t k = 0; k < states.size(); k++)
					states[k] = 0;

				DoubleVec MPE_exhaustive_aug_det;
				maxLL = -inf;
				overflow = false;
				count = 0;
						
				// run over all state combinations, find maximum
				spn->setMaxMode();
				while (!overflow) {
					spn->setVals(Unsigned2DoubleVec(states));
					spn->eval();
					double tmpLL = spn->getTopLayer()->getNode(0)->getLogVal();
					if (tmpLL > maxLL){
						maxLL = tmpLL;
						MPE_exhaustive_aug_det = Unsigned2DoubleVec(states);
					}
					overflow = incStateCounter(states, numStates);
					count++;
				}


				//
				// evaluate likelihoods
				//
				DoubleVec LLorig(3, 0);
				DoubleVec LLaug_uni(3, 0);
				DoubleVec LLaug_det(3, 0);


				// in original SPN
				spn->setSumMode();

				spn->setVals(MPE_exhaustive_orig);
				spn->eval();
				LLorig[0] = spn->getTopLayer()->getNode(0)->getLogVal();
				
				spn->setVals(MPE_aug_uni);
				spn->eval();
				LLorig[1] = spn->getTopLayer()->getNode(0)->getLogVal();

				spn->setVals(MPE_aug_det);
				spn->eval();
				LLorig[2] = spn->getTopLayer()->getNode(0)->getLogVal();
							


				// in augmented SPN, uniform twin-weights
				copyspn->setMaxMode();				

				copyspn->setVals(MPE_exhaustive_aug_uni);
				copyspn->eval();
				LLaug_uni[0] = copyspn->getTopLayer()->getNode(0)->getLogVal();

				copyspn->setVals(MPE_aug_uni);
				copyspn->eval();
				LLaug_uni[1] = copyspn->getTopLayer()->getNode(0)->getLogVal();

				copyspn->setVals(MPE_aug_det);
				copyspn->eval();
				LLaug_uni[2] = copyspn->getTopLayer()->getNode(0)->getLogVal();
								
				

				// in augmented SPN, deterministic twin-weights
				spn->setMaxMode();

				spn->setVals(MPE_exhaustive_aug_det);
				spn->eval();
				LLaug_det[0] = spn->getTopLayer()->getNode(0)->getLogVal();
				
				spn->setVals(MPE_aug_uni);
				spn->eval();
				LLaug_det[1] = spn->getTopLayer()->getNode(0)->getLogVal();

				spn->setVals(MPE_aug_det);
				spn->eval();
				LLaug_det[2] = spn->getTopLayer()->getNode(0)->getLogVal();
								

				outputFile 
					<< LLorig[0] << " "
					<< LLorig[1] << " "
					<< LLorig[2] << " "
					<< LLaug_uni[0] << " "
					<< LLaug_uni[1] << " "
					<< LLaug_uni[2] << " "
					<< LLaug_det[0] << " "
					<< LLaug_det[1] << " "
					<< LLaug_det[2] << endl;
			}

			outputFile.close();
		}
	}

	return 0;
}