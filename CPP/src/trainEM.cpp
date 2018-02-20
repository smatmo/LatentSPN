/*
trainEM.cpp

Loads an SPN and data and runs the EM algorithm

Robert Peharz
October 2016
*/


#include "Network.hpp"
#include "NetworkFactory.hpp"
#include "DataSource.hpp"
#include "Utils.hpp"
#include "FileIO.hpp"
#include "EM.hpp"
#include <iostream>
#include <string>
#include <ctime>


using namespace SPN;
using namespace std;

int main(int argc, char** argv) {

	string useText =
		string("Run EM algorithm for an SPN\n") +
		string("Useage: ") + string(argv[0]) + string(" modelfile datafile outputfile\n") +
		string("\n") +
		string("Parameters: \n") +
		string("\n") +
		string("valDataFile                    = "" \n") +
		string("numIter                        = 100 \n") +
		string("stop_relLikelihoodChange       = 1e-6 \n") +
		string("earlyStoppingK                 = 3 \n") +
		string("updateWeights                  = 1 \n") +
		string("updateMeans                    = 1 \n") +
		string("updateSigmas                   = 1 \n") +
		string("randSeed (0 = no init)         = 0 \n") +
		string("minSigma                       = 1e-6 \n") +
		string("historyFile                    = "" \n") +
		string("PDformat                       = 0 \n") +
		string("width (set >0, if PDformat)    = -1 \n") +
		string("height (set >0, if PDformat)   = -1 \n") +
		string("\n");

	// params:
	//
	// modelfile: file specifiying an SPN,
	// either in our custom format or the format used by Poon & Domingos
	//
	// datafile: file containing data set
	//
	// outputfile: file to store the trained SPN
	//
	//
	// optional params:
	// 
	// valDataFile: data file for validation data -> used for early stopping
	// 
	// numIter: number of iterations for the EM algorithm
	//
	// stop_relLikelihoodChange: minimal relative likelihood change before EM stops 
	// 
	// earlyStoppingK: number of times validation likelihood is allowed to decresae 
	//
	// updateWeights, updateMeans, updateSigmas: indicating which types of parameters 
	// are updated
	//
	// randSeed: seed for random generator. If 0, take stored parameters
	//
	// minSigma: lower bound for Gaussian sigmas
	//
	// historyFile: file where training history (likelihoods) are stored
	//
	// PDformat: use the format of Poon & Domingos?
	//
	// width, height: only relavant for PDformat

	string valDataFile = "";
	size_t numIter = 100;
	double stop_relLikelihoodChange = 1e-6;
	unsigned earlyStoppingK = 3;
	bool updateWeights = true;
	bool updateMeans = true;
	bool updateSigmas = true;
	double minSigma = 1e-6;
	string historyFile = "";
	unsigned randSeed = 0;
	unsigned PDformat = 0;
	int width = -1;
	int height = -1;
	
	map<string, string> optionMap;
	map<string, string>::iterator it;
	stringstream inputStream;

	if (argc < 4) {
		std::cout << useText;
		return -1;
	}

	// parse params
	if (!parseOptions(argc, argv, 4, optionMap)) {
		std::cout << useText;
		return -1;
	}
	else {
		for (it = optionMap.begin(); it != optionMap.end(); it++) {
			if (it->first == "valDataFile") {
				valDataFile = it->second;
			}
			else if (it->first == "numIter") {
				inputStream.clear();
				inputStream.str(it->second);
				inputStream >> numIter;
				if (inputStream.fail()) {
					std::cout << useText;
					return -1;
				}
			}
			else if (it->first == "stop_relLikelihoodChange") {
				inputStream.clear();
				inputStream.str(it->second);
				inputStream >> stop_relLikelihoodChange;
				if (inputStream.fail()) {
					std::cout << useText;
					return -1;
				}
			}
			else if (it->first == "earlyStoppingK") {
				inputStream.clear();
				inputStream.str(it->second);
				inputStream >> earlyStoppingK;
				if (inputStream.fail()) {
					std::cout << useText;
					return -1;
				}
			}
			else if (it->first == "updateWeights") {
				inputStream.clear();
				inputStream.str(it->second);
				inputStream >> updateWeights;
				if (inputStream.fail()) {
					std::cout << useText;
					return -1;
				}
			}
			else if (it->first == "updateMeans") {
				inputStream.clear();
				inputStream.str(it->second);
				inputStream >> updateMeans;
				if (inputStream.fail()) {
					std::cout << useText;
					return -1;
				}
			}
			else if (it->first == "updateSigmas") {
				inputStream.clear();
				inputStream.str(it->second);
				inputStream >> updateSigmas;
				if (inputStream.fail()) {
					std::cout << useText;
					return -1;
				}
			}
			else if (it->first == "randSeed") {
				inputStream.clear();
				inputStream.str(it->second);
				inputStream >> randSeed;
				if (inputStream.fail()) {
					std::cout << useText;
					return -1;
				}
			}
			else if (it->first == "minSigma") {
				inputStream.clear();
				inputStream.str(it->second);
				inputStream >> minSigma;
				if (inputStream.fail()) {
					std::cout << useText;
					return -1;
				}
			}
			else if (it->first == "historyFile") {
				historyFile = it->second;
			}
			else if (it->first == "PDformat") {
				inputStream.clear();
				inputStream.str(it->second);
				inputStream >> PDformat;
				if (inputStream.fail()) {
					std::cout << useText;
					return -1;
				}
			}
			else if (it->first == "width") {
				inputStream.clear();
				inputStream.str(it->second);
				inputStream >> width;
				if (inputStream.fail()) {
					std::cout << useText;
					return -1;
				}
			}
			else if (it->first == "height") {
				inputStream.clear();
				inputStream.str(it->second);
				inputStream >> height;
				if (inputStream.fail()) {
					std::cout << useText;
					return -1;
				}
			}
			else {
				std::cout << "unknown parameter " << it->first << std::endl;
				std::cout << useText;
				return -1;
			}
		}
	}

	if ((PDformat) && ((width <= 0) || (height <= 0))) {
		std::cout << "width and height need to be set for PDformat." << std::endl;
		std::cout << useText;
		return -1;
	}

	// echo params
	cout << "valDataFile                = " << valDataFile << endl;
	cout << "numIter                    = " << numIter << endl;
	cout << "stop_relLikelihoodChange   = " << stop_relLikelihoodChange << endl;
	cout << "earlyStoppingK             = " << earlyStoppingK << endl;
	cout << "updateWeights              = " << updateWeights << endl;
	cout << "updateMeans                = " << updateMeans << endl;
	cout << "updateSigmas               = " << updateSigmas << endl;
	cout << "randSeed                   = " << randSeed << endl;
	cout << "minSigma                   = " << minSigma << endl;
	cout << "historyFile                = " << historyFile << endl;
	cout << "PDformat                   = " << PDformat << endl;
	cout << "width                      = " << width << endl;
	cout << "height                     = " << height << endl;
	cout << endl;

	try {
		NetworkFactory factory;
		Network* spn = factory.generateNetwork();

		// load network
		cout << "Load network" << endl;
		if (PDformat) {			
			loadNetworkPDformat(spn, string(argv[1]), width, height, factory);
		}
		else {
			loadNetwork(spn, argv[1], factory);
		}

		// remove singleton parents
		NodeVec dummNV;
		LayerVec dummLV;
		spn->removeSingletonParents(dummNV, dummLV);
		
		// simple structure check
		spn->performStructureCheck();
		
		if (!spn->isFeedforward()) {			
			std::cout << "not Feedforward " << std::endl;
			throw std::runtime_error("bad structure");
		}

		if (!spn->isComplete()) {
			std::cout << "not Complete " << std::endl;
			throw std::runtime_error("bad structure");
		}

		if (!spn->isDecomposable()) {
			std::cout << "not Decomposable " << std::endl;
			throw std::runtime_error("bad structure");
		}
		
		// load data
		cout << "Load data" << endl;
		DataSource trainData, valData;
		DataSource *valDataPtr = NULL;
		trainData.readFile(argv[2]);
		if (valDataFile != "") {
			valData.readFile(valDataFile);
			valDataPtr = &valData;
		}
		
		EM emlearner;
		
		// set params
		emlearner.maxIter_ = numIter;
		emlearner.stopping_objMinRelChange_ = stop_relLikelihoodChange;
		emlearner.updateGaussianMeans_ = updateMeans;
		emlearner.updateGaussianSigmas_ = updateSigmas;
		emlearner.updateWeights_ = updateWeights;
		emlearner.minSigma_ = minSigma;
		emlearner.maxValLLNotIncreased_ = earlyStoppingK;

		// init
		if (randSeed) {
			srand(randSeed);
			if (updateWeights)
				spn->initWeightsRand();
			if (updateMeans)
				spn->setGaussianMeansRand(-1.0, 1.0);
			if (updateSigmas)
				spn->setGaussianSigmasRand(0.1, 1.0);
		}

		/////////////////////////////////////////////////////
		cout << "start training" << endl;
		
		clock_t start = clock();
		emlearner.learn(spn, &trainData, valDataPtr);
		clock_t end = clock();		
		/////////////////////////////////////////////////////

		double emtime = (double)(end - start) / CLOCKS_PER_SEC * 1000.0;
		
		// save network
		saveNetwork(spn, argv[3]);

		// save history
		if (historyFile != "") {
			DoubleVec2 history;
			history.push_back(DoubleVec(1, emtime));
			history.push_back(emlearner.trainHist_);
			if (valDataFile != "") 
				history.push_back(emlearner.valHist_);
			writeVectorList(historyFile, history);
		}
	}
	catch (std::exception& e) { std::cout << e.what() << '\n'; }

}
