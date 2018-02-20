
#include "FileIO.hpp"

namespace SPN {

	bool findLayerIdx(Network *network, Node *node, size_t &idx) {
		bool foundLayer = false;

		for (unsigned n = 0; n < network->getNumLayers(); n++) {
			if (network->getLayer(n)->containsNode(node)) {
				foundLayer = true;
				idx = n;
				break;
			}
		}

		return foundLayer;
	}

	unsigned getInputIdx(int x, int y, int c, int width, int height, int numChannels, int pixelOrder) {
		unsigned idx;

		if (pixelOrder == 1)
			idx = (unsigned)(y * numChannels * width + c * width + x);
		else
			idx = (unsigned)(x * numChannels * height + c * height + y);

		return idx;
	}

	void saveNetwork(Network *network, std::string fileName) {
		std::ofstream outputFile;
		Node *node, *child;
		size_t layerIdx = 0;
		size_t nodeIdx = 0;

		outputFile.open(fileName.c_str());
		if (!outputFile.is_open())
			throw std::runtime_error("saveNetwork(const Network &network, std::string fileName): unable to open file.");

		outputFile << "SPN_NETWORK" << std::endl;
		outputFile << "NumLayers: " << network->getNumLayers() << std::endl;
		outputFile << "TYPES: ";
		for (unsigned k = 0; k < network->getNumLayers(); k++)
			outputFile << network->getLayer(k)->getType() << " ";
		outputFile << std::endl << std::endl;

		for (unsigned k = 0; k < network->getNumLayers(); k++) {
			outputFile << "LAYER " << k << ": TYPE " << network->getLayer(k)->getType() << std::endl;
			outputFile << "NumNodes: " << network->getLayer(k)->getNumNodes() << std::endl;

			if (network->getLayer(k)->getType() != Layer::VAL_LAYER) {
				for (unsigned l = 0; l < network->getLayer(k)->getNumNodes(); l++) {
					node = network->getLayer(k)->getNode(l);
					outputFile << node->getNumChildren() << ": ";
					for (unsigned m = 0; m < node->getNumChildren(); m++) {
						child = node->getChild(m);

						if (!findLayerIdx(network, child, layerIdx)) {
							outputFile << std::endl << "malformatted network" << std::endl;
							outputFile.close();
							throw std::runtime_error("error in saveNetwork(Network &network, std::string fileName): malformatted network.");
						}

						nodeIdx = network->getLayer(layerIdx)->getNodeIdx(child);
						outputFile << layerIdx << " " << nodeIdx;

						// type dependent extra information 
						if (node->getType() == Node::SUM_NODE)
							outputFile << " " << exp(((SumNode*)node)->getLogWeight(child));
						if (node->getType() == Node::GAUSSIAN_NODE)
							outputFile << " " << ((GaussianNode*)node)->getMean(child) << " " << ((GaussianNode*)node)->getSigma(child);
						if (node->getType() == Node::INDICATOR_NODE)
							outputFile << " " << ((IndicatorNode*)node)->getIndValue(child);

						if (m < node->getNumChildren() - 1)
							outputFile << "|";
					}
					outputFile << std::endl;
				}
			}
			if (k < network->getNumLayers() - 1)
				outputFile << std::endl;
		}
		outputFile.close();
	}

	void loadNetwork(Network *network, std::string fileName, NetworkFactory &factory) {
		std::ifstream inputFile;
		std::string inputLine, tmpString;
		std::stringstream inputStream, tmpStream;
		unsigned numLayers, numNodes, numChildren;
		unsigned layerIdx, nodeIdx;
		int layerType = -1, nodeType = -1;
		size_t pos;
		Layer *layer;
		double weight, mean, sigma, indValue;

		network->clear();

		inputFile.open(fileName.c_str());
		if (!inputFile.is_open())
			throw std::runtime_error("loadNetwork(const Network &network, std::string fileName, NetworkFactory &factory): unable to open file.");

		getline(inputFile, inputLine);
		if (inputLine.find_last_not_of(" \t\f\v\n\r") != std::string::npos)
			inputLine.erase(inputLine.find_last_not_of(" \t\f\v\n\r") + 1);
		if ((inputLine.compare(0, 11, "SPN_NETWORK") != 0) || (inputLine.size() > 11)) {
			inputFile.close();
			throw std::runtime_error(
				"NetworkFactory::loadNetwork(const Network &network, std::string fileName): "
				"Header of file has to be \"SPN_NETWORK\".");
		}

		getline(inputFile, inputLine);
		inputLine.erase(0, 10);
		inputStream.clear();
		inputStream.str(inputLine);
		inputStream >> numLayers;
		if (inputStream.fail()) {
			inputFile.close();
			throw std::runtime_error(
				"NetworkFactory::loadNetwork(const Network &network, std::string fileName, NetworkFactory &factory): cannot read number of layers.");
		}

		// generate layers
		getline(inputFile, inputLine);
		inputLine.erase(0, 6);
		inputStream.clear();
		inputStream.str(inputLine);
		for (unsigned k = 0; k < numLayers; k++) {
			inputStream >> layerType;
			if (inputStream.fail()) {
				inputFile.close();
				throw std::runtime_error(
					"NetworkFactory::loadNetwork(const Network &network, std::string fileName, NetworkFactory &factory): cannot read layer type.");
			}

			layer = factory.generateLayer(layerType, 0);
			network->addLayer(layer);
		}

		for (unsigned k = 0; k < numLayers; k++) {
			// empty line
			getline(inputFile, inputLine);

			// check type of next layer
			getline(inputFile, inputLine);
			pos = inputLine.find_first_of(":");
			inputLine.erase(0, pos + 6);
			inputStream.clear();
			inputStream.str(inputLine);
			inputStream >> layerType;
			if (inputStream.fail()) {
				inputFile.close();
				throw std::runtime_error("NetworkFactory::loadNetwork(const Network &network, std::string fileName): cannot read layer type.");
			}
			if (layerType != network->getLayer(k)->getType()) {
				inputFile.close();
				throw std::runtime_error(
					"NetworkFactory::loadNetwork(const Network &network, std::string fileName): "
					"inconsistent layer types.");
			}

			// find out number of nodes in this layer
			getline(inputFile, inputLine);
			inputLine.erase(0, 9);
			inputStream.clear();
			inputStream.str(inputLine);
			inputStream >> numNodes;
			if (inputStream.fail()) {
				inputFile.close();
				throw std::runtime_error(
					"NetworkFactory::loadNetwork(const Network &network, std::string fileName): "
					"cannot read number of nodes.");
			}

			// generate nodes
			if (layerType == Layer::SUM_LAYER)
				nodeType = Node::SUM_NODE;
			if (layerType == Layer::PROD_LAYER)
				nodeType = Node::PROD_NODE;
			if (layerType == Layer::VAL_LAYER)
				nodeType = Node::VAL_NODE;
			if (layerType == Layer::INDICATOR_LAYER)
				nodeType = Node::INDICATOR_NODE;
			if (layerType == Layer::GAUSSIAN_LAYER)
				nodeType = Node::GAUSSIAN_NODE;

			network->getLayer(k)->insertNodes(factory.generateNodes(nodeType, numNodes));
			if (layerType == Layer::VAL_LAYER)
				continue;

			for (unsigned l = 0; l < numNodes; l++) {
				getline(inputFile, inputLine);
				inputStream.clear();
				inputStream.str(inputLine);				

				// find out num children
				inputStream >> numChildren;
				if (inputStream.fail()) {
					inputFile.close();
					throw std::runtime_error(
						"NetworkFactory::loadNetwork(const Network &network, std::string fileName): "
						"cannot read number of children.");
				}

				pos = inputLine.find_first_of(":");
				inputLine.erase(0, pos + 1);
				inputStream.clear();
				inputStream.str(inputLine);

				for (unsigned m = 0; m < numChildren; m++) {
					inputStream >> layerIdx;
					if (inputStream.fail()) {
						inputFile.close();
						throw std::runtime_error(
							"NetworkFactory::loadNetwork(const Network &network, std::string fileName): "
							"cannot read child layer idx.");
					}
					inputStream >> nodeIdx;
					if (inputStream.fail()) {
						inputFile.close();
						throw std::runtime_error(
							"NetworkFactory::loadNetwork(const Network &network, std::string fileName): "
							"cannot read child node idx.");
					}

					(*network)(k, l)->connectChild((*network)(layerIdx, nodeIdx));

					// type dependent extra information 
					if (nodeType == Node::SUM_NODE) {
						inputStream >> weight;
						if (inputStream.fail()) {
							inputFile.close();
							throw std::runtime_error(
								"NetworkFactory::loadNetwork(const Network &network, std::string fileName): "
								"cannot read weight.");
						}
						((SumNode*)((*network)(k, l)))->setLogWeight((*network)(layerIdx, nodeIdx), log(weight));
					}
					if (nodeType == Node::GAUSSIAN_NODE) {
						inputStream >> mean;
						if (inputStream.fail()) {
							inputFile.close();
							throw std::runtime_error(
								"NetworkFactory::loadNetwork(const Network &network, std::string fileName): "
								"cannot read mean.");
						}
						inputStream >> sigma;
						if (inputStream.fail()) {
							inputFile.close();
							throw std::runtime_error(
								"NetworkFactory::loadNetwork(const Network &network, std::string fileName): "
								"cannot read sigma.");
						}
						((GaussianNode*)((*network)(k, l)))->setMean((*network)(layerIdx, nodeIdx), mean);
						((GaussianNode*)((*network)(k, l)))->setSigma((*network)(layerIdx, nodeIdx), sigma);
					}
					if (nodeType == Node::INDICATOR_NODE) {
						inputStream >> indValue;
						if (inputStream.fail()) {
							inputFile.close();
							throw std::runtime_error(
								"NetworkFactory::loadNetwork(const Network &network, std::string fileName): "
								"cannot read indicator value.");
						}
						((IndicatorNode*)((*network)(k, l)))->setIndValue((*network)(layerIdx, nodeIdx), indValue);
					}

					// prepare string to read next child
					if (m < numChildren - 1) {
						pos = inputLine.find_first_of("|");
						inputLine.erase(0, pos + 1);
						inputStream.clear();
						inputStream.str(inputLine);
					}
				}
			}
		}

		inputFile.close();
	}

	void loadNetworkPDformat(Network *network, std::string fileName, int xDim, int yDim, NetworkFactory &factory) {
		std::ifstream inputFile;
		std::string inputLine, tmpString;
		std::stringstream inputStream;
		int x1, y1, x2, y2;
		int numNodes;

		PDRectangle* rootRectangle;
		std::vector<PDRectangle* > rectangles;
		std::vector<PDPartition* > partitions;
		std::map<PDRectangle, PDRectangle* > rectangleMap;
		std::map<unsigned, std::vector<PDRectangle*> > rectanglesOfLevel;

		PDRectangle tmpR;
		PDRectangle *rPtr, *childR1ptr, *childR2ptr;
		PDPartition *pPtr;
		long int childRegionId1, childRegionId2;
		unsigned childNodeIdx1, childNodeIdx2;

		size_t pos, ePos;

		unsigned inputIdx;

		double mean, weight, weightSum;

		ValLayer *inputLayer = (ValLayer*)factory.generateLayer(Layer::VAL_LAYER, xDim * yDim);
		GaussianLayer *gaussLayer;
		SumLayer *sumLayer;
		ProdLayer *prodLayer;
		GaussianNode *gaussNode;
		SumNode *curSumNode;
		Node *childNode1, *childNode2;
		ProdNode *prodNode;
		std::vector<ProdNode*> childProdNodes;
		DoubleVec weights;

		int maxLevel;

		std::map<unsigned, Node* >::iterator nodeIt;

		rectangles.clear();
		partitions.clear();
		rectangleMap.clear();
		network->clear();
		rootRectangle = NULL;

		inputFile.open(fileName.c_str());
		if (!inputFile.is_open())
			throw std::runtime_error(
			"NetworkFactory::loadNetworkPDformat(Network &network, std::string fileName): "
			"unable to open file.");

		while (!inputFile.eof()) {
			// read the <REGION> line
			getline(inputFile, inputLine);
			if (inputLine.find_last_not_of(" \t\f\v\n\r") == std::string::npos)
				continue;
			else
				inputLine.erase(inputLine.find_last_not_of(" \t\f\v\n\r") + 1);
			if ((inputLine.compare(0, 8, "<REGION>") != 0) || (inputLine.size() > 8)) {
				inputFile.close();
				throw std::runtime_error(
					"NetworkFactory::loadNetworkPDformat(Network &network, std::string fileName): "
					"Expected to read \"<REGION>\".");
			}

			// read coordinates
			getline(inputFile, inputLine);
			inputStream.clear();
			inputStream.str(inputLine);
			inputStream >> x1;
			inputStream >> x2;
			inputStream >> y1;
			inputStream >> y2;
			if (inputStream.fail()) {
				inputFile.close();
				throw std::runtime_error(
					"NetworkFactory::loadNetworkPDformat(Network &network, std::string fileName): "
					"malformatted input (coordinates).");
			}

			tmpR.x_ = x1;
			tmpR.width_ = x2 - x1;
			tmpR.y_ = y1;
			tmpR.height_ = y2 - y1;
			if (rectangleMap.count(tmpR) == 0) {
				rPtr = new PDRectangle;
				*rPtr = tmpR;
				rectangleMap[tmpR] = rPtr;
				rectangles.push_back(rPtr);
			}
			else
				rPtr = rectangleMap[tmpR];

			if ((tmpR.x_ == 0) && (tmpR.y_ == 0) && (tmpR.width_ == (unsigned) xDim) && (tmpR.height_ == (unsigned) yDim)) {
				if (rootRectangle == NULL)
					rootRectangle = rPtr;
				else
					throw std::runtime_error(
					"NetworkFactory::loadNetworkPDformat(Network &network, std::string fileName): "
					"Found 2nd root rectangle.");
			}

			// read the <TYPE> line
			getline(inputFile, inputLine);
			if (inputLine.find_last_not_of(" \t\f\v\n\r") != std::string::npos)
				inputLine.erase(inputLine.find_last_not_of(" \t\f\v\n\r") + 1);
			if ((inputLine.compare(0, 6, "<TYPE>") != 0) || (inputLine.size() > 6)) {
				inputFile.close();
				throw std::runtime_error(
					"NetworkFactory::loadNetworkPDformat(Network &network, std::string fileName): "
					"Expected to read \"<TYPE>\".");
			}

			// read number of "types"
			getline(inputFile, inputLine);
			inputStream.clear();
			inputStream.str(inputLine);
			inputStream >> numNodes;
			if (inputStream.fail()) {
				inputFile.close();
				throw std::runtime_error(
					"NetworkFactory::loadNetworkPDformat(Network &network, std::string fileName): "
					"malformatted input (numTypes).");
			}

			if ((rPtr->width_ == 1) && (rPtr->height_ == 1)) {
				// Pixel: just read the dummy lines
				for (int k = 0; k < numNodes; k++)
					getline(inputFile, inputLine);
			}
			else {
				// larger region
				for (int k = 0; k < numNodes; k++) {
					getline(inputFile, inputLine);
					pos = inputLine.find_first_of("<");
					if (pos == std::string::npos)
						continue;

					if (rPtr->nodePtrs_.count((unsigned)k) == 0) {
						curSumNode = (SumNode*)factory.generateNode(Node::SUM_NODE);
						rPtr->nodePtrs_[(unsigned)k] = curSumNode;
					}
					else
						curSumNode = (SumNode*)rPtr->nodePtrs_[(unsigned)k];

					childProdNodes.clear();
					weights.clear();
					weightSum = 0.0;
					// add children
					do { // while (pos != std::string::npos)
						ePos = inputLine.find_first_of(">");
						tmpString = inputLine.substr(pos + 1, ePos - pos - 1);
						// ...4 15>: -> ePos+2
						inputLine.erase(0, ePos + 2);

						inputStream.clear();
						inputStream.str(tmpString);
						inputStream >> childRegionId1;
						inputStream >> childRegionId2;
						inputStream >> childNodeIdx1;
						inputStream >> childNodeIdx2;
						if (inputStream.fail()) {
							inputFile.close();
							throw std::runtime_error(
								"NetworkFactory::loadNetworkPDformat(Network &network, std::string fileName): "
								"malformatted input (child regions).");
						}

						// get the child rectangles
						// 1
						y2 = childRegionId1 % yDim + 1;
						childRegionId1 = (childRegionId1 - y2 + 1) / yDim;
						y1 = childRegionId1 % yDim;
						childRegionId1 = (childRegionId1 - y1) / yDim;
						x2 = childRegionId1 % xDim + 1;
						x1 = (childRegionId1 - x2 + 1) / xDim;

						tmpR.x_ = x1;
						tmpR.y_ = y1;
						tmpR.width_ = x2 - x1;
						tmpR.height_ = y2 - y1;
						if (rectangleMap.count(tmpR) == 0) {
							childR1ptr = new PDRectangle;
							*childR1ptr = tmpR;
							rectangleMap[tmpR] = childR1ptr;
							rectangles.push_back(childR1ptr);
						}
						else
							childR1ptr = rectangleMap[tmpR];

						// 2
						y2 = childRegionId2 % yDim + 1;
						childRegionId2 = (childRegionId2 - y2 + 1) / yDim;
						y1 = childRegionId2 % yDim;
						childRegionId2 = (childRegionId2 - y1) / yDim;
						x2 = childRegionId2 % xDim + 1;
						x1 = (childRegionId2 - x2 + 1) / xDim;

						tmpR.x_ = x1;
						tmpR.y_ = y1;
						tmpR.width_ = x2 - x1;
						tmpR.height_ = y2 - y1;
						if (rectangleMap.count(tmpR) == 0) {
							childR2ptr = new PDRectangle;
							*childR2ptr = tmpR;
							rectangleMap[tmpR] = childR2ptr;
							rectangles.push_back(childR2ptr);
						}
						else
							childR2ptr = rectangleMap[tmpR];

						// make partition
						pPtr = new PDPartition;
						partitions.push_back(pPtr);
						pPtr->parentRectangle_ = rPtr;
						pPtr->childRectangles_.push_back(childR1ptr);
						pPtr->childRectangles_.push_back(childR2ptr);
						rPtr->childPartitions_.push_back(pPtr);
						childR1ptr->parentPartitions_.push_back(pPtr);
						childR2ptr->parentPartitions_.push_back(pPtr);

						//std::cout << rPtr->x_ << ", " << rPtr->width_ << "   " << rPtr->y_ << ", " << rPtr->height_ << std::endl;
						//std::cout << childR1ptr->x_ << ", " << childR1ptr->width_ << "   " << childR1ptr->y_ << ", " <<childR1ptr->height_ << std::endl;
						//std::cout << childR2ptr->x_ << ", " << childR2ptr->width_ << "   " << childR2ptr->y_ << ", " <<childR2ptr->height_ << std::endl;
						//std::cout << std::endl;

						if (childR1ptr->x_ == childR2ptr->x_) {
							if (childR1ptr->width_ != childR2ptr->width_)
								throw std::runtime_error(
								"NetworkFactory::loadNetworkPDformat(Network &network, std::string fileName): "
								"widths are not correct.");
							if (childR1ptr->height_ + childR2ptr->height_ != rPtr->height_)
								throw std::runtime_error(
								"NetworkFactory::loadNetworkPDformat(Network &network, std::string fileName): "
								"sum heights not correct.");
							if (childR2ptr->y_ != childR1ptr->y_ + childR1ptr->height_)
								throw std::runtime_error(
								"NetworkFactory::loadNetworkPDformat(Network &network, std::string fileName): "
								"y coordinates not correct.");
						}
						else if (childR1ptr->y_ == childR2ptr->y_) {
							if (childR1ptr->height_ != childR2ptr->height_)
								throw std::runtime_error(
								"NetworkFactory::loadNetworkPDformat(Network &network, std::string fileName): "
								"heights are not correct.");
							if (childR1ptr->width_ + childR2ptr->width_ != rPtr->width_)
								throw std::runtime_error(
								"NetworkFactory::loadNetworkPDformat(Network &network, std::string fileName): "
								"sum widths not correct.");
							if (childR2ptr->x_ != childR1ptr->x_ + childR1ptr->width_)
								throw std::runtime_error(
								"NetworkFactory::loadNetworkPDformat(Network &network, std::string fileName): "
								"x coordinates not correct.");
						}
						else
							throw std::runtime_error(
							"NetworkFactory::loadNetworkPDformat(Network &network, std::string fileName): "
							"coordinates are not correct.");

						// get the child nodes
						if (childR1ptr->nodePtrs_.count(childNodeIdx1) == 0) {
							if ((childR1ptr->width_ == 1) && (childR1ptr->height_ == 1))
								childNode1 = factory.generateNode(Node::GAUSSIAN_NODE);
							else
								childNode1 = factory.generateNode(Node::SUM_NODE);
							childR1ptr->nodePtrs_[childNodeIdx1] = childNode1;
						}
						else
							childNode1 = childR1ptr->nodePtrs_[childNodeIdx1];

						if (childR2ptr->nodePtrs_.count(childNodeIdx2) == 0) {
							if ((childR2ptr->width_ == 1) && (childR2ptr->height_ == 1))
								childNode2 = factory.generateNode(Node::GAUSSIAN_NODE);
							else
								childNode2 = factory.generateNode(Node::SUM_NODE);
							childR2ptr->nodePtrs_[childNodeIdx2] = childNode2;
						}
						else
							childNode2 = childR2ptr->nodePtrs_[childNodeIdx2];

						// make a product node and connect
						prodNode = (ProdNode*)factory.generateNode(Node::PROD_NODE);
						curSumNode->connectChild(prodNode);
						prodNode->connectChild(childNode1);
						prodNode->connectChild(childNode2);
						childProdNodes.push_back(prodNode);

						// read weight
						inputStream.clear();
						if (inputLine.find_first_of(":") == std::string::npos)
							inputStream.str(inputLine);
						else
							inputStream.str(inputLine.substr(0, inputLine.find_first_of(":")));
						inputStream >> weight;
						if (inputStream.fail()) {
							inputFile.close();
							throw std::runtime_error(
								"NetworkFactory::loadNetworkPDformat(Network &network, std::string fileName): "
								"malformatted input (weight).");
						}
						weights.push_back(weight);
						weightSum += weight;

						pos = inputLine.find_first_of("<");
					} while (pos != std::string::npos);

					// set weights
					for (unsigned l = 0; l < childProdNodes.size(); l++)
						curSumNode->setLogWeight(childProdNodes[l], log(weights[l] / weightSum));
				}
			}

			// read the </TYPE> line
			getline(inputFile, inputLine);
			if (inputLine.find_last_not_of(" \t\f\v\n\r") != std::string::npos)
				inputLine.erase(inputLine.find_last_not_of(" \t\f\v\n\r") + 1);
			if ((inputLine.compare(0, 7, "</TYPE>") != 0) || (inputLine.size() > 7)) {
				inputFile.close();
				throw std::runtime_error(
					"NetworkFactory::loadNetworkPDformat(Network &network, std::string fileName): "
					"Expected to read \"</TYPE>\".");
			}

			// for pixel region, generate Gaussians
			if ((rPtr->width_ == 1) && (rPtr->height_ == 1)) {
				getline(inputFile, inputLine);
				if (inputLine.compare(0, 6, "<MEAN>") != 0) {
					inputFile.close();
					throw std::runtime_error(
						"NetworkFactory::loadNetworkPDformat(Network &network, std::string fileName): "
						"Expected to read \"<MEAN>\".");
				}
				pos = inputLine.find_first_of(":");
				inputLine.erase(0, pos + 1);
				inputStream.clear();
				inputStream.str(inputLine);
				for (int k = 0; k < numNodes; k++) {
					if (rPtr->nodePtrs_.count((unsigned)k) == 0) {
						gaussNode = (GaussianNode*)factory.generateNode(Node::GAUSSIAN_NODE);
						rPtr->nodePtrs_[(unsigned)k] = gaussNode;
					}
					else
						gaussNode = (GaussianNode*)rPtr->nodePtrs_[(unsigned)k];

					inputIdx = getInputIdx(rPtr->x_, rPtr->y_, 0, xDim, yDim, 1, 1);
					gaussNode->connectChild(inputLayer->getNode(inputIdx));
					inputStream >> mean;
					if (inputStream.fail()) {
						inputFile.close();
						throw std::runtime_error(
							"NetworkFactory::loadNetworkPDformat(Network &network, std::string fileName): "
							"malformatted input (mean).");
					}
					gaussNode->setMean(inputLayer->getNode(inputIdx), mean);
					gaussNode->setSigma(inputLayer->getNode(inputIdx), 1);
				}
				// read the <CNT> line
				getline(inputFile, inputLine);
			}

			// read the </REGION> line
			getline(inputFile, inputLine);
			if (inputLine.find_last_not_of(" \t\f\v\n\r") != std::string::npos)
				inputLine.erase(inputLine.find_last_not_of(" \t\f\v\n\r") + 1);
			if ((inputLine.compare(0, 9, "</REGION>") != 0) || (inputLine.size() > 9)) {
				inputFile.close();
				throw std::runtime_error(
					"NetworkFactory::loadNetworkPDformat(Network &network, std::string fileName): "
					"Expected to read \"</REGION>\".");
			}
		}
		inputFile.close();

		if (rootRectangle == NULL)
			throw std::runtime_error(
			"NetworkFactory::loadNetworkPDformat(Network &network, std::string fileName): "
			"root not found.");
		rootRectangle->setLevelRecursive(0);

		// find the max level of non-pixel rectangles
		maxLevel = 0;
		for (unsigned k = 0; k < rectangles.size(); k++)
			if (((rectangles[k]->width_ > 1) || (rectangles[k]->height_ > 1)) && (rectangles[k]->level_ > maxLevel))
				maxLevel = rectangles[k]->level_;

		maxLevel++;

		// set all pixel rectangles to the maxLevel
		for (unsigned k = 0; k < rectangles.size(); k++)
			if ((rectangles[k]->width_ == 1) && (rectangles[k]->height_ == 1))
				rectangles[k]->level_ = maxLevel;

		// organize level map of rectangles
		for (unsigned k = 0; k < rectangles.size(); k++)
			rectanglesOfLevel[rectangles[k]->level_].push_back(rectangles[k]);

		///////////////////
		// build network //
		///////////////////
		network->addLayer(inputLayer);

		for (unsigned k = 0; k <= (unsigned)maxLevel; k++) {
			unsigned curLevelIdx = (unsigned)maxLevel - k;

			// Gaussian Layer
			if (curLevelIdx == (unsigned)maxLevel) {
				gaussLayer = (GaussianLayer*)factory.generateLayer(Layer::GAUSSIAN_LAYER, 0);
				for (unsigned l = 0; l < rectanglesOfLevel[curLevelIdx].size(); l++) {
					rPtr = rectanglesOfLevel[curLevelIdx][l];
					for (nodeIt = rPtr->nodePtrs_.begin(); nodeIt != rPtr->nodePtrs_.end(); nodeIt++) {
						gaussNode = (GaussianNode*)nodeIt->second;
						// not all Gaussians are used in the PD format
						if (gaussNode->getNumParents() > 0)
							gaussLayer->insertNode(gaussNode);
					}
				}
				network->addLayer(gaussLayer);
				continue;
			}

			// Sum and Product layers
			sumLayer = (SumLayer*)factory.generateLayer(Layer::SUM_LAYER, 0);
			prodLayer = (ProdLayer*)factory.generateLayer(Layer::PROD_LAYER, 0);
			for (unsigned l = 0; l < rectanglesOfLevel[curLevelIdx].size(); l++) {
				rPtr = rectanglesOfLevel[curLevelIdx][l];
				for (nodeIt = rPtr->nodePtrs_.begin(); nodeIt != rPtr->nodePtrs_.end(); nodeIt++) {
					curSumNode = (SumNode*)nodeIt->second;
					sumLayer->insertNode(curSumNode);
					for (unsigned pC = 0; pC < curSumNode->getNumChildren(); pC++)
						prodLayer->insertNode(curSumNode->getChild(pC));
				}
			}
			network->addLayer(prodLayer);
			network->addLayer(sumLayer);
		}

		// clean up
		for (unsigned k = 0; k < rectangles.size(); k++) {
			delete rectangles[k];
		}
		for (unsigned k = 0; k < partitions.size(); k++)
			delete partitions[k];
	}


	void readVectorList(std::string fileName, DoubleVec2 &list) {
		std::ifstream inputFile;
		std::string inputLine;
		std::stringstream inputStream;

		double tmpInput;
		unsigned numVecs;

		// open file
		inputFile.open(fileName.c_str());
		if (!inputFile.is_open())
			throw std::runtime_error("readVectorList(char **filename, "
			"std::vector<std::vector<double> > &list): "
			"unable to open file.");

		// check header
		getline(inputFile, inputLine);
		if (inputLine.find_last_not_of(" \t\f\v\n\r") != std::string::npos)
			inputLine.erase(inputLine.find_last_not_of(" \t\f\v\n\r") + 1);

		if ((inputLine.compare(0, 15, "SPN_VECLIST:txt") != 0) || (inputLine.size() > 15)) {
			inputFile.close();
			throw std::runtime_error("readVectorList(char **filename, "
				"std::vector<std::vector<double> > &list): "
				"Header of file has to be \"SPN_VECLIST:txt\".");
		}

		// get #vectors
		getline(inputFile, inputLine);
		inputStream.clear();
		inputStream.str(inputLine);
		inputStream >> tmpInput;
		if ((inputStream.fail()) || (tmpInput != floor(tmpInput))) {
			inputFile.close();
			throw std::runtime_error("readVectorList(char **filename, "
				"std::vector<std::vector<double> > &list): "
				"cannot read number of vectors.");
		}
		numVecs = (unsigned)tmpInput;

		// read list
		list.resize(numVecs);

		for (unsigned long k = 0; k < numVecs; k++) {
			if (!inputFile.good())
				throw std::runtime_error("readVectorList(char **filename, "
				"std::vector<std::vector<double> > &list): "
				"error reading input file (eof?).");

			getline(inputFile, inputLine);
			if (inputLine.find_first_not_of(" \t\f\v\n\r") == std::string::npos)
				inputLine.clear();
			if (inputLine.find_last_not_of(" \t\f\v\n\r") != std::string::npos)
				inputLine.erase(inputLine.find_last_not_of(" \t\f\v\n\r") + 1);

			if (inputLine.empty()) {
				list[k].clear();
				continue;
			}

			inputStream.clear();
			inputStream.str(inputLine);

			std::vector<double> tmpVector;
			tmpVector.clear();

			while (!inputStream.eof()) {
				std::string tmpString;
				std::stringstream tmpStream;

				inputStream >> tmpString;
				if (
					(tmpString.compare("nan") == 0) ||
					(tmpString.compare("NAN") == 0) ||
					(tmpString.compare("NaN") == 0))
					tmpInput = std::numeric_limits<double>::signaling_NaN();
				else {
					tmpStream.clear();
					tmpStream.str(tmpString);
					tmpStream >> tmpInput;

					if (tmpStream.fail())
						throw std::runtime_error("readVectorList(char **filename, "
						"std::vector<std::vector<double> > &list): "
						"error reading input line (eof?).");
				}
				tmpVector.push_back(tmpInput);
			}

			list[k] = tmpVector;
		}
	}

	void readVectorList(std::string fileName, UnsignedVec2 &list) {
		DoubleVec2 input;
		UnsignedVec2 ret;
		readVectorList(fileName, input);

		ret.resize(input.size());
		for (size_t k = 0; k < input.size(); k++) {
			for (size_t l = 0; l < input[k].size(); l++) {
				double roundin = round(input[k][l]);
				if (fabs(input[k][l] - roundin) > 1e-9)
					throw std::runtime_error("readVectorList: non-integer");
				ret[k].push_back((unsigned) roundin);
			}
		}

		list.swap(ret);
	}


	DoubleVec2 readVectorList(std::string filename) {
		DoubleVec2 ret;
		readVectorList(filename, ret);
		return ret;
	}
	
	void writeVectorList(std::string fileName, DoubleVec2 &list) {
		std::ofstream outputFile;
		
		outputFile.open(fileName.c_str());
		if (!outputFile.is_open())
			throw std::runtime_error("writeVectorList(std::string fileName, DoubleVec2 &list): unable to open file.");

		outputFile << "SPN_VECLIST:txt" << std::endl;
		outputFile << list.size() << std::endl;

		for (size_t k = 0; k < list.size(); k++) {
			for (size_t l = 0; l < list[k].size(); l++)
				outputFile << std::setprecision(20) << list[k][l] << " ";
			outputFile << std::endl;
		}

		outputFile.close();
	}


	void writeDot(Node* subroot, std::string fileName, bool weightLabel, bool indbottom) {

		// sort nodes
		NodeVec sumNodes;
		NodeVec prodNodes;
		NodeVec indNodes;
		NodeVec valNodes;

		std::queue<Node*> queue;
		queue.push(subroot);
		USet<NodeId> ids;

		while (queue.empty() == false) {
			Node* nextNode = queue.front();
			queue.pop();

			if (nextNode->getType() == Node::SUM_NODE)
				sumNodes.push_back(nextNode);
			else if (nextNode->getType() == Node::PROD_NODE)
				prodNodes.push_back(nextNode);
			else if (nextNode->getType() == Node::INDICATOR_NODE)
				indNodes.push_back(nextNode);
			else if (nextNode->getType() == Node::VAL_NODE)
				valNodes.push_back(nextNode);

			for (unsigned k = 0; k < nextNode->getNumChildren(); k++) {
				if (ids.count(nextNode->getChild(k)->getId()) == 0) {
					queue.push(nextNode->getChild(k));
					ids.insert(nextNode->getChild(k)->getId());
				}
			}
		}

		// open file
		std::ofstream outputFile;
		outputFile.open(fileName.c_str());
		if (!outputFile.is_open())
			throw std::runtime_error(
			"selectiveLearning::writeDot(Node* subroot, std::string filename): "
			"unable to open file.");

		outputFile << "digraph spn {" << std::endl;
		outputFile << "splines = true;" << std::endl;
		outputFile << std::endl;

		outputFile << "node [shape=circle, label=\"S\"];"
			<< std::endl;
		for (unsigned k = 0; k < sumNodes.size(); k++) {
			outputFile << "n" << sumNodes[k]->getId() <<
				//              "[label=\"" << sumNodes[k]->getId( ) << "\"]" <<
				";" << std::endl;
		}
		outputFile << std::endl;

		outputFile << "node [shape=circle, label=\"P\"];"
			<< std::endl;
		for (unsigned k = 0; k < prodNodes.size(); k++)
			outputFile << "n" << prodNodes[k]->getId() << ";" << std::endl;
		outputFile << std::endl;

		outputFile << "node [shape=circle, label=\"I\"];"
			<< std::endl;
		for (unsigned k = 0; k < indNodes.size(); k++)
			outputFile << "n" << indNodes[k]->getId() << ";" << std::endl;
		outputFile << std::endl;

		outputFile << "node [shape=circle, label=\" \"];" << std::endl;
		for (unsigned k = 0; k < valNodes.size(); k++) {
			size_t label = ((ValNode*)valNodes[k])->getId();
			outputFile << "n" << valNodes[k]->getId()
				<< " [label=\"" << label << "\"];" << std::endl;
		}
		outputFile << std::endl;

		if (indbottom) {
			outputFile << "subgraph clusterInd" << "{style=\"invis\";";
			for (unsigned k = 0; k < valNodes.size(); k++)
				for (unsigned l = 0; l < valNodes[k]->getNumParents(); l++)
					outputFile << "n" << valNodes[k]->getParent(l)->getId() << ";";
			outputFile << std::endl;
		}

		for (unsigned k = 0; k < valNodes.size(); k++) {
			outputFile << "subgraph cluster" << k << "{style=\"solid\";";
			for (unsigned l = 0; l < valNodes[k]->getNumParents(); l++)
				outputFile << "n" << valNodes[k]->getParent(l)->getId() << ";";
			outputFile << "}" << std::endl;
		}

		if (indbottom) {
			outputFile << "}" << std::endl;
			outputFile << std::endl;
		}

		ids.clear();
		queue.push(subroot);
		while (queue.empty() == false) {
			Node* nextNode = queue.front();
			queue.pop();

			for (unsigned k = 0; k < nextNode->getNumChildren(); k++) {
				if ((nextNode->getType() == Node::SUM_NODE) && (weightLabel)) {
					outputFile << "n" << nextNode->getId()
						<< "->n" << nextNode->getChild(k)->getId() << " [label = \""
						<< exp(((SumNode*)nextNode)->getLogWeight(nextNode->getChild(k)))
						<< "\"];" << std::endl;
				}
				else {
					outputFile << "n" << nextNode->getId()
						<< "->n" << nextNode->getChild(k)->getId() << ";" << std::endl;
				}
				if (ids.count(nextNode->getChild(k)->getId()) == 0) {
					queue.push(nextNode->getChild(k));
					ids.insert(nextNode->getChild(k)->getId());
				}
			}
		}

		outputFile << "}" << std::endl;
		outputFile.close();
	}

}