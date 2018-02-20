
#include "NetworkFactory.hpp"


namespace SPN {

	NetworkFactory::NetworkFactory() {
	}

	NetworkFactory::~NetworkFactory() {
		clear();
	}	

	Network* NetworkFactory::generateNetwork() {
		Network* ptr = new Network;
		generatedNetworkPtrs_[ptr->getId()] = ptr;
		return ptr;
	}

	Layer* NetworkFactory::generateLayer(LayerType type, size_t numNodes) {
		Layer* ptr = NULL;
		NodeVec nodes;
		switch (type) {
		case Layer::SUM_LAYER:
			ptr = new SumLayer;
			nodes = generateNodes(Node::SUM_NODE, numNodes);					
			break;
		case Layer::PROD_LAYER:
			ptr = new ProdLayer;
			nodes = generateNodes(Node::PROD_NODE, numNodes);
			break;
		case Layer::VAL_LAYER:
			ptr = new ValLayer;
			nodes = generateNodes(Node::VAL_NODE, numNodes);
			break;
		case Layer::INDICATOR_LAYER:
			ptr = new IndicatorLayer;
			nodes = generateNodes(Node::INDICATOR_NODE, numNodes);
			break;
		case Layer::GAUSSIAN_LAYER:
			ptr = new GaussianLayer;
			nodes = generateNodes(Node::GAUSSIAN_NODE, numNodes);
			break;
		default:
			throw std::runtime_error(
				"error in NetworkFactory::generateLayer(int type, unsigned numNodes): "
				"unknown type.");
		}
		if (ptr != NULL) {
			generatedLayerPtrs_[ptr->getId()] = ptr;
			for (size_t k = 0; k < numNodes; k++)
				ptr->insertNode(nodes[k]);
		}
		return ptr;
	}

	Node* NetworkFactory::generateNode(NodeType type) {
		Node* ptr = NULL;
		switch (type) {
		case Node::SUM_NODE:
			ptr = new SumNode;
			break;
		case Node::PROD_NODE:
			ptr = new ProdNode;
			break;
		case Node::VAL_NODE:
			ptr = new ValNode;
			break;
		case Node::INDICATOR_NODE:
			ptr = new IndicatorNode;
			break;
		case Node::GAUSSIAN_NODE:
			ptr = new GaussianNode;
			break;
		default:
			throw std::runtime_error("error in NetworkFactory::generateNode(int type): unknown type.");
		}
		if (ptr != NULL)
			generatedNodePtrs_[ptr->getId()] = ptr;
		return ptr;
	}

	NodeVec NetworkFactory::generateNodes(NodeType type, size_t numNodes) {
		NodeVec ret;
		for (size_t k = 0; k < numNodes; k++)
			ret.push_back(generateNode(type));
		return ret;
	}
	
	void NetworkFactory::deleteNetwork(Network *network) {
		NetworkId id = network->getId();
		if (generatedNetworkPtrs_.count(id) == 0)
			throw std::runtime_error("NetworkFactory::deleteNetwork(Network *network): layer not found.");
		delete generatedNetworkPtrs_[id];
		generatedNetworkPtrs_.erase(id);
	}

	void NetworkFactory::deleteNetworkDeep(Network *network) {
		for (size_t k = 0; k < network->getNumLayers(); k++)
			deleteLayerDeep(network->getLayer(k));
		deleteNetwork(network);
	}

	void NetworkFactory::deleteLayer(Layer *layer) {
		LayerId id = layer->getId();
		if (generatedLayerPtrs_.count(id) == 0)
			throw std::runtime_error("NetworkFactory::deleteLayer(Layer *layer): layer not found.");
		delete generatedLayerPtrs_[id];
		generatedLayerPtrs_.erase(id);
	}

	void NetworkFactory::deleteLayerDeep(Layer *layer) {
		for (size_t k = 0; k < layer->getNumNodes(); k++)
			deleteNode(layer->getNode(k));
		deleteLayer(layer);
	}

	void NetworkFactory::deleteNode(Node *node) {
		NodeId id = node->getId();
		if (generatedNodePtrs_.count(id) == 0)
			throw std::runtime_error("NetworkFactory::deleteNode(Node *node): node not found.");
		delete generatedNodePtrs_[id];
		generatedNodePtrs_.erase(id);
	}

	Network* NetworkFactory::copyNetwork(Network *network) {

		Network* newNetwork = generateNetwork();
		UMap<NodeId, Node*> nodeMap;

		for (size_t l = 0; l < network->getNumLayers(); l++) {
			Layer* curLayer = network->getLayer(l);
			Layer* newLayer = generateLayer(curLayer->getType(), curLayer->getNumNodes());
			for (size_t n = 0; n < curLayer->getNumNodes(); n++) {				
				Node* curNode = curLayer->getNode(n);
				Node* newNode = newLayer->getNode(n);
				nodeMap[curNode->getId()] = newNode;
				for (size_t c = 0; c < curNode->getNumChildren(); c++) {
					Node* curChild = curNode->getChild(c);
					Node* newChild = nodeMap.at(curChild->getId());
					newNode->connectChild(newChild);
					if (curNode->getType() == Node::SUM_NODE) {
						double w = ((SumNode*)curNode)->getLogWeight(curChild);
						((SumNode*)newNode)->setLogWeight(newChild, w);
					}
					else if (curNode->getType() == Node::INDICATOR_NODE) {
						double i = ((IndicatorNode*)curNode)->getIndValue(curChild);
						((IndicatorNode*)newNode)->setIndValue(newChild, i);
					}
					else if (curNode->getType() == Node::GAUSSIAN_NODE) {
						double m = ((GaussianNode*)curNode)->getMean(curChild);
						double s = ((GaussianNode*)curNode)->getSigma(curChild);
						((GaussianNode*)newNode)->setMean(newChild, m);
						((GaussianNode*)newNode)->setSigma(newChild, s);
					}
				}
			}
			newNetwork->addLayer(newLayer);
		}
		return newNetwork;
	}

	void NetworkFactory::clear() {
		for (auto it = generatedNetworkPtrs_.begin(); it != generatedNetworkPtrs_.end(); it++)
			delete it->second;
		generatedNetworkPtrs_.clear();

		for (auto it = generatedLayerPtrs_.begin(); it != generatedLayerPtrs_.end(); it++)
			delete it->second;
		generatedLayerPtrs_.clear();

		for (auto it = generatedNodePtrs_.begin(); it != generatedNodePtrs_.end(); it++)
			delete it->second;
		generatedNodePtrs_.clear();
	}
}


