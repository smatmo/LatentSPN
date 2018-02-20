
#include "Network.hpp"
#include <iostream>

namespace SPN {

	NetworkId Network::idCounter_ = 0;

	Network::Network() {
		id_ = idCounter_++;
	}

	Network::~Network() {
	}

	void Network::insertLayer(Layer* layer, size_t idx) {
		if (containsLayer(layer))
			return;

		if (idx >= getNumLayers())
			layers_.push_back(layer);
		else
			layers_.insert(layers_.begin() + idx, layer);
	}

	void Network::addLayer(Layer* layer) {
		insertLayer(layer, getNumLayers());
	}

	void Network::removeLayer(Layer *layer) {
		if (!containsLayer(layer))
			throw std::runtime_error("Network::removeLayer(Layer *layer):tried to remove non-containing layer.");
		layers_.erase(layers_.begin() + getLayerIdx(layer));
	}

	size_t Network::getNumLayers() {
		return layers_.size();
	}

	LayerVec Network::getLayers() {
		return layers_;
	}

	Layer* Network::getLayer(size_t layerIdx) {
		if (layerIdx >= layers_.size())
			throw std::runtime_error("error in Network::getLayer(unsigned idx): idx out of bounds.");
		return layers_[layerIdx];
	}

	Layer* Network::getTopLayer() {
		if (layers_.size() == 0)
			return NULL;
		return layers_[layers_.size()-1];
	}

	size_t Network::getLayerIdx(Layer* layer) {
		unsigned retidx = 0;
		bool found = false;
		LayerId layerId = layer->getId();

		for (unsigned k = 0; k < layers_.size(); k++) {
			if (layers_[k]->getId() == layerId) {
				retidx = k;
				found = true;
				break;
			}
		}

		if (!found)
			throw std::runtime_error("error in Network::getLayerIdx(Layer* layer): network does not contain layer.");

		return retidx;
	}

	Node* Network::getNode(size_t layerIdx, size_t idx) {
		return getLayer(layerIdx)->getNode(idx);
	}

	size_t Network::getLayerIdxOfNode(Node *node) {
		for (unsigned k = 0; k < getNumLayers(); k++)
			if (getLayer(k)->containsNode(node))
				return k;
		throw std::runtime_error("error in Network::getLayerIdxOfNode(Node *node): "
			"network does not contain node.");
	}

	bool Network::containsLayer(Layer *layer) {
		bool found = false;
		LayerId layerId = layer->getId();

		for (unsigned k = 0; k < layers_.size(); k++) {
			if (layers_[k]->getId() == layerId) {
				found = true;
				break;
			}
		}

		return found;
	}

	void Network::setVals(const DoubleVec &values) {
		if (getNumLayers() == 0)
			throw std::runtime_error("error in Network::setVals(const DoubleVec &values): Network has no layers.");

		if (layers_[0]->getType() != Layer::VAL_LAYER)
			throw std::runtime_error("error in Network::setVals(const DoubleVec &values): layer 0 has to be of type Layer::VAL_LAYER.");

		((ValLayer*)layers_[0])->setVals(values);
	}

	DoubleVec Network::getVals() {
		if (getNumLayers() == 0)
			throw std::runtime_error("error in Network::getVals(): Network has no layers.");

		if (layers_[0]->getType() != Layer::VAL_LAYER)
			throw std::runtime_error("error in Network::getVals(): layer 0 has to be of type Layer::VAL_LAYER.");

		return ((ValLayer*)layers_[0])->getVals();
	}

	void Network::eval() {
		for (unsigned k = 0; k < layers_.size(); k++)
			layers_[k]->eval();
	}

	void Network::backprop() {
		for (unsigned k = 0; k < layers_.size(); k++)
			layers_[layers_.size() - k - 1]->backprop();
	}
	
	void Network::MPEinference(size_t rootidx, bool dirty) {
		if (getTopLayer()->getType() != Layer::SUM_LAYER)
			throw std::runtime_error("Network::MPEinference(): last layer has to be a sum layer");

		if (rootidx >= getTopLayer()->getNumNodes()) {
			std::cout << rootidx << " " << getTopLayer()->getNumNodes() << std::endl;
			throw std::runtime_error("Network::MPEinference(): rootidx out of bounds");
		}

		if (dirty) {
			eval();
			backprop();
		}

		std::vector<SumNode* > maxSumNodes;

		maxSumNodes.clear();
		maxSumNodes.push_back((SumNode*)getTopLayer()->getNode(rootidx));

		for (unsigned k = 0; k < maxSumNodes.size(); k++) {
			SumNode* curSumNode = maxSumNodes[k];
			std::vector<ProdNode*> prodNodeOptions;
			prodNodeOptions.clear();
			double maxVal = -inf;

			for (unsigned l = 0; l < curSumNode->getNumChildren(); l++) {
				if (curSumNode->getChild(l)->getType() != Node::PROD_NODE)
					throw std::runtime_error("Network::MPEinference(): children of sum nodes have to be product nodes");

				double curVal = curSumNode->getChild(l)->getLogVal() + curSumNode->getChild(l)->getLogDerivative();
				if ((l == 0) || (curVal > maxVal)) {
					maxVal = curVal;
					prodNodeOptions.clear();
				}
				if (curVal == maxVal)
					prodNodeOptions.push_back((ProdNode*)curSumNode->getChild(l));
			}

			std::queue<ProdNode*> prodQ;
			prodQ.push(prodNodeOptions[rand() % prodNodeOptions.size()]);

			while (!prodQ.empty()) {
				ProdNode* prodNode = prodQ.front();
				prodQ.pop();

				for (unsigned l = 0; l < prodNode->getNumChildren(); l++) {
					if (prodNode->getChild(l)->getType() == Node::SUM_NODE)
						maxSumNodes.push_back((SumNode*)prodNode->getChild(l));
					else if (prodNode->getChild(l)->getType() == Node::GAUSSIAN_NODE) {
						GaussianNode* gaussNode = (GaussianNode*)prodNode->getChild(l);
						for (unsigned m = 0; m < gaussNode->getNumChildren(); m++) {
							if (gaussNode->getChild(m)->getType() != Node::VAL_NODE)
								throw std::runtime_error("error: child of gauss has to be val node.");
							ValNode* valNode = (ValNode*)(gaussNode->getChild(m));
							if (isnan(valNode->getVal()))
								valNode->setVal(gaussNode->getMean(valNode));
						}
					}
					else if (prodNode->getChild(l)->getType() == Node::PROD_NODE)
						prodQ.push((ProdNode*)prodNode->getChild(l));
					else
						throw std::runtime_error("MPEInference: unexpected node type.");
				}
			}
		}
	}
	
	void Network::maxBacktracking() {

		if (getLayer(getNumLayers() - 1)->getNumNodes() != 1)
			throw std::runtime_error("Network::MaxBacktracking(): last layer must contain exactly one node");

		std::queue<Node*> Q;
		Q.push(getLayer(getNumLayers() - 1)->getNode(0));

		while (!Q.empty()) {

			Node* nextNode = Q.front();
			Q.pop();

			if (nextNode->getType() == Node::SUM_NODE) {
				NodeVec nodeOptions;
				nodeOptions.clear();
				double maxVal = -inf;
				for (size_t l = 0; l < nextNode->getNumChildren(); l++) {
					Node* curChild = nextNode->getChild(l);
					double curVal = curChild->getLogVal() + ((SumNode*)nextNode)->getLogWeight(curChild);
					if ((l == 0) || (curVal > maxVal)) {
						maxVal = curVal;
						nodeOptions.clear();
					}
					if (curVal == maxVal)
						nodeOptions.push_back(curChild);
				}
				Q.push(nodeOptions[rand() % nodeOptions.size()]);
			}
			else if (nextNode->getType() == Node::PROD_NODE) {
				for (size_t l = 0; l < nextNode->getNumChildren(); l++)
					Q.push(nextNode->getChild(l));
			}
			else if (nextNode->getType() == Node::INDICATOR_NODE) {
				for (size_t l = 0; l < nextNode->getNumChildren(); l++) {
					Node* curChild = nextNode->getChild(l);
					if (curChild->getType() != Node::VAL_NODE)
						throw std::runtime_error("Network::MaxBacktracking( ): child of indicator is not a valNode.");
					if (isnan(((ValNode*)curChild)->getVal()))
						((ValNode*)curChild)->setVal(((IndicatorNode*)nextNode)->getIndValue(curChild));
				}
			}
			else if (nextNode->getType() == Node::GAUSSIAN_NODE) {
				for (size_t l = 0; l < nextNode->getNumChildren(); l++) {
					Node* curChild = nextNode->getChild(l);
					if (curChild->getType() != Node::VAL_NODE)
						throw std::runtime_error("Network::MaxBacktracking( ): child of Gaussian is not a valNode.");
					if (isnan(((ValNode*)curChild)->getVal()))
						((ValNode*)curChild)->setVal(((GaussianNode*)nextNode)->getMean(curChild));
				}
			}
			else {
				throw std::runtime_error("Network::MaxBacktracking( ): node type not implemented.");
			}
		}
	}


	void Network::setMPEcorrectionWeights() {

		NodeVec desc;
		listDescendants(getLayer(getNumLayers() - 1)->getNode(0), desc);

		for (size_t k = 0; k < desc.size(); k++) {
			if (desc[k]->getType() != Node::SUM_NODE)
				continue;

			SumNode* S = (SumNode*)desc[k];
			

			// find ancestors
			NodeVec anc;
			listAncestors(S, anc);
			USet<NodeId> ancIDs;
			for (size_t l = 0; l < anc.size(); l++)
				ancIDs.insert(anc[l]->getId());

			// find sum ancestors
			NodeVec sumAnc = anc;
			keepNodesOfType(sumAnc, Node::SUM_NODE);

			for (size_t a = 0; a < sumAnc.size(); a++) {
				Node* Sanc = sumAnc[a];
				if (Sanc->getId() == S->getId())
					continue;

				for (unsigned c = 0; c < Sanc->getNumChildren(); c++) {
					if (ancIDs.count(Sanc->getChild(c)->getId()) == 0) {
						double logW = ((SumNode*)Sanc)->getLogWeight(Sanc->getChild(c));
						logW -= log((double)S->getNumChildren());
						((SumNode*)Sanc)->setLogWeight(Sanc->getChild(c), logW);
					}
				}
			}
		}
	}
	
	void Network::setSumMode() {
		for (unsigned k = 0; k < layers_.size(); k++)
			if (layers_[k]->getType() == Layer::SUM_LAYER)
				for (unsigned l = 0; l < layers_[k]->getNumNodes(); l++)
					((SumNode*)layers_[k]->getNode(l))->setSumMode();
	}

	void Network::setMaxMode() {
		for (unsigned k = 0; k < layers_.size(); k++)
			if (layers_[k]->getType() == Layer::SUM_LAYER)
				for (unsigned l = 0; l < layers_[k]->getNumNodes(); l++)
					((SumNode*)layers_[k]->getNode(l))->setMaxMode();
	}

	void Network::setGaussianMeansRand(double minVal, double maxVal) {
		for (unsigned k = 0; k < layers_.size(); k++)
			if (layers_[k]->getType() == Layer::GAUSSIAN_LAYER)
				((GaussianLayer*)layers_[k])->setMeansRand(minVal, maxVal);
	}

	void Network::setGaussianSigmasRand(double minVal, double maxVal) {
		for (unsigned k = 0; k < layers_.size(); k++)
			if (layers_[k]->getType() == Layer::GAUSSIAN_LAYER)
				((GaussianLayer*)layers_[k])->setSigmasRand(minVal, maxVal);
	}

	void Network::setGaussianMeansUniform(double meanVal) {
		for (unsigned k = 0; k < layers_.size(); k++)
			if (layers_[k]->getType() == Layer::GAUSSIAN_LAYER)
				((GaussianLayer*)layers_[k])->setMeansUniform(meanVal);
	}

	void Network::setGaussianSigmasUniform(double sigmaVal) {
		for (unsigned k = 0; k < layers_.size(); k++)
			if (layers_[k]->getType() == Layer::GAUSSIAN_LAYER)
				((GaussianLayer*)layers_[k])->setSigmasUniform(sigmaVal);
	}

	void Network::initWeightsRand() {
		for (unsigned k = 0; k < layers_.size(); k++)
			if (layers_[k]->getType() == Layer::SUM_LAYER)
				((SumLayer*)layers_[k])->initWeightsRand();
	}

	void Network::removeSingletonParents(NodeVec &removedNodes, LayerVec &removedLayers) {
		removedNodes.clear();
		removedLayers.clear();

		for (size_t k = 0; k < getNumLayers(); k++) {
			NodeVec removed;
			for (size_t l = 0; l < getLayer(k)->getNumNodes(); l++) {
				Node* n = getLayer(k)->getNode(l);
				if ((n->getType() != Node::SUM_NODE) && (n->getType() != Node::PROD_NODE))
					continue;
				if (n->getNumChildren() == 1) {
					Node* child = n->getChild(0);
					
					for (size_t p = 0; p < n->getNumParents(); p++) {
						Node* parent = n->getParent(p);
						if (parent->getType() == Node::SUM_NODE) {
							parent->connectChild(child);
							((SumNode*)parent)->setLogWeight(child, ((SumNode*)parent)->getLogWeight(n));
						}
						else if (parent->getType() == Node::PROD_NODE)
							parent->connectChild(child);
						else
							throw std::runtime_error("Network::removeSingletonParents(): unexpected node type");
					}
					n->isolate();
					removed.push_back(n);
				}
			}
			getLayer(k)->removeNodes(removed);
			removedNodes.insert(removedNodes.end(), removed.begin(), removed.end());
		}

		for (size_t k = 0; k < getNumLayers(); k++)
			if (getLayer(k)->getNumNodes() == 0)
				removedLayers.push_back(getLayer(k));

		for (size_t k = 0; k < removedLayers.size(); k++)
			removeLayer(removedLayers[k]);
	}

	void Network::clear() {
		layers_.clear();
	}

	size_t Network::getNumWeights() {
		size_t n = 0;
		for (unsigned k = 0; k < getNumLayers(); k++)
			if (getLayer(k)->getType() == Layer::SUM_LAYER)
				for (unsigned l = 0; l < getLayer(k)->getNumNodes(); l++)
					n += getLayer(k)->getNode(l)->getNumChildren();
		return n;
	}

	size_t Network::getNumGaussians() {
		size_t n = 0;
		for (unsigned k = 0; k < getNumLayers(); k++)
			if (getLayer(k)->getType() == Layer::GAUSSIAN_LAYER)
				n += getLayer(k)->getNumNodes();
		return n;
	}

	size_t Network::getNumSums() {
		size_t n = 0;
		for (unsigned k = 0; k < getNumLayers(); k++)
			if (getLayer(k)->getType() == Layer::SUM_LAYER)
				n += getLayer(k)->getNumNodes();
		return n;
	}

	size_t Network::getNumProds() {
		size_t n = 0;
		for (unsigned k = 0; k < getNumLayers(); k++)
			if (getLayer(k)->getType() == Layer::PROD_LAYER)
				n += getLayer(k)->getNumNodes();
		return n;
	}

	void Network::performStructureCheck() {
		if (getNumLayers() == 0)
			throw std::runtime_error("error in Network::performStructureCheck(): Network has no layers.");
		if (layers_[0]->getType() != Layer::VAL_LAYER)
			throw std::runtime_error("error in Network::performStructureCheck(): layer 0 has to be of type Layer::VAL_LAYER.");
		if (getTopLayer()->getNumNodes() != 1)
			throw std::runtime_error("error in Network::performStructureCheck(): top layer contains not exactly 1 node (currently we assume 1 root)");

		Node* root = getTopLayer()->getNode(0);
		NodeVec desc;
		USet<NodeId> nodeIds;

		// get the nodes reachable from the root
		listDescendants(root, desc);
		for (size_t k = 0; k < desc.size(); k++)
			nodeIds.insert(desc[k]->getId());

		UMap<NodeId, size_t> layerIdx;
		for (size_t k = 0; k < getNumLayers(); k++) {
			for (size_t l = 0; l < getLayer(k)->getNumNodes(); l++) {
				Node* curNode = getLayer(k)->getNode(l);

				if (nodeIds.count(curNode->getId()) == 0) {
					throw std::runtime_error("error in Network::performStructureCheck(): not all nodes reachable from the root");
				}

				// node inserted twice?
				if (layerIdx.count(curNode->getId())) {
					structFeedforward_ = false;
					structComplete_ = false;
					structDecomposable_ = false;
					return;
				}

				layerIdx[curNode->getId()] = k;

				// node gets only input from below?
				for (size_t c = 0; c < curNode->getNumChildren(); c++) {
					Node* child = curNode->getChild(c);
					if ((layerIdx.count(child->getId()) == 0) || (layerIdx[child->getId()] >= k)) {
						structFeedforward_ = false;
						structComplete_ = false;
						structDecomposable_ = false;
						return;
					}
				}
			}
		}

		structFeedforward_ = true;


		// check decomposability and completeness
		UMap<NodeId, Size_tVec > scopes;

		structComplete_ = true;
		structDecomposable_ = true;

		for (size_t k = 1; k < getNumLayers(); k++) {
			if ((!structComplete_) && (!structDecomposable_))
				break;

			for (size_t l = 0; l < getLayer(k)->getNumNodes(); l++) {
				if ((!structComplete_) && (!structDecomposable_))
					break;

				Node* curNode = getNode(k, l);

				// completness -> scopes of all sum children are the same
				if (structComplete_) {
					if (curNode->getType() == Node::SUM_NODE) {
						for (size_t m = 1; m < curNode->getNumChildren(); m++) {
							// it's enough to check all adjacent pairs of children 
							Node* child1 = curNode->getChild(m - 1);
							Node* child2 = curNode->getChild(m);

							if ((scopes.count(child1->getId()) == 0) || (scopes.count(child2->getId()) == 0))
								throw std::runtime_error("error in Network::performStructureCheck(): scope error");

							Size_tVec tmpScope1 = scopes[child1->getId()];
							Size_tVec tmpScope2 = scopes[child2->getId()];

							// if not of same size => not the same scope
							if (tmpScope1.size() != tmpScope2.size()) {
								structComplete_ = false;
								break;
							}
							
							// check element-wise, scope is ordered
							for (size_t n = 0; n < tmpScope1.size(); n++) {
								if (tmpScope1[n] != tmpScope2[n]) {
									structComplete_ = false;
									break;
								}
							}

							if (!structComplete_)
								break;
						}
					}
				}

				// decomposability -> scopes of product children does not overlap
				if (structDecomposable_) {
					if ((curNode->getType() == Node::PROD_NODE) && (curNode->getNumChildren() > 1)) {
						for (size_t m = 0; m + 1 < curNode->getNumChildren(); m++) {
							Node* child1 = curNode->getChild(m);
							Size_tVec tmpScope1 = scopes[child1->getId()];

							for (size_t n = m + 1; n < curNode->getNumChildren(); n++) {
								Node* child2 = curNode->getChild(n);
								Size_tVec tmpScope2 = scopes[child2->getId()];

								// check indices pair-wise
								for (size_t o = 0; o < tmpScope1.size(); o++) {
									for (size_t p = 0; p < tmpScope2.size(); p++) {
										if (tmpScope1[o] == tmpScope2[p]) {
											structDecomposable_ = false;
											break;
										}
									}
									if (!structDecomposable_)
										break;
								}
								if (!structDecomposable_)
									break;
							}
							if (!structDecomposable_)
								break;
						}
					}
				}

				// make scope for current node
				Size_tVec curScope;
				for (size_t m = 0; m < curNode->getNumChildren(); m++) {
					Node* child = curNode->getChild(m);
					if (child->getType() == Node::VAL_NODE) {
						curScope.push_back(getLayer(0)->getNodeIdx(child));
					}
					else {
						if (scopes.count(child->getId()) == 0)
							throw std::runtime_error("error in Network::performStructureCheck(): scope error");
						Size_tVec tmpScope = scopes[child->getId()];
						curScope.insert(curScope.end(), tmpScope.begin(), tmpScope.end());
					}
				}

				std::sort(curScope.begin(), curScope.end());
				Size_tVec::iterator it = std::unique(curScope.begin(), curScope.end());
				curScope.resize(it - curScope.begin());
				scopes[curNode->getId()] = curScope;
			}
		}
	}

	bool Network::isComplete() {
		return structComplete_;
	}

	bool Network::isDecomposable() {
		return structDecomposable_;
	}

	bool Network::isFeedforward() {
		return structFeedforward_;
	}

	Layer* Network::operator()(size_t layerIdx) {
		return getLayer(layerIdx);
	}

	Node* Network::operator()(size_t layerIdx, size_t idx) {
		return getNode(layerIdx, idx);
	}

	NetworkId Network::getId() {
		return id_;
	}

}
