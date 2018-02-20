
#include "Layer.hpp"

namespace SPN {

	LayerId Layer::idCounter_ = 0;

	Layer::Layer(LayerType type) {
		type_ = type;
		id_ = idCounter_++;
		nodes_.clear();
		nodeSet_.clear();
	}

	Layer::~Layer() {
	}

	void Layer::eval() {
		for (unsigned k = 0; k < nodes_.size(); k++)
			nodes_[k]->eval();
	}

	void Layer::backprop() {
		for (unsigned k = 0; k < nodes_.size(); k++)
			nodes_[k]->backprop();
	}

	void Layer::insertNode(Node *node) {
		insertNodes(NodeVec(1, node));
	}

	void Layer::insertNodes(NodeVec nodes) {
		for (size_t k = 0; k < nodes.size(); k++) {
			if (((type_ == Layer::SUM_LAYER) && (nodes[k]->getType() != Node::SUM_NODE)) ||
				((type_ == Layer::PROD_LAYER) && (nodes[k]->getType() != Node::PROD_NODE)) ||
				((type_ == Layer::INDICATOR_LAYER) && (nodes[k]->getType() != Node::INDICATOR_NODE)) ||
				((type_ == Layer::VAL_LAYER) && (nodes[k]->getType() != Node::VAL_NODE)) ||
				((type_ == Layer::GAUSSIAN_LAYER) && (nodes[k]->getType() != Node::GAUSSIAN_NODE)))
				throw std::runtime_error("Layer::insertNodes(NodeVec nodes): wrong type");
			if (nodeSet_.count(nodes[k]->getId()))
				continue;
			nodeSet_.insert(nodes[k]->getId());
			nodes_.push_back(nodes[k]);
		}
	}

	void Layer::removeNode(Node *node) {
		nodes_.erase(nodes_.begin() + getNodeIdx(node));

		if (nodeSet_.count(node->getId()) == 0)
			throw std::runtime_error("Layer::removeNode(Node *node): inconsistent nodeIndMap.");
		nodeSet_.erase(node->getId());
	}

	void Layer::removeNodes(NodeVec nodes) {
		for (size_t k = 0; k < nodes.size(); k++)
			removeNode(nodes[k]);
	}

	size_t Layer::getNumNodes() {
		return nodes_.size();
	}

	NodeVec Layer::getNodes() {
		return nodes_;
	}

	Node* Layer::getNode(size_t idx) {
		if (idx >= nodes_.size())
			throw std::runtime_error("error in Layer::getNode(size_t idx): idx out of bounds.");
		return (nodes_[idx]);
	}

	size_t Layer::getNodeIdx(Node* node) {
		NodeId nodeId = node->getId();
		size_t retidx = 0;
		bool found = false;

		for (unsigned k = 0; k < nodes_.size(); k++) {
			if (nodes_[k]->getId() == nodeId) {
				retidx = k;
				found = true;
				break;
			}
		}

		if (!found)
			throw std::runtime_error("error in Layer::getNodeIdx(Node* node): "
			"layer does not contain node.");

		return retidx;
	}

	bool Layer::containsNode(Node *node) {
		return (nodeSet_.count(node->getId()) != 0);
	}

	int Layer::getType() const {
		return type_;
	}

	LayerId Layer::getId() const {
		return id_;
	}

	Node* Layer::operator()(size_t idx) {
		return getNode(idx);
	}


	/////////////////
	/// Sum Layer ///
	/////////////////
	SumLayer::SumLayer() : Layer(SUM_LAYER) {
	}

	SumLayer::~SumLayer() {
	}

	SumNode* SumLayer::getNode(size_t idx) {
		return (SumNode*)Layer::getNode(idx);
	}

	void SumLayer::initWeightsRand() {
		for (unsigned k = 0; k < getNumNodes(); k++)
			getNode(k)->initWeightsRand();
	}

	SumNode* SumLayer::operator() (size_t idx) {
		return getNode(idx);
	}


	/////////////////////
	/// Product Layer ///
	/////////////////////
	ProdLayer::ProdLayer() : Layer(PROD_LAYER) {
	}

	ProdLayer::~ProdLayer() {
	}

	ProdNode* ProdLayer::getNode(size_t idx) {
		return (ProdNode*)Layer::getNode(idx);
	}

	ProdNode* ProdLayer::operator() (size_t idx) {
		return getNode(idx);
	}


	///////////////////
	/// Value Layer ///
	///////////////////
	ValLayer::ValLayer() : Layer(VAL_LAYER) {
	}

	ValLayer::~ValLayer() {
	}

	ValNode* ValLayer::getNode(size_t idx) {
		return (ValNode*)Layer::getNode(idx);
	}

	ValNode* ValLayer::operator() (size_t idx) {
		return getNode(idx);
	}

	void ValLayer::setVals(DoubleVec values) {
		size_t numNodes = getNumNodes();

		if (numNodes != values.size())
			throw std::runtime_error("error in ValLayer::setVals(DoubleVec values): "
			"values vector must contain getNumNodes() entries.");

		for (unsigned k = 0; k < numNodes; k++)
			((ValNode*)nodes_[k])->setVal(values[k]);
	}

	DoubleVec ValLayer::getVals() {
		size_t numNodes = getNumNodes();
		DoubleVec values;
		values.resize(numNodes);
		for (unsigned k = 0; k < numNodes; k++)
			values[k] = ((ValNode*)nodes_[k])->getVal();
		return values;
	}


	///////////////////////
	/// Indicator Layer ///
	///////////////////////
	IndicatorLayer::IndicatorLayer() : Layer(INDICATOR_LAYER) {
	}

	IndicatorLayer::~IndicatorLayer() {
	}

	IndicatorNode* IndicatorLayer::getNode(size_t idx) {
		return (IndicatorNode*)Layer::getNode(idx);
	}

	IndicatorNode* IndicatorLayer::operator() (size_t idx) {
		return getNode(idx);
	}


	//////////////////////
	/// Gaussian Layer ///
	//////////////////////
	GaussianLayer::GaussianLayer() : Layer(GAUSSIAN_LAYER) {
	}

	GaussianLayer::~GaussianLayer() {
	}

	GaussianNode* GaussianLayer::getNode(size_t idx) {
		return (GaussianNode*)Layer::getNode(idx);
	}

	void GaussianLayer::setMeansUniform(double meanVal) {
		for (unsigned k = 0; k < getNumNodes(); k++) {
			GaussianNode* gaussNode = getNode(k);
			for (unsigned l = 0; l < gaussNode->getNumChildren(); l++) {
				gaussNode->setMean(gaussNode->getChild(l), meanVal);
			}
		}
	}

	void GaussianLayer::setSigmasUniform(double sigmaVal) {
		for (unsigned k = 0; k < getNumNodes(); k++) {
			GaussianNode* gaussNode = getNode(k);
			for (unsigned l = 0; l < gaussNode->getNumChildren(); l++) {
				gaussNode->setSigma(gaussNode->getChild(l), sigmaVal);
			}
		}
	}

	void GaussianLayer::setMeansRand(double minVal, double maxVal) {
		for (unsigned k = 0; k < getNumNodes(); k++)
			getNode(k)->setMeansRand(minVal, maxVal);
	}

	void GaussianLayer::setSigmasRand(double minVal, double maxVal) {
		for (unsigned k = 0; k < getNumNodes(); k++)
			getNode(k)->setSigmasRand(minVal, maxVal);
	}

	GaussianNode* GaussianLayer::operator() (size_t idx) {
		return getNode(idx);
	}

}


