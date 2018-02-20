
#include "Node.hpp"

namespace SPN {

	///////////////////////
	/// Base class Node ///
	///////////////////////
	NodeId Node::idCounter_ = 0;

	Node::Node(NodeType type) {
		if (
			(type != SUM_NODE) &&
			(type != PROD_NODE) &&
			(type != VAL_NODE) &&
			(type != INDICATOR_NODE) &&
			(type != GAUSSIAN_NODE))
			throw std::runtime_error("error in SPN::Node(int type): unknown type.");

		type_ = type;
		id_ = idCounter_++;
		children_.clear();
		parents_.clear();
		childInd_.clear();
		logVal_ = 0;
		logDerivative_ = -inf;
	}

	Node::~Node() {
		isolate();
	}

	NodeId Node::getId() const {
		return id_;
	}

	NodeType Node::getType() const {
		return type_;
	}

	void Node::forceLogVal(double logVal) {
		logVal_ = logVal;
	}

	double Node::getLogVal() const {
		return logVal_;
	}

	double Node::getLogDerivative() const {
		return logDerivative_;
	}

	double Node::getVal() const {
		return exp(logVal_);
	}

	double Node::getDerivative() const {
		return exp(logDerivative_);
	}

	void Node::backprop() {
		size_t numParents = getNumParents();
		if (numParents == 0) {
			logDerivative_ = 0.0;
			return;
		}

		double maxVal = -inf, derivative = 0.0;
		double *computeCache = new double[numParents];

		for (unsigned k = 0; k < numParents; k++) {
			Node *parent = getParent(k);
			computeCache[k] = parent->getLog_Dval_Dchild(this) + parent->getLogDerivative();
			if (computeCache[k] > maxVal)
				maxVal = computeCache[k];
		}

		if (!(maxVal > -inf)) {
			logDerivative_ = -inf;
			delete[] computeCache;
			return;
		}

		for (unsigned k = 0; k < numParents; k++)
			derivative += exp(computeCache[k] - maxVal);

		if (isnan(derivative))
			logDerivative_ = nan;
		else if (derivative > 0)
			logDerivative_ = maxVal + log(derivative);
		else
			logDerivative_ = -inf;

		delete[] computeCache;
	}

	void Node::connectChild(Node *child) {
		if ((isChild(child)) || (child->getId() == getId()))
			return;
		children_.push_back(child);
		child->parents_.push_back(this);
		childInd_.insert(child->getId());
	}

	void Node::connectParent(Node *parent) {
		if (isParent(parent) || (parent->getId() == getId()))
			return;
		parent->connectChild(this);
	}

	void Node::disconnectChild(Node *child) {
		NodeId childId = child->getId();

		bool found = false;
		for (unsigned k = 0; k < children_.size(); k++) {
			if (children_[k]->getId() == childId) {
				found = true;
				childInd_.erase(children_[k]->getId());
				children_.erase(children_.begin() + k);
				break;
			}
		}

		if (!found)
			throw std::runtime_error("Node::disconnectChild(Node &child): tried to disconnect non-child");

		found = false;
		for (unsigned k = 0; k < child->parents_.size(); k++) {
			if (child->parents_[k]->getId() == id_) {
				found = true;
				child->parents_.erase(child->parents_.begin() + k);
				break;
			}
		}

		if (!found)
			throw std::runtime_error("Node::disconnectChild(Node &child): inconsistent parent-child connection");
	}

	void Node::disconnectParent(Node *parent) {
		if (!isParent(parent))
			return;

		parent->disconnectChild(this);
	}

	void Node::isolate() {
		while (getNumParents() != 0)
			disconnectParent(getParent(getNumParents() - 1));

		while (getNumChildren() != 0)
			disconnectChild(getChild(getNumChildren() - 1));
	}

	NodeVec Node::getChildren() {
		return children_;
	}

	NodeVec Node::getParents() {
		return parents_;
	}

	Node* Node::getChild(size_t idx) {
		return (children_[idx]);
	}

	Node* Node::getParent(size_t idx) {
		return (parents_[idx]);
	}

	size_t Node::getNumChildren() {
		return children_.size();
	}

	size_t Node::getNumParents() {
		return parents_.size();
	}

	bool Node::isChild(Node *node) {
		return (childInd_.count(node->getId()) != 0);
	}

	bool Node::isParent(Node *node) {
		return node->isChild(this);
	}

	bool Node::operator==(const Node &rightNode) const {
		return (getId() == rightNode.getId());
	}



	////////////////
	/// Sum Node ///
	////////////////
	SumNode::SumNode() : Node(SUM_NODE) {
		logVal_ = -inf;
		maxLogVal_ = -inf;
		maxIdx_ = 0;
		logWeights_.clear();
		maxMode_ = false;
	}

	SumNode::~SumNode() {
	}

	void SumNode::eval() {
		size_t numChildren = getNumChildren();

		if (numChildren == 0) {
			maxIdx_ = 0;
			maxLogVal_ = -inf;
			logVal_ = -inf;
			return;
		}

		double val = 0.0;
		double *computeCache = new double[logWeights_.size()];

		maxLogVal_ = -inf;
		maxIdx_ = 0;

		for (unsigned k = 0; k < numChildren; k++) {
			Node* child = getChild(k);
			NodeId childId = child->getId();
			computeCache[k] = child->getLogVal() + logWeights_[childId];
			if (computeCache[k] > maxLogVal_) {
				maxLogVal_ = computeCache[k];
				maxIdx_ = k;
			}
		}

		if (maxLogVal_ > -inf) {
			for (unsigned k = 0; k < numChildren; k++)
				val += exp(computeCache[k] - maxLogVal_);
		}
		else
			val = 0.0;

		if (isnan(val))
			logVal_ = nan;
		else if (val > 0)
			logVal_ = maxLogVal_ + log(val);
		else
			logVal_ = -inf;

		delete[] computeCache;
	}

	double SumNode::getLogVal() const {
		if (maxMode_)
			return maxLogVal_;
		else
			return logVal_;
	}

	double SumNode::getVal() const {
		if (maxMode_)
			return exp(maxLogVal_);
		else
			return exp(logVal_);
	}

	unsigned SumNode::getMaxIdx() const {
		return maxIdx_;
	}

	double SumNode::getLog_Dval_Dchild(Node *child) {
		if (!isChild(child))
			throw std::runtime_error("SumNode::getLog_Dval_Dchild(Node *child): child node not found");

		NodeId childId = child->getId();
		if (logWeights_.count(childId) == 0)
			throw std::runtime_error("SumNode::getLog_Dval_Dchild(Node *child): child weight not found");
		return logWeights_[childId];
	}

	void SumNode::connectChild(Node *child) {
		connectChild(child, 0.0);
	}

	void SumNode::connectChild(Node *child, double weight) {
		if ((isChild(child)) || (child->getId() == getId()))
			return;
		if (weight < 0)
			throw std::runtime_error("SumNode::connectChild(Node *child, double weight): weight must be positive");

		Node::connectChild(child);
		setLogWeight(child, log(weight));
	}

	void SumNode::disconnectChild(Node *child) {
		Node::disconnectChild(child);
		NodeId childId = child->getId();
		logWeights_.erase(childId);
	}

	void SumNode::setLogWeight(Node *child, double logWeight) {
		if (!isChild(child))
			throw std::runtime_error("error in SumNode::setLogWeight(Node *child, double weight): not a child.");

		size_t childId = child->getId();
		logWeights_[childId] = logWeight;
	}

	void SumNode::setLogWeights(const DoubleVec &logWeights) {
		size_t numChildren = getNumChildren();

		if (numChildren != logWeights.size())
			throw std::runtime_error("error in SumNode::setLogWeights(DoubleVec logWeights): "
			"weights vector must contain getNumChildren() entries.");

		for (unsigned k = 0; k < numChildren; k++)
			setLogWeight(getChild(k), logWeights[k]);
	}

	double SumNode::getLogWeight(Node *child) {
		if (!isChild(child))
			throw std::runtime_error("SumNode::getLogWeight(Node *child): not a child.");
		return logWeights_[child->getId()];
	}

	// stupid way to draw from Dirichlet distribution with uniform alpha=1
	void SumNode::initWeightsRand() {
		double sum = 0.0;
		DoubleVec r;

		r.resize(getNumChildren());
		for (unsigned k = 0; k < getNumChildren(); k++) {

			double gsample1 = 0.0;
			for (unsigned l = 0; l < 48; l++)
				gsample1 += ((double)rand() / double(RAND_MAX)) - 0.5;
			gsample1 = gsample1 * 0.5;

			double gsample2 = 0.0;
			for (unsigned l = 0; l < 48; l++)
				gsample2 += ((double)rand() / double(RAND_MAX)) - 0.5;
			gsample2 = gsample2 * 0.5;

			r[k] = gsample1 * gsample1 + gsample2 * gsample2;
			sum += r[k];
		}
		for (unsigned k = 0; k < logWeights_.size(); k++)
			r[k] = log(r[k] / sum);

		setLogWeights(r);
	}

	void SumNode::setSumMode() {
		maxMode_ = false;
	}

	void SumNode::setMaxMode() {
		maxMode_ = true;
	}


	////////////////////
	/// Product Node ///
	////////////////////
	ProdNode::ProdNode() : Node(PROD_NODE) {
		logVal_ = 0.0;
		zeroChildFlag_ = 0;
	}

	ProdNode::~ProdNode() {
	}

	void ProdNode::eval() {
		logVal_ = 0.0;
		zeroChildFlag_ = 0;

		for (unsigned k = 0; k < getNumChildren(); k++) {
			double childLogVal = getChild(k)->getLogVal();
			if (isnan(childLogVal))
				logVal_ = nan;
			else if (childLogVal > -inf)
				logVal_ += childLogVal;
			else {
				zeroChildFlag_++;
				if (zeroChildFlag_ == 2) {
					logVal_ = -inf;
					break;
				}
			}
		}
	}

	double ProdNode::getLogVal() const {
		if (zeroChildFlag_ == 0)
			return logVal_;
		else
			return -inf;
	}

	double ProdNode::getVal() const {
		if (zeroChildFlag_ == 0)
			return exp(logVal_);
		else
			return 0;
	}

	double ProdNode::getLog_Dval_Dchild(Node *child) {
		if (!isChild(child))
			throw std::runtime_error(" ProdNode::getLog_Dval_Dchild(Node *child): child node not found");

		if (zeroChildFlag_ == 0)
			return logVal_ - child->getLogVal();
		else if ((zeroChildFlag_ == 1) && (!(child->getLogVal() > -inf)))
			return logVal_;
		else
			return -inf;
	}


	//////////////////
	/// Value Node ///
	//////////////////
	ValNode::ValNode(double val) : Node(VAL_NODE) {
		logVal_ = log(val);
	}

	ValNode::ValNode() : Node(VAL_NODE) {
		logVal_ = -inf;
	}

	ValNode::~ValNode() {
	}

	void ValNode::eval() {
	}

	double ValNode::getLog_Dval_Dchild(Node *child) {
		return -inf;
	}

	void ValNode::setVal(double val) {
		val_ = val;
	}

	double ValNode::getVal() const {
		return val_;
	}

	void ValNode::setLogVal(double logVal) {
		val_ = exp(logVal);
	}

	double ValNode::getLogVal() const {
		return log(val_);
	}


	//////////////////////
	/// Indicator Node ///
	//////////////////////
	double IndicatorNode::tol_ = 1e-12;

	IndicatorNode::IndicatorNode() : Node(INDICATOR_NODE) {
		logVal_ = 0.0;
		indValues_.clear();
	}

	IndicatorNode::~IndicatorNode() {
	}

	void IndicatorNode::eval() {
		logVal_ = 0.0;
		size_t numChildren = getNumChildren();

		for (unsigned k = 0; k < numChildren; k++) {
			Node* child = getChild(k);

			// when child is a val node, whose value is set to nan -> interpret as missing value
			if ((child->getType() == VAL_NODE) && (isnan(child->getVal())))
				continue;

			if (fabs(indValues_[child->getId()] - child->getVal()) > tol_) {
				logVal_ = -inf;
				break;
			}
		}
	}

	double IndicatorNode::getLog_Dval_Dchild(Node *child) {
		return logVal_;
	}

	void IndicatorNode::connectChild(Node *child) {
		connectChild(child, 0.0);
	}

	void IndicatorNode::connectChild(Node *child, double indValue) {
		if ((isChild(child)) || (child->getId() == getId()))
			return;

		Node::connectChild(child);
		setIndValue(child, indValue);
	}

	void IndicatorNode::disconnectChild(Node *child) {
		Node::disconnectChild(child);
		indValues_.erase(child->getId());
	}

	double IndicatorNode::getLogVal() const {
		return logVal_;
	}

	double IndicatorNode::getVal() const {
		if (logVal_ == -inf)
			return 0.0;
		else if (logVal_ == 0.0)
			return 1.0;
		else
			throw std::runtime_error("IndicatorNode::getVal(): inconsistent value");
	}

	void IndicatorNode::setIndValue(Node *child, double indValue) {
		if (!isChild(child))
			throw std::runtime_error("error in IndicatorNode::setIndValue(Node *child, double indValue): not a child.");

		indValues_[child->getId()] = indValue;
	}

	double IndicatorNode::getIndValue(Node *child) {
		if (!isChild(child))
			throw std::runtime_error("error in IndicatorNode::getIndValue(Node *child): not a child.");

		return indValues_[child->getId()];
	}

	void IndicatorNode::setTolerance(double tol) {
		tol_ = tol;
	}

	double IndicatorNode::getTolerance() {
		return tol_;
	}


	/////////////////////
	/// Gaussian Node ///
	/////////////////////
	const double GaussianNode::logSqrt2Pi = log(sqrt(2 * 3.141592653589793238462643383279));

	GaussianNode::GaussianNode() : Node(GAUSSIAN_NODE) {
		means_.clear();
		sigmas_.clear();
		mode_ = PDF_MODE;
	}

	GaussianNode::~GaussianNode() {
	}

	void GaussianNode::eval() {
		logVal_ = 0.0;
		size_t numChildren = getNumChildren();
		if (numChildren == 0)
			return;

		if (mode_ == PDF_MODE) {
			for (unsigned k = 0; k < numChildren; k++) {
				Node* child = getChild(k);

				// when child is a val node, whose value is set to nan -> interpret as missing value
				if ((child->getType() == VAL_NODE) && (isnan(child->getVal())))
					continue;

				NodeId childId = child->getId();
				logVal_ += (-log(sigmas_[childId]) - logSqrt2Pi) -
					(0.5 * pow(((child->getVal() - means_[childId]) / sigmas_[childId]), 2));

			}
		}
		else if (mode_ == CDF_MODE) {
			for (unsigned k = 0; k < numChildren; k++) {
				Node* child = getChild(k);

				// when child is a val node, whose value is set to nan -> interpret as missing value
				if ((child->getType() == VAL_NODE) && (isnan(child->getVal())))
					continue;

				NodeId childId = child->getId();
				double tmp = 0.5 * (1 + erf((child->getVal() - means_[childId]) / sqrt(2 * pow(sigmas_[childId], 2))));
				if (tmp > 0)
					logVal_ += log(tmp);
				else {
					logVal_ = -inf;
					break;
				}
			}
		}
		else
			throw std::runtime_error("GaussianNode::eval(): unknown mode");
	}

	double GaussianNode::getLog_Dval_Dchild(Node *child) {
		if (!isChild(child))
			throw std::runtime_error("GaussianNode::getLog_Dval_Dchild(Node *child): child node not found");

		if (isnan(child->getVal()))
			return nan;

		double logDval;
		NodeId childId = child->getId();

		if (mode_ == PDF_MODE)
			logDval = (-log(sigmas_[childId]) - logSqrt2Pi) -
			(0.5 * pow(((child->getVal() - means_[childId]) / sigmas_[childId]), 2));
		else {
			logDval = 0.5 * (1 + erf((child->getVal() - means_[childId]) / sqrt(2 * pow(sigmas_[childId], 2))));
			if (logDval > 0)
				logDval = log(logDval);
			else
				logDval = -inf;
		}
		return logDval;
	}

	void GaussianNode::connectChild(Node *child) {
		connectChild(child, 0.0, 1.0);
	}

	void GaussianNode::connectChild(Node *child, double mean, double sigma) {
		if ((isChild(child)) || (child->getId() == getId()))
			return;

		Node::connectChild(child);
		setMean(child, mean);
		setSigma(child, sigma);
	}

	void GaussianNode::disconnectChild(Node *child) {
		Node::disconnectChild(child);
		NodeId childId = child->getId();
		means_.erase(childId);
		sigmas_.erase(childId);
	}

	void GaussianNode::setMean(Node *child, double mean) {
		if (!isChild(child))
			throw std::runtime_error("error in GaussianNode::setMean(Node *child, double mean): not a child.");
		means_[child->getId()] = mean;
	}

	double GaussianNode::getMean(Node *child) {
		if (!isChild(child))
			throw std::runtime_error("error in GaussianNode::getMean(Node *child): not a child.");
		return means_[child->getId()];
	}

	void GaussianNode::setMeansRand(double minVal, double maxVal) {
		for (unsigned k = 0; k < getNumChildren(); k++) {
			setMean(getChild(k), minVal + (maxVal - minVal) * (((double)rand()) / ((double)RAND_MAX)));
		}
	}

	void GaussianNode::setSigma(Node *child, double sigma) {
		if (!isChild(child))
			throw std::runtime_error("error in GaussianNode::setSigma(Node *child, double sigma): not a child.");
		if (sigma <= 0)
			throw std::runtime_error("error in GaussianNode::setSigma(Node *child, double sigma): "
			"sigma has to be larger than 0.");
		sigmas_[child->getId()] = sigma;
	}

	double GaussianNode::getSigma(Node *child) {
		if (!isChild(child))
			throw std::runtime_error("error in GaussianNode::getSigma(Node *child): not a child.");
		return sigmas_[child->getId()];
	}

	void GaussianNode::setSigmasRand(double minVal, double maxVal) {
		for (unsigned k = 0; k < getNumChildren(); k++) {
			setSigma(getChild(k), minVal + (maxVal - minVal) * (((double)rand()) / ((double)RAND_MAX)));
		}
	}
	
	void GaussianNode::setMode(int mode) {
		mode_ = mode;
	}

}
