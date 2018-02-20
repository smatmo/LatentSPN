/*
Node.hpp, Node.cpp, SumNode.cpp, ProdNode.cpp, GaussianNode.cpp, IndicatorNode.cpp, ValNode.cpp

NOTE: this code generally operates with POINTERS of objects.
For example, nodes are passed to functions simply via their pointers.
In that way, we effectively have Python-like objects.

The Node class is an base class of the nodes which are the atomic units in SPNs. 
We identify node objects uniquely in memory. Therefore, the copy constructor should not be 
called and is private. There is no public constructor, i.e. the Node class is abstract.
There is the protected Node(int type) constructor which has to be called
from derived node classes.

There are currently 5 types of nodes:
SumNode.......... representing weighted sums
ProdNode......... representing products
ValNode.......... representing 'values' or actually random variables
IndicatorNode.... representing indicator variables 
GaussianNode..... representing Gaussian distributions

The connectChild and connectParent allow you to connect nodes with each other.
The eval() method represents the functionality of the node, and has to be
implemented in derived classes.
The backprop() method performs the backpropagation algorithm at (*this) node
object. This method is already implemented.
However, the backprop method needs the log of the derivative of each parent
nodes after (*this) node. Therefore, the getLog_Dval_Dchild(Node *child)
method has to be implemented in derived classes.

Robert Peharz
October 2016
*/


#ifndef NODE_HPP
#define	NODE_HPP

#include <vector>
#include "TypeDefs.hpp"
#include "Utils.hpp"


namespace SPN {

	class Node;
	typedef std::vector<Node*> NodeVec;
	typedef size_t NodeId;
	typedef int NodeType;


	///////////////////////
	/// Base class Node ///
	///////////////////////
	class Node {

	private:
		// this stores the type of the node (sum, product, ...)
		NodeType type_;

		// unique id of the node
		static NodeId idCounter_;
		NodeId id_;

		// connections are implemented as double-linked lists
		NodeVec children_;
		NodeVec parents_;

		// for quick check is some node is a child
		USet<NodeId> childInd_;

		// we don't allow to directly copy nodes
		Node(const Node& node);

	protected:
		Node(NodeType type);

		// we represent values and derivatives in the log domain
		double logVal_;
		double logDerivative_;


	public:
		static const NodeType SUM_NODE = 0x01;
		static const NodeType PROD_NODE = 0x02;
		static const NodeType VAL_NODE = 0x03;
		static const NodeType INDICATOR_NODE = 0x04;
		static const NodeType GAUSSIAN_NODE = 0x05;

		virtual ~Node();

		virtual void eval() {
		}

		virtual double getLog_Dval_Dchild(Node *child) {
			throw std::runtime_error("double getLog_Dval_Dchild(Node *child): must not call");
			return 0.0;
		}

		// perform backprob at this node, assuming that the derivatives of the parents are already computed.
		// requires the local derivatives which are obtained from the parents by calling getLog_Dval_Dchild().
		virtual void backprop();

		// get value/derivative
		virtual double getLogVal() const;
		virtual double getLogDerivative() const;
		virtual double getVal() const;
		virtual double getDerivative() const;

		// useful for learning
		void forceLogVal(double logVal);

		// connect/disconnect
		virtual void connectChild(Node *child);
		virtual void connectParent(Node *parent);
		virtual void disconnectChild(Node *child);
		virtual void disconnectParent(Node *parent);

		// disconnect all parents/children
		void isolate();

		// number of children/parents
		size_t getNumChildren();
		size_t getNumParents();

		// get child or parent node with certain idx -- idx assigned in the order of connection
		Node* getChild(size_t idx);
		Node* getParent(size_t idx);

		// get children/parents
		NodeVec getChildren();
		NodeVec getParents();

		// checks if other node is a child/parent
		bool isChild(Node* node);
		bool isParent(Node* node);

		// check equality
		bool operator==(const Node &rigthNode) const;

		// get node id
		NodeId getId() const;

		// get node type
		int getType() const;
	};


	////////////////
	/// Sum Node ///
	////////////////
	class SumNode : public Node {
	private:
		// weights of outgoing edges, represented in the log domain
		UMap<NodeId, double> logWeights_;

	protected:
		// we also store the maximum in addition to the sum
		double maxLogVal_;
		unsigned maxIdx_;
		bool maxMode_;

	public:
		SumNode();
		virtual ~SumNode();

		virtual void eval();

		virtual double getLog_Dval_Dchild(Node *childId);
		virtual void connectChild(Node *child);
		virtual void connectChild(Node *child, double weight);
		virtual void disconnectChild(Node *child);

		virtual double getLogVal() const;
		virtual double getVal() const;

		virtual unsigned getMaxIdx() const;

		// set/get/init weights
		virtual void setLogWeight(Node *child, double logWeight);
		virtual void setLogWeights(const DoubleVec &logWeights);
		virtual double getLogWeight(Node *child);
		virtual void initWeightsRand();
		
		// sum nodes can be switched to be max nodes
		virtual void setSumMode();
		virtual void setMaxMode();
	};


	////////////////////
	/// Product Node ///
	////////////////////
	class ProdNode : public Node {
	protected:
		// this is for the efficient back-prop trick in products, described in
		// A. Darwiche, "A Differential Approach to Inference in Bayesian Networks", ACM, 2003.
		int zeroChildFlag_;

	public:
		ProdNode();
		virtual ~ProdNode();
		virtual void eval();
		virtual double getLogVal() const;
		virtual double getVal() const;
		virtual double getLog_Dval_Dchild(Node *child);
	};


	//////////////////
	/// Value Node ///
	//////////////////
	class ValNode : public Node {
	private:
		// A valNode represents inputs, which can also be negative
		// We use val_ instead of logVal_ here.
		double val_;

	public:
		ValNode();
		ValNode(double logVal);
		virtual ~ValNode();
		virtual void eval();
		virtual double getLog_Dval_Dchild(Node *child);

		virtual double getLogVal() const;
		virtual double getVal() const;
		virtual void setLogVal(double logVal);
		virtual void setVal(double val);
	};


	//////////////////////
	/// Indicator Node ///
	//////////////////////
	class IndicatorNode : public Node {
	protected:
		// for storing the indicator values of the child nodes
		UMap<NodeId, double> indValues_;
		// tolerance, child values can deviate from their indValue
		static double tol_;

	public:
		IndicatorNode();
		virtual ~IndicatorNode();
		virtual void eval();
		virtual double getLog_Dval_Dchild(Node *child);

		virtual void connectChild(Node *child);
		virtual void connectChild(Node *child, double indValue);
		virtual void disconnectChild(Node *child);

		virtual double getLogVal() const;
		virtual double getVal() const;

		virtual void setIndValue(Node *child, double indValue);
		virtual double getIndValue(Node *child);
		
		static void setTolerance(double tol);
		static double getTolerance();
	};


	/////////////////////
	/// Gaussian Node ///
	/////////////////////
	class GaussianNode : public Node {
	protected:
		// means/standard deviations o fthe child nodes
		UMap<NodeId, double> means_;
		UMap<NodeId, double> sigmas_;
		int mode_;

	public:
		// log(sqrt(2*pi))
		static const double logSqrt2Pi;
		static const int PDF_MODE = 1;
		static const int CDF_MODE = 2;

		GaussianNode();
		virtual ~GaussianNode();
		virtual void eval();
		virtual double getLog_Dval_Dchild(Node *child);

		virtual void connectChild(Node *child);
		virtual void connectChild(Node *child, double mean, double sigma);
		virtual void disconnectChild(Node *child);

		virtual void setMean(Node *child, double mean);
		virtual double getMean(Node *child);
		virtual void setMeansRand(double minVal, double maxVal);
		virtual void setSigma(Node *child, double sigma);
		virtual double getSigma(Node *child);
		virtual void setSigmasRand(double minVal, double maxVal);

		virtual void setMode(int mode);
	};


}

#endif	/* NODE_HPP */

