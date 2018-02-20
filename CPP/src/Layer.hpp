/*
Layer.hpp, Layer.cpp

NOTE: this code generally operates with POINTERS of objects.
For example, nodes are passed to functions simply via their pointers.
In that way, we effectively have Python-like objects.

Implements layers of nodes, in order to interpret SPNs as layered neural networks.
Like nodes, they are of a certain type. Currently there are

SumLayer
ProdLayer
ValLayer
IndicatorLayer
GaussianLayer

mirroring the corresponding node types. There are some manipulation routines for 
inserting/removing nodes. There are also array-like functions, which call a certain 
function for each contained node, e.g. eval() calls eval() for each contained node.

Robert Peharz
October 2016
*/


#ifndef LAYER_HPP
#define	LAYER_HPP

#include "TypeDefs.hpp"
#include "Node.hpp"

namespace SPN {

	class Layer;

	typedef std::vector<Layer*> LayerVec;
	typedef size_t LayerId;
	typedef int LayerType;
	
	////////////////////////
	/// Base class Layer ///
	////////////////////////
	class Layer {
	private:
		LayerType type_;
		LayerId id_;
		static LayerId idCounter_;
		
		Layer(const Layer&);

	protected:
		// list of node ptrs in this layer
		NodeVec nodes_;

		// for quick check, if node is contained in layer
		USet<NodeId> nodeSet_;

		Layer(LayerType type);

	public:
		static const LayerType SUM_LAYER = 0x01;
		static const LayerType PROD_LAYER = 0x02;
		static const LayerType VAL_LAYER = 0x03;
		static const LayerType INDICATOR_LAYER = 0x04;
		static const LayerType GAUSSIAN_LAYER = 0x05;

		virtual ~Layer();

		// call eval or backprop for all nodes in this layer
		void eval();
		void backprop();

		// insert and remove nodes
		void insertNode(Node *node);
		void insertNodes(NodeVec nodes);
		void removeNode(Node *node);
		void removeNodes(NodeVec nodes);

		// get methods
		virtual size_t getNumNodes();
		virtual NodeVec getNodes();
		virtual Node* getNode(size_t idx);
		virtual size_t getNodeIdx(Node* node);
		virtual bool containsNode(Node* node);
				
		Node* operator()(size_t idx);

		LayerType getType() const;
		LayerId getId() const;
	};


	/////////////////
	/// Sum Layer ///
	/////////////////
	class SumLayer : public Layer {
	
	public:
		SumLayer();
		virtual ~SumLayer();

		virtual SumNode* getNode(size_t idx);
		virtual void initWeightsRand();

		SumNode* operator() (size_t idx);
	};


	/////////////////////
	/// Product Layer ///
	/////////////////////
	class ProdLayer : public Layer {
	private:

	protected:

	public:
		ProdLayer();
		virtual ~ProdLayer();
		virtual ProdNode* getNode(size_t idx);
		ProdNode* operator() (size_t idx);
	};


	///////////////////
	/// Value Layer ///
	///////////////////
	class ValLayer : public Layer {
	private:

	protected:

	public:
		ValLayer();
		virtual ~ValLayer();

		virtual ValNode* getNode(size_t idx);
		virtual void setVals(DoubleVec values);
		virtual DoubleVec getVals();

		ValNode* operator() (size_t idx);
	};


	///////////////////////
	/// Indicator Layer ///
	///////////////////////
	class IndicatorLayer : public Layer {
	private:

	protected:

	public:
		IndicatorLayer();
		virtual ~IndicatorLayer();

		virtual IndicatorNode* getNode(size_t idx);

		IndicatorNode* operator() (size_t idx);
	};


	//////////////////////
	/// Gaussian Layer ///
	//////////////////////
	class GaussianLayer : public Layer {
	private:

	protected:

	public:
		GaussianLayer();
		virtual ~GaussianLayer();

		virtual GaussianNode* getNode(size_t idx);

		virtual void setMeansUniform(double meanVal);
		virtual void setSigmasUniform(double sigmaVal);
		virtual void setMeansRand(double minVal, double maxVal);
		virtual void setSigmasRand(double minVal, double maxVal);		

		GaussianNode* operator() (size_t idx);
	};
}


#endif	/* LAYER_HPP */

