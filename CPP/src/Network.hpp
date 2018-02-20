/*
Network.hpp, Network.cpp

NOTE: this code generally operates with POINTERS of objects.
For example, nodes are passed to functions simply via their pointers.
In that way, we effectively have Python-like objects.

Implements feedforward sum-product networks, organized as a stack of layers.
The input layer (ValLayer) is assumed to be layer 0, the output layer is the top layer.
There are some manipulation routines like inserting/removing layers. There are also
array-like functions, which call a certain function for each contained layer, 
e.g. eval() calls eval() for each contained layer, bottom-to-top, while backprop()
calls backprop() for each contained layer, top-to-bottom.

Robert Peharz
October 2016
*/


#ifndef NETWORK_HPP
#define	NETWORK_HPP

#include "TypeDefs.hpp"
#include "Node.hpp"
#include "Layer.hpp"
#include "GraphTools.hpp"
#include <queue>

namespace SPN {

	class Network;
	typedef std::vector<Network*> NetworkVec;
	typedef size_t NetworkId;

	class Network {
	private:
		NetworkId id_;
		static NetworkId idCounter_;

	protected:
		// list of layers
		LayerVec layers_;

		// result flags of performStructureCheck()
		bool structComplete_;
		bool structDecomposable_;
		bool structFeedforward_;

	public:
		Network();
		virtual ~Network();

		// insert layer at a certain position
		virtual void insertLayer(Layer* layer, size_t idx);
		virtual void removeLayer(Layer* layer);
		// insert layer on top
		virtual void addLayer(Layer* layer);

		// get methods
		virtual size_t getNumLayers();
		virtual LayerVec getLayers();
		virtual Layer* getLayer(size_t layerIdx);
		virtual Layer* getTopLayer();
		virtual size_t getLayerIdx(Layer* layer);
		virtual Node* getNode(size_t layerIdx, size_t idx);
		virtual bool containsLayer(Layer *layer);
		virtual size_t getLayerIdxOfNode(Node *node);

		// set values in the bottom layer
		virtual void setVals(const DoubleVec &values);
		virtual DoubleVec getVals();

		// eval the network from bottom to top
		virtual void eval();

		// backprop from top to bottom
		virtual void backprop();

		// inference 
		virtual void MPEinference(size_t rootidx = 0, bool dirty = true);
		
		virtual void maxBacktracking();
		//virtual std::vector<double> MaxBacktrackingAssignSums();
		virtual void setMPEcorrectionWeights();

		// for all sum nodes in the network, switch between sum and max mode
		virtual void setSumMode();
		virtual void setMaxMode();

		// parameter initialization
		virtual void initWeightsRand();
		virtual void setGaussianMeansRand(double minVal, double maxVal);
		virtual void setGaussianSigmasRand(double minVal, double maxVal);
		virtual void setGaussianMeansUniform(double meanVal);
		virtual void setGaussianSigmasUniform(double sigmaVal);

		// remove all layers
		virtual void clear();

		// remove parents with only one child
		virtual void removeSingletonParents(NodeVec &removedNodes, LayerVec &removedLayers);
		
		// renormalize
		//virtual void normalizeWeights();
		
		// count number of weights, Gaussians, etc.
		virtual size_t getNumWeights();
		virtual size_t getNumGaussians();
		virtual size_t getNumSums();
		virtual size_t getNumProds();

		// for structure check: first call performStructureCheck(), 
		// then call isComplete(), isDecomposable(), or isFeedforward()
		virtual void performStructureCheck();
		virtual bool isComplete();
		virtual bool isDecomposable();
		virtual bool isFeedforward();

		Layer* operator()(size_t layerIdx);
		Node* operator()(size_t layerIdx, size_t idx);

		NetworkId getId();
	};

}


#endif	/* NETWORK_HPP */

