/*
NetworkFactory.hpp, NetworkFactory.cpp

A factory class to generate parts of SPNs, i.e. Node, Layer and Network objects.
This is in fact just a memory manager, as it generate objects on demand and takes 
care that they are properly deleted as soon as the NetworkFactory is deleted.


Robert Peharz
October 2016
*/


#ifndef NETWORKFACTORY_HPP
#define	NETWORKFACTORY_HPP

#include "Network.hpp"

namespace SPN {

	class NetworkFactory {
	private:
		UMap<NetworkId, Network*> generatedNetworkPtrs_;
		UMap<LayerId, Layer*> generatedLayerPtrs_;
		UMap<NodeId, Node* > generatedNodePtrs_;
		
	public:
		NetworkFactory();
		virtual ~NetworkFactory();
		
		Network* generateNetwork();
		void deleteNetwork(Network *network);
		void deleteNetworkDeep(Network *network);
		Layer* generateLayer(LayerType type, size_t numNodes);
		void deleteLayer(Layer *layer);
		void deleteLayerDeep(Layer *layer);
		NodeVec generateNodes(NodeType type, size_t numNodes);
		Node* generateNode(NodeType type);
		void deleteNode(Node *node);

		Network* copyNetwork(Network *network);		
		
		void clear();
	};

}

#endif	/* NETWORKFACTORY_HPP */

