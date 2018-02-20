/*
GraphTools.hpp, GraphTools.cpp

Some useful routines operating on acylic directed graphs.

Robert Peharz
October 2016
*/


#ifndef GRAPHTOOLS_HPP
#define	GRAPHTOOLS_HPP

#include "Node.hpp"
#include "TypeDefs.hpp"
#include <queue>

namespace SPN {

	// takes a node and returns a list of its descendants in desc,
	// including the node itself
	void listDescendants(Node* subroot, NodeVec &desc);

	// takes a node and returns a list of its ancestors in anc,
	// including the node itself
	void listAncestors(Node* subroot, NodeVec &anc);

	// filters nodes of a certain type from a list
	void keepNodesOfType(NodeVec &nodes, NodeType type);

	// find a topological order of a node list
	bool topOrder(NodeVec &nodes, NodeVec &order, Node* root = NULL);

	// check if a node list is topologically ordered
	bool checkOrder(NodeVec & desc);
	
}

#endif // GRAPHTOOLS_HPP