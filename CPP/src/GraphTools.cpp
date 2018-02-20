
#include "GraphTools.hpp"

namespace SPN {

	void listDescendants(Node* subroot, NodeVec &desc) {

		USet<NodeId> seen;
		std::queue<Node*> queue;

		desc.clear();

		queue.push(subroot);
		seen.insert(subroot->getId());

		while (queue.empty() == false) {
			Node* nextNode = queue.front();
			queue.pop();
			desc.push_back(nextNode);

			for (unsigned k = 0; k < nextNode->getNumChildren(); k++) {
				Node* child = nextNode->getChild(k);
				if (seen.count(child->getId()) == 0) {
					queue.push(child);
					seen.insert(child->getId());
				}
			}
		}

	}
	

	void listAncestors(Node* subroot, NodeVec &anc) {
		USet<NodeId> seen;
		std::queue<Node*> queue;

		anc.clear();

		queue.push(subroot);
		seen.insert(subroot->getId());

		while (queue.empty() == false) {
			Node* nextNode = queue.front();
			queue.pop();
			anc.push_back(nextNode);

			for (unsigned k = 0; k < nextNode->getNumParents(); k++) {
				Node* parent = nextNode->getParent(k);
				if (seen.count(parent->getId()) == 0) {
					queue.push(parent);
					seen.insert(parent->getId());
				}
			}
		}
	}

	void keepNodesOfType(NodeVec &nodes, NodeType type) {
		NodeVec ret;
		for (size_t k = 0; k < nodes.size(); k++)
			if (nodes[k]->getType() == type)
				ret.push_back(nodes[k]);
		ret.swap(nodes);
	}

	bool topOrder(NodeVec &nodes, NodeVec &order, Node* root) {

		order.clear();
		order.reserve(nodes.size());

		if (nodes.size() == 0)
			return true;

		USet<NodeId> innodes;
		UMap<NodeId, unsigned> numPar;

		for (unsigned k = 0; k < nodes.size(); k++)
			innodes.insert(nodes[k]->getId());

		for (unsigned k = 0; k < nodes.size(); k++)
			for (unsigned p = 0; p < nodes[k]->getNumParents(); p++)
				if (innodes.count(nodes[k]->getParent(p)->getId()))
					numPar[nodes[k]->getId()]++;

		if (root != NULL) {
			if (innodes.count(root->getId()) == 0)
				throw std::runtime_error("topOrder(NodeVec &nodes, bool isSubgraph, Node* root): provided root is not in nodes");
			if (numPar[root->getId()] != 0)
				throw std::runtime_error("topOrder(NodeVec &nodes, bool isSubgraph, Node* root): provided node (root) is not a root");
		}
		else {
			bool foundroot = false;
			for (unsigned k = 0; k < nodes.size(); k++) {
				if (numPar[nodes[k]->getId()] == 0) {
					root = nodes[k];
					foundroot = true;
					break;
				}
			}
			if (!foundroot) {
				order.clear();
				return false;
			}
		}

		std::queue<Node*> queue;
		queue.push(root);
		UMap<NodeId, unsigned> parcount;

		while (!queue.empty()) {
			Node* nextNode = queue.front();
			queue.pop();
			order.push_back(nextNode);
			for (unsigned k = 0; k < nextNode->getNumChildren(); k++) {
				Node* child = nextNode->getChild(k);
				parcount[child->getId()]++;
				if (parcount[child->getId()] == numPar[child->getId()])
					queue.push(child);
			}
		}

		if (nodes.size() != order.size()) {
			order.clear();
			return false;
		}

		return true;
	}

	bool checkOrder(NodeVec & desc) {
		USet<NodeId> seen;
		USet<NodeId> suspect;

		for (unsigned k = 0; k < desc.size(); k++) {
			if (suspect.count(desc[k]->getId()))
				return false;
			for (unsigned p = 0; p < desc[k]->getNumParents(); p++) {
				if (seen.count(desc[k]->getParent(p)->getId()) == 0)
					suspect.insert(desc[k]->getParent(p)->getId());
			}
			seen.insert(desc[k]->getId());
		}
		return true;
	}
	
}