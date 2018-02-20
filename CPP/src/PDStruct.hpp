/*
PDStruct.hpp, PDStruct.cpp

Some classes to represent the PD architecture.

PDRectangle
Represents a rectangle of an image (or more generally, a rectangular data array, e.g. a spectrogram).
A PDRectangle has child partitions and parent partitions, represented via class PDPartition (see below).
A PDRectangle has the coordinates of the rectangle it represents, x, y, width and height.
Furthermore, it contains a list of nodes (nodePtrs_), which have the rectangle as scope (sums or input 
distributions), and a list of product nodes (prodNodePtrs_) which are children of the nodes in nodePtrs_ 
(if they are sums).

PDPartition
Represents a partitions of a PDRectangle into several PDRectangles (here always two).
A PDPartition has a parent rectangle and a list of child rectangles. The child rectangles must be 
non-overlapping and their union must be the parent rectangle.

In summary: PDRectangles and PDPartitions form a bi-partitite doubly connect list.

PDProdNodeID:
Represents a virtual id of product nodes: In the PD architecture, products are uniquely identified by 
the two sub-rectangles they link and a node index within these sub-rectangles. There can be a huge number 
of (potential) product nodes. Representing them virtually via their id makes learning tractable. In the 
orginal code by Poon & Domingos, this id is encoded by an integer. As noted by Poon & Domingos, this will 
overflow for larger structures. Therefore, we use a separate class here.


Robert Peharz
October 2016
*/


#ifndef PDSTRUCT_HPP
#define	PDSTRUCT_HPP

#include "Node.hpp"
#include "TypeDefs.hpp"
#include "Utils.hpp"

namespace SPN {

	class PDRectangle;
	class PDPartition;

	typedef std::vector<PDRectangle*> PDRectangleVec;
	typedef std::vector<PDPartition*> PDPartitionVec;

	class PDRectangle {
	public:
		unsigned x_, y_, width_, height_;
		int level_;
		NodeVec prodNodePtrs_;
		std::map<unsigned, Node* > nodePtrs_;
		PDPartitionVec parentPartitions_;
		PDPartitionVec childPartitions_;

		unsigned MAPNodeIdx_;
		double MAPNodeLogVal_;

		PDRectangle();
		void setLevelRecursive(int level);

		bool operator<(const PDRectangle &rhs) const;
		bool operator==(const PDRectangle &rhs) const;
		bool operator!=(const PDRectangle &rhs) const;
	};

	class PDPartition {
	public:
		PDRectangle* parentRectangle_;
		PDRectangleVec childRectangles_;

		PDPartition();
	};

	class PDProdNodeId {
		friend class PDLearning;
	private:
		unsigned x1_;
		unsigned y1_;
		unsigned width1_;
		unsigned height1_;
		unsigned x2_;
		unsigned y2_;
		unsigned width2_;
		unsigned height2_;
		unsigned idx1_;
		unsigned idx2_;
	public:
		PDRectangle getRectangle1();
		unsigned getIdx1();
		PDRectangle getRectangle2();
		unsigned getIdx2();
		void set(PDRectangle &r1, PDRectangle &r2, unsigned idx1, unsigned idx2);
		PDProdNodeId(PDRectangle &r1, PDRectangle &r2, unsigned idx1, unsigned idx2);
		PDProdNodeId();
		bool operator<(const PDProdNodeId &rhs) const;
		bool operator==(const PDProdNodeId &rhs) const;
		bool operator!=(const PDProdNodeId &rhs) const;
	};

}

#endif	/* PDSTRUCT_HPP */

