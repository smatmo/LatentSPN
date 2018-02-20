
#include "PDStruct.hpp"

namespace SPN {

	// Partition
	PDPartition::PDPartition() {
		parentRectangle_ = NULL;
		childRectangles_.clear();
	}


	// Rectangle
	PDRectangle::PDRectangle() {
		nodePtrs_.clear();
		prodNodePtrs_.clear();
		parentPartitions_.clear();
		childPartitions_.clear();
		x_ = -1;
		y_ = -1;
		width_ = -1;
		height_ = -1;
		level_ = -1;
	}

	bool PDRectangle::operator<(const PDRectangle &rhs) const {
		if (x_ < rhs.x_)
			return true;
		else if (x_ > rhs.x_)
			return false;
		if (y_ < rhs.y_)
			return true;
		else if (y_ > rhs.y_)
			return false;
		if (width_ < rhs.width_)
			return true;
		else if (width_ > rhs.width_)
			return false;
		if (height_ < rhs.height_)
			return true;
		else if (height_ > rhs.height_)
			return false;
		return false;
	}

	bool PDRectangle::operator==(const PDRectangle &rhs) const {
		return ((x_ == rhs.x_) && (y_ == rhs.y_) && (width_ == rhs.width_) && (height_ == rhs.height_));
	}

	bool PDRectangle::operator!=(const PDRectangle &rhs) const {
		return ((x_ != rhs.x_) || (y_ != rhs.y_) || (width_ != rhs.width_) || (height_ != rhs.height_));
	}

	void PDRectangle::setLevelRecursive(int level) {
		if (level_ >= level)
			return;
		level_ = level;
		for (unsigned k = 0; k < childPartitions_.size(); k++) {
			for (unsigned l = 0; l < childPartitions_[k]->childRectangles_.size(); l++) {
				childPartitions_[k]->childRectangles_[l]->setLevelRecursive(level + 1);
			}
		}
	}


	// Product node ID
	PDProdNodeId::PDProdNodeId() {

	}

	PDProdNodeId::PDProdNodeId(PDRectangle &r1, PDRectangle &r2, unsigned idx1, unsigned idx2) {
		set(r1, r2, idx1, idx2);
	}

	PDRectangle PDProdNodeId::getRectangle1() {
		PDRectangle r;
		r.x_ = x1_;
		r.y_ = y1_;
		r.width_ = width1_;
		r.height_ = height1_;
		return r;
	}

	unsigned PDProdNodeId::getIdx1() {
		return idx1_;
	}

	PDRectangle PDProdNodeId::getRectangle2() {
		PDRectangle r;
		r.x_ = x2_;
		r.y_ = y2_;
		r.width_ = width2_;
		r.height_ = height2_;
		return r;
	}

	unsigned PDProdNodeId::getIdx2() {
		return idx2_;
	}

	void PDProdNodeId::set(PDRectangle &r1, PDRectangle &r2, unsigned idx1, unsigned idx2) {
		if (r1 < r2) {
			x1_ = r1.x_;
			y1_ = r1.y_;
			width1_ = r1.width_;
			height1_ = r1.height_;
			idx1_ = idx1;
			x2_ = r2.x_;
			y2_ = r2.y_;
			width2_ = r2.width_;
			height2_ = r2.height_;
			idx2_ = idx2;
		}
		else {
			x1_ = r2.x_;
			y1_ = r2.y_;
			width1_ = r2.width_;
			height1_ = r2.height_;
			idx1_ = idx2;
			x2_ = r1.x_;
			y2_ = r1.y_;
			width2_ = r1.width_;
			height2_ = r1.height_;
			idx2_ = idx1;
		}
	}

	bool PDProdNodeId::operator<(const PDProdNodeId &rhs) const {
		if (x1_ < rhs.x1_)
			return true;
		else if (x1_ > rhs.x1_)
			return false;
		if (x2_ < rhs.x2_)
			return true;
		else if (x2_ > rhs.x2_)
			return false;
		if (y1_ < rhs.y1_)
			return true;
		else if (y1_ > rhs.y1_)
			return false;
		if (y2_ < rhs.y2_)
			return true;
		else if (y2_ > rhs.y2_)
			return false;
		if (width1_ < rhs.width1_)
			return true;
		else if (width1_ > rhs.width1_)
			return false;
		if (width2_ < rhs.width2_)
			return true;
		else if (width2_ > rhs.width2_)
			return false;
		if (height1_ < rhs.height1_)
			return true;
		else if (height1_ > rhs.height1_)
			return false;
		if (height2_ < rhs.height2_)
			return true;
		else if (height2_ > rhs.height2_)
			return false;
		if (idx1_ < rhs.idx1_)
			return true;
		else if (idx1_ > rhs.idx1_)
			return false;
		if (idx2_ < rhs.idx2_)
			return true;
		else if (idx2_ > rhs.idx2_)
			return false;
		return false;
	}

	bool PDProdNodeId::operator==(const PDProdNodeId &rhs) const {
		if ((x1_ == rhs.x1_) && (x2_ == rhs.x2_) && (y1_ == rhs.y1_) && (y2_ == rhs.y2_) &&
			(width1_ == rhs.width1_) && (width2_ == rhs.width2_) &&
			(height1_ == rhs.height1_) && (height2_ == rhs.height2_) &&
			(idx1_ == rhs.idx1_) && (idx2_ == rhs.idx2_))
			return true;
		else
			return false;
	}

	bool PDProdNodeId::operator!=(const PDProdNodeId &rhs) const {
		return !(*this == rhs);
	}

}


