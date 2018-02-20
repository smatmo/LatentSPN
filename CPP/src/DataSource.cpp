
#include "DataSource.hpp"

namespace SPN {

	DataSource::DataSource() {
		dataRef_ = NULL;
		dataRefAllocated_ = false;
	}

	DataSource::DataSource(DoubleVec2 &dataRef) {
		feed(dataRef);
	}

	DataSource::~DataSource() {
		if (dataRefAllocated_)
			delete dataRef_;
	}

	void DataSource::feed(DoubleVec2 &dataRef) {
		if (dataRefAllocated_)
			delete dataRef_;
		dataRef_ = &dataRef;
		dataRefAllocated_ = false;
	}

	void DataSource::readFile(std::string fileName) {
		if (dataRefAllocated_)
			delete dataRef_;
		dataRef_ = new DoubleVec2;
		dataRefAllocated_ = true;
		readVectorList(fileName, *dataRef_);
		if (!isMat(*dataRef_))
			throw std::runtime_error("DataSource::readFile(std::string fileName): not a matrix");
	}

	size_t DataSource::getNumSamples() {
		return dataRef_->size();
	}

	DoubleVec &DataSource::getSample(size_t sampleIdx) {
		return (*dataRef_)[sampleIdx % dataRef_->size()];
	}

	DoubleVec &DataSource::operator() (size_t sampleIdx) {
		return getSample(sampleIdx);
	}

	double DataSource::operator() (size_t sampleIdx, size_t dimIdx) {
		return getSample(sampleIdx)[dimIdx];
	}

	void DataSource::makeDataCopy(DoubleVec2 &C) {
		C = (*dataRef_);
	}

}
