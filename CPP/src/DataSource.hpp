/*
DataSource.hpp, DataSource.cpp

Represents a data source and is in fact just a wrapper for DoubleVec2

Robert Peharz
October 2016
*/

#ifndef DATASOURCE_HPP
#define	DATASOURCE_HPP

#include "TypeDefs.hpp"
#include "FileIO.hpp"
#include "Utils.hpp"

namespace SPN {

	class DataSource {
	private:

	protected:
		bool dataRefAllocated_;
		DoubleVec2 *dataRef_;

	public:
		DataSource();
		DataSource(DoubleVec2 &dataRef);
		virtual ~DataSource();

		void feed(DoubleVec2 &dataRef);
		void readFile(std::string fileName);

		virtual size_t getNumSamples();
		virtual DoubleVec &getSample(size_t idx);
		void makeDataCopy(DoubleVec2 &C);

		DoubleVec & operator() (size_t sampleIdx);
		double operator() (size_t sampleIdx, size_t dimIdx);
	};

}

#endif	/* DATASOURCE_HPP */

