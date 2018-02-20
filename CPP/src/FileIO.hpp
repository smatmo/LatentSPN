/*
FileIO.hpp, FileIO.cpp

Input/Output routines 

Robert Peharz
October 2016
*/

#ifndef FILEIO
#define	FILEIO

#include "Network.hpp"
#include "NetworkFactory.hpp"
#include "PDStruct.hpp"
#include <string>
#include <iomanip>

namespace SPN {
	
	// save/load SPNs as text files
	void saveNetwork(Network *network, std::string fileName);
	void loadNetwork(Network *network, std::string fileName, NetworkFactory &factory);

	// load SPNs from files used by Poon & Domingos, UAI 2011.
	void loadNetworkPDformat(Network *network, std::string fileName, int xDim, int yDim, NetworkFactory &factory);

	// read a list of vectors (double, unsigned)
	void readVectorList(std::string fileName, DoubleVec2 &list);
	void readVectorList(std::string fileName, UnsignedVec2 &list);
	DoubleVec2 readVectorList(std::string fileName);

	// write a lit of vectors (double)
	void writeVectorList(std::string fileName, DoubleVec2 &list);

	// compile spn into a dot file, which can then be comiled into a PDF depicting the SPN.
	// for big SPNs, this becomes quickly crowded
	void writeDot(Node* subroot, std::string fileName, bool weightLabel, bool indbottom);
	
}

#endif //FILEIO