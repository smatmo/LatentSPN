/*
Typedefs.hpp

Some convenient typedefs.


Robert Peharz
October 2016
*/


#ifndef TYPEDEFS
#define	TYPEDEFS

#include <cstddef>
#include <limits>
#include <vector>
#include <unordered_set>
#include <unordered_map>

#define UMap std::unordered_map
#define USet std::unordered_set

namespace SPN {
	typedef std::vector<double> DoubleVec;
	typedef std::vector<unsigned> UnsignedVec;
	typedef std::vector<size_t> Size_tVec;

	typedef std::vector<std::vector<double>> DoubleVec2;
	typedef std::vector<std::vector<unsigned>> UnsignedVec2;
	typedef std::vector<std::vector<size_t>> Size_tVec2;

	typedef std::vector<std::vector<std::vector<double>>> DoubleVec3;
	typedef std::vector<std::vector<std::vector<unsigned>>> UnsignedVec3;
	typedef std::vector<std::vector<std::vector<size_t>>> Size_tVec3;

	const double inf = std::numeric_limits<double>::infinity();
	const double nan = std::numeric_limits<double>::signaling_NaN();
	const double numericZero = 1e-300;
}

#endif // TYPEDEFS