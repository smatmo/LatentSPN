/*
Utils.hpp, Utils.cpp

Some utility classes (not very organized).


Robert Peharz
October 2016
*/


#ifndef UTILS_HPP
#define	UTILS_HPP

#include <map>
#include <vector>
#include <string>
#include <sstream>
#include <stdexcept>
#include <limits>
#include <math.h>
#include <algorithm>
#include <fstream>
#include "TypeDefs.hpp"

namespace SPN {
	
	const int SORTASCEND = 0;
	const int SORTDESCEND = 1;

	std::string int2str(int i);
	std::string unsigned2str(unsigned i);
	std::string double2str(double i);
	bool parseOptions(int argc, char **argv, int offset, std::map<std::string, std::string> &optionMap);

	bool doubleUnsignedPairSmaller___(std::pair<double, unsigned> x, std::pair<double, unsigned> y);
	bool doubleUnsignedPairLarger___(std::pair<double, unsigned> x, std::pair<double, unsigned> y);

	void sortWithIdx(std::vector<double> &X, std::vector<unsigned> &idx, int asc_des_flag = SORTASCEND);

	std::vector<double> mean(const std::vector<std::vector<double> > &mat, unsigned dim);

	std::vector<double> std(const std::vector<std::vector<double> > &mat, unsigned dim);

	void cov(const std::vector<std::vector<double> > &data,
		std::vector<std::vector<double> > &C);

	void corrcoef(const std::vector<std::vector<double> > &data,
		std::vector<std::vector<double> > &C);

	std::string vec2string(const std::vector<double> &vec);
	std::string vec2string(const std::vector<unsigned> &vec);
	std::string vec2string(const std::vector<int> &vec);
	std::string vec2string(const std::vector<bool> &vec);

	std::vector<unsigned> string2uintvec(std::string s);

	template <class T>
	bool isMat(const std::vector<std::vector<T> > &mat) {
		if (mat.empty())
			return true;

		if (mat[0].empty())
			return false;

		bool ret = true;
		for (unsigned k = 1; k < mat.size(); k++) {
			if (mat[k].size() != mat[0].size()) {
				ret = false;
				break;
			}
		}
		return ret;
	}

	std::string mat2string(std::vector<std::vector<double> > mat);
	std::string mat2string(std::vector<std::vector<unsigned> > vec);
	std::string mat2string(std::vector<std::vector<int> > vec);

	void vecAbsInPlace(std::vector<double> &vec);
	std::vector<double> vecAbs(const std::vector<double> &vec);

	void matAbsInPlace(std::vector<std::vector<double> > &mat);
	std::vector<std::vector<double> > matAbs(const std::vector<std::vector<double> > &mat);

	double logSumExp(double v1, double v2);
	double logSumExp(const std::vector<double> &vec);
	std::vector<double> logSumExp(const std::vector<std::vector<double> > &mat, unsigned dim);

	template <class T>
	std::vector<T> vecUnique(const std::vector<T> &vec) {
		std::vector<T> ret = vec;
		typename std::vector<T>::iterator it;

		std::sort(ret.begin(), ret.end());
		it = std::unique(ret.begin(), ret.end());
		ret.resize(std::distance(ret.begin(), it));

		return ret;
	}

	template <class T>
	std::vector<T> vecConcatenate(const std::vector<T> &vec1, const std::vector<T> &vec2) {
		std::vector<T> ret;

		for (unsigned k = 0; k < vec1.size(); k++)
			ret.push_back(vec1[k]);

		for (unsigned k = 0; k < vec2.size(); k++)
			ret.push_back(vec2[k]);

		return ret;
	}

	template <class T>
	std::vector<unsigned> vecEqual(const std::vector<T> &vec, const std::vector<T> &compVec) {
		std::vector<unsigned> ret;
		ret.resize(vec.size());
		bool hit;

		for (unsigned k = 0; k < vec.size(); k++) {
			hit = false;
			for (unsigned l = 0; l < compVec.size(); l++) {
				if (vec[k] == compVec[l]) {
					hit = true;
					break;
				}
			}
			if (hit)
				ret[k] = 1;
			else
				ret[k] = 0;
		}
		return ret;
	}

	template <class T>
	std::vector<unsigned> vecEqual(const std::vector<T> &vec, T compVal) {
		return vecEqual(vec, std::vector<T >(1, compVal));
	}

	template <class T>
	std::vector<unsigned> vecFind(const std::vector<T> &vec) {
		std::vector<unsigned> ret;
		ret.clear();

		for (unsigned k = 0; k < vec.size(); k++)
			if (vec[k] != 0)
				ret.push_back(k);

		return ret;
	}

	template <class T>
	void vecFilter(std::vector<T> &vec, const std::vector<bool> &keepFlags) {
		if (vec.size() != keepFlags.size())
			throw std::runtime_error("vecErase: vec.size() != keepFlags.size()");

		size_t ridx = 0;
		size_t widx = 0;

		while (ridx < vec.size()) {
			if (keepFlags[ridx]) {
				vec[widx] = vec[ridx];
				widx++;
			}
			ridx++;
		}
		vec.resize(widx);
	}

	template <class T>
	size_t findFirstInSortedVec_geq(const std::vector<T> &vec, T val) {
		if (vec.empty())
			return 0;
		if (vec.back() < val)
			return vec.size();
		if (vec.front() >= val)
			return 0;

		size_t l = 0;
		size_t r = vec.size() - 1;
		size_t ret;

		while (r > l) {
			ret = (r + l) / 2;
			if (vec[ret] >= val) {
				if ((ret == 0) || (vec[ret - 1] < val))
					return ret;
				r = ret - 1;
			}
			else {
				if ((ret == vec.size() - 1) || (vec[ret + 1] >= val))
					return ret + 1;
				l = ret + 1;
			}
		}

		if (r < l)
			throw std::logic_error("r < l");

		return r;
	}

	template <class T>
	std::vector<T> Vec_Ind(const std::vector<T> &vec, std::vector<unsigned> ind) {
		if (vec.size() != ind.size())
			throw std::runtime_error("Vec_Ind: sizes inconsistent");

		std::vector<T> ret;
		ret.clear();

		for (unsigned k = 0; k < vec.size(); k++)
			if (ind[k])
				ret.push_back(vec[k]);

		return ret;
	}

	std::vector<std::vector<double> > Mat_colInd(
		const std::vector<std::vector<double> > &mat,
		std::vector<unsigned> ind);
	std::vector<std::vector<double> > Mat_rowInd(
		const std::vector<std::vector<double> > &mat,
		std::vector<unsigned> ind);

	template <class T>
	std::vector<std::vector<T> > Mat_transpose(const std::vector<std::vector<T> > &mat) {
		if (!isMat(mat))
			throw std::runtime_error("Mat_transpose: mat is not a matrix");

		std::vector<std::vector<T> > ret;
		if (mat.empty()) {
			ret.clear();
			return ret;
		}

		ret.resize(mat[0].size());
		for (unsigned k = 0; k < ret.size(); k++) {
			ret[k].resize(mat.size());
			for (unsigned l = 0; l < mat.size(); l++)
				ret[k][l] = mat[l][k];
		}

		return ret;
	}

	bool parseOptions(int argc, char **argv, int offset, std::map<std::string, std::string> &optionMap);

}

#endif	/* UTILS_HPP */

