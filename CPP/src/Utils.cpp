#include "Utils.hpp"
#include <iostream>

namespace SPN {

	std::string int2str(int i) {
		std::stringstream convert;
		convert << i;
		if (!convert.good())
			throw std::runtime_error("std::string int2str(int i): error while converting");
		return convert.str();
	}

	std::string unsigned2str(unsigned i) {
		std::stringstream convert;
		convert << i;
		if (!convert.good())
			throw std::runtime_error("std::string unsigned2str(unsigned i): error while converting");
		return convert.str();
	}

	std::string double2str(double i) {
		std::stringstream convert;
		convert << i;
		if (!convert.good())
			throw std::runtime_error("std::string double2str(double i): error while converting");
		return convert.str();
	}

	std::string vec2string(const std::vector<double> &vec) {
		std::string ret;
		if (vec.size() == 0) {
			ret = "[]";
			return ret;
		}
		ret = "[";
		for (unsigned k = 0; k < vec.size() - 1; k++)
			ret += double2str(vec[k]) + ", ";
		ret += double2str(vec[vec.size() - 1]) + "]";
		return ret;
	}

	std::string vec2string(const std::vector<unsigned> &vec) {
		std::string ret;
		if (vec.size() == 0) {
			ret = "[]";
			return ret;
		}
		ret = "[";
		for (unsigned k = 0; k < vec.size() - 1; k++)
			ret += unsigned2str(vec[k]) + ", ";
		ret += unsigned2str(vec[vec.size() - 1]) + "]";
		return ret;
	}

	std::string vec2string(const std::vector<int> &vec) {
		std::string ret;
		if (vec.size() == 0) {
			ret = "[]";
			return ret;
		}
		ret = "[";
		for (unsigned k = 0; k < vec.size() - 1; k++)
			ret += int2str(vec[k]) + ", ";
		ret += int2str(vec[vec.size() - 1]) + "]";
		return ret;
	}

	std::string vec2string(const std::vector<bool> &vec) {
		std::string ret;
		if (vec.size() == 0) {
			ret = "[]";
			return ret;
		}
		ret = "[";
		for (unsigned k = 0; k < vec.size() - 1; k++) {
			if (vec[k])
				ret += "1, ";
			else
				ret += "0, ";
		}
		if (vec[vec.size() - 1])
			ret += "1]";
		else
			ret += "0]";
		return ret;
	}

	std::vector<unsigned> string2uintvec(std::string s) {
		std::stringstream inputStream;
		std::vector<unsigned> ret;
		ret.clear();

		if (
			(s.length() < 2) ||
			(s[0] != '[') ||
			(s[s.length() - 1] != ']') ||
			(s.find_first_not_of(" \t[],0123456789") != std::string::npos))
			throw std::runtime_error(
			"std::vector<unsigned> string2uintvec(std::string string): "
			"malformatted string");

		s = s.substr(1, s.length() - 2);
		while (!s.empty()) {
			unsigned val;
			size_t pos = s.find_first_of(",");
			if (pos == std::string::npos) {
				inputStream.clear();
				inputStream.str(s);
				inputStream >> val;
				if (inputStream.fail())
					throw std::runtime_error(
					"std::vector<unsigned> string2uintvec(std::string string): "
					"malformatted string");
				ret.push_back(val);
				break;
			}
			else {
				inputStream.clear();
				inputStream.str(s.substr(0, pos));
				inputStream >> val;
				if (inputStream.fail())
					throw std::runtime_error(
					"std::vector<unsigned> string2uintvec(std::string string): "
					"malformatted string");
				ret.push_back(val);
				s = s.substr(pos + 1, std::string::npos);
			}
		}

		return ret;
	}

	std::string mat2string(std::vector<std::vector<double> > mat) {
		if (!isMat(mat))
			throw std::runtime_error("mat2string: input is not a correct matrix");

		std::string ret;
		if (mat.size() == 0) {
			ret = "[]";
			return ret;
		}

		ret = "[";
		for (unsigned k = 0; k < mat.size(); k++) {
			for (unsigned l = 0; l < mat[k].size(); l++) {
				ret += double2str(mat[k][l]);
				if (l + 1 < mat[k].size())
					ret += " ";
			}
			if (k + 1 < mat.size())
				ret += "\n";
		}
		ret += "]";

		return ret;
	}

	std::string mat2string(std::vector<std::vector<unsigned> > mat) {
		if (!isMat(mat))
			throw std::runtime_error("mat2string: input is not a correct matrix");

		std::string ret;
		if (mat.size() == 0) {
			ret = "[]";
			return ret;
		}

		ret = "[";
		for (unsigned k = 0; k < mat.size(); k++) {
			for (unsigned l = 0; l < mat[k].size(); l++) {
				ret += unsigned2str(mat[k][l]);
				if (l + 1 < mat[k].size())
					ret += " ";
			}
			if (k + 1 < mat.size())
				ret += "\n";
		}
		ret += "]";

		return ret;
	}

	std::string mat2string(std::vector<std::vector<int> > mat) {
		if (!isMat(mat))
			throw std::runtime_error("mat2string: input is not a correct matrix");

		std::string ret;
		if (mat.size() == 0) {
			ret = "[]";
			return ret;
		}

		ret = "[";
		for (unsigned k = 0; k < mat.size(); k++) {
			for (unsigned l = 0; l < mat[k].size(); l++) {
				ret += int2str(mat[k][l]);
				if (l + 1 < mat[k].size())
					ret += " ";
			}
			if (k + 1 < mat.size())
				ret += "\n";
		}
		ret += "]";

		return ret;
	}

	bool doubleUnsignedPairSmaller___(std::pair<double, unsigned> x, std::pair<double, unsigned> y) {
		return (x.first < y.first);
	}

	bool doubleUnsignedPairLarger___(std::pair<double, unsigned> x, std::pair<double, unsigned> y) {
		return (x.first > y.first);
	}

	void sortWithIdx(std::vector<double> &X, std::vector<unsigned> &idx, int asc_des_flag) {
		std::vector<std::pair<double, unsigned> > pairVec;
		pairVec.resize(X.size());
		for (unsigned k = 0; k < X.size(); k++) {
			pairVec[k].first = X[k];
			pairVec[k].second = k;
		}

		if (asc_des_flag == SORTASCEND)
			sort(pairVec.begin(), pairVec.end(), doubleUnsignedPairSmaller___);
		else
			sort(pairVec.begin(), pairVec.end(), doubleUnsignedPairLarger___);

		idx.resize(X.size());
		for (unsigned k = 0; k < X.size(); k++)
			idx[k] = 0;

		for (unsigned k = 0; k < X.size(); k++) {
			X[k] = pairVec[k].first;
			idx[k] = pairVec[k].second;
		}
	}

	std::vector<double> mean(const std::vector<std::vector<double> > &mat, unsigned dim) {
		if ((dim != 1) && (dim != 2))
			throw std::runtime_error("mean: dim has to be 1 or 2.");
		if (!isMat(mat))
			throw std::runtime_error("mean: mat is not a matrix.");

		std::vector<double> ret;

		if (mat.empty()) {
			ret.clear();
			return ret;
		}

		size_t N = mat.size();
		size_t D = mat[0].size();

		if (dim == 1) {
			ret.resize(D);
			for (unsigned d = 0; d < D; d++)
				ret[d] = 0.0;

			for (unsigned n = 0; n < N; n++)
				for (unsigned d = 0; d < D; d++)
					ret[d] += mat[n][d];

			for (unsigned d = 0; d < D; d++)
				ret[d] /= N;
		}
		else {
			if (D == 0) {
				ret.clear();
				return ret;
			}

			ret.resize(N);
			for (unsigned n = 0; n < N; n++)
				ret[n] = 0.0;

			for (unsigned n = 0; n < N; n++)
				for (unsigned d = 0; d < D; d++)
					ret[n] += mat[n][d];

			for (unsigned n = 0; n < N; n++)
				ret[n] /= D;
		}

		return ret;
	}

	std::vector<double> std(const std::vector<std::vector<double> > &mat, unsigned dim) {
		if ((dim != 1) && (dim != 2))
			throw std::runtime_error("mean: dim has to be 1 or 2.");
		if (!isMat(mat))
			throw std::runtime_error("mean: mat is not a matrix.");

		std::vector<double> ret;

		if (mat.empty()) {
			ret.clear();
			return ret;
		}

		std::vector<double> mu = mean(mat, dim);

		size_t N = mat.size();
		size_t D = mat[0].size();

		if (dim == 1) {
			ret.resize(D);
			for (unsigned d = 0; d < D; d++)
				ret[d] = 0.0;

			for (unsigned n = 0; n < N; n++)
				for (unsigned d = 0; d < D; d++)
					ret[d] += pow(mat[n][d] - mu[d], 2.0);

			for (unsigned d = 0; d < D; d++) {
				if (N > 1) {
					ret[d] /= (N - 1);
					ret[d] = sqrt(ret[d]);
				}
				else
					ret[d] = 0.0;
			}
		}
		else {
			if (D == 0) {
				ret.clear();
				return ret;
			}

			ret.resize(N);
			for (unsigned n = 0; n < N; n++)
				ret[n] = 0.0;

			for (unsigned n = 0; n < N; n++)
				for (unsigned d = 0; d < D; d++)
					ret[n] += pow(mat[n][d] - mu[n], 2.0);

			for (unsigned n = 0; n < N; n++) {
				if (N > 1) {
					ret[n] /= (D - 1);
					ret[n] = sqrt(ret[n]);
				}
				else
					ret[n] = 0.0;
			}
		}

		return ret;
	}

	void cov(const std::vector<std::vector<double> > &data,
		std::vector<std::vector<double> > &C) {

		if (data.size() < 1)
			throw std::logic_error("SPN::cov: data must contain at least 1 sample");

		size_t N = data.size();
		size_t D = data[0].size();

		for (unsigned n = 0; n < data.size(); n++)
			if (data[n].size() != D)
				throw std::logic_error("SPN::cov: data vectors must be of same length");

		if (D == 1) {
			C.resize(1);
			C[0].resize(1);
			C[0][0] = 1.0;
			return;
		}

		// means
		std::vector<double> means;
		means.resize(D);
		for (unsigned n = 0; n < data.size(); n++)
			for (unsigned d = 0; d < D; d++)
				means[d] += data[n][d];

		for (unsigned d = 0; d < D; d++)
			means[d] /= (double)N;

		// init C
		C.resize(D);
		for (unsigned d1 = 0; d1 < D; d1++) {
			C[d1].resize(D);
			for (unsigned d2 = 0; d2 < D; d2++)
				C[d1][d2] = 0.0;
		}

		std::vector<double> zeroMeanSample;
		zeroMeanSample.resize(D);
		for (unsigned n = 0; n < N; n++) {
			for (unsigned d = 0; d < D; d++)
				zeroMeanSample[d] = data[n][d] - means[d];

			for (unsigned d1 = 0; d1 < D; d1++) {
				for (unsigned d2 = d1; d2 < D; d2++) {
					C[d1][d2] += zeroMeanSample[d1] * zeroMeanSample[d2];
					if (d1 != d2)
						C[d2][d1] = C[d1][d2];
				}
			}
		}

		if (N > 1) {
			for (unsigned d1 = 0; d1 < D; d1++)
				for (unsigned d2 = 0; d2 < D; d2++)
					C[d1][d2] /= (double)(N - 1);
		}

		return;
	}

	void corrcoef(const std::vector<std::vector<double> > &data,
		std::vector<std::vector<double> > &C) {

		cov(data, C);

		std::vector<double> s;
		s.resize(C.size());
		for (unsigned k = 0; k < s.size(); k++)
			s[k] = sqrt(C[k][k]);

		for (unsigned d = 0; d < C.size(); d++) {
			for (unsigned d2 = 0; d2 < C[d].size(); d2++) {
				if ((s[d] == 0.0) || (s[d2] == 0.0))
					C[d][d2] = nan;
				else {
					if (d != d2)
						C[d][d2] /= (s[d] * s[d2]);
					else
						C[d][d2] = 1.0;
				}
			}
		}
	}

	/*
	void corrcoef(
	const std::vector<std::vector<double> > &Y,
	const std::vector<double> &meanY,
	const std::vector<double> &stdY,
	const std::vector<double> &X,
	double meanX,
	double stdX,
	std::vector<double> &C) {

	if (X.size() != Y.size())
	throw std::runtime_error("SPN::corrcoef: X and Y must contain the same number of samples");

	if (Y.size() < 1)
	throw std::logic_error("SPN::corrcoef: data must contain at least 1 sample");

	unsigned N = Y.size();
	unsigned D = Y[0].size();
	for (unsigned n = 0; n < N; n++)
	if (Y[n].size() != D)
	throw std::logic_error("SPN::corrcoef: vectors in Y must all be of same length");

	if (meanY.size() != D)
	throw std::logic_error("SPN::corrcoef: meanY != D");

	if (stdY.size() != D)
	throw std::logic_error("SPN::corrcoef: stdY != D");

	// initialize C
	C.resize(D);
	for (unsigned d = 0; d < D; d++)
	C[d] = 0.0;

	for (unsigned n = 0; n < N; n++)
	for (unsigned d = 0; d < D; d++)
	C[d] += (Y[n][d] - meanY[d]) * (X[n] - meanX);

	for (unsigned d = 0; d < D; d++)
	C[d] /= stdX * stdY[d];

	return;
	}
	*/

	void vecAbsInPlace(std::vector<double> &vec) {
		for (unsigned k = 0; k < vec.size(); k++)
			vec[k] = fabs(vec[k]);
	}

	std::vector<double> vecAbs(const std::vector<double> &vec) {
		std::vector<double> ret = vec;
		vecAbsInPlace(ret);
		return ret;
	}

	void matAbsInPlace(std::vector<std::vector<double> > &mat) {
		for (unsigned k = 0; k < mat.size(); k++)
			for (unsigned l = 0; l < mat[k].size(); l++)
				mat[k][l] = fabs(mat[k][l]);
	}

	std::vector<std::vector<double> > matAbs(const std::vector<std::vector<double> > &mat) {
		std::vector<std::vector<double> > ret = mat;
		matAbsInPlace(ret);
		return ret;
	}

	double logSumExp(double v1, double v2) {
		if ((v1 == -inf) && (v2 == -inf))
			return -inf;
		if ((v1 == inf) && (v2 == inf))
			return inf;

		if (v1 < v2)
			return log(exp(v1 - v2) + 1) + v2;
		else
			return log(exp(v2 - v1) + 1) + v1;
	}

	double logSumExp(const std::vector<double> &vec) {
		double ret = 0.0;
		double max = 0.0;

		for (unsigned k = 0; k < vec.size(); k++)
			if ((k == 0) || vec[k] > max)
				max = vec[k];

		if (max == inf)
			return inf;
		else if (max == -inf)
			return -inf;
		else {
			for (unsigned k = 0; k < vec.size(); k++)
				ret += exp(vec[k] - max);
		}

		if (ret == 0.0)
			return -inf;

		return log(ret) + max;
	}

	std::vector<double> logSumExp(const std::vector<std::vector<double> > &mat, unsigned dim) {
		if ((dim != 1) && (dim != 2))
			throw std::runtime_error("logSumExp: dim has to be 1 or 2.");
		if (!isMat(mat))
			throw std::runtime_error("logSumExp: mat is not a matrix.");

		std::vector<double> ret;
		if (mat.empty()) {
			ret.clear();
			return ret;
		}

		if (dim == 1) {
			std::vector<double> max;
			max.resize(mat[0].size(), -inf);
			ret.resize(mat[0].size(), 0.0);

			for (unsigned k = 0; k < mat.size(); k++)
				for (unsigned l = 0; l < mat[k].size(); l++)
					if (mat[k][l] > max[l])
						max[l] = mat[k][l];
			for (unsigned k = 0; k < mat.size(); k++)
				for (unsigned l = 0; l < mat[k].size(); l++)
					ret[l] += exp(mat[k][l] - max[l]);

			for (unsigned l = 0; l < ret.size(); l++)
				ret[l] = log(ret[l]) + max[l];
		}
		else {
			std::vector<double> max;
			max.resize(mat.size(), -inf);
			ret.resize(mat.size(), 0.0);

			for (unsigned k = 0; k < mat.size(); k++)
				for (unsigned l = 0; l < mat[k].size(); l++)
					if (mat[k][l] > max[k])
						max[k] = mat[k][l];
			for (unsigned k = 0; k < mat.size(); k++)
				for (unsigned l = 0; l < mat[k].size(); l++)
					ret[k] += exp(mat[k][l] - max[k]);

			for (unsigned k = 0; k < ret.size(); k++)
				ret[k] = log(ret[k]) + max[k];
		}

		return ret;
	}

	std::vector<std::vector<double> > Mat_colInd(
		const std::vector<std::vector<double> > &mat,
		std::vector<unsigned> ind) {

		std::vector<std::vector<double> > ret;

		if (mat.empty()) {
			if (!ind.empty())
				throw std::runtime_error("Mat_colInd: ind inconsistent");
			ret.clear();
			return ret;
		}

		for (unsigned k = 1; k < mat.size(); k++)
			if (mat[k].size() != mat[0].size())
				throw std::runtime_error("Mat_colInd: mat is not a matrix");

		if (mat[0].size() != ind.size())
			throw std::runtime_error("Mat_colInd: ind inconsistent");

		for (unsigned k = 0; k < ind.size(); k++)
			if ((ind[k] != 0) && (ind[k] != 1))
				throw std::runtime_error("Mat_colInd: ind must contain only 0 or 1");

		std::vector<unsigned> idx = vecFind(ind);
		ret.resize(mat.size());
		for (unsigned k = 0; k < mat.size(); k++) {
			ret[k].resize(idx.size());
			for (unsigned l = 0; l < idx.size(); l++)
				ret[k][l] = mat[k][idx[l]];
		}

		return ret;
	}

	std::vector<std::vector<double> > Mat_rowInd(
		const std::vector<std::vector<double> > &mat,
		std::vector<unsigned> ind) {

		std::vector<std::vector<double> > ret;

		if (mat.empty()) {
			if (!ind.empty())
				throw std::runtime_error("Mat_colInd: ind inconsistent");
			ret.clear();
			return ret;
		}

		for (unsigned k = 1; k < mat.size(); k++)
			if (mat[k].size() != mat[0].size())
				throw std::runtime_error("Mat_colInd: mat is not a matrix");

		if (mat.size() != ind.size())
			throw std::runtime_error("Mat_colInd: ind inconsistent");

		for (unsigned k = 0; k < ind.size(); k++)
			if ((ind[k] != 0) && (ind[k] != 1))
				throw std::runtime_error("Mat_colInd: ind must contain only 0 or 1");

		std::vector<unsigned> idx = vecFind(ind);
		ret.clear();
		for (unsigned k = 0; k < idx.size(); k++)
			ret.push_back(mat[idx[k]]);

		return ret;
	}

	bool parseOptions(int argc, char **argv, int offset,
		std::map<std::string, std::string> &optionMap) {
		std::string option;
		std::string value;
		size_t pos;

		optionMap.clear();

		for (int k = offset; k < argc; k++) {
			option = std::string(argv[k]);
			pos = option.find_first_of("=");

			if (pos == std::string::npos) {
				if (argc < k + 2)
					return false;

				value = std::string(argv[++k]);
				pos = value.find_first_of("=");
				if (pos + 1 == value.size()) {
					if (argc < k + 2)
						return false;
					value = std::string(argv[++k]);
				}
				else
					value = value.substr(pos + 1, value.size() - pos - 1);
			}
			else {
				if (pos + 1 == option.size()) {
					if (argc < k + 2)
						return false;
					value = std::string(argv[++k]);
				}
				else
					value = option.substr(pos + 1, option.size() - pos - 1);
				option = option.substr(0, pos);
			}
			optionMap[option] = value;
		}

		return true;
	}


}
