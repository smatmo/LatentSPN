/*
EM.hpp, EM.cpp

Implements the EM algorithm for SPNs.

Robert Peharz
October 2016
*/

#ifndef EM_HPP
#define	EM_HPP

#include "GraphTools.hpp"
#include "Network.hpp"
#include "DataSource.hpp"
#include <iostream>
#include <iomanip>

namespace SPN {

	class EM {
	private:

		// stores statistics for sum and Gauss nodes
		// first index is a "slot" index
		// 2nd, 3rd indices represent parent and child nodes
		UMap<unsigned, UMap<NodeId, UMap<NodeId, double>>> store_logweights_;
		UMap<unsigned, UMap<NodeId, UMap<NodeId, double>>> store_means_;
		UMap<unsigned, UMap<NodeId, UMap<NodeId, double>>> store_sigmas_;

		// store/load parameters in "slots"
		void storeparams(Network* spn, unsigned idx);
		void loadparams(Network* spn, unsigned idx);


	protected:

	public:
		
		EM();

		// max. number of iterations
		size_t maxIter_;

		// EM stops of improvment in log-likelihood is below this threshold
		double stopping_objMinRelChange_;

		// update weights/means/sigmas?
		bool updateGaussianMeans_;
		bool updateGaussianSigmas_;
		bool updateWeights_;

		// minimal value for sigmas
		double minSigma_;

		// stores the log-likelihood over iterations
		DoubleVec trainHist_;
		DoubleVec valHist_;

		// counter for early stopping
		unsigned maxValLLNotIncreased_;		
		
		// main routine
		void learn(Network* spn, DataSource *trainData, DataSource *valData);
	};

}

#endif	/* EM_HPP */

