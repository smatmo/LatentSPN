
#include "EM.hpp"

namespace SPN {

	EM::EM() {
		maxIter_ = 20;
		stopping_objMinRelChange_ = 1e-6;
		updateGaussianMeans_ = true;
		updateGaussianSigmas_ = false;
		updateWeights_ = true;
		minSigma_ = 1e-6;
		maxValLLNotIncreased_ = 5;
	}

	// store current params in "slot" idx
	void EM::storeparams(Network* spn, unsigned idx) {
		store_logweights_[idx].clear();
		store_means_[idx].clear();
		store_sigmas_[idx].clear();

		for (size_t k = 0; k < spn->getNumLayers(); k++) {
			for (size_t l = 0; l < spn->getLayer(k)->getNumNodes(); l++) {
				Node* N = spn->getLayer(k)->getNode(l);

				if (N->getType() == Node::PROD_NODE) {
					// 
				}
				else if (N->getType() == Node::SUM_NODE) {
					for (size_t cc = 0; cc < N->getNumChildren(); cc++) {
						Node* C = N->getChild(cc);
						store_logweights_[idx][N->getId()][C->getId()] = ((SumNode*)N)->getLogWeight(C);
					}
				}
				else if (N->getType() == Node::GAUSSIAN_NODE) {
					for (size_t cc = 0; cc < N->getNumChildren(); cc++) {
						Node* C = N->getChild(cc);
						store_means_[idx][N->getId()][C->getId()] = ((GaussianNode*)N)->getMean(C);
						store_sigmas_[idx][N->getId()][C->getId()] = ((GaussianNode*)N)->getSigma(C);
					}
				}
				else if (N->getType() == Node::VAL_NODE) {
					// 
				}
				else {
					throw std::runtime_error("EM::storeparams: Node type not implemented");
				}
			}
		}
	}

	// load current params from "slot" idx
	void EM::loadparams(Network* spn, unsigned idx) {
		if ((store_logweights_.count(idx) == 0) || (store_means_.count(idx) == 0) || (store_sigmas_.count(idx) == 0))
			throw std::runtime_error("EM::loadparams: index error");

		size_t numSums = 0;
		size_t numGauss = 0;

		for (size_t k = 0; k < spn->getNumLayers(); k++) {
			for (size_t l = 0; l < spn->getLayer(k)->getNumNodes(); l++) {
				Node* N = spn->getLayer(k)->getNode(l);

				if (N->getType() == Node::PROD_NODE) {
					//  
				}
				else if (N->getType() == Node::SUM_NODE) {
					if (store_logweights_.at(idx).at(N->getId()).size() != N->getNumChildren())
						throw std::runtime_error("EM::loadparams: inconsistent stash (num sum-children)");

					numSums++;

					for (size_t cc = 0; cc < N->getNumChildren(); cc++) {
						Node* C = N->getChild(cc);
						double lw = store_logweights_.at(idx).at(N->getId()).at(C->getId());
						((SumNode*)N)->setLogWeight(C, lw);
					}
				}
				else if (N->getType() == Node::GAUSSIAN_NODE) {
					if (store_means_.at(idx).at(N->getId()).size() != N->getNumChildren())
						throw std::runtime_error("EM::loadparams: inconsistent stash (num gauss-children)");
					if (store_sigmas_.at(idx).at(N->getId()).size() != N->getNumChildren())
						throw std::runtime_error("EM::loadparams: inconsistent stash (num gauss-children)");

					numGauss++;

					for (size_t cc = 0; cc < N->getNumChildren(); cc++) {
						Node* C = N->getChild(cc);

						double mu = store_means_.at(idx).at(N->getId()).at(C->getId());
						double sigma = store_sigmas_.at(idx).at(N->getId()).at(C->getId());
						((GaussianNode*)N)->setMean(C, mu);
						((GaussianNode*)N)->setSigma(C, sigma);
					}
				}
				else if (N->getType() == Node::VAL_NODE) {
					// 
				}
				else {
					throw std::runtime_error("EM::loadparams: Node type not implemented");
				}
			}
		}

		if (numSums != store_logweights_.at(idx).size())
			throw std::runtime_error("EM::loadparams: inconsistent slot (num sums)");

		if (numGauss != store_means_.at(idx).size())
			throw std::runtime_error("EM::loadparams: inconsistent slot (num Gauss)");

		if (numGauss != store_sigmas_.at(idx).size())
			throw std::runtime_error("EM::loadparams: inconsistent slot (num Gauss)");
	}

	//
	// acutal learning routine
	//
	void EM::learn(Network* spn, DataSource *trainData, DataSource *valData) {

		if (spn == NULL)
			throw std::runtime_error("EM::learn: no reference set to SPN object.");
		if (trainData == NULL)
			throw std::runtime_error("EM::learn: no reference set to DATA source.");
		if (spn->getTopLayer()->getNumNodes() != 1)
			throw std::runtime_error("EM::learn: non-unique root");

		Layer* valLayer = spn->getLayer(0);

		if (valLayer->getType() != Layer::VAL_LAYER)
			throw std::runtime_error("EM: first layer is not a val layer");

		// sumCounts represent the soft counts in the log domain, for numerical stability.
		// logMax keeps track of the maximum log-soft count.
		// This is in essence the log-sum-exp trick.
		UMap<NodeId, UMap<NodeId, double> > sumCounts;
		UMap<NodeId, UMap<NodeId, double> > logMax;

		// vector of all sum nodes in the SPN
		NodeVec sums;

		// vector of vector of Gaussian nodes: gaussians[i] are the Gauss nodes for input dimension i.
		std::vector<NodeVec> gaussians;

		// this represents the posterior of the Gaussians' distribution selector, conditioned 
		// on evidence.
		UMap<NodeId, DoubleVec> logGaussPost;



		// collect sum nodes
		for (size_t lc = 0; lc < spn->getNumLayers(); lc++) {
			for (size_t nc = 0; nc < spn->getLayer(lc)->getNumNodes(); nc++) {
				Node* N = spn->getLayer(lc)->getNode(nc);
				if (N->getType() == Node::SUM_NODE)
					sums.push_back(N);
			}
		}

		// collect Gauss nodes
		for (size_t k = 0; k < valLayer->getNumNodes(); k++) {
			NodeVec tmpGauss;
			for (size_t p = 0; p < valLayer->getNode(k)->getNumParents(); p++) {
				Node* parent = valLayer->getNode(k)->getParent(p);
				if (parent->getType() != Node::GAUSSIAN_NODE)
					throw std::runtime_error("EM: parent of val node is not a Gaussian");
				if (parent->getNumChildren() != 1)
					throw std::runtime_error("EM: Gaussian has more than 1 child -- sorry, currently not implemented for more children");
				tmpGauss.push_back(parent);
				logGaussPost[parent->getId()].resize(trainData->getNumSamples());
			}
			gaussians.push_back(tmpGauss);
		}


		//
		// likelihoods
		//
		double LL = 0.0;
		double valLL = 0.0;
		double bestValLL = -inf;
		unsigned numValLLNotIncreased = 0;

		// train set
		for (size_t sampleC = 0; sampleC < trainData->getNumSamples(); sampleC++) {
			spn->setVals(trainData->getSample(sampleC));
			spn->eval();
			LL += spn->getTopLayer()->getNode(0)->getLogVal();
		}
		trainHist_.push_back(LL);

		// validation set
		if (valData != NULL) {
			for (size_t sampleC = 0; sampleC < valData->getNumSamples(); sampleC++) {
				spn->setVals(valData->getSample(sampleC));
				spn->eval();
				valLL += spn->getTopLayer()->getNode(0)->getLogVal();
			}
			valHist_.push_back(valLL);
			bestValLL = valLL;
		}

		// display
		std::cout << std::endl;
		std::cout << "LL train: " << std::setprecision(4) << std::fixed << LL << "   " << "LL val: " << valLL << std::endl;

		storeparams(spn, 0);

		//
		// EM main loop
		//
		for (size_t iterC = 0; iterC < maxIter_; iterC++) {

			for (size_t sampleC = 0; sampleC < trainData->getNumSamples(); sampleC++) {

				// evaluate SPN
				spn->setVals(trainData->getSample(sampleC));
				spn->eval();
				spn->backprop();

				// root value
				double logS = spn->getTopLayer()->getNode(0)->getLogVal();

				//
				// sums
				//
				for (size_t s = 0; s < sums.size(); s++) {
					SumNode* sum = (SumNode*)sums[s];
					NodeId sumId = sum->getId();
					double logDS = sum->getLogDerivative();
					for (size_t c = 0; c < sum->getNumChildren(); c++) {
						Node* child = sum->getChild(c);
						NodeId childId = child->getId();

						// posterior S(Z, Y | e), where Z is the latent variable associated with sum and
						// Y is the switching parent
						double inc = logDS + child->getLogVal() + sum->getLogWeight(child) - logS;

						// this accumulates inc in the log-domain (log-sum-exp trick)
						if (sampleC == 0) {
							sumCounts[sumId][childId] = 1.0;
							logMax[sumId][childId] = inc;
						}
						else {
							if (inc > logMax[sumId][childId]) {
								double sumCount = log(sumCounts[sumId][childId]) + logMax[sumId][childId];
								sumCounts[sumId][childId] = exp(sumCount - inc) + 1.0;
								logMax[sumId][childId] = inc;
							}
							else {
								if (logMax[sumId][childId] > -inf)
									sumCounts[sumId][childId] += exp(inc - logMax[sumId][childId]);
							}
						}
					}
				}

				//
				// Gaussians
				//
				for (size_t k = 0; k < gaussians.size(); k++) {
					for (size_t g = 0; g < gaussians[k].size(); g++)
						logGaussPost[gaussians[k][g]->getId()][sampleC] = gaussians[k][g]->getLogDerivative() + gaussians[k][g]->getLogVal() - logS;
				}
			}


			//
			// set sum weights
			//
			if (updateWeights_) {
				for (size_t s = 0; s < sums.size(); s++) {
					SumNode* sum = (SumNode*)sums[s];
					NodeId sumId = sum->getId();
					DoubleVec w;
					for (size_t c = 0; c < sum->getNumChildren(); c++) {
						Node* child = sum->getChild(c);
						NodeId childId = child->getId();
						w.push_back(log(sumCounts[sumId][childId]) + logMax[sumId][childId]);
					}

					double lse = logSumExp(w);

					if (lse > -inf)
						for (size_t c = 0; c < sum->getNumChildren(); c++)
							sum->setLogWeight(sum->getChild(c), w[c] - lse);
					else
						for (size_t c = 0; c < sum->getNumChildren(); c++)
							sum->setLogWeight(sum->getChild(c), -log((double)sum->getNumChildren()));
				}
			}

			// 
			// set Gaussians
			//
			if ((updateGaussianMeans_) || (updateGaussianSigmas_)) {
				for (size_t k = 0; k < gaussians.size(); k++) {
					DoubleVec dataCol;
					for (size_t n = 0; n < trainData->getNumSamples(); n++)
						dataCol.push_back(trainData->getSample(n)[k]);

					for (size_t g = 0; g < gaussians[k].size(); g++) {
						double lse = logSumExp(logGaussPost[gaussians[k][g]->getId()]);
						double mu = 0.0;
						double sumW = 0.0;

						// update means
						if (updateGaussianMeans_) {
							for (size_t n = 0; n < dataCol.size(); n++) {
								double w = exp(logGaussPost[gaussians[k][g]->getId()][n] - lse);
								sumW += w;
								double x = 0.0;

								// missing value? 
								// -> take the expected value of the first natural parameter 'x'.
								// this is simply the current mean, since we are integrating over the whole real line
								//
								// else, take observed value
								if (isnan(dataCol[n]))
									x = ((GaussianNode*)gaussians[k][g])->getMean(gaussians[k][g]->getChild(0));
								else
									x = dataCol[n];

								// accumulate
								mu += w * x;
							}
							mu = mu / sumW;

							((GaussianNode*)gaussians[k][g])->setMean(gaussians[k][g]->getChild(0), mu);
						}

						// update sigmas
						if (updateGaussianSigmas_) {
							if (!updateGaussianMeans_)
								mu = ((GaussianNode*)gaussians[k][g])->getMean(gaussians[k][g]->getChild(0));
							double sigma = 0.0;
							sumW = 0.0;

							for (size_t n = 0; n < dataCol.size(); n++) {
								double w = exp(logGaussPost[gaussians[k][g]->getId()][n] - lse);
								sumW += w;
								double d = 0.0;

								// missing value? 
								// -> take the expected value of the second natural parameter 'x^2'.
								// this is simply the current variance (plus acutally mean^2), 
								// since we are integrating over the whole real line
								//
								// else, take observed value
								if (isnan(dataCol[n]))
									d = ((GaussianNode*)gaussians[k][g])->getSigma(gaussians[k][g]->getChild(0));
								else
									d = dataCol[n] - mu;

								// accumulate
								sigma += w * pow(d, 2);
							}
							sigma = sqrt(sigma / sumW);
							if (sigma < minSigma_)
								sigma = minSigma_;

							((GaussianNode*)gaussians[k][g])->setSigma(gaussians[k][g]->getChild(0), sigma);
						}
					}
				}
			}

			//
			// likelihoods
			//
			LL = 0.0;
			for (size_t sampleC = 0; sampleC < trainData->getNumSamples(); sampleC++) {
				spn->setVals(trainData->getSample(sampleC));
				spn->eval();
				LL += spn->getTopLayer()->getNode(0)->getLogVal();
			}
			trainHist_.push_back(LL);

			valLL = 0.0;
			if (valData != NULL) {
				for (size_t sampleC = 0; sampleC < valData->getNumSamples(); sampleC++) {
					spn->setVals(valData->getSample(sampleC));
					spn->eval();
					valLL += spn->getTopLayer()->getNode(0)->getLogVal();
				}
				valHist_.push_back(valLL);
			}

			//
			// Display likelihoods, check early stopping
			//
			std::cout << "LL train: " << std::setprecision(4) << std::fixed << LL;
			if (valData != NULL) {
				std::cout << "   " << "LL val: " << std::setprecision(4) << std::fixed << valLL;
				if (valLL > bestValLL) {
					bestValLL = valLL;
					numValLLNotIncreased = 0;
					storeparams(spn, 0);
					std::cout << " (* store parameters)";
				}
				else {
					numValLLNotIncreased++;
					std::cout << " (" << numValLLNotIncreased << "/" << maxValLLNotIncreased_ << ")";
					if (numValLLNotIncreased == maxValLLNotIncreased_) {
						std::cout << std::endl << "Val Loglikelihood has not increased for " << numValLLNotIncreased << " times, break" << std::endl;
						break;
					}
				}
			}
			std::cout << std::endl;

			// the likelihood must never decrease in a proper EM derivation
			// accept 1e-6 tolerance due to numerical issues
			if (LL < (trainHist_[iterC] - 1e-6)) {
				std::cout << "***************************" << std::endl;
				std::cout << "Likelihood decreased!" << std::endl;
				std::cout << "***************************" << std::endl;
			}
		}

		loadparams(spn, 0);
	}
}


