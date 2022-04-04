// ***************************************************************************
// CloneCaller.h (c) 2021 Zhenhua Yu <qasim0208@163.com>
// Health Informatics Lab, Ningxia University
// All rights reserved.

#ifndef _CLONECALLER_H
#define _CLONECALLER_H

#include <vector>
#include <map>
#include <pthread.h>

#include "Matrix.h"

class Solution {
	public:
		Solution() {}
		Solution(int num_cluster, double alpha, double beta, Matrix<double>& pie, Matrix<int>& states) :
			num_cluster(num_cluster), alpha(alpha), beta(beta), pie(pie), states(states) {}
		
		int num_cluster; //number of clusters
		double ll; //log-likelihood
		double bic;
		double acc; //accuracy
		double alpha; //false positive rate
		double beta; //false negative rate
		Matrix<int> states; //mutation states of clusters
		Matrix<double> pie; //occurrence probability of each clone
		Matrix<double> post_probs; //posterior probability
};

class CloneCaller {
	private:
		Matrix<int> obsData;
		Matrix<int> uniqueData;
		Matrix<int> realData;
		
		vector<int> doublet_indxs;
		double missing_rate;
		
		Solution best_s;
		
		void preProcess();
		
		static void* inferClones(const void *arg);
		void fineTune(Solution& s);
		
		void obsProbs(Matrix<int>& states, Matrix<long double>& probs, double alpha, double beta);
		double getObsProb(int gt, int observed, double alpha, double beta);
		
		void printSolution(Solution& s);
		
		void loadMutationData();
		void loadRealData();
		
		void evalAccuracy(Solution& s);
		void refineSolution(Solution& s);
		
		void saveResults();
		
	public:
		CloneCaller();
	
		Matrix<int>& getObsData() {return obsData;}
		Matrix<int>& getUniqueData() {return uniqueData;}
		
		void loadData();
		void call();
};

#endif
