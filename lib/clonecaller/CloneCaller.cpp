// ***************************************************************************
// CloneCaller.cpp (c) 2021 zhenhua yu <qasim0208@163.com>
// Health Informatics Lab, Ningxia University
// All rights reserved.

#include <iostream>
#include <cstdio>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <algorithm>
#include <string>
#include <limits>
#include <array>
#include <tuple>

#include "CloneCaller.h"
#include "split.h"
#include "MyDefine.h"

using namespace std;

CloneCaller::CloneCaller() {
}

void CloneCaller::loadData() {
	loadMutationData();
	loadRealData();
}

void CloneCaller::loadMutationData() {
	string inputFile = config.getStringPara("input");
	ifstream ifs;
	ifs.open(inputFile.c_str());
	if(!ifs.is_open()) {
		cerr << "Error: cannot open file " << inputFile << endl;
		exit(-1);
	}
	
	FILE *fp;
	char buf[500];
	string cmd = "cat "+inputFile+" | wc -l";
	fp = popen(cmd.c_str(), "r");
	if(!fp) {
		cerr << "Error: cannot open file " << inputFile << endl;
		exit(-1);
	}
	char *p = fgets(buf, 500, fp);
	fclose(fp);
	int num_cell = atoi(buf);
	
	int i, j, k, indx;
	long line_num = 0;
	string line;
	int num_muta = 0;
	while(getline(ifs, line)) {
		line_num++;
		vector<string> fields = split(line, '\t');
		if(num_muta == 0) {
			num_muta = fields.size();
			obsData.resize(num_cell, num_muta, false);
		}
		if(num_muta != fields.size()) {
			cerr << "Error: malformed input file " << inputFile <<
                    ", inconsistent number of fields @line " << line_num << endl;
            cerr << buf << endl;
			exit(1);
		}
		for(i = 0; i < num_muta; i++) {
			indx = (line_num-1)*num_muta+i;
			k = atoi(fields[i].c_str());
			if(k != 0 && k != 1 && k != 3) {
				cerr << "unsupported genotype value " << k << " at line " << line_num << " in file " << inputFile << endl;
				exit(1);
			}
			obsData[indx] = k;
		}
	}
	ifs.close();
	assert(line_num == num_cell);
	
	cerr << "Total " << num_cell << " cells and " << num_muta << " mutations were loaded from file " << inputFile << endl;
}

void CloneCaller::loadRealData() {
	string inputFile = config.getStringPara("rinput");
	if(inputFile.empty()) {
		return;
	}
	
	FILE *fp;
	char buf[500];
	string cmd = "cat "+inputFile+" | wc -l";
	fp = popen(cmd.c_str(), "r");
	if(!fp) {
		cerr << "Error: cannot open file " << inputFile << endl;
		exit(-1);
	}
	char *p = fgets(buf, 500, fp);
	fclose(fp);
	int num_cell = atoi(buf)-1;
	
	int i, j, k, indx;
	long line_num = 0;
	string line;
	int num_muta = 0;
	
	ifstream ifs;
	ifs.open(inputFile.c_str());
	if(!ifs.is_open()) {
		cerr << "Error: cannot open file " << inputFile << endl;
		exit(-1);
	}
	getline(ifs, line);
	vector<string> fields = split(line, ':');
	if(fields.size() > 1) {
		fields = split(fields[1], '\t');
		for(i = 0; i < fields.size(); i++) {
			doublet_indxs.push_back(atoi(fields[i].c_str())-1);
		}
	}
	
	while(getline(ifs, line)) {
		line_num++;
		vector<string> fields = split(line, '\t');
		if(num_muta == 0) {
			num_muta = fields.size();
			realData.resize(num_cell, num_muta, false);
		}
		if(num_muta != fields.size()) {
			cerr << "Error: malformed input file " << inputFile <<
                    ", inconsistent number of fields @line " << line_num << endl;
            cerr << buf << endl;
			exit(1);
		}
		for(i = 0; i < num_muta; i++) {
			indx = (line_num-1)*num_muta+i;
			k = atoi(fields[i].c_str());
			realData[indx] = k;
		}
	}
	ifs.close();
	assert(line_num == num_cell);
	
	cerr << "real mutation data were loaded from file " << inputFile << endl;
}

void CloneCaller::preProcess() {
	int num_muta = obsData.getCOLS();
	int num_cell = obsData.getROWS();
	int i, j, k;
	
	Matrix<int> tmp = obsData;
	uniqueData = obsData;
	
	missing_rate = 0;
	for(i = 0; i < num_cell; i++) {
		for(j = 0; j < num_muta; j++) {
			k = obsData[i*num_muta+j];
			if(k == 3) {
				tmp[i*num_muta+j] = 0;
				missing_rate++;
			}
		}
	}
	missing_rate /= (num_cell*num_muta);
	
	int rindx = 0;
	for(i = 0; i < num_cell; i++) {
		int flag = 0;
		for(k = 0; k < rindx; k++) {
			for(j = 0; j < num_muta; j++) {
				if(uniqueData[k*num_muta+j] != tmp[i*num_muta+j]) {
					break;
				}
			}
			if(j == num_muta) {
				flag = 1;
				break;
			}
		}
		if(flag == 0) {
			for(j = 0; j < num_muta; j++) {
				uniqueData[rindx*num_muta+j] = tmp[i*num_muta+j];
			}
			rindx++;
		}
	}
	uniqueData.resize(rindx, num_muta, true);
	
}

void CloneCaller::call() {
	int i, j, k, n;
	
	preProcess();
	
	int num_muta = obsData.getCOLS();
	int num_cell = obsData.getROWS();
	int num_cell_u = uniqueData.getROWS();
	
	/*** find the optimal number of clusters based on BIC ***/
	int maxc = config.getIntPara("maxc");
	if(maxc <= 0) {
		maxc = min(500, num_cell/5);
	}
	else {
		maxc = min(maxc, num_cell/5);
	}
	
	int epoches = 100;
	
	//vector<Solution> solus;
	Solution* solus = new Solution[epoches];
	for(n = 0; n < epoches; n++) {
		/*** init paras ***/
		Matrix<double> pie(1, maxc);
		Matrix<int> states(maxc, num_muta);
		
		vector<int> indxs;
		for(i = 0; i < maxc; i++) {
			k = randomInteger(0, num_cell_u);
			if(indxs.empty()) {
				indxs.push_back(k);
			}
			else {
				vector<int>::iterator it = find(indxs.begin(), indxs.end(), k);
				while(it != indxs.end()) {
					k = randomInteger(0, num_cell_u);
					it = find(indxs.begin(), indxs.end(), k);
				}
				indxs.push_back(k);
			}
		}
		
		//vector<int> indxs = dkm::details::random_plusplus_2(data_for_cluster, maxc, n);
		uniqueData.Rows(indxs, states);
		
		pie.set(1.0/maxc);
		
		Solution s0(maxc, 0.01, 0.01, pie, states);
		s0.post_probs.resize(num_cell, maxc, false);
		solus[n] = s0;
		tp->add_work(&CloneCaller::inferClones, &solus[n], n);
	}
	tp->wait();
	
	//sort the solutions by log-likelihood
	for(i = 0; i < epoches-1; i++) {
		k = i;
		for(j = i+1; j < epoches; j++) {
			if(solus[j].ll > solus[k].ll) {
				k = j;
			}
		}
		if(k != i) {
			Solution tmp = solus[k];
			solus[k] = solus[i];
			solus[i] = tmp;
		}
	}
	
	//compare the top-10 solutions
	int best_s_indx = 0;
	for(i = 0; i < 10; i++) {
		cerr << "--------------- screening report -----------------" << endl;
		printSolution(solus[i]);
		cerr << "--------------- screening report -----------------" << endl;
		/*
		if(solus[i].score > solus[best_s_indx].score) {
			best_s_indx = i;
		}
		*/
		if(solus[i].bic < solus[best_s_indx].bic) {
			best_s_indx = i;
		}
	}
	best_s = solus[best_s_indx];
	delete[] solus;
	
	/*** best solution ***/
	cerr << "--------------- best solution -----------------" << endl;
	printSolution(best_s);
	cerr << "-----------------------------------------------" << endl;
	
	/*** save results ***/
	saveResults();
}

void* CloneCaller::inferClones(const void *arg) {
	Solution& s = *((Solution*) arg);
	
	int i, j, k, n;
	Matrix<int>& obsData = clonecaller.getObsData();
	int num_muta = obsData.getCOLS();
	int num_cell = obsData.getROWS();
	
	int num_cluster = s.num_cluster;
	double& alpha = s.alpha;
	double&	beta = s.beta;
	Matrix<int>& states = s.states;
	Matrix<double>& pie = s.pie;
	Matrix<double>& post_probs = s.post_probs;
	post_probs.resize(num_cell, num_cluster, false);
	
	long double p1, p0;
	double p, ll;
	Matrix<long double> probs(num_cell, num_cluster);
	Matrix<double> post_probs_pre(num_cell, num_cluster), post_probs_a(num_cell, num_cluster);
	Matrix<int> deltas(1, num_cell), deltas_1(1, num_cell);
	Matrix<double> states_p1(num_cluster, num_muta), states_p1_pre(num_cluster, num_muta);
	Matrix<double> pie_a(1, num_cluster);
	
	double eps = numeric_limits<double>::epsilon();
	double max_alpha = 0.1, max_beta = 0.99;
	long iter = 0, max_iter = config.getIntPara("max_iter");
	long iter_max;
	while(iter < max_iter) {
		clonecaller.obsProbs(states, probs, alpha, beta);
		
		// posterior probability calculation 
		ll = 0;
		for(i = 0; i < num_cell; i++) {
			long double sum = 0;
			for(k = 0; k < num_cluster; k++) {
				sum += probs[i*num_cluster+k]*pie[k];
			}
			for(k = 0; k < num_cluster; k++) {
				post_probs[i*num_cluster+k] = probs[i*num_cluster+k]*pie[k]/sum;
			}
			ll += logl(sum);
		}
		
		// sample Delta
		post_probs.cumsum(post_probs_a);
		map<int, vector<int>> cell_assigns;
		for(i = 0; i < num_cell; i++) {
			deltas[i] = randIndx(&post_probs_a[i*num_cluster], num_cluster);
			cell_assigns[deltas[i]].push_back(i);
		}
		
		// sample genotypes
		for(k = 0; k < num_cluster; k++) {
			vector<int>& tmp = cell_assigns[k];
			for(j = 0; j < num_muta; j++) {
				p1 = p0 = 1;
				if(tmp.size() > 0) {
					for(n = 0; n < tmp.size(); n++) {
						i = tmp[n];
						p1 *= clonecaller.getObsProb(1, obsData[i*num_muta+j], alpha, beta);
						p0 *= clonecaller.getObsProb(0, obsData[i*num_muta+j], alpha, beta);
					}
					double prob = p1/(p1+p0);
					states_p1[k*num_muta+j] = prob;
					p = randomDouble(0, 1);
					states[k*num_muta+j] = (p <= prob)? 1:0;
				}
			}
		}
		
		// estimate beta and alpha
		double tmp1 = 0, tmp2 = 0, tmp3 = 0, tmp4 = 0;
		for(i = 0; i < num_cell; i++) {
			k = deltas[i];
			for(j = 0; j < num_muta; j++) {
				if(obsData[i*num_muta+j] == 3) {
					continue;
				}
				int s = obsData[i*num_muta+j];
				tmp1 += (1-states[k*num_muta+j])*s;
				tmp2 += (1-states[k*num_muta+j]);
				tmp3 += states[k*num_muta+j]*(1-s);
				tmp4 += states[k*num_muta+j];
			}
		}
		alpha = tmp1/(tmp2+eps);
		beta = tmp3/(tmp4+eps);
		
		alpha = min(max(1e-5, alpha), max_alpha);
		beta = min(max(1e-5, beta), max_beta);
		
		
		// estimate pie
		pie.set(0.0);
		for(i = 0; i < num_cell; i++) {
			k = deltas[i];
			pie[k] += 1;
		}
		for(k = 0; k < num_cluster; k++) {
			pie[k] /= num_cell;
		}
		
		// check if the posterior distribution of z converges
		if(iter > 0) {
			double diff = 0;
			for(i = 0; i < num_cell; i++) {
				for(k = 0; k < num_cluster; k++) {
					diff += fabs(post_probs[i*num_cluster+k]-post_probs_pre[i*num_cluster+k]);
				}
			}
			diff /= (num_cell*num_cluster);
			if(diff < 1e-6) {
				break;
			}
		}
		post_probs_pre = post_probs;
		
		iter++;
	}
	
	// estimate genotypes
	for(k = 0; k < num_cluster; k++) {
		for(j = 0; j < num_muta; j++) {
			states[k*num_muta+j] = (states_p1[k*num_muta+j] > 0.5)? 1:0;
		}
	}
	
	// fine-tune the paras
	clonecaller.fineTune(s);
	clonecaller.refineSolution(s);
	clonecaller.fineTune(s);

	clonecaller.evalAccuracy(s);
	
	double lambda = min(0.15, 50.0/sqrt(num_muta*num_cell));
	s.bic = -2*s.ll+lambda*log(num_cell)*(s.num_cluster*(num_muta+1)+1);
	
	return NULL;
}

void CloneCaller::fineTune(Solution& s) {
	int i, j, k, n;
	int num_muta = obsData.getCOLS();
	int num_cell = obsData.getROWS();
	
	int num_cluster = s.num_cluster;
	double& alpha = s.alpha;
	double&	beta = s.beta;
	Matrix<int>& states = s.states;
	Matrix<double>& pie = s.pie;
	Matrix<double>& post_probs = s.post_probs;
	
	Matrix<long double> probs(num_cell, num_cluster);
	double eps = numeric_limits<double>::epsilon();
	double ll, pre_ll = numeric_limits<long>::min();
	while(1) {
		obsProbs(states, probs, alpha, beta);
		ll = 0;
		for(i = 0; i < num_cell; i++) {
			long double sum = 0;
			for(k = 0; k < num_cluster; k++) {
				sum += probs[i*num_cluster+k]*pie[k];
			}
			for(k = 0; k < num_cluster; k++) {
				post_probs[i*num_cluster+k] = probs[i*num_cluster+k]*pie[k]/sum;
			}
			ll += logl(sum);
		}
		s.ll = ll;
		double tmp1 = 0, tmp2 = 0, tmp3 = 0, tmp4 = 0;
		for(i = 0; i < num_cell; i++) {
			for(k = 0; k < num_cluster; k++) {
				double gamma = post_probs[i*num_cluster+k];
				for(j = 0; j < num_muta; j++) {
					if(obsData[i*num_muta+j] == 3) {
						continue;
					}
					int s = obsData[i*num_muta+j];
					tmp1 += gamma*(1-states[k*num_muta+j])*s;
					tmp2 += gamma*(1-states[k*num_muta+j]);
					tmp3 += gamma*states[k*num_muta+j]*(1-s);
					tmp4 += gamma*states[k*num_muta+j];
				}
			}
		}
		alpha = tmp1/(tmp2+eps);
		beta = tmp3/(tmp4+eps);
		for(k = 0; k < num_cluster; k++) {
			double sum = 0;
			for(i = 0; i < num_cell; i++) {
				sum += post_probs[i*num_cluster+k];
			}
			pie[k] = sum/num_cell;
		}
		if(fabs(pre_ll-ll) < 1e-4) {
			break;
		}
		else {
			pre_ll = ll;
		}
	}
}

void CloneCaller::printSolution(Solution& s) {
	cerr << "#clusters: " << s.num_cluster << endl;
	cerr << "LL: " << s.ll << endl;
	cerr << "BIC: " << s.bic << endl;
	
	if(s.acc >= 0) {
		cerr << "Acc: " << s.acc << endl;
	}
	cerr << "alpha: " << s.alpha << endl;
	cerr << "beta: " << s.beta << endl;
}

void CloneCaller::evalAccuracy(Solution& s) {
	if(realData.getROWS() == 0) {
		s.acc = -1;
		return;
	}
	int i, j, k;
	int num_cell = obsData.getROWS();
	int num_muta = obsData.getCOLS();
	
	Matrix<int>& states = s.states;
	Matrix<double>& post_probs = s.post_probs;
	Matrix<int> s_indxs;
	post_probs.max(2, s_indxs);
	Matrix<double> results(num_cell, num_muta, 0.0);
	int num_cluster = s.num_cluster;
	for(i = 0; i < num_cell; i++) {
		k = s_indxs[i];
		for(j = 0; j < num_muta; j++) {
			results[i*num_muta+j] = states[k*num_muta+j];
		}
	}
	
	double acc = 0;
	int n = 0;
	for(i = 0; i < num_cell; i++) {
		for(j = 0; j < doublet_indxs.size(); j++) {
			if(doublet_indxs[j] == i) {
				break;
			}
		}
		if(j < doublet_indxs.size()) {
			continue;
		}
		n++;
		for(j = 0; j < num_muta; j++) {
			k = (results[i*num_muta+j] > 0.5)? 1:0;
			if(k == realData[i*num_muta+j])	acc++;
		}
	}
	s.acc = acc/(n*num_muta);
}

void CloneCaller::refineSolution(Solution& s) {
	int i, j, k;
	int num_cell = obsData.getROWS();
	int num_muta = obsData.getCOLS();
	
	int& num_cluster = s.num_cluster;
	Matrix<int>& states = s.states;
	Matrix<double>& pie = s.pie;
	Matrix<double>& post_probs = s.post_probs;
	Matrix<int> s_indxs;
	post_probs.max(2, s_indxs);
	Matrix<int> flags(1, num_cluster, 0);
	vector<int> indxs;
	for(i = 0; i < num_cell; i++) {
		k = s_indxs[i];
		if(flags[k] == 0) {
			indxs.push_back(k);
		}
		flags[k] = 1;
	}
	
	num_cluster = indxs.size();
	states = states.Rows(indxs);
	post_probs = post_probs.Cols(indxs);
	post_probs.normalize(1);
	pie = pie.Cols(indxs);
	pie.normalize(1);
}

double CloneCaller::getObsProb(int gt, int observed, double alpha, double beta) {
	double prob = 1, eps = numeric_limits<double>::epsilon();
	if(gt == 1) {
		if(observed == 0) {
			prob = beta+eps;
		}
		else if(observed == 1) {
			prob = 1-beta+eps;
		}
	}
	else {
		if(observed == 0) {
			prob = 1-alpha+eps;
		}
		else if(observed == 1) {
			prob = alpha+eps;
		}
	}
	return prob;
}

void CloneCaller::obsProbs(Matrix<int>& states, Matrix<long double>& probs, double alpha, double beta) {
	int i, j, k;
	int num_muta = obsData.getCOLS();
	int num_cell = obsData.getROWS();
	int num_cluster = states.getROWS();
	long double eps = numeric_limits<long double>::epsilon();
	
	probs.resize(num_cell, num_cluster, false);
	
	for(i = 0; i < num_cell; i++) {
		for(k = 0; k < num_cluster; k++) {
			long double prob = 1;
			for(j = 0; j < num_muta; j++) {
				if(states[k*num_muta+j] == 1) {
					if(obsData[i*num_muta+j] == 0) {
						prob *= beta;
					}
					else if(obsData[i*num_muta+j] == 1) {
						prob *= 1-beta;
					}
				}
				else {
					if(obsData[i*num_muta+j] == 0) {
						prob *= 1-alpha;
					}
					else if(obsData[i*num_muta+j] == 1) {
						prob *= alpha;
					}
				}
			}
			probs[i*num_cluster+k] = prob;
		}
	}
}

void CloneCaller::saveResults() {
	int i, j, k;
	ofstream ofs;
	string fn, outputPrefix = config.getStringPara("output");
	
	int num_cell = obsData.getROWS();
	int num_muta = obsData.getCOLS();
	
	int num_cluster = best_s.num_cluster;
	Matrix<int>& states = best_s.states;
	fn = outputPrefix+".clones.txt";
	ofs.open(fn.c_str());
	if(!ofs.is_open()) {
		cerr << "Error: cannot open file " << fn << endl;
		exit(-1);
	}
	for(i = 0; i < num_cluster; i++) {
		for(j = 0; j < num_muta; j++) {
			if(j < num_muta-1) {
				ofs << states[i*num_muta+j] << '\t';
			}
			else {
				ofs << states[i*num_muta+j] << endl;
			}
		}
	}
	ofs.close();
	
	Matrix<double>& post_probs = best_s.post_probs;
	Matrix<int> s_indxs;
	post_probs.max(2, s_indxs);
	fn = outputPrefix+".assignment.txt";
	ofs.open(fn.c_str());
	if(!ofs.is_open()) {
		cerr << "Error: cannot open file " << fn << endl;
		exit(-1);
	}
	for(i = 0; i < num_cell-1; i++) {
		ofs << s_indxs[i] << " ";
	}
	ofs << s_indxs[num_cell-1] << endl;
	ofs.close();
}


