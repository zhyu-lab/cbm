// ***************************************************************************
// MyDefine.cpp (c) 2021 Zhenhua Yu <qasim0208@163.com>
// Health Informatics Lab, Ningxia University
// All rights reserved.

#include <cstdio>
#include <fstream>
#include <vector>
#include <cmath>
#include <limits>

#include "MyDefine.h"

using namespace std;

MyDefine::MyDefine() {
}

MyDefine::~MyDefine() {
}

string current_version = "1.0";

/*** definition of global vars ***/
Config config;
InputParser parser;
CloneCaller clonecaller;
ThreadPool *tp;
/*** end of definition ***/

/*** definition of general functions ***/

double normrand(double mu, double sigma) {
	double eps = numeric_limits<double>::epsilon();
    double pi = 3.14159265358979323846;

    static double z0, z1;
    static bool flag = true;

    if(!flag) {
       return z1*sigma+mu;
	}
	flag = !flag;
	
    double u1, u2;
    do {
       u1 = rand()*(1.0/RAND_MAX);
       u2 = rand()*(1.0/RAND_MAX);
	}while(u1 <= eps);

    z0 = sqrt(-2.0*log(u1))*cos(2*pi*u2);
    z1 = sqrt(-2.0*log(u1))*sin(2*pi*u2);
	
    return z0*sigma+mu;
}

double normrand_n(double mu, double sigma) {
	static default_random_engine generator;
	static normal_distribution<double> distribution(0, 1);
    double z = distribution(generator);
	
    return z*sigma+mu;
}



//****** produce a matrix containing random index according to the probability distribution ******//
void randIndx(int m, int n, Matrix<int>& ret, Matrix<double>& prob) {
	if(m <= 0 || n <= 0) {
		return;
	}
	//assert(prob.getROWS() == 1 && prob.getCOLS() >= 1);
	
	ret.resize(m, n, false);
	int i, j;
	for(i = 0; i < m; i++) {
		for(j = 0; j < n; j++) {
			ret[i*n+j] = randIndx(prob);
		}
	}
}

int randIndx(Matrix<double>& prob) {
	//assert(prob.getROWS() == 1 && prob.getCOLS() >= 1);
	int ac = prob.getCOLS();
	
	double r = randomDouble(0, 1);
	for(int k = 0; k < ac; k++) {
		if(r <= prob[k]) {
			return k;
		}
	}
	return ac-1;
}

int randIndx(double* prob, int n) {
	double r = randomDouble(0, 1);
	for(int k = 0; k < n; k++) {
		if(r <= prob[k]) {
			return k;
		}
	}
	return n-1;
}

//****** produce random double ******//
double randomDouble(double start, double end) {
	//return start+(end-start)*(rand()/(RAND_MAX+1.0));
	return tp->randomDouble(start, end);
}

//****** produce random int ******//
long randomInteger(long start, long end) {
	//return start+(end-start)*(rand()/(RAND_MAX+1.0));
	return tp->randomInteger(start, end);
}

