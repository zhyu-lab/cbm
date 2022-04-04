// ***************************************************************************
// MyDefine.h (c) 2021 Zhenhua Yu <qasim0208@163.com>
// Health Informatics Lab, Ningxia University
// All rights reserved.

#ifndef _MYDEFINE_H
#define _MYDEFINE_H

#include <iostream>
#include <cstring>
#include <cstdlib>
#include <random>

#include "Config.h"
#include "InputParser.h"
#include "CloneCaller.h"
#include "Matrix.h"
#include "ThreadPool.h"

class MyDefine {
	public:
		MyDefine();
		~MyDefine();
};

/*** declaration of global vars ***/
extern string current_version;
extern Config config;
extern InputParser parser;
extern CloneCaller clonecaller;
extern ThreadPool *tp;
/*** end of declaration ***/

/*** declaration of general functions ***/

double normrand(double mu, double sigma);
double normrand_n(double mu, double sigma);
void randIndx(int m, int n, Matrix<int>& ret, Matrix<double>& prob);
int randIndx(Matrix<double>& prob);
int randIndx(double* prob, int n);
double randomDouble(double start, double end);
long randomInteger(long start, long end);

/*** end of declaration ***/


#endif

