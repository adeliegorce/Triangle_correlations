#ifndef _CC_H_
#define _CC_H_
#pragma once

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex>
#include <string>
#include <omp.h>
#include <fftw3.h>
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <chrono>
using namespace std;

static const double pi = 3.141592653589793;
static const double s3 = sqrt(3.0);

static const int nthreads = 5; //number of threads for parallelisation
static const int N = 50; //sampling number
static const double L = 40; //length of the realspace box in Mpc
static const string filename_box1 = "Field_10binary_bubbles_nooverlap=False_radius=2_xhII0.05_N50_2D";
static const string filename_box2 = "Field_10binary_bubbles_nooverlap=False_radius=2_xhII0.05_N50_2D";

static const int nbins = 10;//number of bins for correlation scales
static const double rmin = 0.5;
static const double rmax = 5.; //min and max value of correlation scales probed in Mpc 

#endif /* CC_H */