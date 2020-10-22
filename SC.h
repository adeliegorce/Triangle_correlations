#ifndef _SC_H_
#define _SC_H_
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

static const int nthreads = 20; //number of threads for parallelisation
static const int N = 512; //sampling number
static const double L = 400; //length of the realspace box in Mpc
static const string filename_box = "Field_20binary_bubbles_nooverlap=False_radius=20_xhII0.092_N512_2D"; //"Field_625binary_bubbles_nooverlap=False_radius=8_xhII0.075_N256_3D";

static const int nbins = 100;//number of bins for correlation scales
static const double rmin = 0.5;
static const double rmax = 30.; //min and max value of correlation scales probed in Mpc 

// #define pi (double) 3.141592653589793
// #define s3 (double) sqrt(3.0)

// #define nthreads (int) 20 //number of threads for parallelisation
// #define N (int) 512 //sampling number
// #define L (double) 400 //length of the realspace box
// #define NDIM (int) 2 //dimension of the box
// #define filename_box (str) "Field_20binary_bubbles_nooverlap=False_radius=20_xhII0.100_N512_2D"

// #define nbins (int) 100 //number of bins for correlation scales
// #define rmin (double) 1.
// #define rmax (double) 30.  //min and max value of correlation scales probed 

#endif /* SC_H */