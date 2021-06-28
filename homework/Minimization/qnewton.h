#include<stdio.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<float.h>
#ifndef DBL_EPSILON
#define DBL_EPSILON 2.22045e-16
#endif

void numeric_gradient(double F(gsl_vector*), gsl_vector*x, gsl_vector*grad);

int qnewton(double F(gsl_vector*), gsl_vector*x, double acc);
