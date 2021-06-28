#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include "GramSchmidt.h"

void least_square(int numOfPoints, int numOfFuncs, double* x, double* y, double* unc_y, double f(int,double), gsl_vector* coefficients, gsl_matrix* covMatrix); 
