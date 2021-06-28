#include <float.h>
#include <assert.h>
#include <math.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include "GramSchmidt.h"

void newtonMethod(void f(gsl_vector* x, gsl_vector* f_values), gsl_vector* start, double eps);
