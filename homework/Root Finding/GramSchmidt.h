#include <assert.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <math.h>

void GS_decomp( gsl_matrix* A, gsl_matrix* R ); 
void GS_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x);
void GS_inverse(gsl_matrix* Q, gsl_matrix* R, gsl_matrix* B);
