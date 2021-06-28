#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include "GramSchmidt.h"


void least_square(int numOfPoints, int numOfFuncs, double* x, double* y, double* unc_y, double f(int,double), gsl_vector* coefficients, gsl_matrix* covMatrix){
	
	gsl_matrix* A = gsl_matrix_alloc(numOfPoints, numOfFuncs);
	gsl_matrix* R = gsl_matrix_alloc(numOfFuncs, numOfFuncs);
    
	for(int i = 0; i < numOfPoints; i++){
		for(int j = 0; j < numOfFuncs; j++)
		{
            gsl_matrix_set(A, i, j, (double) f(j, x[i])/unc_y[i]);
		}
    }

	gsl_vector* b = gsl_vector_alloc(numOfPoints);
    
	for(int i = 0; i < b->size; i++){
		gsl_vector_set(b,i, (double) y[i]/unc_y[i]);
    }
	GS_decomp(A, R);
	GS_solve(A, R, b, coefficients);

	gsl_matrix* R_inverse = gsl_matrix_alloc(numOfFuncs, numOfFuncs);
    gsl_matrix* I = gsl_matrix_alloc(numOfFuncs, numOfFuncs);
    gsl_matrix_set_identity(I);
	GS_inverse(I, R, R_inverse);
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1, R_inverse, R_inverse, 0, covMatrix);

	gsl_matrix_free(A);
    gsl_matrix_free(R);
    gsl_vector_free(b);
    gsl_matrix_free(I);
    gsl_matrix_free(R_inverse);
} 
