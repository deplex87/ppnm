#ifndef EIGENVALUES_JACOBI_H
#define EIGENVALUES_JACOBI_H

#include <math.h>
#include <gsl/gsl_matrix.h>

void timesJ(gsl_matrix* matrix, int i, int j, double angle);
void Jtimes(gsl_matrix* matrix, int i, int j, double angle);
void jacobi_diag(gsl_matrix* OriginalMatrix, gsl_matrix* EigenVectorMatrix);

#endif
