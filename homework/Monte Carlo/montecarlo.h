#include <math.h>
#include <complex.h>
#include <assert.h>
#include <stdlib.h>
#define RND (double)rand()/RAND_MAX

complex PlainMonteCarlo(int dim, double f(int dim, double* x), double* a, double* b, int N);

void halton(int n, int dim, double *a, double *b, double *x);

void halton2(int n, int dim, double *a, double *b, double *x);

double corput(int n, int base);

complex plain_quasi_MonteCarlo(int dim, double f(int dim, double* x), double* a, double* b, int N);
