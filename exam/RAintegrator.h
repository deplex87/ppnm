#include <math.h>
#include <assert.h>
#include <complex.h>

double complex adapt24(double complex f(double complex), double complex a, double complex b, double delta, double eps, double complex f2, double complex f3, int numOfEvals);

double complex integrate(double complex f(double complex), double complex a, double complex b, double delta, double eps);

double complex OpenQuadCC_integrate(double complex f(double complex), double complex a, double complex b, double delta, double eps);
