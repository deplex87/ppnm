#include<math.h>
#include<assert.h>

double adapt24(double f(double), double a, double b, double delta, double eps, double f2, double f3, int numOfEvals, double* err);

double integrate(double f(double), double a, double b, double delta, double eps, double* err);

double OpenQuadCC_integrate(double f(double), double a, double b, double delta, double eps, double* err);
