#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <math.h>

void rkstep12(
	void f(double t, gsl_vector* y,gsl_vector* dydt), /* the f from dy/dt=f(t,y) */
	double t,              /* the current value of the variable */
	gsl_vector* y_val,            /* the current value y(t) of the sought function */
	double step,              /* the step to be taken */
	gsl_vector* y_step_val,             /* output: y(t+h) */
	gsl_vector* err             /* output: error estimate */
);

void driver(
	void f(double,gsl_vector*,gsl_vector*), /* right-hand-side of dy/dt=f(t,y) */
	double start_point,                     /* the start-point a */
	gsl_vector* ya,                     /* y(a) */
	double end_point,                     /* the end-point of the integration */
	gsl_vector* yb,                     /* y(b) to be calculated */
	double step,                     /* initial step-size */
	double acc,                   /* absolute accuracy goal */
	double eps,                    /* relative accuracy goal */
    FILE* filepath
);
