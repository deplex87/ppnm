#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#include<math.h>


int binsearch(int n, double* x, double z);

typedef struct {int n; double *x,*y,*b,*c,*d;} cspline;					// Der skal igen dannes en ny strucktur med et parameter mere, men ellers gør vi meget det samme som ved qspline

cspline* cspline_alloc(int n, double *x, double *y);
double cspline_eval(cspline *s,double z);

double cub_integral(cspline *s, double z);

double cub_diff(cspline *s, double z);

void cspline_free(cspline *s);
