#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#include<math.h>

typedef struct {int n ; double *x , *y , *b , *c;} qspline; 	// Der definieres en struct (qspline) hvor alt vores data skal opsamles i

qspline* qspline_alloc ( int n , double* x , double* y );


double qspline_eval(qspline *s, double z);

double qua_integral(qspline *s, double z);

double qua_diff(qspline *s, double z);

void qspline_free(qspline *s);
