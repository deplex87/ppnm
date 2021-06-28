#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#include<math.h>
#include"Binary_Search.h"

		// Qudratic spline //

typedef struct {int n ; double *x , *y , *b , *c;} qspline; 	// Der definieres en struct (qspline) hvor alt vores data skal opsamles i

qspline* qspline_alloc ( int n , double* x , double* y ){	//qspline bygges her
	qspline  *s = (qspline*) malloc (sizeof(qspline));		// Først laves der plads til alt data i hukommelsen
	s->b = (double*) malloc((n-1)*sizeof(double));
	s->c = (double*) malloc((n-1)*sizeof(double));
	s->x = (double*) malloc(n*sizeof(double)) ;
	s->y = (double*) malloc(n*sizeof(double)) ;
	s->n = n;
	int i; // Vi gidder ikke skrive int i alle vores forlykker
	for (i = 0; i<n; i++){s->x[i]=x[i]; s->y[i]=y[i];}		//Der henvises til vores data
	double a[n-1] , dx[n-1];
	for (i = 0; i<n-1; i++){
		dx[i]=x[i+1]-x[i];
		a[i]=(y[i+1]-y[i])/dx[i];}
	//Nu er alt sat op til at udregne en quadratic spline, som gøres på to måder, da starten ikke kendes
	s->c[0]=0;
	for( i=0; i<n-2; i++){
		s->c[i+1]=(a[i+1]-a[i]-s->c[i]*dx[i])/dx[i+1];}		//Sætning 13 fra kapitel 1
	s->c[n-2]/=2;
	for (i=n-3; i >=0;i--){
		s->c[i]=(a[i+1]-a[i]-s->c[i+1]*dx[i+1])/dx[i];}		//sætning 14 fra kapitel 1
	for (i =0; i<n-1; i++){
		s->b[i]=a[i]-s->c[i]*dx[i];}				//Sætning 15 fra kapitel 1
return s;
}


double qspline_eval(qspline *s, double z){			//Functionen som via qspline udregner s(z)
	int i=binsearch(s->n,s->x,z);
	double dx=z-s->x[i];
	double f_z=s->y[i]+dx*(s->b[i]+dx*s->c[i]);
return f_z;
}

double qua_integral(qspline *s, double z){
	int I=binsearch(s->n,s->x,z);
	double S=0;

	for(int i=0;i<I;i++){
	double dx=s->x[i+1]-s->x[i]; assert(dx>0);
		S+=dx*(s->y[i]+dx*(s->b[i]/2+dx*s->c[i]/3));
	}

	double dx=z-s->x[I]; assert(dx>0);
	double F_z=dx*(s->y[I]+dx*(s->b[I]/2+dx*s->c[I]/3));
	S+=F_z;
return S;
}

double qua_diff(qspline *s, double z){
	int i=binsearch(s->n,s->x,z);
	double dx=z-s->x[i];
	double df_dz=s->b[i]+2*dx*s->c[i];
return df_dz;
}

void qspline_free(qspline *s){
free(s->x);free(s->y);free(s->b);free(s->c);free(s);
}
