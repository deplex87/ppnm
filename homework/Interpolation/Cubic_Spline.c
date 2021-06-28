#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#include<math.h>
#include"Binary_Search.h"
 
	// Cubic spline //

typedef struct {int n; double *x,*y,*b,*c,*d;} cspline;					// Der skal igen dannes en ny strucktur med et parameter mere, men ellers gør vi meget det samme som ved qspline

cspline* cspline_alloc(int n, double *x, double *y){
	cspline* s = (cspline*)malloc(sizeof(cspline));
	s->x= (double*)malloc(n*sizeof(double));
	s->y= (double*)malloc(n*sizeof(double));
	s->b= (double*)malloc(n*sizeof(double));
	s->c= (double*)malloc((n-1)*sizeof(double));
	s->d= (double*)malloc((n-1)*sizeof(double));
	s->n= n;

	int i;
	for(i=0;i<n;i++){s->x[i]=x[i];s->y[i]=y[i];}

	double a[n-1] , dx[n-1];
	for (i = 0; i<n-1; i++){
		dx[i]=x[i+1]-x[i];
		a[i]=(y[i+1]-y[i])/dx[i];}
	// Her laves systemmet beskrevet ved sætning 21 fra kapitel 1
	double D[n], Q[n-1], B[n];
	D[0]=2; D[n-1]=2; Q[0]=1;
	for(i=0;i<n-2;i++)D[i+1]=2*dx[i]/dx[i+1]+2;
	for(i=0;i<n-2;i++)Q[i+1]=dx[i]/dx[i+1];
	for(i=0;i<n-2;i++)B[i+1]=3*(a[i]+a[i+1]*dx[i]/dx[i+1]);
	//solving the system:
	B[0]=3*a[0]; B[n-1]=3*a[n-2];
	for(i=1;i<n;i++){
		D[i]-=Q[i-1]/D[i-1];
		B[i]-=B[i-1]/D[i-1];}
	s->b[n-1]=B[n-1]/D[n-1];
	for(i=n-2;i>=0;i--)s->b[i]=(B[i]-Q[i]*s->b[i+1])/D[i];
	for(i=0;i<n-1;i++){
		s->c[i]=(-2*s->b[i]-s->b[i+1]+3*a[i])/dx[i];
		s->d[i]=(s->b[i]+s->b[i+1]-2*a[i])/dx[i]/dx[i];}
return s;
}

double cspline_eval(cspline *s,double z){
	int i=binsearch(s->n,s->x,z);
	double	dx=z-s->x[i];
	double f_z=s->y[i]+dx*(s->b[i]+dx*(s->c[i]+dx*s->d[i]));
return f_z;
}

double cub_integral(cspline *s, double z){
	int I=binsearch(s->n,s->x,z);
	double S=0;

	for(int i=0;i<I;i++){
	double dx=s->x[i+1]-s->x[i]; assert(dx>0);
		S+=dx*(s->y[i]+dx*(s->b[i]/2+dx*(s->c[i]/3+dx*s->d[i]/4)));
	}

	double dx=z-s->x[I]; assert(dx>0);
	double F_z=dx*(s->y[I]+dx*(s->b[I]/2+dx*(s->c[I]/3+dx*s->d[I]/4)));
	S+=F_z;
return S;
}

double cub_diff(cspline *s, double z){
	int i=binsearch(s->n,s->x,z);
	double dx=z-s->x[i];
	double df_dz=s->b[i]+dx*(2*s->c[i]+3*dx*s->d[i]);
return df_dz;
}
	
	
void cspline_free(cspline *s){
free(s->x);free(s->y);free(s->b);free(s->c);free(s->d);free(s);
}
