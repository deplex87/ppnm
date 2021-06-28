#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#include<math.h>
#include"Binary_Search.h"

		// linear spline //

double linterp(int n, double x[], double y[], double z){
	int i=binsearch(n,x,z);

	double dy=y[i+1]-y[i];
	double dx=x[i+1]-x[i]; assert(dx>0);
	double f_z=y[i]+dy/dx*(z-x[i]);
return f_z;
}

double linterp_integral(int n, double x[], double y[], double z){
	int I=binsearch(n,x,z);
	double s=0;

	for(int i=0;i<I;i++){
	double dy=y[i+1]-y[i];
	double dx=x[i+1]-x[i]; assert(dx>0);
	double ai=dy/dx;
		s+=y[i]*(x[i+1]-x[i])+ai*pow((x[i+1]-x[i]),2)/2;
	}

	double dy=y[I+1]-y[I];
	double dx=x[I+1]-x[I]; assert(dx>0); //
	double az=dy/dx;
	double F_z=y[I]*(z-x[I])+az*pow((z-x[I]),2)/2;
	s+=F_z;
return s;
} 
