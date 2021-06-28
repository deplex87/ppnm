#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<math.h>
#include<stdio.h>

#include "qnewton.h"

double *Energy;
double *CrossSection;
double *Error;

double BreitWigner(double E,gsl_vector*v){
	return gsl_vector_get(v,2)/(pow((E-gsl_vector_get(v,0)),2)+pow(gsl_vector_get(v,1)/2,2));
}

double deviation(gsl_vector* v) {
    double Energy[] =  {101, 103, 105, 107, 109, 111, 113, 115, 117, 119, 121, 123, 125, 127, 129, 131, 133, 135, 137, 139, 141, 143, 145, 147, 149, 151, 153, 155, 157, 159};
	double CrossSection[] = {-0.25,-0.30,-0.15,-1.71,0.81,0.65,-0.91,0.91,0.96,-2.52,-1.01,2.01,4.83,4.58,1.26,1.01,-1.26,0.45,0.15,-0.91,-0.81,-1.41,1.36,0.50,-0.45,1.61,-2.21,-1.86,1.76,-0.50};
	double Error[] = {2, 2, 1.9, 1.9, 1.9, 1.9, 1.9, 1.9, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 0.9, 0.9, 0.9};
    double sum = 0;
    for (int i = 0; i < 30; i++){
	sum += (BreitWigner(Energy[i],v) - CrossSection[i]) * (BreitWigner(Energy[i],v) - CrossSection[i]) / (Error[i] * Error[i]);
	}
	return sum;
}

double rosenbrockValley(gsl_vector* values) {
    double x = gsl_vector_get(values, 0);
    double y = gsl_vector_get(values, 1);

    double f_value = (1-x)*(1-x) + 100*(y-x*x)*(y-x*x);
    return f_value;
}

double himmelblau(gsl_vector* values) {
    double x = gsl_vector_get(values, 0);
    double y = gsl_vector_get(values, 1);

    double f_value = (x*x+y-11)*(x*x + y - 11) + (x + y*y - 7)*(x + y*y - 7);
    return f_value;
}

int main(){
	printf( "Task A\n\n");

    printf("We begin by testing the routine by finding the minimum of the Rosenbrock's valley function f(x,y) = (1 - x)^2 + 100*(y - x^2)^2.\n");
    printf("We guess that the minimum is found at x = 0 and y = 0.\n");
    
    gsl_vector* minimum = gsl_vector_alloc(2);

    gsl_vector_set(minimum, 0, 0.0);
    gsl_vector_set(minimum, 1, 0.0);

    int nsteps = qnewton(rosenbrockValley, minimum, 1e-5);
    printf("The minimum is found to be at x = %g and y = %g\n", gsl_vector_get(minimum, 0), gsl_vector_get(minimum, 1));
    printf("The real minimum is at x = 1 and y = 1.\n");
    printf("This calculation required in total %i steps.\n\n", nsteps);

    printf("Now we continue and try to minimize the Himmelblaus function, f(x,y) = (x^2 + y - 11)^2 + (x + y^2 - 7)^2.\n");
    printf("We guess that the minimum is found at x = 2.5 and y = 1.5.\n");

    gsl_vector_set(minimum, 0, 2.5);
    gsl_vector_set(minimum, 1, 1.5);

    nsteps = qnewton(himmelblau, minimum, 1e-5);
    printf("The minimum is found to be at x = %g and y = %g\n", gsl_vector_get(minimum, 0), gsl_vector_get(minimum, 1));
    printf("The actual minima of the function are (x, y) = (3, 2), (-3.78, -3.28), (-2.80, 3.13) and (3.58, -1.84)\n");
    printf("This calculation required in total %i steps.\n\n", nsteps);


    printf("\nTask B\n\n");
    gsl_vector* guess     = gsl_vector_alloc(30);
    gsl_vector_set(guess, 0, 125);
    gsl_vector_set(guess, 1, 3);
    gsl_vector_set(guess, 2, 8);
    double eps = 1e-6;
    qnewton(deviation,guess,eps);
    printf("Using the minimization routine to fit the data we find the following results:\n");
    printf("Mass = %f GeV\nWidth = %f\nA = %f \n",fabs(gsl_vector_get(guess,0)),fabs(gsl_vector_get(guess,1)), gsl_vector_get(guess,2));
    
    FILE* plot = fopen("higgs_fit.txt","w");
    for(double e=90;e<170;e++)fprintf(plot,"%6g %6g\n",e,BreitWigner(e,guess));
    fclose(plot);
    
return 0;
}
