#include <stdio.h>
#include <math.h>
#include <complex.h>
#include "montecarlo.h"

double f(int dim, double* x){
	return 1/sqrt(x[0] + x[1]);
}

double g(int dim,double* x){
    return x[1]/exp(x[0]*x[0]);
}

double h(int dim, double * x){
	return 1/((1 - cos(x[0]) * cos(x[1]) * cos(x[2])) * (M_PI*M_PI*M_PI));
}


int main(){
    
    printf("Task A\n\n");

    printf("We start by trying to calculate the integral of 1/sqrt(x + y) from 0 to 1 in both x and y:\n");
	double a_f[] = {0, 0};
	double b_f[] = {1, 1};
	complex Integral_f = PlainMonteCarlo(2, f, a_f, b_f, 1e7);
	printf("Numerical result = %g +/- %g\n", creal(Integral_f),cimag(Integral_f));
    printf("Analytical result = %g\n\n", (double) (8.0/3.0)*(sqrt(2.0) - 1.0));

    printf("Next we integrate y/e^(x^2) from 0 to 1 in both x and y:\n");
	complex Integral_g = PlainMonteCarlo(2, g, a_f, b_f, 1e7);
	printf("Numerical result = %g +/- %g\n", creal(Integral_g),cimag(Integral_g));
    printf("Analytical result = %g\n\n", (double) (1.0/4.0)*sqrt(M_PI)*erf(1.0));

    printf("Lastly we calculate the integral given to us in the description of the task:\n");
	double a_h[] = {0,0,0};
	double b_h[] = {M_PI, M_PI, M_PI};
	complex Integral_h = PlainMonteCarlo(3, h, a_h, b_h, 1e7);
	printf ("Numerical result = %g +/- %g \n", creal(Integral_h), cimag(Integral_h));
	printf ("Analytical result = %g \n\n\n", pow(tgamma(1./4), 4)/(4*M_PI*M_PI*M_PI));


    printf("Task B\n\n");
    
    printf("Now we test the integrator with quasi-random sequences. \nFirst we do the calculation for the integral 1/sqrt(x + y) from 0 to 1 in both x and y:\n");
	complex result_quasi_f = plain_quasi_MonteCarlo(2, f, a_f, b_f, 1e7);
	printf("Numerical result = %g +/- %g\n", creal(result_quasi_f), cimag(result_quasi_f));
    printf("Analytical result = %g\n\n", (double) (8.0/3.0)*(sqrt(2.0) - 1.0));  
    
    printf("Now we do the same thing for the integral y/e^(x^2) from 0 to 1 in both x and y:\n");
	complex result_quasi_g = plain_quasi_MonteCarlo(2, g, a_f, b_f, 1e7);
	printf("Numerical result = %g +/- %g\n",creal(result_quasi_g), cimag(result_quasi_g));
    printf("Analytical result = %g\n\n", (double) (1.0/4.0)*sqrt(M_PI)*erf(1.0));    
    
    printf("Lastly we do it again for the integral given in task A:\n");
	complex result_quasi_h = plain_quasi_MonteCarlo(3, h, a_h, b_h, 1e7);
	printf("Numerical result = %g +/- %g\n",creal(result_quasi_h),cimag(result_quasi_h));
	printf ("Analytical result = %g \n\n\n", pow(tgamma(1./4), 4)/(4*M_PI*M_PI*M_PI));    


	FILE* error_data = fopen("error_data.txt","w");
	for(int N = 1000; N < 100000; N += 1000){
        complex result_plain = PlainMonteCarlo(2, g, a_f, b_f, N);
		complex result_quasi = plain_quasi_MonteCarlo(2, g, a_f, b_f, N);
		fprintf(error_data, "%20d %20g %20g\n", N, cimag(result_plain), cimag(result_quasi));
	}

	fclose(error_data);

	return 0;
} 
