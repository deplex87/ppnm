#include <math.h>
#include <stdio.h>
#include <complex.h>
#include "RAintegrator.h"

int main(){
    
    printf("Exam Project\n\n");
    printf("We test our adaptive integrator on different integrals:\n");
    
    int numOfEvals = 0;

    double complex func1(double complex z){
        numOfEvals += 1;
        return 1/(2*csqrt(z));
    }
    
    double complex func2(double complex z){
        numOfEvals += 1;
        return clog(z)/csqrt(z);
    }
    
    double complex b = 3.0 + 5*I;
    double complex IntegralResult = integrate(func1, (double complex) 0, b, 0.001, 0.001);
    
    printf("Using the recursive adaptive routine on the complex function 1/(2*sqrt(z)) from 0 to 3.0 + 5i gives the following results:\n");
    printf("Numerical result = %g + %gi\n", creal(IntegralResult), cimag(IntegralResult));
    printf("Exact result = sqrt(3 + 5i)\n");
    printf("Exact w. decimals = 2.101303 + 1.189737i\n");
    printf("numOfEvals = %i\n\n", numOfEvals);
    
    numOfEvals = 0;
    
    b = 3.0 + 5*I;
    IntegralResult = OpenQuadCC_integrate(func1, (double complex) 0, b, 0.001, 0.001, NULL);
    
    printf("Using open quadrature with Clenshaw–Curtis variable transformation on the complex function 1/(2*sqrt(z)) from 0 to 3.0 + 5i gives the following results:\n");
    printf("Numerical result = %g + %gi\n", creal(IntegralResult), cimag(IntegralResult));
    printf("Exact result = sqrt(3 + 5i)\n");
    printf("Exact w. decimals = 2.101303 + 1.189737i\n");
    printf("numOfEvals = %i\n", numOfEvals);
    printf("Here we see that it takes a lot less evaluations using this instead of recursive adaptive integrator.\n\n");
    
    numOfEvals = 0;
    
    b = 2.0 + 1.5*I;
    FILE* line = fopen("line.txt", "w");
    IntegralResult = OpenQuadCC_integrate(func2, (double complex) -1.5 - 1.8*I, b, 0.001, 0.001, line);
    
    printf("Using open quadrature with Clenshaw–Curtis variable transformation on the complex function log(z)/sqrt(z) from 0 to 4.0 + 8i gives the following results:\n");
    printf("Numerical result = %g + %gi\n", creal(IntegralResult), cimag(IntegralResult));
    printf("Exact result = 4 sqrt(1 + 2 i) (-2 + log(4 + 8 i))\n");
    printf("Exact w. decimals = -2.509655 + 6.233921i\n");
    printf("numOfEvals = %i\n\n", numOfEvals);
    
    FILE* plane = fopen("plane.txt", "w");
	for(int i = -60; i < 60; i++){
		for(int j = -60; j < 60; j++){
			fprintf(plane,"%g\t%g\t%g\n", (double) i/20, (double) j/20, pow(pow(creal((func2(CMPLX((double) i/20, (double) j/20)))),2)+pow(cimag((func2(CMPLX((double) i/20, (double) j/20)))),2),0.5)); 
			}
		}
	fclose(plane);
    
    return 0;
}
