#include <math.h>
#include <stdio.h>
#include <complex.h>
#include "RAintegrator.h"

int main(){
    
    printf("Task A\n\n");
    printf("We test our adaptive integrator on different integral:\n");
    
    int numOfEvals = 0;

    double complex func1(double complex z){
        numOfEvals += 1;
        return 1/(csqrt(z));
    }
    
    double complex b = 1.0 + 2*I;
    double complex IntegralResult = integrate(func1, (double complex) 0, b, 0.001, 0.001);
    
    printf("Our numerical integration of the complex function 1/(2*sqrt(z)) from 0 to 1.0 + 2i gives the following results:\n");
    printf("Numerical result = %g + %gi\n", creal(IntegralResult), cimag(IntegralResult));
    printf("numOfEvals = %i\n\n", numOfEvals);
    
    
    return 0;
}
