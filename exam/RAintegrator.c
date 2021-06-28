#include <math.h>
#include <assert.h>
#include <complex.h>

double complex adapt24(double complex f(double complex), double complex a, double complex b, double delta, double eps, double complex f2, double complex f3, int numOfEvals){
    
    assert(numOfEvals < 100000);
    
    double complex f1 = f(a + (b - a) / 6);
    double complex f4 = f(a + 5 * (b - a) / 6);
    
    double complex Q = (2 * f1 + f2 + f3 + 2 * f4) / 6 * (b - a);
    double complex q = (f1 + f4 + f2 + f3) / 4 * (b - a);
    
    double complex_norm(double complex z){
        return sqrt(creal(z)*creal(z) + cimag(z)*cimag(z));
    }
    
    double tol = delta + eps*complex_norm(Q);
    double error = complex_norm(Q - q);
    
    if(error < tol){
        return Q;
    }else{
        double complex new_Q1 = adapt24(f, a, (a + b) / 2, delta/sqrt(2.), eps, f1, f2, numOfEvals + 1);
        double complex new_Q2 = adapt24(f, (a + b) / 2, b, delta/sqrt(2.), eps, f3, f4, numOfEvals + 1);
        return new_Q1 + new_Q2;
    }
}

double complex integrate(double complex f(double complex), double complex a, double complex b, double delta, double eps){
    int numOfEvals = 0;
    double complex f2 = f(a + 2 * (b - a) / 6);
    double complex f3 = f(a + 4 * (b - a) / 6);
    return adapt24 (f, a, b, delta, eps, f2, f3, numOfEvals);
}

double complex OpenQuadCC_integrate(double complex f(double complex), double complex a, double complex b, double delta, double eps){
    int numOfEvals = 0;
    double complex f2 = f(a + 2 * (b - a) / 6);
    double complex f3 = f(a + 4 * (b - a) / 6);
    
    double complex IntervalTransformation (double complex theta){
        return f((a + b) / 2 + (b - a) * cos(theta) / 2 )*(b - a) * sin(theta) / 2;
    }

    return adapt24(IntervalTransformation, 0, M_PI, delta, eps, f2, f3, numOfEvals);
}
