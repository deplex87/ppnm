#include <math.h>
#include <complex.h>
#include <stdlib.h>
#include <assert.h>
#define RND (double)rand()/RAND_MAX

complex PlainMonteCarlo(int dim, double f(int dim, double* x), double* a, double* b, int N){
    double volume = 1;
    
    for(int i = 0; i < dim; i++){
        volume *= b[i] - a[i];
    }
    
    double sum = 0;
    double sumSquared = 0;
    double x[dim];
    double fval;
    for(int j = 0; j < N; j++){
        for(int i = 0; i < dim; i++){
            x[i] = a[i] + RND * (b[i] - a[i]);
        }
        fval = f(dim, x);
        sum += fval;
        sumSquared += pow(fval, 2);
    }
    
    double mean = sum/N;
    double sigma = sqrt(sumSquared/N - pow(mean, 2));
    complex result = mean*volume + I*sigma*volume/sqrt(N);
    return result;
}

double corput(int n, int base){
    double CorputNumber = 0;
    double bk = (double) 1 / base;
    while(n > 0){
        CorputNumber += (n % base)*bk;
        n /= base;
        bk /= base;
    }
    return CorputNumber;
}

void halton(int n, int dim, double *a, double *b, double *x){
    int base[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73};
    int maxd = sizeof(base) / sizeof(int);  assert(dim <= maxd);
    for(int i = 0; i < dim; i++){
        x[i] = a[i] + corput(n + 1, base[i]) * (b[i] - a[i]);
    }
}

void halton2(int n, int dim, double *a, double *b, double *x){
    int base[] = {5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67};
    int maxd = sizeof(base) / sizeof(int);  assert(dim <= maxd);
    for(int i = 0; i < dim; i++){
        x[i] = a[i] + corput(n + 1, base[i]) * (b[i] - a[i]);
    }
}

complex plain_quasi_MonteCarlo(int dim, double f(int dim, double* x), double* a, double* b, int N){
    double volume = 1;
    
    for(int i = 0; i < dim; i++){
        volume *= b[i] - a[i];
    }
    
    double sum = 0;
    double sum2 = 0;
    double x[dim];
    
    for(int i = 0; i < N/2; i++){
        halton(i, dim, a, b, x);
        sum += f(dim, x);
    }
    for(int i = 0; i < N/2; i++){
        halton2(i, dim, a, b, x);
        sum2 += f(dim, x);
    }
    double integral = (sum + sum2) / N * volume;
	double error = fabs(sum - sum2) / N * volume;
	return integral+I*error;
}
