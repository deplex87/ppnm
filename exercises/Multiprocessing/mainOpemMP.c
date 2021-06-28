#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<omp.h>

int InsideOrOutside(double x, double y){
    if(sqrt(pow(x, 2) + pow(y, 2)) < 1){
        return 1;
    } else {
        return 0;
    }
}

int main(){
    int N = 100000;
    double x[N];
    double y[N];
    
    unsigned int *seeds = malloc(2);
    seeds[1] = 42;
    seeds[2] = 1337;
    
#pragma omp parallel sections
    {
#pragma omp section
        {
            for (int i = 0; i < N; i++) x[i]= (double)rand_r(&seeds[1]) / (double)RAND_MAX;
        }
#pragma omp section
        {
            for (int i = 0; i < N; i++) y[i]= (double)rand_r(&seeds[2]) / (double)RAND_MAX;
        }
    }
    
    int N_inside = 0;
    
    for(int i = 0; i < N; i++){
        N_inside += InsideOrOutside(x[i], y[i]);
    }
    
    double pi_value = 4*(double)N_inside/(double)N;
    
    printf("Pi is approximately equal to %g\n", pi_value);  
}
