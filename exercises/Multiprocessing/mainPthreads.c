#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<pthread.h>

struct params {int n_points, points_inside; unsigned int seed;};

int InsideOrOutside(double x, double y){
    if(sqrt(pow(x, 2) + pow(y, 2)) < 1){
        return 1;
    } else {
        return 0;
    }
}

void* generate_points(void* args){
    struct params * p = (struct params*)args;
    double x[p->n_points];
    double y[p->n_points];
    for (int i = 0; i < p->n_points; i++){
        x[i]= (double)rand_r(&(p->seed)) / (double)RAND_MAX;
        y[i]= (double)rand_r(&(p->seed)) / (double)RAND_MAX;
        p->points_inside += InsideOrOutside(x[i], y[i]);
    }
    return NULL;
}

int main(){
    int N = 1000000;
    
    struct params param1 = {.n_points = N/2, .points_inside = 0, .seed = 42};
    struct params param2 = {.n_points = N/2, .points_inside = 0, .seed = 1337};
    
    pthread_t t1, t2;
    pthread_create(&t1, NULL, generate_points, (void*)&param1);
    pthread_create(&t2, NULL, generate_points, (void*)&param2);
    pthread_join(t1, NULL);
    pthread_join(t2, NULL);
    
    int N_inside = param1.points_inside + param2.points_inside;
    int N_total = param1.n_points + param2.n_points;
    double pi_value = 4*(double)N_inside/N_total;
    printf("Pi is approximately equal to %g\n", pi_value);
}
