#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <math.h>

void rkstep12(
	void f(double t, gsl_vector* y,gsl_vector* dydt), /* the f from dy/dt=f(t,y) */
	double t,              /* the current value of the variable */
	gsl_vector* y_val,            /* the current value y(t) of the sought function */
	double step,              /* the step to be taken */
	gsl_vector* y_step_val,             /* output: y(t+h) */
	gsl_vector* err             /* output: error estimate */
){
    
    int order = y_val->size;
    
    gsl_vector* k0 = gsl_vector_alloc(order);
    gsl_vector* k12 = gsl_vector_alloc(order);
    gsl_vector* y_val_temp = gsl_vector_alloc(order);
    
    f(t, y_val, k0);
    
    for(int i = 0; i < order; i++){
        gsl_vector_set(y_val_temp, i, gsl_vector_get(y_val, i) + 0.5*step*gsl_vector_get(k0, i));
    }
    
    f(t + 0.5*step, y_val_temp, k12);
    
    for(int i = 0; i < order; i++){
        gsl_vector_set(y_step_val, i, gsl_vector_get(y_val, i) + step*gsl_vector_get(k12, i));
    }
    
    for(int i = 0; i < order; i++){
        gsl_vector_set(err, i, (gsl_vector_get(k0, i) - gsl_vector_get(k12, i))*step/2);
    }
    
    gsl_vector_free(k0);
    gsl_vector_free(k12);
    gsl_vector_free(y_val_temp);
}

void driver(
	void f(double,gsl_vector*,gsl_vector*), /* right-hand-side of dy/dt=f(t,y) */
	double start_point,                     /* the start-point a */
	gsl_vector* ya,                     /* y(a) */
	double end_point,                     /* the end-point of the integration */
	gsl_vector* yb,                     /* y(b) to be calculated */
	double step,                     /* initial step-size */
	double acc,                   /* absolute accuracy goal */
	double eps,                    /* relative accuracy goal */
    FILE* filepath
){
    
    int order = ya->size;
    double err, tol, norm_y_val;
    gsl_vector* y_step_val = gsl_vector_alloc(order);
    gsl_vector* y_val_err = gsl_vector_alloc(order);
    gsl_vector* current_y_val = gsl_vector_alloc(order);
    gsl_vector_memcpy(current_y_val, ya);
    
    double current_point = start_point;
    
    while(current_point < end_point){
        
        fprintf(filepath, "%.6g\t", current_point);
        
        for(int i = 0; i < order; ++i){
            fprintf(filepath,"%.6g\t", gsl_vector_get(current_y_val, i));
        }
        
        fprintf(filepath, "\n");
        
        double nextStep;
        
        if(current_point + step > end_point){
            step = end_point - current_point;
        }
        
        do{
        rkstep12(f, current_point, current_y_val, step, y_step_val, y_val_err);
        
        err = gsl_blas_dnrm2(y_val_err);
        norm_y_val = gsl_blas_dnrm2(y_step_val);
        tol = (norm_y_val * eps + acc) * sqrt(step / (end_point - start_point));
        nextStep = step;
        if(err > 0) step *= pow(tol / err, 0.25)*0.95; else step *= 2;
        }while(err > tol);
        
        current_point += nextStep;
        gsl_vector_memcpy(current_y_val, y_step_val);
    }
    gsl_vector_memcpy(yb, y_step_val);
    
    gsl_vector_free(y_step_val);
    gsl_vector_free(y_val_err);
}
