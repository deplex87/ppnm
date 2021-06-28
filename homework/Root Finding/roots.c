#include <float.h>
#include <assert.h>
#include <math.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include "GramSchmidt.h"

void newtonMethod(void f(gsl_vector* x, gsl_vector* f_values), gsl_vector* start, double eps){
    double dx = sqrt(eps);
    int dim = start -> size;

    // We define the jacobi matrix and the triangular matrix we need later when we solve the linear equations
    gsl_matrix* j_matrix = gsl_matrix_alloc(dim, dim);
    gsl_matrix* t_matrix = gsl_matrix_alloc(dim, dim);
    
    gsl_vector* next_x = gsl_vector_alloc(dim);
    gsl_vector* f_value = gsl_vector_alloc(dim);
    gsl_vector* f_value_copy = gsl_vector_alloc(dim);
    gsl_vector* next_f_value = gsl_vector_alloc(dim);
    
    f(start, f_value); // We get the function value at the start
    
    gsl_vector* solution = gsl_vector_alloc(dim);
    gsl_vector* scaled_solution = gsl_vector_alloc(dim);

    
    int count = 0;
    while (gsl_blas_dnrm2(f_value) > eps) {
        count++;
        assert(count < 1e5);

        for (int i = 0; i < dim; i++){
            gsl_vector_set(start, i, gsl_vector_get(start, i) + dx);
            f(start, f_value_copy);

            for (int j = 0; j < dim; j++){
                gsl_matrix_set(j_matrix, j, i, (double) (gsl_vector_get(f_value_copy, j) - gsl_vector_get(f_value, j)) / dx);
            }
            gsl_vector_set(start, i, gsl_vector_get(start, i) - dx);
        }
        gsl_vector_scale(f_value, -1.0);
        GS_decomp(j_matrix, t_matrix);
        GS_solve(j_matrix, t_matrix, f_value, solution);
        gsl_vector_scale(f_value, -1.0);

        double scale = 2;
        while ((gsl_blas_dnrm2(next_f_value) >= (1 - scale / 2) * gsl_blas_dnrm2(f_value) ) && scale >= 0.02 ){
            scale /= 2;
            gsl_vector_memcpy(scaled_solution, solution);
            gsl_vector_scale(scaled_solution, scale);
            gsl_vector_memcpy(next_x, scaled_solution);
            gsl_vector_add(next_x, start);

            f(next_x, next_f_value);
        }
        gsl_vector_memcpy(start, next_x);
        gsl_vector_memcpy(f_value, next_f_value);

        if (gsl_blas_dnrm2(solution) < dx){
            break;
        }
    }

    gsl_matrix_free(t_matrix);
    gsl_matrix_free(j_matrix);
    gsl_vector_free(f_value);
    gsl_vector_free(f_value_copy);
    gsl_vector_free(next_x);
    gsl_vector_free(next_f_value);
    gsl_vector_free(solution);
    gsl_vector_free(scaled_solution);
}

