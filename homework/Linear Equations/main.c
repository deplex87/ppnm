#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <math.h>
#include "GramSchmidt.h"

void vector_print(char s[], gsl_vector* v){
    printf("%s\n",s);
    for(int i = 0; i < v->size; i++){
        printf("%10g \n",gsl_vector_get(v,i));
    }
    printf("\n");
}

void matrix_print(char s[], gsl_matrix* A){
    int n = A->size1;
    int m = A->size2;
    
    for(int i = 0; i < n; i++){
        for(int j = 0; j < m; j++){
            if(fabs(gsl_matrix_get(A,i,j))<10e-7)gsl_matrix_set(A,i,j,0);
        }
    }
    printf("%s\n",s);
    
    for(int i = 0; i < n; i++){
        for(int j = 0; j < m; j++){
            printf("%10g ",gsl_matrix_get(A,i,j));
        }
        printf("\n");
    }
    printf("\n");
}

void generate_random_matrix(gsl_matrix* A){
    for(int i = 0; i < A->size1; i++){
        for(int j = 0; j < A->size2; j++){
            gsl_matrix_set(A, i, j, (double)rand()/RAND_MAX);
        }
    }
}

void generate_random_vector(gsl_vector* b){
    for(int i = 0; i < b->size; i++){
        gsl_vector_set(b, i, (double)rand()/RAND_MAX);
    }
}

int main(){
    int n = 8; int m = 7;
    
    gsl_matrix* A = gsl_matrix_alloc(n, m);
    gsl_matrix* A_copy = gsl_matrix_alloc(n, m);
    gsl_matrix* R = gsl_matrix_alloc(m, m);
    gsl_matrix* QR = gsl_matrix_alloc(n, m);
    gsl_matrix* QTQ = gsl_matrix_alloc(m, m);
    
    generate_random_matrix(A);
    gsl_matrix_memcpy(A_copy, A);
    GS_decomp(A, R);
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, A, A, 0, QTQ);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, A, R, 0, QR);
    
    printf("Task A part 1\n\n");
    matrix_print("We start by generating a random matrix of size (8, 7):", A_copy);
    matrix_print("Then we factorize the matrix into the matrices Q", A);
    matrix_print("and R, which we see is upper triangular:", R);
    matrix_print("We furthermore see that Q^TQ equals the identity matrix:", QTQ);
    matrix_print("And lastly we see that QR equals our original matrix:", QR);
    
    free(A);
    free(A_copy);
    free(R);
    free(QR);
    free(QTQ);
    
    gsl_matrix* B = gsl_matrix_alloc(n, n);
    gsl_matrix* B_copy = gsl_matrix_alloc(n, n);
    gsl_vector* b = gsl_vector_alloc(n);
    gsl_matrix* P = gsl_matrix_alloc(n, n);
    gsl_vector* x = gsl_vector_alloc(n);
    gsl_vector* Ax = gsl_vector_alloc(n);
    
    generate_random_matrix(B);
    gsl_matrix_memcpy(B_copy, B);
    generate_random_vector(b);
    GS_decomp(B, P);
    GS_solve(B, P, b, x);
    gsl_blas_dgemv(CblasNoTrans, 1, B_copy, x, 0, Ax);
    
    printf("\n\nTask A part 2\n\n");
    matrix_print("We start by generating a random square matrix of size 8:", B_copy);
    vector_print("And a random vector of the same size:", b);
    matrix_print("Then we factorize the matrix into the matrices Q", B);
    matrix_print("and R, which we see is upper triangular:", P);
    vector_print("Then we can solve the equation QRx = b for x:", x);
    vector_print("And finally see that Ax = b:", Ax);

    free(B);
    free(B_copy);
    free(P);
    free(b);
    free(x);
    free(Ax);
    
    gsl_matrix* C = gsl_matrix_alloc(n, n);
    gsl_matrix* C_copy = gsl_matrix_alloc(n, n);
    gsl_matrix* S = gsl_matrix_alloc(n, n);
    gsl_matrix* D = gsl_matrix_alloc(n, n);
    gsl_matrix* AB = gsl_matrix_alloc(n, n);
    gsl_matrix* BA = gsl_matrix_alloc(n, n);
    
    generate_random_matrix(C);
    gsl_matrix_memcpy(C_copy, C);
    GS_decomp(C, S);
    GS_inverse(C, S, D);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, C_copy, D, 0, AB);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, D, C_copy, 0, BA);
    
    printf("\n\nTask B\n\n");
    matrix_print("We start by generating a random square matrix of size 8:", C_copy);
    matrix_print("Then we factorize the matrix into the matrices Q", C);
    matrix_print("and R, which we see is upper triangular:", S);
    matrix_print("Then we calculate the inverse matrix:", D);
    matrix_print("Which we can verify by observing that AB =", AB);
    matrix_print("And BA =", BA);
    
    free(C);
    free(C_copy);
    free(S);
    free(D);
    free(AB);
    free(BA);
    
    return 0;
}
