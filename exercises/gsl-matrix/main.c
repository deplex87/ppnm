#include<stdio.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h>

void vector_print(char s[], gsl_vector* v){
	printf("%s\n",s);
	for(int i=0;i<v->size;i++){
		printf("%10g ",gsl_vector_get(v,i));
	}
	printf("\n");
}

int main(){
	int n=3;
	gsl_matrix* A=gsl_matrix_alloc(n,n);
	gsl_matrix* Acopy=gsl_matrix_alloc(n,n);
	gsl_vector* b=gsl_vector_alloc(n);
	gsl_vector* x=gsl_vector_alloc(n);
	gsl_vector* y=gsl_vector_calloc(n);

	gsl_matrix_set(A,0,0,6.13);
    gsl_matrix_set(A,0,1,-2.90);
    gsl_matrix_set(A,0,2,5.86);
    gsl_matrix_set(A,1,0,8.08);
    gsl_matrix_set(A,1,1,-6.31);
    gsl_matrix_set(A,1,2,-3.89);
    gsl_matrix_set(A,2,0,-4.36);
    gsl_matrix_set(A,2,1,1.00);
    gsl_matrix_set(A,2,2,0.19);

	gsl_vector_set(b,0,6.23);
	gsl_vector_set(b,1,5.37);
	gsl_vector_set(b,2,2.29);

	gsl_matrix_memcpy(Acopy,A);

	gsl_linalg_HH_solve(Acopy,b,x);
	gsl_blas_dgemv(CblasNoTrans,1,A,x,0,y);
	vector_print("We have the left side b =",b);
	vector_print("Using GSL's Householding solver to solve for x in Ax = b, we should now expect that b = Ax =",y);
	printf("Which we see is true, so the solver most work.\n The resulting x vector is:\n");
	vector_print("x =",x);
	gsl_matrix_free(A);
	gsl_matrix_free(Acopy);
	gsl_vector_free(b);
	gsl_vector_free(x);
	gsl_vector_free(y);
return 0;
}
