#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include "LeastSquares.h"

double fun(int i, double x){
	switch(i){
		case  0: return 1  ; break;
		case  1: return -x  ; break;
		default: return NAN;
		}
	}

	
void vector_print(char s[], gsl_vector* v){
	printf("%s\n",s);
	for(int i=0;i< v->size ;i++)printf("%10g \n",gsl_vector_get(v,i));
	printf("\n");
}

void matrix_print(char s[], gsl_matrix* A){
	printf("%s\n",s);
	for(int i=0;i< A->size1 ;i++){							// Note til selv size1=vertical, size2=horisontal
		for(int j=0;j< A->size2 ;j++)printf("%10g ",gsl_matrix_get(A,i,j));
		printf("\n");}
	printf("\n");
}
	
	
int main(){
    
    double t[] = {1, 2, 3, 4, 6, 9, 10, 13, 15};
    double y[] = {117, 100, 88, 72, 53, 29.5, 25.2, 15.2, 11.1};
    
	int numOfPoints = 9;
    int numOfFuncs = 2;
    
    gsl_vector* coefficients = gsl_vector_alloc(numOfFuncs);
	gsl_matrix* covMatrix = gsl_matrix_alloc(numOfFuncs, numOfFuncs);
    
	double unc_y[numOfPoints], ln_y[numOfPoints], unc_ln_y[numOfPoints];
    
	for(int i = 0; i < numOfPoints; i++){
        unc_y[i] = y[i]/20.0;
        ln_y[i] = log(y[i]);
        unc_ln_y[i] = unc_y[i]/y[i];
    }
    
	double (*f)(int,double);
    f=fun;
	least_square(numOfPoints, numOfFuncs, t, ln_y, unc_ln_y, f, coefficients, covMatrix);
    
    
	printf("#index -1: Task A and B\n\n");
	printf("The coefficients found from doing the linear fit is:\n");
    printf("a      = %lg +/- %lg \n", gsl_vector_get(coefficients,0), sqrt(gsl_matrix_get(covMatrix, 0, 0)));
    printf("lambda = %lg +/- %lg \n", gsl_vector_get(coefficients,1), sqrt(gsl_matrix_get(covMatrix, 1, 1)));
	printf("From here we find that the half-life time of 224Ra is:\n");
    printf("t_1/2 = %lg +/- %lg days \n", log(2)/gsl_vector_get(coefficients,1), sqrt(gsl_matrix_get(covMatrix,1,1))/gsl_vector_get(coefficients,1)/gsl_vector_get(coefficients,1));
	matrix_print("which is calculated from the covariance matrix of the coefficients given by:", covMatrix);
	printf("The result is not consistent with the modern half-life time of 3.66 days, which could be due to the presence of other radium isoptopes in the source. For example 223Ra have a half life time of 11.43 days and 226Ra have a half life time of 1600 years. Both 223Ra and 226Ra come from alpha decay of radon isotopes, just like 224Ra does, so this could definitely be an explanation.\n\n");
    
    FILE* outFileStream  =  fopen("data.txt", "w");
    
    fprintf(outFileStream, "#index 0: t y log(y) unc_y unc_log(y)\n");
	for(int i = 0; i < numOfPoints; i++){
        fprintf(outFileStream, "%7g %7g %7g %7g %7g\n", t[i], y[i], ln_y[i], unc_y[i], unc_ln_y[i]);
    }

	fprintf(outFileStream, "\n\n#index 1: t_linspace, fit, fit + uncertainty, fit - uncertainty\n");
	int N=250;
    double t_linspace[N], fit[N], fitplus[N], fitminus[N];
	for(int i = 0; i < N; i++){
        t_linspace[i]=(double)(i*16)/N;
        
        fit[i] = gsl_vector_get(coefficients, 0) - t_linspace[i]*gsl_vector_get(coefficients, 1);
        
        fitplus[i] = gsl_vector_get(coefficients, 0) + sqrt(gsl_matrix_get(covMatrix, 0, 0)) - t_linspace[i]*(gsl_vector_get(coefficients, 1) + sqrt(gsl_matrix_get(covMatrix, 1, 1)));
        
        fitminus[i] = gsl_vector_get(coefficients, 0) - sqrt(gsl_matrix_get(covMatrix, 0, 0)) - t_linspace[i]*(gsl_vector_get(coefficients, 1) - sqrt(gsl_matrix_get(covMatrix, 1, 1)));
        
        fprintf(outFileStream, "%7g %7g %7g %7g\n", t_linspace[i], fit[i], fitplus[i], fitminus[i]);
    }
        fclose(outFileStream);


return 0;
} 
