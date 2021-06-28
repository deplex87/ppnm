#include <assert.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <math.h>

void GS_decomp( gsl_matrix* A, gsl_matrix* R ){


  int n = (int) A->size1; int m = (int) A->size2;
  

  assert(n >= m);  // Tjekker om kravet fra opgaven om at n>=m er overholdt


  for (int i = 0; i < m; i++){  // Algoritmen kører igennem alle søjler i = 1 til m


    double Norm_Squared = 0;
    
    for (int j = 0; j < n; j++){
        Norm_Squared += gsl_matrix_get(A, j, i)*gsl_matrix_get(A, j, i);
    }
    
    double Norm = sqrt(Norm_Squared);                        // Vi beregner normkvadratet på den i'de søjle

    gsl_matrix_set(R, i, i, Norm);                                  // Vi sætter nu diagonalelementet i R matricen lig normkvadratet


    for (int j = 0; j < n; j++){
        double a_ji = gsl_matrix_get(A, j, i);
        if (Norm > 1e-15){
            gsl_matrix_set(A, j, i, a_ji/Norm);
        }else{
            gsl_matrix_set(A, j, i, 0);
            printf("The matrix we are working with are singular.\n");
        }
    }            


    for (int j = i + 1; j < m; j++ ){
        
        double R_ij = 0;
        
        for (int k = 0; k < n; k++){
            R_ij += gsl_matrix_get(A, k, j)*gsl_matrix_get(A, k, i);
        }
        gsl_matrix_set(R, i, j, R_ij);
        
        for (int k = 0; k < n; k++){
            gsl_matrix_set(A,k,j,gsl_matrix_get(A,k,j)-gsl_matrix_get(A,k,i)*R_ij);
        }
    }
  }
}

void GS_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x){
    
	gsl_blas_dgemv(CblasTrans,1,Q,b,0,x); //Vi transformerer fra QRx = b til Rx=Q^Tb
    
    for(int i = x->size-1; i >= 0; i--){  // Derefter laver vi bare baglæns substitution
        
		double s = gsl_vector_get(x,i);
        
		for(int k = i+1; k < R->size1; k++){
            s -= gsl_matrix_get(R,i,k)*gsl_vector_get(x,k);
        }
		gsl_vector_set(x, i, s/gsl_matrix_get(R,i,i));
	}
}

void GS_inverse(gsl_matrix* Q, gsl_matrix* R, gsl_matrix* B){

	gsl_vector* b = gsl_vector_alloc(R -> size1);
	gsl_vector* x = gsl_vector_alloc(R -> size1);
    
	for(int i = 0; i < R -> size1; i++){

		gsl_vector_set(b, i, 1);

		GS_solve(Q, R, b, x);

		for(int k = 0; k < R -> size1; k++){
			gsl_matrix_set(B, k, i, gsl_vector_get(x,k));
		}
		gsl_vector_set(b, i, 0);
	}	
}
