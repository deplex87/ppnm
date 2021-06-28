#include<math.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include "Jacobi.h"

void generate_symmetric_matrix(gsl_matrix* A){
	for(int i = 0; i < A->size1; i++){
		gsl_matrix_set(A, i, i, (double)rand()/RAND_MAX);
		for(int j = i + 1; j < A->size2; j++){
			double number = (double)rand()/RAND_MAX;
			gsl_matrix_set(A, i, j, number);
			gsl_matrix_set(A, j, i, number);
        }
	}
}

void matrix_print(char s[], gsl_matrix* A){
	int n=A->size1, m=A->size2;
	for(int i=0;i<n;i++){
		for(int j=0; j<m;j++){
			if(fabs(gsl_matrix_get(A,i,j))<10e-6)gsl_matrix_set(A,i,j,0);
		}
	}
	printf("%s\n",s);
	for(int i=0;i< n ;i++){
		for(int j=0;j< m ;j++)printf("%10g ",gsl_matrix_get(A,i,j));
		printf("\n");}
	printf("\n");
}

int main(){

    int n = 10;
    gsl_matrix* Matrix = gsl_matrix_alloc(n,n);
    gsl_matrix* Matrix_copy = gsl_matrix_alloc(n,n);
    generate_symmetric_matrix(Matrix);
    gsl_matrix_memcpy(Matrix_copy, Matrix);
    
    gsl_matrix* EigenVectorMatrix = gsl_matrix_alloc(n,n);
    jacobi_diag(Matrix, EigenVectorMatrix);
    
    gsl_matrix* AV = gsl_matrix_alloc(n,n);
    gsl_matrix* VTAV = gsl_matrix_alloc(n,n);
    gsl_matrix* VD = gsl_matrix_alloc(n,n);
    gsl_matrix* VDVT = gsl_matrix_alloc(n,n);
    gsl_matrix* VTV = gsl_matrix_alloc(n,n);
    
    printf("Task A\n\n");
    
    matrix_print("To test that the implementation works, we start by generating a random symmetric matrix A", Matrix_copy);
    matrix_print("Then we run the Jacobi eigenvalue algorithm and see that the we have the eigenvalues in a matrix D given by", Matrix);
    matrix_print("with corresponding eigen vectors in a matrix V given by", EigenVectorMatrix);

    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, EigenVectorMatrix, Matrix_copy, 0, AV);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, AV, EigenVectorMatrix, 0, VTAV);
    matrix_print("Now we can calculate V^(T)AV and see that we obtain D again:", VTAV);

    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, EigenVectorMatrix, Matrix, 0, VD);
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1, VD, EigenVectorMatrix, 0, VDVT);
    matrix_print("Furthermore we calculate VDV^(T) and see that we get back our original matrix", VDVT);

    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, EigenVectorMatrix, EigenVectorMatrix, 0, VTV);
    matrix_print("Lastly we can calculate that V^(T)V is equal to", VTV);
    printf("and hence we can conclude that the implementation works as we would expect.\n\n");
    
    printf("\n\nTask B\n\n");
    

    FILE* myOutputFilestream = fopen("EigenEnergies.txt", "w");

    int N = 60;
    double s = 1.0 / (N + 1);
    gsl_matrix* hamiltonian = gsl_matrix_alloc(N, N);
    for(int i = 0; i < N - 1; i++){
        gsl_matrix_set(hamiltonian, i, i, -2);
        gsl_matrix_set(hamiltonian, i, i + 1, 1);
        gsl_matrix_set(hamiltonian, i + 1, i, 1);
    }
    gsl_matrix_set(hamiltonian, N - 1, N - 1, -2);
    gsl_matrix_scale(hamiltonian, -1 / s / s);
    matrix_print("We start by bulding the hamiltonian matrix:", hamiltonian);

    gsl_matrix* eigStates = gsl_matrix_alloc(N, N);
    jacobi_diag(hamiltonian, eigStates);
    matrix_print("Then we diagonalize our hamiltonian with our algorithm and get the eigenvalues:", hamiltonian);

    printf("Now we check that the eigen-energies are correct:\n");
    printf("#\tCalculated\tExact\n");
    for (int k = 0; k < N / 3; k++){
        double exact = M_PI*M_PI*(k + 1)*(k + 1);
        double calculated = gsl_matrix_get(hamiltonian, k, k);
        printf("%i\t%.5g\t%.5g\n", k, calculated, exact);
    }

    for(int k = 0; k < 3; k++) {
        fprintf(myOutputFilestream, "0\t0\t");
    }
    fprintf(myOutputFilestream, "\n");
    for(int i = 0; i < N; i++){
        fprintf(myOutputFilestream, "%.5g\t", (i + 1.0) / (N + 1));
        for(int j = 0; j < 3; j++) {
            fprintf(myOutputFilestream, "%.5g\t", gsl_matrix_get(eigStates, i, j));
        }
            double k=(i+1.0)/(N+1);
            fprintf(myOutputFilestream, "%.5g\t%.5g\t%.5g\t", sqrt(1/M_PI)*sin(k*M_PI), -sqrt(1/M_PI)*sin(k*2*M_PI)+1.5, sqrt(1/M_PI)*sin(3*k*M_PI)+3);
        fprintf(myOutputFilestream, "\n");
    }
    fprintf(myOutputFilestream, "1\t");
    for(int k = 0; k < 5; k++) {
        fprintf(myOutputFilestream, "0\t");
    }
    fprintf(myOutputFilestream, "\n");

    return 0;
} 
