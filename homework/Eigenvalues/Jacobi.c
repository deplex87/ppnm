#include <math.h>
#include <gsl/gsl_matrix.h>

void timesJ(gsl_matrix* matrix, int i, int j, double angle){
    double c = cos(angle);
    double s = sin(angle);
    for(int k = 0; k < matrix->size1; k++){
        double new_aip = c * gsl_matrix_get(matrix, k, i)  -  s * gsl_matrix_get(matrix, k, j);
        double new_aiq = s * gsl_matrix_get(matrix, k, i)  +  c * gsl_matrix_get(matrix, k, j);

        gsl_matrix_set(matrix, k, i,  new_aip);
        gsl_matrix_set(matrix, k, j, new_aiq);
    }
}


void Jtimes(gsl_matrix* matrix, int i, int j, double angle){
    double c = cos(angle);
    double s = sin(angle);
    for(int k = 0; k < matrix->size2; k++){
        double new_apj =   c * gsl_matrix_get(matrix, i, k)  +  s * gsl_matrix_get(matrix, j, k);
        double new_aqj = - s * gsl_matrix_get(matrix, i, k)  +  c * gsl_matrix_get(matrix, j, k);
        gsl_matrix_set(matrix, i,  k, new_apj);
        gsl_matrix_set(matrix, j, k, new_aqj);
    }
}


void jacobi_diag(gsl_matrix* OriginalMatrix, gsl_matrix* EigenVectorMatrix){
    int changed;
    gsl_matrix_set_identity(EigenVectorMatrix);
    do{
        changed=0;
        for(int i = 0; i < OriginalMatrix->size1 - 1; i++){
            for(int j = i + 1; j < OriginalMatrix->size2; j++){
                
                double apq = gsl_matrix_get(OriginalMatrix, i, j);
                double app = gsl_matrix_get(OriginalMatrix, i, i);
                double aqq = gsl_matrix_get(OriginalMatrix, j, j);
                
                double angle = 0.5 * atan2(2 * apq, aqq - app);
                
                double c = cos(angle);
                double s = sin(angle);
                double new_app = c*c*app - 2*s*c*apq + s*s*aqq;
                double new_aqq = s*s*app + 2*s*c*apq + c*c*aqq;
                if(new_app!=app || new_aqq!=aqq){
                    changed=1;
                    timesJ(OriginalMatrix, i, j, angle);
                    Jtimes(OriginalMatrix, i, j, -angle);
                    timesJ(EigenVectorMatrix, i, j, angle);
                }
            }
        }
    }while(changed!=0);
}
