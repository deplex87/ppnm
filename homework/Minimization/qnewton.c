#include<stdio.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<float.h>

void numeric_gradient(double f(gsl_vector*), gsl_vector*x, gsl_vector*grad){
    
    long double StepSize = sqrt(2.22045e-16);
    
	double f_value =f(x);
    
	for(int i = 0; i < x->size; i++){
        
		double step, xi = gsl_vector_get(x,i);
        
		if(fabs(xi) < sqrt(StepSize)){
            step = StepSize;
        }else{
            step = fabs(xi) * StepSize;
        }
        
		gsl_vector_set(x, i, xi + step);
		gsl_vector_set(grad, i, (f(x) - f_value) / step);
		gsl_vector_set(x, i, xi);
	}
}

int qnewton(double F(gsl_vector*), gsl_vector*x, double acc) {
    
    long double StepSize = sqrt(2.22045e-16);
    
	int dim = x->size;
	int numOfSteps = 0;
    
	gsl_matrix* inverse_HMatrix = gsl_matrix_alloc(dim,dim);
	gsl_vector* gradient = gsl_vector_alloc(dim);
	gsl_vector* Dx = gsl_vector_alloc(dim);
	gsl_vector* nextMinimum = gsl_vector_alloc(dim);
	gsl_vector* nextGradient = gsl_vector_alloc(dim);
	gsl_vector* solution = gsl_vector_alloc(dim);
	gsl_vector* Dy = gsl_vector_alloc(dim);
    
	gsl_matrix_set_identity(inverse_HMatrix);
	numeric_gradient(F,x,gradient);
    
	double f_value = F(x);
    double next_f_value;
    
	while(numOfSteps < 1000){
		numOfSteps++;
        
		gsl_blas_dgemv(CblasNoTrans, -1, inverse_HMatrix, gradient, 0, Dx);
        
		if(gsl_blas_dnrm2(Dx) < StepSize * gsl_blas_dnrm2(x)) break;
		if(gsl_blas_dnrm2(gradient) < acc) break;
        
		double lambda=1;
		while(1){
            
			gsl_vector_memcpy(nextMinimum, x);
			gsl_vector_add(nextMinimum, Dx);
            
			next_f_value = F(nextMinimum);
			double sTg; gsl_blas_ddot(Dx, gradient, &sTg);
            
			if(next_f_value < f_value + 0.01*sTg) break;
			if(lambda < StepSize){ gsl_matrix_set_identity(inverse_HMatrix); break; }
			
			lambda *= 0.5;
			gsl_vector_scale(Dx, 0.5);
		}
		
		numeric_gradient(F,nextMinimum, nextGradient);
		
        gsl_vector_memcpy(solution, nextGradient);
		gsl_blas_daxpy(-1,gradient, solution); 
		gsl_vector_memcpy(Dy, Dx);
		gsl_blas_dgemv(CblasNoTrans, -1, inverse_HMatrix, solution, 1, Dy);
        
		double sTy, uTy;
		gsl_blas_ddot(Dx, solution, &sTy);
        
		if(fabs(sTy)>1e-12){
			gsl_blas_ddot(Dy, solution, &uTy);
            
			double gamma = uTy / 2 / sTy;
            
			gsl_blas_daxpy(-gamma, Dx, Dy);
			gsl_blas_dger(1.0/sTy, Dy, Dx, inverse_HMatrix);
			gsl_blas_dger(1.0/sTy, Dx, Dy, inverse_HMatrix);
		}
		gsl_vector_memcpy(x, nextMinimum);
		gsl_vector_memcpy(gradient, nextGradient);
		f_value = next_f_value;
	}
gsl_matrix_free(inverse_HMatrix);
gsl_vector_free(gradient);
gsl_vector_free(Dx);
gsl_vector_free(nextMinimum);
gsl_vector_free(nextGradient);
gsl_vector_free(solution);
gsl_vector_free(Dy);

return numOfSteps;
}
