#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include "RKODE.h"

void diffeq(double t, gsl_vector* y, gsl_vector* diffs){
    double dydt = gsl_vector_get(y, 1);
    double d2ydt2 = -gsl_vector_get(y, 0);
    gsl_vector_set(diffs, 0, dydt);
    gsl_vector_set(diffs, 1, d2ydt2);
}

void SIRmodel(double t, gsl_vector* y, gsl_vector* diffs){
    double population = 5843347;
    double Tr = 10;
    double Tc = 2;
    
    double diffSusceptiblePop = -(gsl_vector_get(y, 1) * gsl_vector_get(y, 0)) / (population * Tc);
    double diffRemovedPop = gsl_vector_get(y, 1) / Tr;
	double diffInfectiousPop = -diffSusceptiblePop - diffRemovedPop;
	gsl_vector_set(diffs,0, diffSusceptiblePop);
	gsl_vector_set(diffs,1, diffInfectiousPop);
	gsl_vector_set(diffs,2, diffRemovedPop);
}

int main(){
    
    // Solving the differential equation y'' = -y
    
    gsl_vector* FunctionValuesStart = gsl_vector_alloc(2);
    gsl_vector* FunctionValuesEnd = gsl_vector_alloc(2);
    gsl_vector_set(FunctionValuesStart, 1, 1);     
    
    FILE* diffeqpath = fopen("DifferentialEquation.txt", "w");
    
    driver(diffeq, 0.0, FunctionValuesStart, 2*M_PI, FunctionValuesEnd, 2*M_PI/100, 0.001, 0.001, diffeqpath);
    fclose(diffeqpath);
    
    FILE* sinefile = fopen("sine.txt", "w");
    for(int i = 0; i < 101; i++){
                double pos = 2*M_PI*i/100;
                fprintf(sinefile,"%.5g\t%.5g\n", pos, sin(pos));
            }
    fclose(sinefile);
    
    gsl_vector_free(FunctionValuesEnd);
    gsl_vector_free(FunctionValuesStart);
    
    // Solving the SIR model
    
    double start   =   0.0;
    double end  =   100;
    gsl_vector* SIR_end_val  =  gsl_vector_alloc(3);
    gsl_vector* SIR_start_val   =  gsl_vector_calloc(3);

    double population = 5843347;
    double recovered = 85000;
    double Infected = 92000;
    double stillInfected = Infected - recovered;
    double vaccinated = 10000;
    double removed = recovered + vaccinated;

    // Set the initial values for the population
    gsl_vector_set(SIR_start_val, 0, population - stillInfected - removed);
    gsl_vector_set(SIR_start_val, 1, stillInfected);
    gsl_vector_set(SIR_start_val, 2, removed);

    FILE* SIRpath = fopen("SIRmodel.txt", "w");
    driver(SIRmodel, start, SIR_start_val, end, SIR_end_val, 1, 0.001, 0.001, SIRpath);
    fclose(SIRpath);// Close file stream
    
    return 0;
}
