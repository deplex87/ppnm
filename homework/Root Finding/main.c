#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "RKODE.h"
#include "GramSchmidt.h"
#include "roots.h"

double energy_solution = 0.0;

void TestFunction(gsl_vector* vals, gsl_vector* f_values){
    double x = gsl_vector_get(vals, 0);
    double y = gsl_vector_get(vals, 1);

    gsl_vector_set(f_values, 0, 2*x*(x - 5));
    gsl_vector_set(f_values, 1, 2*y*(y - 7));

}

void Rosenbrocks_Valley_Gradient(gsl_vector* points, gsl_vector* f_values) {
    double x = gsl_vector_get(points, 0);
    double y = gsl_vector_get(points, 1);

    gsl_vector_set(f_values, 0, (-1)*2*(1 - x) + (-2*x)*2*100*(y - x*x));
    gsl_vector_set(f_values, 1, 2*100*(y - x*x));
}

void Schrodingers_Equation(double var, gsl_vector* f_values, gsl_vector* derivatives){
    double f_value = gsl_vector_get(f_values, 0);
    double FirstDerivative = gsl_vector_get(f_values, 1);
    double SecondDerivative = (-2)*(1.0/var + energy_solution)*f_value;

    gsl_vector_set(derivatives, 0, FirstDerivative);
    gsl_vector_set(derivatives, 1, SecondDerivative);
}

void WaveFunction(gsl_vector* vals, gsl_vector* f_values){

    gsl_vector* final_f_value  =  gsl_vector_alloc(2);
    gsl_vector* start_f_value   =  gsl_vector_calloc(2);

    energy_solution              =   gsl_vector_get(vals, 0);

    gsl_vector_set(start_f_value, 0, (1e-3 - 1e-3*1e-3));
    gsl_vector_set(start_f_value, 1, (1 - 2*1e-3));

    driver(Schrodingers_Equation, 1e-3, start_f_value, 10.0, final_f_value, (10.0 - 1e-3) / 10, 1e-3, 1e-3, NULL);
    gsl_vector_set(f_values, 0, gsl_vector_get(final_f_value,0) );
}

int main(){

    printf("Task A\n\n");

    printf("I test my root finding routine with the function f(x) = 1 + (x - 5)^2 + (y - 7)^2 and get the following results:\n");

    gsl_vector* minimum_guess = gsl_vector_alloc(2);
    gsl_vector_set(minimum_guess, 0, 3);
    gsl_vector_set(minimum_guess, 1, 12);

    printf("My guess of the root: x = %g and y = %g\n", gsl_vector_get(minimum_guess, 0), gsl_vector_get(minimum_guess, 1));
    newtonMethod(TestFunction, minimum_guess, 1e-3);
    printf("The root found by routine: x = %g and y = %g\n", gsl_vector_get(minimum_guess, 0), gsl_vector_get(minimum_guess, 1));
    printf("The real root is: x = 5 and y = 7\n");


    printf("\nNow we test the routine on the Rosenbrock Valley function f(x,y) = (1-x)^2+100(y-x^2)^2:\n");

    gsl_vector_set(minimum_guess, 0, 0.5);
    gsl_vector_set(minimum_guess, 1, 0.5);

    printf("My guess of the root: x = %g and y = %g\n", gsl_vector_get(minimum_guess, 0), gsl_vector_get(minimum_guess, 1));
    newtonMethod(Rosenbrocks_Valley_Gradient, minimum_guess, 1e-5);
    printf("The root found by routine: x = %g and y = %g\n", gsl_vector_get(minimum_guess, 0), gsl_vector_get(minimum_guess, 1));
    printf("The real root is: x = 1 and y = 1\n");


    printf("\nTask B\n\n");
    printf("Now we want to find the lowest energy solution to the hydrogen atom by using the root routine:\n");
    
    gsl_vector* minimum = gsl_vector_alloc(1);
    gsl_vector_set(minimum, 0, -3);
    
    newtonMethod(WaveFunction, minimum, (double) 1e-5);

    printf("Lowest energy state found by root finding = %g\n Furthermore see hydrogen.png", energy_solution);

    gsl_vector* final_f_value  =  gsl_vector_alloc(2);
    gsl_vector* start_f_value   =  gsl_vector_calloc(2);
    double start = 1e-5;
    double end = 8.0;

    gsl_vector_set(start_f_value, 0, (start - start*start));
    gsl_vector_set(start_f_value, 1, (1 - 2.*start));

    FILE* solutionfile = fopen("hydrogen.txt", "a");
    
    double AA = 1e-5;
    double RA = 1e-5;
    double step = (end - start) / 10;

    driver(Schrodingers_Equation, start, start_f_value, end, final_f_value, step, AA, RA, solutionfile);
    fclose(solutionfile);


    gsl_vector_free(minimum);
    gsl_vector_free(minimum_guess);


    return 0;
}
