#include <math.h>
#include <stdio.h>
#include <gsl/gsl_integration.h>

#include "RAintegrator.h"

int main(){
    
    printf("Task A\n\n");
    printf("We test our adaptive integrator on different integral:\n");
    
    int numOfEvals = 0;

    double func1(double t){
        numOfEvals += 1;
        return sqrt(t);
    }
    
    double IntegrationError = 0;
    double IntegralResult = integrate(func1, 0, 1, 0.001, 0.001, &IntegrationError);
    
    printf("Our numerical integration of the integral of sqrt(x) from 0 to 1 gives the following results:\n");
    printf("Analytical result = %.10g\n", (double) 2/3);
    printf("Numerical result = %.10g\n", IntegralResult);
    printf("Estimated error = %.10g\n", IntegrationError);
    printf("Actual error = %.10g\n", fabs(IntegralResult - 2.0/3.0));
    printf("numOfEvals = %i\n\n", numOfEvals);
    
    numOfEvals = 0;
    IntegrationError = 0;
    
    double func2(double t){
        numOfEvals += 1;
        return 4 * sqrt(1 - t * t);
    }
    
    IntegralResult = integrate(func2, 0, 1, 0.001, 0.001, &IntegrationError);
    
    printf("Our numerical integration of the integral of 4 * sqrt(1 - x^2) from 0 to 1 gives the following results:\n");
    printf("Analytical result = %.10g\n", M_PI);
    printf("Numerical result = %.10g\n", IntegralResult);
    printf("Estimated error = %.10g\n", IntegrationError);
    printf("Actual error = %.10g\n", fabs(IntegralResult - M_PI));
    printf("numOfEvals = %i\n\n", numOfEvals);
    
    printf("Task B\n\n");
    
    numOfEvals = 0;
    IntegrationError = 0;
    
    double func3(double t){
        numOfEvals += 1;
        return 1 / sqrt(t);
    }
    
    IntegralResult = OpenQuadCC_integrate(func3, 0, 1, 0.001, 0.001, &IntegrationError);
    
    printf("Our numerical integration of the integral of 1 / sqrt(x) from 0 to 1 using open quadrature with Clenshaw–Curtis variable transformation gives the following results:\n");
    printf("Analytical result = %.d\n", 2);
    printf("Numerical result = %.10g\n", IntegralResult);
    printf("Estimated error = %.10g\n", IntegrationError);
    printf("Actual error = %.10g\n", fabs(IntegralResult - 2));
    printf("numOfEvals = %i\n\n", numOfEvals);
    
    size_t gslnumOfEvals;
    double result;
    double gslerror;
    gsl_integration_cquad_workspace* workspace = gsl_integration_cquad_workspace_alloc(1000);
    gsl_function function1_gsl;
    double gslfunction1(double t, void* params){
        params = NULL;
        return func3(t);
    }
    function1_gsl.function = gslfunction1;
    function1_gsl.params = NULL;
    gsl_integration_cquad(&function1_gsl, 0, 1, 0.001, 0.001, workspace, &result, &gslerror, &gslnumOfEvals);
    printf("GSL numerical integration of the integral of 1 / sqrt(x) from 0 to 1 using gsl_integration_cquad the following results:\n");
    printf("Analytical result = %.d\n", 2);
    printf("Numerical result = %.10g\n", result);
    printf("Estimated error = %.10g\n", gslerror);
    printf("Actual error = %.10g\n", fabs(result - 2));
    printf("numOfEvals = %li\n\n", gslnumOfEvals);
    
    numOfEvals = 0;
    IntegrationError = 0;
    
    double func4(double t){
        numOfEvals += 1;
        return log(t) / sqrt(t);
    }
    
    IntegralResult = OpenQuadCC_integrate(func4, 0, 1, 0.001, 0.001, &IntegrationError);
    
    printf("Our numerical integration of the integral of ln(x) / sqrt(x) from 0 to 1 using open quadrature with Clenshaw–Curtis variable transformation gives the following results:\n");
    printf("Analytical result = %.d\n", -4);
    printf("Numerical result = %.10g\n", IntegralResult);
    printf("Estimated error = %.10g\n", IntegrationError);
    printf("Actual error = %.10g\n", fabs(IntegralResult + 4));
    printf("numOfEvals = %i\n\n", numOfEvals);
    
    gsl_function function2_gsl;
    double gslfunction2(double t, void* params){
        params = NULL;
        return func4(t);
    }
    function2_gsl.function = gslfunction2;
    function2_gsl.params = NULL;
    gsl_integration_cquad(&function2_gsl, 0, 1, 0.001, 0.001, workspace, &result, &gslerror, &gslnumOfEvals);
    printf("GSL numerical integration of the integral of ln(x) / sqrt(x) from 0 to 1 using gsl_integration_cquad the following results:\n");
    printf("Analytical result = %.d\n", -4);
    printf("Numerical result = %.10g\n", result);
    printf("Estimated error = %.10g\n", gslerror);
    printf("Actual error = %.10g\n", fabs(result + 4));
    printf("numOfEvals = %li\n\n", gslnumOfEvals);
    
    numOfEvals = 0;
    IntegrationError = 0;
    
    IntegralResult = OpenQuadCC_integrate(func2, 0, 1, 0.000001, 0.000001, &IntegrationError);
    
    printf("Our numerical integration of the integral of 4 * sqrt(1 - x^2) from 0 to 1 using open quadrature with Clenshaw–Curtis variable transformation gives the following results:\n");
    printf("Analytical result = %.30g\n", M_PI);
    printf("Numerical result = %.30g\n", IntegralResult);
    printf("Estimated error = %.30g\n", IntegrationError);
    printf("Actual error = %.30g\n", fabs(IntegralResult - M_PI));
    printf("numOfEvals = %i\n\n", numOfEvals);
    
        gsl_function function3_gsl;
    double gslfunction3(double t, void* params){
        params = NULL;
        return func2(t);
    }
    function3_gsl.function = gslfunction3;
    function3_gsl.params = NULL;
    gsl_integration_cquad(&function3_gsl, 0, 1, 0.0000000000001, 0.0000000000001, workspace, &result, &gslerror, &gslnumOfEvals);
    printf("GSL numerical integration of the integral of 4 * sqrt(1 - x^2) from 0 to 1 using gsl_integration_cquad the following results:\n");
    printf("Analytical result = %.25g\n", M_PI);
    printf("Numerical result = %.25g\n", result);
    printf("Estimated error = %.25g\n", gslerror);
    printf("Actual error = %.25g\n", fabs(result - M_PI));
    printf("numOfEvals = %li\n\n", gslnumOfEvals);

    printf("Task C\n\n");
    
//     numOfEvals = 0;
//     IntegrationError = 0;
//     
//     double func5(double t){
//         numOfEvals += 1;
//         return 1 / sqrt(1 + t * t);
//     }
//     
//     IntegralResult = integrate(func5, 0, INFINITY, 0.0001, 0.0001, &IntegrationError);
//     
//     printf("Our numerical integration of the integral of 1 / sqrt(1 + x^2) from -inf to inf gives the following results:\n");
//     printf("Analytical result = %.10g\n", M_PI/2);
//     printf("Numerical result = %.10g\n", IntegralResult);
//     printf("Estimated error = %.10g\n", IntegrationError);
//     printf("Actual error = %.10g\n", fabs(IntegralResult - M_PI/2));
//     printf("numOfEvals = %i\n\n", numOfEvals);
//     
//     double func6(double t){
//         numOfEvals++;
//         return exp(-t*t);
//     }
//     
//     numOfEvals = 0;
//     IntegrationError = 0;
//     IntegralResult = integrate(func6, -INFINITY, INFINITY, 0.001, 0.001, &IntegrationError);
//     printf("Our numerical integration of the integral of 1 / sqrt(1 + x^2) from -inf to inf gives the following results:\n");
//     printf("Analytical result = %.10g\n", M_PI/2);
//     printf("Numerical result = %.10g\n", IntegralResult);
//     printf("Estimated error = %.10g\n", IntegrationError);
//     printf("Actual error = %.10g\n", fabs(IntegralResult - M_PI/2));
//     printf("numOfEvals = %i\n\n", numOfEvals);
    
    return 0;
}
