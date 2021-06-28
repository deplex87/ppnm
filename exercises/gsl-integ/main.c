#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>

double f (double x, void * params) {
  double f = log(x) / sqrt(x);
  return f;
}

int main ()
{
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

    double result, error;

    gsl_function F;
    F.function = &f;

  gsl_integration_qags (&F, 0, 1, 0, 1e-7, 1000, w, &result, &error);

  printf("Integration of the function log(x)/sqrt(x) from 0 to 1:\n");
  printf("Exact result     = % .20f\n", -4.0);
  printf("Numerical result = % .20f\n", result);
  printf("Estimated error  = % .20f\n", error);
  printf("Actual error     = % .20f\n", result + 4);

  gsl_integration_workspace_free (w);

  return 0;
}
