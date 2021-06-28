My studentnumber is 201707128, which means that I've done examproject number 6.

----------------6----------------

Adaptive integration of complex-valued functions

Implement an adaptive integrator which calculates the integral of a complex-valued function f(z) of a complex variable z along a straight line between two points in the complex plane.

----------------------------------------------------------------------------------------------------------

In this examproject I've created a routine that can take a complex function that takes a complex variable
as an argument and integrate the function along a straight line in the complex plane.

I've done this by updating the algorithm I've already made in homework 6 about adaptive integration.
I updated the routine by using the headerfile complex.h to take care of all the calculations with
complex numbers. Changing the variables in the code I created for homework 6 from 'double' to
'double complex' and using the complex norm to calculate the error, with the trapezium rule as the higher order
rule and the rectangle rule as the lower order rule, as this was what I did in homework 6.


