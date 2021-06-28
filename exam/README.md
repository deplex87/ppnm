My studentnumber is 201707128, which means that I've done examproject number 6.

-----------------------------------------------------------------6---------------------------------------------------------------

Adaptive integration of complex-valued functions

Implement an adaptive integrator which calculates the integral of a complex-valued function f(z) of a complex variable z along a straight line between two points in the complex plane.

----------------------------------------------------------------------------------------------------------------------------------

In this examproject I've created a routine that can take a complex function that takes a complex variable
as an argument and integrate the function along a straight line in the complex plane.

I've done this by updating the algorithm I've already made in homework 6 about adaptive integration.
I updated the routine by using the headerfile complex.h to take care of all the calculations with
complex numbers. Changing the variables in the code I created for homework 6 from 'double' to
'double complex' and using the complex norm to calculate the error, with the trapezium rule as the higher order
rule and the rectangle rule as the lower order rule, as this was what I did in homework 6.

I also implemented open quadrature with Clenshaw–Curtis variable transformation like in homework 6 for complex
functions. In the file output.txt are both results with recursive adaptive integration and variable transformation
shown. Furthermore, I've plotted one of the integrations, which also shows the first division of the straight line
in the algorithm.

----------------------------------------------------------Files in folder---------------------------------------------------------
IntegralPlot.gpi    :   This file holds all the pyxplot code that generates the plot IntegralPlot.png

IntegralPlot.png    :   This image shows the line which i integrate along for one of the functions that
                        I integrate with my algorithm
                        
line.txt            :   Contains points for the line on the plot

main.c              :   In this file we test the integration algorithms with different complex integrals

out.txt             :   Contains all the results of the different integrals integrated with the algorithm

plane.txt           :   Contains all the points for the complex plane on the plot

points.txt          :   Contains the points for the star points on the plot

RAintegrator.c      :   Contains all the functions that calculate the integrals



