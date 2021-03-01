#include<stdio.h>
#include<limits.h>
#include<float.h>
#include "functions.h"

int main()
{
    printf("1.i\n");
    int i = 1;
    while(i + 1 > i){i++;}
    printf("While loop\n");
    printf("My biggest int = %i \n", i);
    printf("INT_MAX = %i \n", INT_MAX);
    
    i = 1;
    for (; i < i + 1; i++){}
    printf("For loop\n");
    printf("My biggest int = %i \n", i);
    printf("INT_MAX = %i \n", INT_MAX);

    i = 1;
    do{i++;} while (i + 1 > i);
    printf("Do-While loop\n");
    printf("My biggest int = %i \n", i);
    printf("INT_MAX = %i \n", INT_MAX);
    
    printf("\n 1.ii\n");
    i = -1;
    while(i - 1 < i){i--;}
    printf("While loop\n");
    printf("My smallest int = %i \n", i);
    printf("INT_MIN = %i \n", INT_MIN);
    
    i = -1;
    for (; i > i - 1; i--){}
    printf("For loop\n");
    printf("My smallest int = %i \n", i);
    printf("INT_MIN = %i \n", INT_MIN);

    i = -1;
    do{i--;} while (i - 1 < i);
    printf("Do-While loop\n");
    printf("My smallest int = %i \n", i);
    printf("INT_MIN = %i \n", INT_MIN);
    
    printf("\n 1.iii\n");
    
    float my_float_epsilon;
    float x=1; while( (float) 1+x != 1){my_float_epsilon=x; x/=2;}
    printf("While loop\n");
    printf("My float epsilon = %.10f \n", my_float_epsilon);
    printf("FLT_EPSILON = %.10f \n", FLT_EPSILON);
    
    double my_double_epsilon;
    double y=1; while( (double) 1+y != 1){my_double_epsilon=y; y/=2;}
    printf("My double epsilon = %.10g \n", my_double_epsilon);
    printf("DBL_EPSILON = %.10g \n", DBL_EPSILON);
    
    long double my_long_double_epsilon;
    long double z=1.0L; while( (long double) 1+z != 1){my_long_double_epsilon=z; z/=2;}
    printf("My long double epsilon = %.10Lg \n", my_long_double_epsilon);
    printf("LDBL_EPSILON = %.10Lg \n", LDBL_EPSILON);
    
    printf("\nFor loop\n");
    for(x=1; 1+x != 1; x/=2){my_float_epsilon = x;}
    printf("My float epsilon = %.10f \n", my_float_epsilon);
    printf("FLT_EPSILON = %.10f \n", FLT_EPSILON);
    
    for(y=1; 1+y != 1; y/=2){my_double_epsilon = y;}
    printf("My double epsilon = %.10g \n", my_double_epsilon);
    printf("DBL_EPSILON = %.10g \n", DBL_EPSILON);
    
    for(z=1; 1+z != 1; z/=2){my_long_double_epsilon = z;}
    printf("My long double epsilon = %.10Lg \n", my_long_double_epsilon);
    printf("LDBL_EPSILON = %.10Lg \n", LDBL_EPSILON);
    
    printf("\nDo-while loop\n");
    do{my_float_epsilon = x*2; x /= 2;} while (1+x != 1);
    printf("My float epsilon = %.10f \n", my_float_epsilon);
    printf("FLT_EPSILON = %.10f \n", FLT_EPSILON);
    
    do{my_double_epsilon = y*2; y /= 2;} while (1+y != 1);
    printf("My double epsilon = %.10g \n", my_double_epsilon);
    printf("DBL_EPSILON = %.10g \n", DBL_EPSILON);
    
    do{my_long_double_epsilon = z*2; z /= 2;} while (1+z != 1);
    printf("My long double epsilon = %.10Lg \n", my_long_double_epsilon);
    printf("LDBL_EPSILON = %.10Lg \n", LDBL_EPSILON);
    
    printf("\n 2.i\n");
    
    int max = INT_MAX/2;
    float sum_up_float = 0.0;
    int k = 1;
    while(k < max+1){ sum_up_float = sum_up_float + 1.0f/k; k++;}
    
    float sum_down_float = 0.0;
    k = max;
    while(k > 0){ sum_down_float = sum_down_float + 1.0f/k; k--;}

    
    printf("sum_up_float = %.10f \n", sum_up_float);
    printf("sum_down_float = %.10f \n", sum_down_float);
     
    printf("\n 2.ii\n");
    printf("The difference between sum_up_float and sum_down_float is due to rounding errors. When we calculate sum_up_float we start with a large number and add ever so smaller numbers to it. At some point the numbers we add get smaller then the precision of the number we started with and the sum will terminate. On the other hand when we calculate sum_down_float, we start with a smaller number and add increasingly larger numbers to the sum. Doing this, the sum will not terminate because no term is smaller than the precision of the sum. This results in a smaller rounding error and thus we get a larger number as we see.");
    
    printf("\n 2.iii\n");
    printf("No I don't think the sum converges. Only the first sum where the last terms vanish due to rounding.");
    
    printf("\n 2.iv\n");
    max = INT_MAX/2;
    double sum_up_double = 0.0;
    k = 1;
    while(k < max+1){ sum_up_double = sum_up_double + 1.0f/k; k++;}
    
    double sum_down_double = 0.0;
    k = max;
    while(k > 0){ sum_down_double = sum_down_double + 1.0f/k; k--;}
    
    printf("sum_up_double = %.10g \n", sum_up_double);
    printf("sum_down_double = %.10g \n", sum_down_double);
    printf("Here we don't get the rounding errors we did with float, since double has a higher precision than float.");
    
    printf("\n 3.i\n");
    printf("Function implemented\n");
}
