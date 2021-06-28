#include<math.h>
#include<assert.h>

double adapt24(double f(double), double a, double b, double delta, double eps, double f2, double f3, int numOfEvals, double* err){
    
    assert(numOfEvals < 100000);
    
    double f1 = f(a + (b - a) / 6);
    double f4 = f(a + 5 * (b - a) / 6);
    
    double Q = (2 * f1 + f2 + f3 + 2 * f4) / 6 * (b - a);
    double q = (f1 + f4 + f2 + f3) / 4 * (b - a);
    
    double tol = delta + eps*fabs(Q);
    double error = fabs(Q - q);
    
    if(error < tol){
        *err += error;
        return Q;
    }else{
        double new_Q1 = adapt24(f, a, (a + b) / 2, delta/sqrt(2.), eps, f1, f2, numOfEvals + 1, err);
        double new_Q2 = adapt24(f, (a + b) / 2, b, delta/sqrt(2.), eps, f3, f4, numOfEvals + 1, err);
        return new_Q1 + new_Q2;
    }
}

double integrate(double f(double), double a, double b, double delta, double eps, double* err){
    int numOfEvals = 0;
    double f2 = f(a + 2 * (b - a) / 6);
    double f3 = f(a + 4 * (b - a) / 6);
    return adapt24 (f, a, b, delta, eps, f2, f3, numOfEvals, err);
    
//     if(isinf(-a)){          // Check if we integrate from -infinity
//         if(isinf(b)){       // Check if we integrate to infinity while integrating from -infinity
// 
//             double IntervalTransformation (double t){
//                 return f(t / (1 - t * t)) * (1 + t * t) / pow(1 - t * t, 2);
//             }
//             
//             return adapt24 (IntervalTransformation, -1, 1, delta, eps, f2, f3, numOfEvals, err);
//             
//         }else{              // Check if we integrate to a finite number while integrating from -infinity
//             
//             double IntervalTransformation (double t){
//                 return f(b - (1 - t) / t) / (t * t);
//             }
//             
//             return adapt24 (IntervalTransformation, 0, 1, delta, eps, f2, f3, numOfEvals, err);
//         }
//     }else if(isinf(b)){      // Check if we integrate to infinity while integrating from a finite number
//         
//         double IntervalTransformation (double t){
//             return f(a + t / (1 - t)) / (1 - t * t);
//         }
//         
//         return adapt24 (IntervalTransformation, 0, 1, delta, eps, f2, f3, numOfEvals, err);
//     }else{
//         
//     }
}

double OpenQuadCC_integrate(double f(double), double a, double b, double delta, double eps, double* err){
    int numOfEvals = 0;
    double f2 = f(a + 2 * (b - a) / 6);
    double f3 = f(a + 4 * (b - a) / 6);
    
    double IntervalTransformation (double theta){
        return f((a + b) / 2 + (b - a) * cos(theta) / 2 )*(b - a) * sin(theta) / 2;
    }

    return adapt24(IntervalTransformation, 0, M_PI, delta, eps, f2, f3, numOfEvals, err);
}
