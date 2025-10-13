/**
 * @file            newtonsMethod.c
 * @description     This program approximates a solution to x^2 - 5 = 0
 *                  using Newton's Method: x_{n+1} = x_n - f(x_n)/f'(x_n)
 *                  For f(x) = x^2 - 5, we get: x_{n+1} = x_n - (x_n^2 - 5)/(2*x_n)
 * @course          Physics 177
 * @assignment      Homework 1 Problem 2
 * @date            4/7/2011
 * @author          Matthew Maldonado
 * @version         2.0
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

const int N = 10;   /* number of iterations */

/**
 * Newton's Method iteration function for solving x^2 - 5 = 0
 * Formula: x_{n+1} = x_n - f(x_n)/f'(x_n) = x_n - (x_n^2 - 5)/(2*x_n)
 * @param initialPoint is the current approximation
 * @pre   initialPoint is initialized by the calling program
 * @post  returns the next approximation in the iteration
 * @return the new value x_{n+1}
 */
double newtonIteration(double initialPoint);

/**
 * @return zero if application executes successfully
 */
int main()
{
    double initialPoint = 1.0;      /* initial guess */
    double newPoint = initialPoint; /* x_n */
    double previousPoint;           /* x_{n-1} */
    double convergenceError;        /* |x_n - sqrt(5)| */
    int i;                          /* iteration counter */
    
    printf("#Newton's Method for solving x^2 - 5 = 0\n");
    printf("#n\tX_{n-1}\t\tX_n\t\t|x_n - sqrt(5)|\n");
    
    for(i = 0; i < N; i++){
        previousPoint = newPoint;
        newPoint = newtonIteration(newPoint);
        convergenceError = fabs(newPoint - sqrt(5));
        printf("%d\t%.10g\t%.10g\t%.10g\n", i+1, previousPoint, newPoint, convergenceError);
    }
    
    printf("\n# Exact value: sqrt(5) = %.15f\n", sqrt(5));
    printf("# Final approximation: %.15f\n", newPoint);
    printf("# Final error: %.15e\n", convergenceError);
    
    return 0;
}

double newtonIteration(double initialPoint){
    double x = initialPoint;
    
    /* Check if we're too close to zero (division by zero protection) */
    if (fabs(x) < 1e-15) {
        fprintf(stderr, "Error: x value too close to zero, cannot compute Newton iteration.\n");
        exit(1);
    }
    
    /* Newton's method: x_{n+1} = x - f(x)/f'(x) = x - (x^2-5)/(2x) */
    return (x - ((x*x) - 5.0)/(2.0*x));
}