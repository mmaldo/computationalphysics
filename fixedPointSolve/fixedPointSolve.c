/**
 * @file            fixedPointSolve.c
 * @description     This program approximates a solution to x^2 - 5 = 0
 *                  using fixed point iteration: x_{n+1} = x_n - A(x_n^2 - 5)
 *                  where A is a convergence parameter.
 *                  For convergence near sqrt(5), we need |1 - 2A*sqrt(5)| < 1
 *                  which gives 0 < A < 2/sqrt(5) ≈ 0.894
 * @course          Physics 177
 * @assignment      Homework 1 Problem 2
 * @date            4/7/2011
 * @author          Matthew Maldonado
 * @version         2.0
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

const int N = 10;        /* number of iterations */
const double A = 0.1;    /* convergence parameter - must be < 2/sqrt(5) ≈ 0.894 for convergence */

/**
 * Fixed point iteration function for solving x^2 - 5 = 0
 * Formula: x_{n+1} = x_n - A(x_n^2 - 5)
 * @param initialPoint is the current approximation
 * @param alpha is the convergence parameter
 * @pre   initialPoint is initialized by the calling program
 * @pre   alpha must satisfy convergence criteria: 0 < alpha < 2/sqrt(5)
 * @post  returns the next approximation in the iteration
 * @return the new value x_{n+1}
 */
double fixedPointIteration(double initialPoint, double alpha);

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
    
    printf("#Fixed Point Iteration for solving x^2 - 5 = 0\n");
    printf("#Using iteration: x_{n+1} = x_n - A(x_n^2 - 5) with A = %g\n", A);
    printf("#n\tX_{n-1}\t\tX_n\t\t|x_n - sqrt(5)|\n");
    
    for(i = 0; i < N; i++){
        previousPoint = newPoint;
        newPoint = fixedPointIteration(newPoint, A);
        convergenceError = fabs(newPoint - sqrt(5));
        printf("%d\t%.10g\t%.10g\t%.10g\n", i+1, previousPoint, newPoint, convergenceError);
    }
    
    printf("\n# Exact value: sqrt(5) = %.15f\n", sqrt(5));
    printf("# Final approximation: %.15f\n", newPoint);
    printf("# Final error: %.15e\n", convergenceError);
    printf("\n# Note: Fixed point iteration converges linearly (slower than Newton's method)\n");
    printf("# Convergence rate depends on |1 - 2A*sqrt(5)| = |1 - 2*%g*%.6f| = %.6f\n", 
           A, sqrt(5), fabs(1.0 - 2.0*A*sqrt(5)));
    
    return 0;
}

double fixedPointIteration(double initialPoint, double alpha){
    double x = initialPoint;
    
    /* Fixed point iteration: x_{n+1} = x - alpha*(x^2 - 5) */
    return (x - alpha*((x*x) - 5.0));
}