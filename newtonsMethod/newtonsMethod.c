#include <stdio.h>
#include <math.h>

#define EPSILON 1e-15

// Version 2.0
void newtonIteration(double initialPoint) {
    // Implementation of Newton's Method
    double x0 = initialPoint;
    double x1;
    int iterations = 0;
    double tolerance = 1e-6;
    
    // Loop until convergence
    do {
        x1 = x0 - fun(x0) / (derivative(x0) + EPSILON); // Improved division by zero check
        iterations++;
        
        printf("Iteration %d: x = %lf\n", iterations, x1);
        if (fabs(x1 - x0) < tolerance) {
            break;
        }
        x0 = x1;
    } while (1);

    printf("Converged to %lf after %d iterations\n", x1, iterations);
}

// Example function
double fun(double x) {
    // Your function implementation
}

double derivative(double x) {
    // Derivative implementation
}

int main() {
    double initialPoint = 0.0; // Example initial point
    newtonIteration(initialPoint);
    return 0;
}