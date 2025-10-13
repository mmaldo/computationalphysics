/**
 * @file            logmap.c
 * @description     Logistic map bifurcation diagram with proper thermalization
 *                  x_{n+1} = r * x_n * (1 - x_n)
 * @author          Matthew Maldonado
 * @version         2.0
 */

#include <stdio.h>

int main() {
    double r, x_n;
    int i, j;
    const int THERM_STEPS = 100;  /* thermalization iterations */
    const int DATA_STEPS = 100;   /* data collection iterations */
    
    printf("#x_n\tr\n");  /* Fixed: added \n */
    
    for (r = 2.5; r <= 4.0; r += 0.001) {
        x_n = 0.25;  /* Reset initial condition for each r */
        
        /* Thermalization: let transients die out */
        for (i = 0; i < THERM_STEPS; i++) {
            x_n = r * x_n * (1.0 - x_n);
        }
        
        /* Data collection: record points on attractor */
        for (j = 0; j < DATA_STEPS; j++) {
            x_n = r * x_n * (1.0 - x_n);
            printf("%g\t%g\n", x_n, r);  /* Fixed: removed extra space */
        }
    }
    
    return 0;
}