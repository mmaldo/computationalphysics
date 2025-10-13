/**
 * @file            test.c
 * @description     Driven damped nonlinear pendulum using RK4 method
 *                  Equations: theta' = v
 *                            v' = -v/Q - sin(theta) + A*cos(w*t)
 * @author          Matthew Maldonado
 * @version         2.0
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct {
    double tmax;    /* total evolution time */
    double w;       /* driving frequency */
    double A;       /* driving amplitude */
    double q;       /* quality factor (damping) */
    double x0;      /* initial angle (theta) */
    double v0;      /* initial angular velocity */
    int n;          /* number of steps */
    double h;       /* time step size */
    double energ0;  /* initial energy */
} params;

/* FIXED: 2.0/3.0 instead of 2/3.0 to avoid integer division */
params p = {100.0, 2.0/3.0, 1.5, 1.372, 0.0, 0.0, 100, 0.0, 0.0};
params p1 = {100.0, 2.0/3.0, 1.5, 1.372, 1.0e-6, 1.0e-6, 100, 0.0, 0.0};

/**
 * Calculate acceleration for driven damped pendulum
 */
double acceleration(double v, double x, double t, params par) {
    return (-v/par.q) - sin(x) + par.A*cos(par.w*t);
}

/**
 * Calculate energy for the driven damped pendulum
 */
double energy(double x, double v) {
    return 0.5*v*v + (1.0 - cos(x));
}

/**
 * Normalize angle to [-π, π]
 */
double normalize_angle(double angle) {
    while (angle > M_PI) angle -= 2.0*M_PI;
    while (angle < -M_PI) angle += 2.0*M_PI;
    return angle;
}

/**
 * 4th order Runge-Kutta integrator
 */
void rk4() {
    double t = 0.0;
    double x = p.x0;
    double v = p.v0;
    double x1 = p1.x0;
    double v1 = p1.v0;
    
    p.h = p.tmax / (double)p.n;
    
    double xk1, xk2, xk3, xk4, vk1, vk2, vk3, vk4;
    double x1k1, x1k2, x1k3, x1k4, v1k1, v1k2, v1k3, v1k4;
    double x_temp, v_temp, x1_temp, v1_temp;
    
    p.energ0 = energy(p.x0, p.v0);
    p1.energ0 = energy(p1.x0, p1.v0);
    
    printf("#Driven Damped Pendulum - RK4 Integration\n");
    printf("#Parameters: A=%.3f, w=%.3f, Q=%.3f, dt=%.6f\n", p.A, p.w, p.q, p.h);
    printf("#time\ttheta1\tv1\tE1\ttheta2\tv2\tE2\t|Δtheta|\t|Δv|\tlog10|Δtheta|\tlog10|Δv|\n");
    
    double diff_x = fabs(x - x1);
    double diff_v = fabs(v - v1);
    double energ = energy(x, v);
    double energ1 = energy(x1, v1);
    printf("%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
           t, x, v, energ, x1, v1, energ1, diff_x, diff_v,
           log10(diff_x), log10(diff_v));
    
    for (int i = 0; i < p.n; i++) {
        /* RK4 for first trajectory - PROPERLY implemented */
        xk1 = v;
        vk1 = acceleration(v, x, t, p);
        
        x_temp = x + 0.5*p.h*xk1;
        v_temp = v + 0.5*p.h*vk1;
        xk2 = v_temp;
        vk2 = acceleration(v_temp, x_temp, t + 0.5*p.h, p);
        
        x_temp = x + 0.5*p.h*xk2;
        v_temp = v + 0.5*p.h*vk2;
        xk3 = v_temp;
        vk3 = acceleration(v_temp, x_temp, t + 0.5*p.h, p);
        
        x_temp = x + p.h*xk3;
        v_temp = v + p.h*vk3;
        xk4 = v_temp;
        vk4 = acceleration(v_temp, x_temp, t + p.h, p);
        
        /* RK4 for second trajectory */
        x1k1 = v1;
        v1k1 = acceleration(v1, x1, t, p1);
        
        x1_temp = x1 + 0.5*p.h*x1k1;
        v1_temp = v1 + 0.5*p.h*v1k1;
        x1k2 = v1_temp;
        v1k2 = acceleration(v1_temp, x1_temp, t + 0.5*p.h, p1);
        
        x1_temp = x1 + 0.5*p.h*x1k2;
        v1_temp = v1 + 0.5*p.h*v1k2;
        x1k3 = v1_temp;
        v1k3 = acceleration(v1_temp, x1_temp, t + 0.5*p.h, p1);
        
        x1_temp = x1 + p.h*x1k3;
        v1_temp = v1 + p.h*v1k3;
        x1k4 = v1_temp;
        v1k4 = acceleration(v1_temp, x1_temp, t + p.h, p1);
        
        /* Update positions and velocities */
        x = x + (p.h/6.0)*(xk1 + 2.0*xk2 + 2.0*xk3 + xk4);
        v = v + (p.h/6.0)*(vk1 + 2.0*vk2 + 2.0*vk3 + vk4);
        
        x1 = x1 + (p.h/6.0)*(x1k1 + 2.0*x1k2 + 2.0*x1k3 + x1k4);
        v1 = v1 + (p.h/6.0)*(v1k1 + 2.0*v1k2 + 2.0*v1k3 + v1k4);
        
        /* Normalize angles */
        x = normalize_angle(x);
        x1 = normalize_angle(x1);
        
        t += p.h;
        
        energ = energy(x, v);
        energ1 = energy(x1, v1);
        diff_x = fabs(x - x1);
        diff_v = fabs(v - v1);
        
        printf("%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
               t, x, v, energ, x1, v1, energ1, diff_x, diff_v,
               log10(diff_x), log10(diff_v));
    }
}

int main() {
    rk4();
    return 0;
}
