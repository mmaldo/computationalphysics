/**
 * @file            finalproject.c
 * @description     Solves 2D Laplace equation for electrostatic potential between
 *                  rectangular and circular electrodes using Gauss-Seidel/SOR method
 * @course          Computational Physics Final Project
 * @author          Matthew Maldonado
 * @version         2.0 - Fixed C/C++ mixing, memory leaks, and RNG inefficiency
 */

//#define DEBUG 1 
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define A(i,j) v1[(p.Nx+1)*(i)+(j)]

typedef struct {
  double hx; /* coordinate step */
  double Xmax; /* real-space dimensions */
  double Ymax;
  int Nx;    /* how many steps in x dimension */
  int Ny;    /* how many steps in y dimension */
  int Nt;    /* how many time steps */
  double eps; /* convergence threshold */
  double alpha; /* convergence parameter <=2 */
  
  int Nb;    /* how many boundaries */
  const double * Vext; /* potentials on the boundaries */
  double R;   /* circle radius */
  double Xc;  /* center X */
  double Yc;  /* center Y */
} params;

/* Global RNG variables - initialized once for efficiency */
gsl_rng * r;   
const gsl_rng_type * T;

const double VV[]={1.0, 2.0};
params p = { 0.1, 3.0, 3.0, 30, 30, 100, 1.0e-10, 1.0,
             2, VV ,
             1., 1.5, 1.5}; 

double * v1; /* pointer to the potential array */

/* Global RNG for efficient random number generation */
gsl_rng * rng_generator = NULL;

/* Initialize the random number generator once at program start */
void init_rng(void){
  const gsl_rng_type * T;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  rng_generator = gsl_rng_alloc(T);
}

/* Generate a random number using the global RNG */
double randomNumber(void){
  return gsl_rng_uniform(rng_generator);
}
     


/**
 * Determine which region a grid point belongs to.
 * @param i Grid index in x direction
 * @param j Grid index in y direction
 * @return 0 for outer electrode, 1 for inner circular electrode, -1 for interior
 */
int in_region(int i, int j){
  double x=i*p.hx-p.Xc;
  double y=j*p.hx-p.Yc;
  if((i<=0)||(j<=0)||(i>=p.Nx)||(j>=p.Ny))
    return 0; /* electrode 0 */
  else if (x*x+y*y<p.R*p.R)
    return 1; /* electrode 1 */
  return (-1); /* inside the region */
}

/**
 * Initialize simulation parameters and allocate memory.
 * Sets up the global RNG, parses command line arguments, and initializes the potential grid.
 * @param argc Number of command line arguments
 * @param argv Array of command line argument strings
 */
void initialize(int argc, char * argv[]){
  int i,j,reg;
  
  /* Initialize GSL random number generator once for efficiency */
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
  
  /* Parse command line arguments */
  for(i=1;i<argc;i++){
    if(sscanf(argv[i],"hx=%lg",& p.hx)==1){			   
      printf("# read hx=%g\n",p.hx);				   
    }
    else if (sscanf(argv[i],"alpha=%lg",& p.alpha)==1){      
      printf("# read alpha=%g\n",p.alpha);				   
    }
    else {
      fprintf(stderr,"error in argv[%d]=%s\n" 
	      "\tusage: %s [hx=#] [alpha=#]\n"
	      "\twhere 0<=alpha<2 is the convergence factor\n"
	      ,i,argv[i],argv[0]);
      exit(-1);
    }
  }
  
  /* Validate parameters */
  if(p.hx<=0){
    printf("ERROR: step size hx=%g should be positive\n",p.hx);
    exit(-1); /*error */
  }
  if((p.alpha<0)||(p.alpha>=2)){
    printf("ERROR: convergence factor alpha=%g"
           " should satisfy 0<=alpha<2\n",p.alpha);
    exit(-1); /*error */
  }
  
  /* Calculate grid dimensions */
  p.Nx=ceil(p.Xmax/p.hx);
  p.Ny=ceil(p.Ymax/p.hx);
  
  printf("# %s: %s relaxation for Laplace equation in 2D\n"
         "# params: hx=%g Nx=%d Ny=%d Nt=%d eps=%g alpha=%g",
         argv[0], p.alpha==1?"Gauss-Seidel":"SOR",
         p.hx, p.Nx, p.Ny, p.Nt, p.eps, p.alpha);
  
  for(i=0;i<p.Nb;i++)
    printf(" Vext[%d]=%g",i,p.Vext[i]);
  printf("\n");
    v1=calloc((p.Nx+1)*(p.Ny+1),sizeof(double)); 
  /* allocate mem for the potential */
  if(v1==NULL){
    fprintf(stderr,"failed to allocate mem for v1\n");
    exit(-1);
  }
  
  /* Set boundary values */
  for(i=0;i<=p.Nx;i++)
    for(j=0;j<=p.Ny;j++){
      reg=in_region(i,j);
      if (reg>=0)
        A(i,j)=p.Vext[reg];  /* set boundary values */      
    }  
}

/**
 * Clean up and deallocate memory.
 * Frees the potential grid and the global RNG.
 */
void finalize(void){
  free(v1);
  if(rng_generator != NULL)
    gsl_rng_free(rng_generator);  /* free the RNG */
}


/**
 * Output the potential matrix for debugging.
 * Prints the entire potential grid in a human-readable format.
 * @param v1 Pointer to the potential array
 */
void output_matrix(double *v1){
  int i,j;
  for(i=0;i<=p.Nx;i++){
    printf("#");
    for(j=0;j<=p.Ny;j++)
      printf(" %4.2f",A(i,j));
    printf("\n");
  }
  printf("#\n");
}

/**
 * Calculate the charge on a specific electrode.
 * Uses Gauss's law by summing the electric field at the electrode boundary.
 * @param v1 Pointer to the potential array
 * @param n Electrode number
 * @return Charge on electrode n
 */
double do_charge(double *v1, int n){
  int i,j;
  double charg=0;
  const double pi=M_PI;
  for(i=0;i<p.Nx;i++)
    for(j=0;j<p.Ny;j++)
      if(in_region(i,j)==n){ /* bonds beginning on n; field added */
        if(in_region(i+1,j)==-1)
          charg+=A(i+1,j)-A(i,j);
        if(in_region(i,j+1)==-1)
          charg+=A(i,j+1)-A(i,j);
      }
      else if(in_region(i,j)==-1){ /* bonds ending on n; field subtracted */
        if(in_region(i+1,j)==n)
          charg-=A(i+1,j)-A(i,j);
        if(in_region(i,j+1)==n)
          charg-=A(i,j+1)-A(i,j);
      }
  return charg/(2.0*pi);
}

/**
 * Output the charges on all electrodes.
 * @param v Pointer to the potential array
 */
void output_charge( double *v){
  int i;
  printf("# charges: ");
  for(i=0;i<p.Nb;i++)
    printf(" %g",do_charge(v,i));
  printf("\n");  
}

/**
 * Output the potential grid in column format for gnuplot's splot.
 * Outputs x, y, and potential values in a format suitable for 3D plotting.
 * @param v1 Pointer to the potential array
 */
void output_col(double *v1){
  int i,j;
  for(i=0;i<=p.Nx;i++){
    for(j=0;j<=p.Ny;j++)
      printf("%g \t %g \t %g\n",p.hx*i,p.hx*j,A(i,j));
    printf("\n"); 
  }
  printf("\n");
}

/* Does Gauss-Seidel method to solve differential equation
 * Returns the number of iterations performed (positive if converged, negative if max iterations reached)
 */
int do_gs(int *iterations_out){
  int j,i,k;
  double error, diff, val;
  for(k=0;k<p.Nt;k++){
    error=0;           /* maximum difference from previous iteration */
    for(i=1;i<p.Nx;i++){
      for(j=1;j<p.Ny;j++){ /* update inside values */
        if(in_region(i,j)<0){
          val=(1-p.alpha)*A(i,j)+
            p.alpha*0.25*(A(i-1,j)+A(i+1,j)+A(i,j-1)+A(i,j+1));
          diff=fabs(val-A(i,j)); /* difference */
          A(i,j)=val;                  
          if(error<diff)
            error=diff;  /* update the error */
        }
      }
    }
#if DEBUG
    printf("# this is v1, k=%d\n",k);
    output_matrix(v1);  
#endif /* DEBUG */
    if(error<p.eps){
      *iterations_out = k+1;  /* converged */
      return 1;
    }
  }
  *iterations_out = p.Nt;  /* max iterations reached */
  return 0;
}
/* The main function of the program. This is where the program actually runs. */
int main(int argc, char * argv[]){
  int iter=0;
  int gs_iterations=0;
  int converged;
  
  init_rng();  /* Initialize the random number generator */
  initialize(argc, argv);

  do{
    output_charge(v1);
#if DEBUG
    output_matrix(v1);
#endif /* DEBUG */
    
    output_col(v1);
    fflush(stdout);
    iter++;    
  }
  while ((converged = !do_gs(&gs_iterations)) && iter<p.Nt); 

  output_col(v1); /* final output */
  output_charge(v1);
    
  finalize();
  
  fprintf(stderr,"hx=%5g alpha=%5g iterations=%d %s %g seconds\n",
          p.hx, p.alpha, gs_iterations,
          converged ? "(max reached)" : "(converged)",
          (double)(clock())/CLOCKS_PER_SEC);
  
  return 0;
}
