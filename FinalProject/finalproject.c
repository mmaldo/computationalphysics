/*
 * 2D Laplace Equation Solver using Gauss-Seidel and SOR Methods
 * 
 * This program solves the 2D Laplace equation (∇²φ = 0) for electrostatic 
 * potential in a capacitor-like geometry:
 * - Outer rectangular boundary at potential V₀ = 1.0V
 * - Inner circular electrode at potential V₁ = 2.0V
 * 
 * The Gauss-Seidel method (alpha=1) or Successive Over-Relaxation (SOR, alpha>1)
 * iteratively updates each grid point to the weighted average of its four neighbors
 * until the solution converges to the specified tolerance.
 * 
 * The charge calculation uses Gauss's law: Q = ε₀∮E·dA, computed by integrating
 * the electric field (gradient of potential) normal to each electrode surface.
 * 
 * SOR parameter alpha:
 * - alpha = 1: Standard Gauss-Seidel
 * - 1 < alpha < 2: Over-relaxation for faster convergence
 * - Optimal alpha ≈ 2/(1 + π/N) for N×N grid
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
  double eps;
  double alpha; /* convergence parameter <=2 */
  
  int Nb;    /* how many boundaries */
  const double * Vext; /* potentials on the boundaries */
  double R;   /* circle radius */
  double Xc;  /* center X */
  double Yc;  /* center Y */
} params;

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
     


//What region are you stepping through
int in_region(int i, int j){
  double x=i*p.hx-p.Xc;
  double y=j*p.hx-p.Yc;
  if((i<=0)||(j<=0)||(i>=p.Nx)||(j>=p.Ny))
    return 0; /* electrode 0 */
  else if (x*x+y*y<p.R*p.R)
    return 1; /* electrode 1 */
  return (-1); /* inside the region */
}

void initialize(int argc, char * argv[]){
  int i,j,reg;
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
  if(p.hx<=0){
    printf("ERROR: step size hx=%g should be positive\n",p.hx);
    exit(-1); /*error */
  }
  if((p.alpha<0)||(p.alpha>=2)){
    printf("ERROR: convergence factor alpha=%g"
           " should satisfy 0<=alpha<2\n",p.alpha);
    exit(-1); /*error */
  }
  p.Nx=ceil(p.Xmax/p.hx); //Max steps base on total dived by steps for your matrix
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
  for(i=0;i<=p.Nx;i++)
    for(j=0;j<=p.Ny;j++){
      reg=in_region(i,j);
      if (reg>=0)
        A(i,j)=p.Vext[reg];  /* set boundary values */      
    }  
}

void finalize(void){
  free(v1);
  if(rng_generator != NULL)
    gsl_rng_free(rng_generator);  /* free the RNG */
}


void output_matrix(double *v1){ /* convenient function for debugging */
  int i,j;
  for(i=0;i<=p.Nx;i++){
    printf("#");
    for(j=0;j<=p.Ny;j++)
      printf(" %4.2f",A(i,j));
    printf("\n");
  }
  printf("#\n");
}

double do_charge(double *v1, int n){ /* charge on electrode n */
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

void output_charge( double *v){
  int i;  /* output charge on every electrode */
  printf("# charges: ");
  for(i=0;i<p.Nb;i++)
    printf(" %g",do_charge(v,i));
  printf("\n");  
}

void output_col(double *v1){ /* output for gnuplot's splot  */
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
