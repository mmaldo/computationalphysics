/**
 * ============================================================================
 * 2D Laplace Equation Solver with Thermal Effects
 * ============================================================================
 * 
 * PHYSICS PROBLEM:
 * This program solves a modified 2D Laplace equation with thermal effects
 * using the Gauss-Seidel/SOR (Successive Over-Relaxation) iterative method.
 * 
 * STANDARD LAPLACE EQUATION:
 *   ∇²φ = 0
 *   Discretized: A(i,j) = 0.25*(A(i+1,j) + A(i-1,j) + A(i,j+1) + A(i,j-1))
 * 
 * THERMAL MODIFICATION:
 * The thermal factor exp(-V/T) models temperature-dependent relaxation:
 *   A(i,j) = 0.25 * Σ(neighbors) A(neighbor) * exp(-(V(i,j)+V(neighbor))/(2*T))
 * 
 * where:
 *   - T is the temperature parameter (controls thermal fluctuations)
 *   - V(i,j) = A(i,j) + F*i is the total potential including external field
 *   - F is the external field strength (creates linear potential gradient)
 * 
 * PHYSICAL INTERPRETATION:
 *   - At T → 0: reduces to standard Laplace equation (no thermal effects)
 *   - At T > 0: thermal fluctuations affect the relaxation dynamics
 *   - The Boltzmann-like factor exp(-V/T) weights updates by energy differences
 *   - Can model thermal disorder, random electrode configurations, or
 *     temperature-dependent conductivity in electrostatic systems
 * 
 * GAUSSIAN CONVOLUTION KERNEL:
 * The program uses convolution with a Gaussian kernel to create spatially
 * correlated random initial conditions:
 *   - Random field B(m,n) is generated uniformly
 *   - Gaussian kernel C(k,l) = exp(-(x² + y²)/(2*L²)) with correlation length L
 *   - Convolution A = B ⊗ C creates smooth, correlated disorder
 * 
 * KEY PARAMETERS:
 *   - hx: Grid spacing in x and y directions
 *   - Nx, Ny: Number of grid points in x and y
 *   - Nt: Maximum number of iterations
 *   - eps: Convergence criterion (max change between iterations)
 *   - alpha: SOR parameter (1 = Gauss-Seidel, >1 = over-relaxation, must be <2)
 *   - T: Temperature parameter for thermal effects
 *   - F: External field strength
 *   - L: Correlation length for Gaussian kernel (set to hx*10)
 * 
 * APPLICATIONS:
 *   - Classical diffusion with strong correlated disorder
 *   - Thermal fluctuations in electrostatic systems
 *   - Random electrode potentials with spatial correlations
 *   - Temperature-dependent conductivity problems
 * ============================================================================
 */

//#define DEBUG 1 
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
/* Grid dimensions for array allocation */
#define GRID_SIZE_A 61  /* Size for potential array A and auxiliary array V */
#define GRID_SIZE_C 51  /* Size for Gaussian kernel array C (must be odd for symmetry) */
#define KERNEL_SIZE 31  /* Size of convolution kernel (smaller than GRID_SIZE_C) */

/* Macro accessors for 2D arrays stored as 1D */
#define A(i,j) v1[(GRID_SIZE_A)*(i)+(j)]  /* Potential field φ(i,j) */
#define B(m,n) v2[(p.Nx)*(m)+(n)]         /* Random number field for initial conditions */
#define C(k,l) v3[(GRID_SIZE_C)*(k)+(l)]  /* Gaussian convolution kernel */
#define V(i,j) v4[(GRID_SIZE_A)*(i)+(j)]  /* Total potential V = A + F*i (includes external field) */

/**
 * Parameters structure for the simulation
 */
typedef struct {
  double hx;    /* Grid spacing (coordinate step) */
  double Xmax;  /* Domain size in x direction */
  double Ymax;  /* Domain size in y direction */
  int Nx;       /* Number of grid points in x dimension */
  int Ny;       /* Number of grid points in y dimension */
  int Nt;       /* Maximum number of iterations */
  int l;        /* Computed grid dimension for kernel */
  int k;        /* Computed grid dimension for kernel */
  double eps;   /* Convergence criterion (tolerance) */
  double alpha; /* SOR parameter: 1=Gauss-Seidel, 1<alpha<2 for over-relaxation */
  
  int Nb;       /* Number of boundary electrodes */
  const double * Vext; /* Array of potentials on the boundaries */
  double R;     /* Inner circular electrode radius */
  double Xc;    /* Center X coordinate of circular electrode */
  double Yc;    /* Center Y coordinate of circular electrode */
  double T;     /* Temperature parameter for thermal effects (Boltzmann factor) */
  double F;     /* External field strength (creates linear potential gradient) */
} params;

/* GSL random number generator */
gsl_rng * r;                    /* RNG instance */
const gsl_rng_type * T;         /* RNG type */

/* Boundary electrode potentials */
const double VV[]={1.0, 2.0};   /* Outer and inner electrode potentials */

/* Default simulation parameters */
params p = { 0.1, 6.0,6.0, 30, 30, 30,30 ,100, 1.0e-10, 1.0,
             2, VV ,
             1., 1.5, 1.5,2.0,2.0}; 

/* Dynamic arrays for field storage */
double * v1; /* Potential field A(i,j) */
double * v2; /* Random number field B(m,n) for initial conditions */
double * v3; /* Gaussian kernel C(k,l) for convolution */
double * v4; /* Total potential V(i,j) = A(i,j) + F*i */

/**
 * Generate uniform random number in [0,1)
 * @return Random number from GSL uniform distribution
 */
double randomNumber(){
  double u;
  u = gsl_rng_uniform (r);
  return u;
}    

/**
 * Determine which region a grid point belongs to
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
 * Initialize simulation parameters and allocate memory
 * - Parses command line arguments
 * - Allocates memory for all arrays
 * - Generates random field B(m,n)
 * - Creates Gaussian convolution kernel C(k,l)
 * - Sets boundary conditions on electrodes
 * @param argc Number of command line arguments
 * @param argv Array of command line argument strings
 */
void initialize(int argc, char * argv[]){
  int i,k,l,m,n;
  int reg;  /* Region identifier for boundary conditions */
  p.l = ceil((p.Xmax+1)/p.hx);
 p.k = ceil((p.Ymax+1)/p.hx);
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
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
  p.Nx=ceil(p.Xmax/p.hx);
  p.Ny=ceil(p.Ymax/p.hx);
  
  printf("# %s: %s relaxation for Laplace equation in 2D\n"
         "# params: hx=%g Nx=%d Ny=%d Nt=%d eps=%g alpha=%g",
         argv[0], p.alpha==1?"Gauss-Seidel":"SOR",
         p.hx, p.Nx, p.Ny, p.Nt, p.eps, p.alpha);
  
  /* Allocate memory for random number field B(m,n) */
  v2=calloc((p.Nx+1)*(p.Ny),sizeof(double)); 
  if(v2==NULL){
    fprintf(stderr,"failed to allocate mem for v2 (random field)\n");
    exit(-1);
  }
    
  /* Allocate memory for Gaussian convolution kernel C(k,l) */
  v3=calloc((GRID_SIZE_C)*(GRID_SIZE_C),sizeof(double)); 
  if(v3==NULL){
    fprintf(stderr,"failed to allocate mem for v3 (Gaussian kernel)\n");
    exit(-1);
  }
  
  /* Allocate memory for total potential V(i,j) = A(i,j) + F*i */
  v4=calloc((GRID_SIZE_A)*(GRID_SIZE_A),sizeof(double)); 
  if(v4==NULL){
    fprintf(stderr,"failed to allocate mem for v4 (total potential)\n");
    exit(-1);
  }
  
  for(i=0;i<p.Nb;i++)
    printf(" Vext[%d]=%g",i,p.Vext[i]);
  printf("\n");
  
  /* Allocate memory for potential field A(i,j) */
  v1=calloc((GRID_SIZE_A)*(GRID_SIZE_A),sizeof(double)); 
  if(v1==NULL){
    fprintf(stderr,"failed to allocate mem for v1 (potential field)\n");
    exit(-1);
  }
  
  /* Generate uniform random field B(m,n) for initial conditions */
  for(m=0;m <p.Nx;m++){
    for(n = 0;n<p.Ny;n++){
      B(m,n)=randomNumber();
    }
  }
  
  /* Create Gaussian convolution kernel C(k,l) 
   * This creates spatial correlations in the initial random field
   * Kernel: C(k,l) = exp(-(x² + y²)/(2*L²))
   * where L is the correlation length
   */
  double L =p.hx*10;  /* Correlation length for Gaussian kernel */
  double x;
  double y;
  
  for(k=0;k<KERNEL_SIZE;k++){
    for(l=0;l<KERNEL_SIZE;l++){
      x=p.hx*k-1.5;
      y=p.hx*l-1.5;
      C(k,l)=exp(-1*((x*x)+(y*y))/(2*(L*L)));
    }
  }

  /* Set boundary conditions on electrodes 
   * - Outer boundary (edges): V = Vext[0]
   * - Inner circular electrode: V = Vext[1]
   * - Interior points will be solved iteratively
   */
  for(i=0;i<=p.Nx;i++){
    for(int j=0;j<=p.Ny;j++){ 
      reg=in_region(i,j);
      if (reg>=0){
        A(i,j)=p.Vext[reg];  /* Set boundary electrode potentials */
      }
    }
  }  
}

/**
 * Clean up and deallocate memory
 * Frees all dynamically allocated arrays and the RNG
 */
void finalize(void){
  free(v1);  /* Free potential field */
  free(v2);  /* Free random field */
  free(v3);  /* Free Gaussian kernel */
  free(v4);  /* Free total potential array */
  gsl_rng_free (r);  /* Free random number generator */
}

/**
 * Perform 2D convolution of random field with Gaussian kernel
 * Creates spatially correlated initial conditions: A = B ⊗ C
 * 
 * The convolution creates smooth, correlated disorder from uniform random field:
 *   A(i,j) = Σ_k Σ_l B(m,n) * C(k,l)
 * where m = (kernel_size-1) - k + i, n = (kernel_size-1) - l + j
 * 
 * This models:
 * - Spatial correlations in random electrode potentials
 * - Smooth disorder in dielectric materials
 * - Correlated thermal fluctuations
 * 
 * @param v1 Pointer to potential array A (output)
 * @param v2 Pointer to random field B (input)
 * @param v3 Pointer to Gaussian kernel C (filter)
 */
void convolute(double *v1,double *v2,double *v3){
  int i,j,k,l,m,n;
  int a = KERNEL_SIZE;  /* Size of convolution kernel */
  int b = p.Nx;         /* Size of random field */
  
  /* Perform 2D convolution */
  for(i=0;i < (b-a)+1;i++){ 
    for(j=0;j< (b-a)+1;j++){ 
      for(k=0;k < a;k++){ 
	for(l=0;l < a;l++){ 
          m=(a-1)-k+i;  /* Flipped kernel index for convolution */
          n=(a-1)-l+j;
	  A(i,j) += B(m,n)*C(k,l);  /* Accumulate convolution sum */
	}
      }
    }
  }
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
/**
 * Output Gaussian kernel in gnuplot format
 * @param v3 Pointer to kernel array C
 */
void output_colg(double *v3){ 
  int k,l;
  
  for(k=0;k<KERNEL_SIZE;k++){
    for(l=0;l<KERNEL_SIZE;l++)
      printf("%g \t %g \t %g\n", p.hx*k-1.5, p.hx*l-1.5, C(k,l));
    printf("\n"); 
  }
  printf("\n");
}

/**
 * Gauss-Seidel/SOR iteration with thermal effects
 * 
 * STANDARD LAPLACE: A(i,j) = 0.25*(A(i+1,j) + A(i-1,j) + A(i,j+1) + A(i,j-1))
 * 
 * THERMAL MODIFICATION: Each neighbor contribution is weighted by a Boltzmann factor
 *   A(i,j) = 0.25 * Σ A(neighbor) * exp(-(V(i,j) + V(neighbor))/(2*T))
 * 
 * where:
 *   V(i,j) = A(i,j) + F*i  (total potential including external field)
 *   T = temperature parameter
 *   F = external field strength
 * 
 * PHYSICAL INTERPRETATION:
 * - The factor exp(-(V(i,j)+V(neighbor))/(2*T)) models thermal relaxation
 * - At high temperature (large T): factor ≈ 1, standard averaging
 * - At low temperature (small T): exponentially suppresses high-energy configurations
 * - External field F creates a linear potential gradient across the domain
 * 
 * The method iterates until convergence (error < eps) or max iterations reached.
 * 
 * @return 1 if converged, 0 if max iterations reached without convergence
 */
int do_gs(){
  int j,i,k;
  double error, diff, val;
  
  /* Iterate until convergence or max iterations */
  for(k=0;k<p.Nt;k++){
    error=0;  /* Track maximum change in this iteration */
    
    /* Update interior grid points (boundaries are fixed) */
    for(i=1;i<p.Nx;i++){
      for(j=1;j<p.Ny;j++){
        /* Skip boundary electrode points - only update interior */
        if(in_region(i,j)<0){
          /* Calculate total potential V(i,j) = A(i,j) + F*i
           * This includes the external field contribution */
          V(i,j)=A(i,j)+p.F*i;
          
          /* Also compute V for all four neighbors */
          V(i+1,j)=A(i+1,j)+p.F*(i+1);
          V(i-1,j)=A(i-1,j)+p.F*(i-1);
          V(i,j+1)=A(i,j+1)+p.F*i;
          V(i,j-1)=A(i,j-1)+p.F*i;
          
          /* Thermal-weighted average of four neighbors
           * Each term: A(neighbor) * exp(-(V(center)+V(neighbor))/(2*T))
           * This is the core thermal modification to Laplace equation */
          val=0.25*(A(i+1,j)*exp(-(V(i,j)+V(i+1,j))/(2*p.T)) +
                    A(i-1,j)*exp(-(V(i,j)+V(i-1,j))/(2*p.T)) +
                    A(i,j+1)*exp(-(V(i,j)+V(i,j+1))/(2*p.T)) +
                    A(i,j-1)*exp(-(V(i,j)+V(i,j-1))/(2*p.T)));
          
          diff=fabs(val-A(i,j));  /* Change from previous iteration */
          A(i,j)=val;              /* Update with new value (Gauss-Seidel) */
          
          if(error<diff)
            error=diff;  /* Track maximum change across all points */
        }
      }
    }
    
#if DEBUG
    printf("# Iteration k=%d, error=%g\n",k,error);
    output_matrix(v1);  
#endif /* DEBUG */
    
    /* Check for convergence */
    if(error<p.eps){
      fprintf(stderr,"# Converged after %d iterations (error=%g < eps=%g)\n",
              k+1, error, p.eps);
      return 1; /* Successfully converged */
    }
  }
  
  fprintf(stderr,"# Max iterations reached (%d), error=%g >= eps=%g\n",
          p.Nt, error, p.eps);
  return 0; /* Did not converge within max iterations */
}

int main(int argc, char * argv[]){
  int iter=0;
  
  initialize(argc, argv);
  convolute(v1,v2,v3);

  do{
    //  output_charge(v1);
#if DEBUG
    output_matrix(v1);
    output_matrix(v2);
#endif /* DEBUG */
    
    output_col(v1);
    fflush(stdout);
    iter++;    
  }
  while ((! do_gs()) && iter<p.Nt); 
  
  //output_col(v1); /* final output */
  //output_charge(v1);
  //output_matrix(v3);  
  // output_colg(v3);
  //output_matrix(v1);
   output_col(v1);
  finalize();
  fprintf(stderr,"hx=%5g alpha=%5g %g seconds\n",p.hx,
	  p.alpha,(double) (clock())/CLOCKS_PER_SEC);
  
  return 0;
}