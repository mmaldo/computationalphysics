//#define DEBUG 1 
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define A(i,j) v1[(30)*(i)+(j)]  // potnetial
#define B(m,n) v2[(p.Nx)*(m)+(n)]  //RandomNumbers
#define C(k,l) v3[(31)*(k)+(l)]  //Kernal
#define V(i,j) v4[(31)*(i)+(j)]

typedef struct {
  double hx; /* coordinate step */
  double Xmax; /* real-space dimensions */
  double Ymax;
  int Nx;    /* how many steps in x dimension */
  int Ny;    /* how many steps in y dimension */
  int Nt;    /* how many time steps */
  int l;
  int k;
  double eps;
  double alpha; /* convergence parameter <=2 */
  
  int Nb;    /* how many boundaries */
  const double * Vext; /* potentials on the boundaries */
  double R;   /* circle radius */
  double Xc;  /* center X */
  double Yc;  /* center Y */
  double T;
  double F; 
} params;

gsl_rng * r;   
const gsl_rng_type * T;        /* global generator */
const double VV[]={1.0, 2.0};
params p = { 0.1, 6.0,6.0, 30, 30, 30,30 ,100, 1.0e-10, 1.0,
             2, VV ,
             1., 1.5, 1.5,2.0,2.0}; 

double * v1; /* pointer to a vector */
double * v2;
double * v3;
double * v4;


double randomNumber(){
  double u;
  u = gsl_rng_uniform (r);
  return u;
}    
     
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
  int i,k,l,m,n;
  // int reg;
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
	      "\tusage: %s [hx=#] [alhpa=#]\n"
	      "\twhere 0<=alpha<2 is theconvergence factor\n"
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
  
  v2=calloc((p.Nx+1)*(p.Ny),sizeof(double)); 
  /* allocate mem for the A */
    if(v2==NULL){
       fprintf(stderr,"failed to allocate mem for v2\n");
       exit(-1);
     }
    
    v3=calloc((51)*(51),sizeof(double)); 
  /* allocate mem for the B */
   if(v3==NULL){
   fprintf(stderr,"failed to allocate mem for v3\n");
   exit(-1);
   }
   v4=calloc((61)*(61),sizeof(double)); 
  /* allocate mem for the B */
   if(v4==NULL){
   fprintf(stderr,"failed to allocate mem for v3\n");
   exit(-1);
   }
  
  for(i=0;i<p.Nb;i++)
    printf(" Vext[%d]=%g",i,p.Vext[i]);
    printf("\n");
    v1=calloc((61)*(61),sizeof(double)); 
  /* allocate mem for the potential */
  if(v1==NULL){
    fprintf(stderr,"failed to allocate mem for v1\n");
    exit(-1);
  }
  
 for(m=0;m <p.Nx;m++){
    for(n = 0;n<p.Ny;n++){
          
      B(m,n)=randomNumber();
    }
  }
 double L =p.hx*10;
 double x;
 double y;
  
 for(k=0;k<31;k++){
   for(l=0;l<31;l++){
     x=p.hx*k-1.5, y=p.hx*l-1.5;
     C(k,l)=exp(-1*((x*x)+(y*y))/(2*(L*L)));
   }
 }

 //  for(i=0;i<=p.Nx;i++){
 // for(j=0;j<=p.Ny;j++){
 //   reg=in_region(i,j);
 //   if (reg>=0){
 //     A(i,j)=p.Vext[reg];
 //     }  /* set boundary values */      
  // }
 // }  
}

void finalize(void){
  free(v1);
  free(v2);
  free(v3);   /* deallocate the arrays */
  gsl_rng_free (r);
}

void convolute(double *v1,double *v2,double *v3){
  int i,j,k,l,m,n;
  int a =31;
  int b = p.Nx;
  for(i=0;i < (b-a)+1;i++){
    for(j=0;j< (b-a)+1;j++){
      for(k=0;k < a;k++){
	for(l=0;l < a;l++){
          m=(a-1)-k+i;
          n=(a-1)-l+j;
	  A(i,j)= B(m,n)*C(k,l);
	  
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
void output_colg(double *v3){ /* output for gnuplot's splot  */
  int k,l;
  
  for(k=0;k<31;k++){
    for(l=0;l<31;l++)
      printf("%g \t %g \t %g\n", p.hx*k-1.5, p.hx*l-1.5, C(k,l));
    printf("\n"); 
  }
  printf("\n");
}

//A(i,j)=0.25*(A(i+1,j)*exp(-1*(v(i,j)+v(i+1,j))/2T)+(A(i-1,j)*exp(-1*(v(i,j)+v(i-1,j))/2T)+(A(i,j+1)*exp(-1*(v(i,j)+v(i,j+1))/2T)+(A(i,j-1)*exp(-1*(v(i,j)+v(i,j-1))/2T)))));
//V=A(i,j)+F*i;
int do_gs(){
  int j,i,k;
  double error, diff, val;
  for(k=0;k<p.Nt;k++){
    error=0;           /* maximum difference from previous iteration */
    for(i=1;i<p.Nx;i++){
      for(j=1;j<p.Ny;j++){ /* update inside values */
	// if(in_region(i,j)<0){
	  //  val=(1-p.alpha)*A(i,j)+
          //  p.alpha*0.25*(A(i-1,j)+A(i+1,j)+A(i,j-1)+A(i,j+1));
           V(i,j)=A(i,j)+p.F*i;
          val=0.25*(A(i+1,j)*exp(-1*(V(i,j)+V(i+1,j))/(2*p.T))+(A(i-1,j)*exp(-1*(V(i,j)+V(i-1,j))/(2*p.T))+(A(i,j+1)*exp(-1*(V(i,j)+V(i,j+1))/(2*p.T))+(A(i,j-1)*exp(-1*(V(i,j)+V(i,j-1))/(2*p.T))))));
          diff=fabs(val-A(i,j)); /* difference */
          A(i,j)=val;                  
          if(error<diff)
            error=diff;  /* update the error */
	  //  }
      }
    }
#if DEBUG
    printf("# this is v1, k=%d\n",k);
    output_matrix(v1);  
#endif /* DEBUG */
    if(error<p.eps) return 1; /* done */
  }
  return 0;
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
  // output_matrix(v1);
   output_col(v1);
  finalize();
  fprintf(stderr,"hx=%5g alpha=%5g %g seconds\n",p.hx,
	  p.alpha,(double) (clock())/CLOCKS_PER_SEC);
  
  return 0;
}
