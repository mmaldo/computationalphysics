/**
This is a dampened harmonic osciliattor
RK2
theta= v
v = -v/Q-Theta + Acoswt
Author: Matthew Maldoando
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


typedef struct {  /* this type will be used to hold parameters */
  double tmax;      /* evolution time, seconds */
  double w;        /* frequency in cos */
  double A;        /* Amplitude of the cos */
  double q;	    /*quality factor*/
  double x0;        /* initial position */
  double v0;        /* initial velocity */
  int n;            /* how many steps */
  double h;         /* step size */
  double energ0;    /* initial energy */
} params;
params p = {100,  2/3.0 , 1.5 , 1.372 , 0, 0 , 100 , 0,0};
params p1 = {100, 2/3.0,1.5, 1.372,1/1000000.0,1/1000000.0,100,0,0};



void rk4(){
double t = 0;
double v = p.v0;
double v1 = p1.v0;
double x = p.x0;
double x1 = p1.x0;
p.h = (2*M_PI/p.w)/p.tmax;
double xk1, xk2, xk3, xk4, vk1,vk2,vk3,vk4;
double x1k1, x1k2, x1k3, x1k4, v1k1,v1k2,v1k3,v1k4;
double tcount;
double difference=fabs(x-x1);
double differencev=fabs(v-v1);
double ld=log(difference);
double lv=log(differencev);
p.energ0=p.w*(p.x0*p.x0*0.5*(1+0.05*p.x0*p.x0))+0.5*p.v0*p.v0; /* initial energy */
p1.energ0=p1.w*(p1.x0*p1.x0*0.5*(1+0.05*p1.x0*p1.x0))+0.5*p1.v0*p1.v0;
printf("#time\t x\t v\t energy\t x1\t v1\t energy1\t x-x1\t v-v1\t logx\t logv\n");
 printf("%g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\n",t,x/(M_PI),v/(M_PI),p.energ0, x1/(M_PI),v1/(M_PI), p1.energ0, difference, differencev,ld,lv);

for(tcount=0; tcount < (p.tmax); tcount= t + p.h){
	xk1 = v;
	x1k1 = v1;
	vk1 =(-v/p.q) - sin(x) + p.A*cos(p.w*t);
	v1k1 =(-v1/p1.q) - sin(x1) + p1.A*cos(p1.w*t); 
	
	t += 0.5*p.h;
	xk2 = v + ((0.5)*p.h*xk1);
	x1k2 = v1 + ((0.5)*p.h*x1k1);
	vk2=(-v/p.q)- sin(x) + p.A*cos(p.w*t) + ((0.5)*p.h*vk1);
	v1k2=(-v1/p1.q)-sin(x1) + p1.A*cos(p1.w*t) + ((0.5)*p.h*v1k1);
	
	xk3 = v + ((0.5)*p.h*xk2);
	x1k3 = v1 + ((0.5)*p.h*x1k2);
	vk3= (-v/p.q)- sin(x) + p.A*cos(p.w*t) + ((0.5)*p.h*vk2);
	v1k3= (-v1/p1.q)- sin(x1) + p1.A*cos(p1.w*t) + (0.5)* p.h *v1k2;
	t += (0.5*p.h);
	
	xk4 = v + (p.h*xk3);
	x1k4 = v1 + (p.h*x1k3);
	vk4=(-v/p.q)- sin(x) + p.A*cos(p.w*t) + ((0.5)*p.h*vk3);
	v1k4=(-v1/p1.q)- sin(x1) + p1.A*cos(p1.w*t) + ((0.5)*p.h*v1k3);
	
	x = x + ((1.0/6.0)*p.h)*(xk1 + (2*xk2) + (2*xk3) + xk4);
	x1 = x1 + ((1.0/6.0)*p.h)*(x1k1 + (2*x1k2) + (2*x1k3) + x1k4);
	v = v + ((1.0/6.0)*p.h)*(vk1 + (2*vk2) + (2*vk3) + vk4);
	v1 = v1 + ((1.0/6.0)*p.h)*(v1k1 + (2*v1k2) + (2*v1k3) + v1k4);
      
      double energ=p.w*(x*x*0.5+x*x*x*x*0.025)+0.5*v*v;
      double energ1=p1.w*(x1*x1*0.5+x1*x1*x1*x1*0.025)+0.5*v1*v1;
	
	
	//if(x > M_PI){
	//x-=(2*M_PI);
	//}else if(x < -M_PI){
	//x+=(2*M_PI);
	//}else if(x1 > M_PI){
	//x1-=(2*M_PI);
	//}else if(x1 < -M_PI){
	//x1+=(2*M_PI);
	//}
        
	difference = fabs(x-x1);
	differencev = fabs(v-v1);
        ld =log10(difference);
	lv = log10(differencev);
	printf("%g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\n", t , x/(M_PI),v/(M_PI), energ, x1/(M_PI),v1/(M_PI), energ1, difference, differencev,ld,lv);
       }

}



int main(){
rk4();

return 0;

}
