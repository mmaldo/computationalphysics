#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(){
int i;
int max = 1000;
int thermalization = 500;  /* discard first iterations to reach attractor */
double r;
double x_n1;
double x_n = 0.25;
printf("#x_n1     r\n");
for(r = 2.5; r <=4; r+=0.001){
	x_n = 0.25;  /* reset initial condition for each r */
	for(i = 0; i < max; i++){
		x_n1 = r*x_n*(1-x_n);
		x_n= x_n1;
		/* only output after thermalization period */
		if(i >= thermalization){
			printf("%g\t %g\n", x_n1,r);
		}
	}
}
return 0;
}

