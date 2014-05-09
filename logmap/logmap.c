#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(){
int i;
int max = 1000;
double r;
double x_n1;
double x_n = 0.25;
printf("#x_n1     r");
for(r = 2.5; r <=4; r+=0.001){
	//printf("%g\n",r);
	for(i = 0; i < max; i++){
		x_n1 = r*x_n*(1-x_n);
		x_n= x_n1;
		printf("%g\t %g\n", x_n1,r);
	}
}
return 0;
}

