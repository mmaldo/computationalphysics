/**
 * @file            fixedPointsolve.cpp
 * @description     This program aprroximates a solution to a algebraic equation
 *                  using fixed point iteration.
 * @course          Physics 177
 * @assignment      Homework 1 Problem 2
 * @date            4/7/2011
 * @author          Matthew Maldoando
 * @version         1.0
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

const int N = 10;   //number of iterations
double A = 0.1;  //constant number
/**
 * Function used for approximation f(x)=x-A(x^2-5)
 * @param initialPoint is the initial guess
 * @param A is a constant
 * @pre   initialPoint is initialized by the calling program
 * @post  initialPoint is sent into the approixmation function
 *        to get a new vslue
 * @return the the new value of the function
 */
double fun(double intialPoint, double A);

/**
 * @retunrn zero if application executes succesfully
 */
int main()
{
    double initialPoint = 1;        //initial guess
    double newPoint = initialPoint; //variable to store newpoint
    double previousPoint;           //variable to store x_n-1
    double convergenceCompare;        //|x_n-sqrt(5)|
    int i;                        //counter variabel
   // FILE *fp;                       //file to write result to
    //enter file to write to if you want
    //fp = fopen("C:/Users/mmaldonado/Documents/My Dropbox/UCR/Spring Quarter/Physics 177/Assignments/gnuplot/fixedPointsolve.data", "w");
    printf("#Fixed Point Solve\n#n   X_n-1   X_n   |x_n-sqrt(5)|\n");
    for(i = 0; i < N; i++){
        previousPoint = newPoint;
        newPoint = fun(newPoint,A);
        convergenceCompare = fabs(newPoint - sqrt(5));
        //fprintf(fp,"%d      %g      %g      %g\n",i+1, previousPoint, newPoint, convergenceCompare);
        printf("%d      %g      %g     %g\n",i+1, previousPoint, newPoint, convergenceCompare);
    }
    //fclose(fp);  undo commentif writing to file
    return 0;
}

double fun(double intialPoint, double A){
    double x = intialPoint;
    double a = A;
return (x - a*((x*x)-5));

}




