/**
 * @file            NewtonsMethod.cpp
 * @description     This program aprroximates a solution to x^2-5
 *                  using Newtons Method.
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

/**
 * Function used for approximation of x^2-5 is f(x)=x-A(x^2-5)
 * @param initialPoint is the initial guess
 * @pre   initialPoint is initialized by the calling program
 * @post  initialPoint is sent into the approixmation function
 *        to get a new vslue
 * @return the the new value of the function, also new x value
 */
double fun(double intialPoint);

/**
 * @retunrn zero if application executes succesfully
 */
int main()
{
    double initialPoint = 1;        //initial guess
    double newPoint = initialPoint; // x_n
    double previousPoint;           //x_n-1
    double convergenceCompare;        //|x_n-sqrt(5)|
    int i;                          //counter
   // FILE *fp;                       //file to write to
        //enter file to write to
    //fp = fopen("C:/Users/mmaldonado/Documents/My Dropbox/UCR/Spring Quarter/Physics 177/Assignments/gnuplot/newtonsMethod.data", "w");
    printf("#Newton's Method\n#n   X_n-1   X_n   |x_n-sqrt(5)|\n");
    for(i = 0; i < N; i++){
        previousPoint = newPoint;
        newPoint = fun(newPoint);
        convergenceCompare = fabs(newPoint - sqrt(5));
       // fprintf(fp,"%d      %g      %g      %g\n",i+1, previousPoint, newPoint, convergenceCompare);
        printf("%d      %g      %g     %g\n",i+1, previousPoint, newPoint, convergenceCompare);
    }
    //fclose(fp); undo comment to write to file
    return 0;
}

double fun(double intialPoint){
    double x = intialPoint;

    return (x - ((x*x)-5)/(2*x));

}

