#ifndef _calgauss20_h_
#define _calgauss20_h_
#include <iostream>
#include <math.h>
#include "g20.h"
#include "w20.h"

double calgauss20(double (*f)(double),double x1,double x2)
{double A=(x2-x1)/2;
double B=(x2+x1)/2;
double wynik=0;
for (int j=0; j<=19; j++)
wynik+=f(A*g20[j]+B)*w20[j];
return wynik*A;
}
double calgauss40(double (*f)(double),double x1,double x2)
{return calgauss20(f,x1,(x1+x2)/2)+calgauss20(f,(x1+x2)/2,x2);}

double calgauss60(double (*f)(double),double x1,double x2)
{return calgauss20(f,x1,(x2+2*x1)/3)+calgauss20(f,(x2+2*x1)/3,(x1+2*x2)/3)+
calgauss20(f,(x1+2*x2)/3,x2);}

double calgauss80(double (*f)(double),double x1,double x2)
{return calgauss20(f,x1,(x2+3*x1)/4)+calgauss20(f,(x2+3*x1)/4,(x2+x1)/2)
+calgauss20(f,(x1+x2)/2,(x1+3*x2)/4)+calgauss20(f,(x1+3*x2)/4,x2);}

double calgauss100(double (*f)(double),double x1,double x2)
{return 
calgauss20(f,x1,(x2+4*x1)/5)+
calgauss20(f,(x2+4*x1)/5,(2*x2+3*x1)/5)+
calgauss20(f,(3*x1+2*x2)/5,(2*x1+3*x2)/5)+
calgauss20(f,(2*x1+3*x2)/5,(x1+4*x2)/5)+
calgauss20(f,(x1+4*x2)/5,x2);}


#endif


