#ifndef _calgauss20_h_
#define _calgauss20_h_
#include <iostream>
#include <math.h>

double w20[20]=
{
    0.017614, 0.0406014, 0.062672, 0.0832767, 
    0.10193,  0.118195, 0.131689, 0.142096, 
    0.149173, 0.152753, 0.152753, 0.149173, 
    0.142096, 0.131689, 0.118195, 0.10193, 
    0.0832767, 0.062672, 0.0406014,0.017614
};

double g20[20]=
{
    -0.9931285991850298,-0.9639719272747169,-0.9122344282502259,-0.839116971821869,
    -0.7463319064598453,-0.6360536807265005,-0.5108670019508232,-0.3737060887154194,
    -0.227785851141645,-0.07652652113349734,0.07652652113349734,0.227785851141645,
    0.3737060887154201,0.5108670019508297,0.6360536807264445,0.7463319064602532,
    0.8391169718218283,0.9122344282514332,0.9639719272771623,0.9931285991846996
};

double calgauss20(double (*f)(double),double x1,double x2)
{
    double A=(x2-x1)/2;
    double B=(x2+x1)/2;
    double wynik=0;
    for (int j=0; j<=19; j++)
        wynik += f(A*g20[j]+B)*w20[j];
    return wynik*A;
}

double calgauss40(double (*f)(double),double x1,double x2)
{
    return calgauss20(f,x1,(x1+x2)/2)+calgauss20(f,(x1+x2)/2,x2);
}

double calgauss60(double (*f)(double),double x1,double x2)
{
    return  calgauss20(f,x1,(x2+2*x1)/3)+
            calgauss20(f,(x2+2*x1)/3,(x1+2*x2)/3)+
            calgauss20(f,(x1+2*x2)/3,x2);
}

double calgauss80(double (*f)(double),double x1,double x2)
{
    return calgauss20(f,x1,(x2+3*x1)/4)
            +calgauss20(f,(x2+3*x1)/4,(x2+x1)/2)
            +calgauss20(f,(x1+x2)/2,(x1+3*x2)/4)+calgauss20(f,(x1+3*x2)/4,x2);
}

double calgauss100(double (*f)(double),double x1,double x2)
{
    return
            calgauss20(f,x1,(x2+4*x1)/5)+
            calgauss20(f,(x2+4*x1)/5,(2*x2+3*x1)/5)+
            calgauss20(f,(3*x1+2*x2)/5,(2*x1+3*x2)/5)+
            calgauss20(f,(2*x1+3*x2)/5,(x1+4*x2)/5)+
            calgauss20(f,(x1+4*x2)/5,x2);
}


#endif
