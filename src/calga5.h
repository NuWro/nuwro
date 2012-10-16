#ifndef _calga_h_
#define _calga_h_
#include <math.h>
#include <errno.h>


template <class F>
double calg5(F f,double x1,double x2,int ile=20)
{//double static const sqrt2=sqrt(2);
if(ile==1)
{
 static double const a=0.53846931010568309105; //sqrt(245-14*sqrt(70))/21
 static double const b=0.90617984593866399282; //sqrt(245+14*sqrt(70))/21
 static double const f0=128.0/225;
 static double const fb=(1.0/3 -a*a*(1-f0/2))/(b*b-a*a); // 0.23692689
 static double const fa=1-f0/2-fb; // 0.47862867
 double h=(x2-x1)/2;
 double xs=x1+h;
 double   y1m=f(xs-h*a);
 double   y1p=f(xs+h*a);
 double   y2m=f(xs-h*b);
 double   y2p=f(xs+h*b);
 double   ys=f(xs);
 //double r1=(y1m+y1p)/2-ys;
 //double r2=(y2m+y2p)/2-ys;
 return h*(f0*ys+fa*(y1m+y1p)+fb*(y2m+y2p));  
}
else
{double d=(x2-x1)/ile;
 double suma=0;
 while(ile-->0)
 {suma+=calg5(f,x1+ile*d,x1+ile*d+d,1);
 }
 return suma;
} 
}



template <class F>
double calg20(F f,double x1,double x2, int n=3)
{
static const double w20[10]={
  0.017614007139152118311861962351853,
  0.040601429800386941331039952274932,
  0.062672048334109063569506535187042,
  0.083276741576704748724758143222046,
  0.101930119817240435036750135480350,
  0.118194531961518417312377377711382,
  0.131688638449176626898494499748163,
  0.142096109318382051329298325067165,
  0.149172986472603746787828737001969,
  0.152753387130725850698084331955098
};

static const double g20[10]={
  0.993128599185094924786122388471320,
  0.963971927277913791267666131197277,
  0.912234428251325905867752441203298,
  0.839116971822218823394529061701521,
  0.746331906460150792614305070355642,
  0.636053680726515025452836696226286,
  0.510867001950827098004364050955251,
  0.373706088715419560672548177024927,
  0.227785851141645078080496195368575,
  0.076526521133497333754640409398838
};

//if (n>1) return calg20(f,x1,x1+(x2-x1)*(n/2)/n,n/2)+calg20(f,x1+(x2-x1)*(n/2)/n,x2,n-n/2);

double wynik=0;
double A=(x2-x1)/2/n;
 for (int i=0;i<n;i++)
 {
 double B=(x1+x2)/2+(2*i+1-n)*A;
 for (int j=0; j<10; j++)
   wynik+=(f(B+A*g20[j])+f(B-A*g20[j]))*w20[j];
 }
 return wynik*A;
}

template <class F>
double calg40(F f,double x1,double x2)
{ 
  return calg20(f,x1,(x1+x2)/2)+calg20(f,(x1+x2)/2,x2);
}

template <class F>
double calg60(F f,double x1,double x2)
{ double a=(x2+2*x1)/3, b=(x1+2*x2)/3;
 return calg20(f,x1,a)
 +calg20(f,a,b)+
  calg20(f,b,x2);
}

template <class F>
double calg80(F f,double x1,double x2)
{ double b=(x1+x2)/2,a=(x1+b)/2,c=(b+x2)/2;
  return  calg20(f,x1,a)
         +calg20(f,a,b)
         +calg20(f,b,c)
	 +calg20(f,c,x2);
}

template <class F>
double calg100(F f,double x1,double x2)
{return 
calg20(f,x1,(x2+4*x1)/5)+
calg20(f,(x2+4*x1)/5,(2*x2+3*x1)/5)+
calg20(f,(3*x1+2*x2)/5,(2*x1+3*x2)/5)+
calg20(f,(2*x1+3*x2)/5,(x1+4*x2)/5)+
calg20(f,(x1+4*x2)/5,x2);
}


#endif
