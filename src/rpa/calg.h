#ifndef _calg_h_
#define _calg_h_
#include <cmath>
#include <cerrno>

namespace rpa
{
double defdokl=1e-12;
int lcalg;
int lcalg3;
int lcalg5;
int lcalk;

double calk(double (*f)(double),double x1,double x2,double y1=0,double y2=0,int depth=0)
{
	if(x1==x2) 
	{ 
		return 0;
	}
	else
	{
		if(depth==0)
			lcalk=0;
		
		if(y1==0)
		{  
			y1=(*f)(x1);lcalk++;
		}
		
		if(y2==0)
		{     
			y2=(*f)(x2);lcalk++;
		}
		
		double xs=(x1+x2)/2;
		double ys=(*f)(xs);lcalk++;
		
		if(errno) 
			return errno=0;
		
		if( (fabs(2*ys-y1-y2)<=0.0000001 *(fabs(y1)+fabs(y2)+fabs(ys)) || depth>15) && depth>6)
			return (4*ys+y1+y2)*(x2-x1)/6; //metoda Simpsona (dopasowanie parabol±)
		else
			return calk(f,x1,xs,y1,ys,depth+1)+calk(f,xs,x2,ys,y2,depth+1);    
	}
}


//metoda Gaussa G2
double calg(double (*f)(double),double x1,double x2,int depth=0)
{
	double static const sqrt3=sqrt(3);
	if(x1==x2) 
	{ 
		return 0;
	}
	else
	{ 
		if(depth==0)
    		lcalg=0;
		double h=(x2-x1)/2;
		double xs=x1+h;
		 //double ys=(*f)(xs);
		if(errno) 
			return errno=0;
		 if(depth>5)
		 {
			double   y1=(*f)(xs+h/sqrt3);
			double   y2=(*f)(xs-h/sqrt3);
			lcalg+=2;
			return h*(y1+y2); 
		 } 
		 else
			return calg(f,x1,xs,depth+1)+calg(f,xs,x2,depth+1);    
	}
}

double calg3(double (*f)(double),double x1,double x2,double dokl=defdokl)
{
	double static const s6=sqrt(0.6);
	double static const sqrt2=sqrt(2);
	if(dokl==defdokl)
		lcalg3=0;
	double h=(x2-x1)/2;
	double xs=x1+h;
	double   y1=(*f)(xs+h*s6);
	double   y2=(*f)(xs-h*s6);
	double   ys=(*f)(xs);
	double   sigma=fabs(y1+y2-2*ys);

	lcalg3+=3;
	if(sigma*2*h<dokl)
		return h*(5*(y1+y2)+8*ys)/9; 
	else 
		return calg3(f,x1,xs,dokl/sqrt2) +calg3(f,xs,x2,dokl/sqrt2); 
}

double calg3dokl(double (*f)(double,double),double x1,double x2,double dokl=1e-10)
{
	double static const s6=sqrt(0.6);
	double static const sqrt2=sqrt(2);
	double h=(x2-x1)/2;
	double xs=x1+h;
	double   y1=(*f)(xs+h*s6,dokl/h/3);
	double   y2=(*f)(xs-h*s6,dokl/h/3);
	double   ys=(*f)(xs,dokl/h/3);
	double   sigma=fabs(y1+y2-2*ys);
	if(sigma*2*h<dokl)
		return h*(5*(y1+y2)+8*ys)/9; 
	else 
		return calg3dokl(f,x1,xs,dokl/sqrt2) +calg3dokl(f,xs,x2,dokl/sqrt2); 
}

double calg5(double (*f)(double),double x1,double x2,double dokl=defdokl,int ile=1)
{//double static const sqrt2=sqrt(2);
	if(ile==1)
	{
		static double const a=0.53846931010568309105; //sqrt(245-14*sqrt(70))/21
		static double const b=0.90617984593866399282; //sqrt(245+14*sqrt(70))/21
		static double const f0=128.0/225;
		static double const fb=(1.0/3 -a*a*(1-f0/2))/(b*b-a*a); // 0.23692689
		static double const fa=1-f0/2-fb; // 0.47862867
		double h=(x2-x1)/2;
		if(dokl==defdokl)
			lcalg5=0;
		double xs=x1+h;
		double   y1m=(*f)(xs-h*a);
		double   y1p=(*f)(xs+h*a);
		double   y2m=(*f)(xs-h*b);
		double   y2p=(*f)(xs+h*b);
		double   ys=(*f)(xs);
		double r1=(y1m+y1p)/2-ys;
		double r2=(y2m+y2p)/2-ys;
		lcalg5+=5;
		if( fabs((r2*a*a-r1*b*b)*10*h)<dokl)     
			return h*(f0*ys+fa*(y1m+y1p)+fb*(y2m+y2p));  
		else
		{
			dokl/=2;h/=2;
			return calg5(f,x1,x1+h,dokl)+calg5(f,x1+h,xs,dokl)
				  +calg5(f,xs,xs+h,dokl)+calg5(f,xs+h,x2,dokl);    
		}  
	}
	else
	{
		double d=(x2-x1)/ile;
		dokl/=sqrt(ile);
		double suma=0;
		while(ile-->0)
		{
			suma+=calg5(f,x1+ile*d,x1+ile*d+d,dokl,1);
		}
		return suma;
	}
}

double calg5_int(double (*f)(int,double),int z, double x1,double x2,double dokl=defdokl,int ile=1)
{//double static const sqrt2=sqrt(2);
	if(ile==1)
	{
		static double const a=0.53846931010568309105; //sqrt(245-14*sqrt(70))/21
		static double const b=0.90617984593866399282; //sqrt(245+14*sqrt(70))/21
		static double const f0=128.0/225;
		static double const fb=(1.0/3 -a*a*(1-f0/2))/(b*b-a*a); // 0.23692689
		static double const fa=1-f0/2-fb; // 0.47862867
		double h=(x2-x1)/2;
		if(dokl==defdokl)
			lcalg5=0;
		double xs=x1+h;
		double   y1m=(*f)(z,xs-h*a);
		double   y1p=(*f)(z,xs+h*a);
		double   y2m=(*f)(z,xs-h*b);
		double   y2p=(*f)(z,xs+h*b);
		double   ys=(*f)(z,xs);
		double r1=(y1m+y1p)/2-ys;
		double r2=(y2m+y2p)/2-ys;
		lcalg5+=5;
		if( fabs((r2*a*a-r1*b*b)*10*h)<dokl)     
			return h*(f0*ys+fa*(y1m+y1p)+fb*(y2m+y2p));  
		else
		{
			dokl/=2;h/=2;
			return calg5_int(f,z,x1,x1+h,dokl)+calg5_int(f,z,x1+h,xs,dokl)
				  +calg5_int(f,z,xs,xs+h,dokl)+calg5_int(f,z,xs+h,x2,dokl);    
		}  
	}
	else
	{
		double d=(x2-x1)/ile;
		 dokl/=sqrt(ile);
		 double suma=0;
		 while(ile-->0)
		 {
		 	suma+=calg5_int(f,z,x1+ile*d,x1+ile*d+d,dokl,1);
		 }
		 return suma;
	} 
}

double calg5dokl(double (*f)(double,double),double x1,double x2,double dokl=defdokl,int ile=1)
{//double static const sqrt2=sqrt(2);
	if(x1==x2) 
	{ 
		return 0;
	}
	else
	{
		if(ile==1)
		{
			static double const a=0.53846931010568309105; //sqrt(245-14*sqrt(70))/21
			static double const b=0.90617984593866399282; //sqrt(245+14*sqrt(70))/21
			static double const f0=128.0/225;
			static double const fb=(1.0/3 -a*a*(1-f0/2))/(b*b-a*a); // 0.23692689
			static double const fa=1-f0/2-fb; // 0.47862867
			double h=(x2-x1)/2;
			if(dokl==defdokl)
			lcalg5=0;
			double xs=x1+h;
			double   y1m=(*f)(xs-h*a,dokl/h/5);
			double   y1p=(*f)(xs+h*a,dokl/h/5);
			double   y2m=(*f)(xs-h*b,dokl/h/5);
			double   y2p=(*f)(xs+h*b,dokl/h/5);
			double   ys=(*f)(xs,dokl/h/5);
			double r1=(y1m+y1p)/2-ys;
			double r2=(y2m+y2p)/2-ys;
			lcalg5+=5;
			if( fabs((r2*a*a-r1*b*b)*10*h)<dokl)     
				return h*(f0*ys+fa*(y1m+y1p)+fb*(y2m+y2p));  
			else
			{
				dokl/=2;h/=2;
				return calg5dokl(f,x1,x1+h,dokl)+calg5dokl(f,x1+h,xs,dokl)
					  +calg5dokl(f,xs,xs+h,dokl)+calg5dokl(f,xs+h,x2,dokl);    
			}  
		}
		else
		{
			 double d=(x2-x1)/ile;
			 dokl/=sqrt(ile);
			 double suma=0;
			 while(ile-->0)
			 {
				suma+=calg5dokl(f,x1+ile*d,x1+ile*d+d,dokl,1);
			 }
			 return suma;

		}
	}
}


double calg5x(double (*f)(double,double),double x,double x1,double x2,double dokl=defdokl,int ile=1)
{//double static const sqrt2=sqrt(2);
	if(ile==1)
	{ 
		static double const a=0.53846931010568309105; //sqrt(245-14*sqrt(70))/21
		static double const b=0.90617984593866399282; //sqrt(245+14*sqrt(70))/21
		static double const f0=128.0/225;
		static double const fb=(1.0/3 -a*a*(1-f0/2))/(b*b-a*a); // 0.23692689
		static double const fa=1-f0/2-fb; // 0.47862867
		double h=(x2-x1)/2;
		if(dokl==defdokl)
			lcalg5=0;
		double xs=x1+h;
		double   y1m=(*f)(xs-h*a,x);
		double   y1p=(*f)(xs+h*a,x);
		double   y2m=(*f)(xs-h*b,x);
		double   y2p=(*f)(xs+h*b,x);
		double   ys=(*f)(xs,x);
		double r1=(y1m+y1p)/2-ys;
		double r2=(y2m+y2p)/2-ys;
		lcalg5+=5;
		if( fabs((r2*a*a-r1*b*b)*10*h)<dokl)     
			return h*(f0*ys+fa*(y1m+y1p)+fb*(y2m+y2p));  
		else
		{
			dokl/=2;h/=2;
			return calg5x(f,x,x1,x1+h,dokl)+calg5x(f,x,x1+h,xs,dokl)
				  +calg5x(f,x,xs,xs+h,dokl)+calg5x(f,x,xs+h,x2,dokl);    
		}  
	}
	else
	{
		double d=(x2-x1)/ile;
		dokl/=sqrt(ile);
		double suma=0;
		while(ile-->0)
		{
			suma+=calg5x(f,x,x1+ile*d,x1+ile*d+d,dokl,1);
		}
		return suma;
	}
}

double calg5x_3int(double (*f)(double,double,int),double x,int z,double x1,double x2,double dokl=defdokl,int ile=1)
{//double static const sqrt2=sqrt(2);
	if(ile==1)
	{ 
		static double const a=0.53846931010568309105; //sqrt(245-14*sqrt(70))/21
		static double const b=0.90617984593866399282; //sqrt(245+14*sqrt(70))/21
		static double const f0=128.0/225;
		static double const fb=(1.0/3 -a*a*(1-f0/2))/(b*b-a*a); // 0.23692689
		static double const fa=1-f0/2-fb; // 0.47862867
		double h=(x2-x1)/2;
		if(dokl==defdokl)
			lcalg5=0;
		double xs=x1+h;
		double   y1m=(*f)(xs-h*a,x,z);
		double   y1p=(*f)(xs+h*a,x,z);
		double   y2m=(*f)(xs-h*b,x,z);
		double   y2p=(*f)(xs+h*b,x,z);
		double   ys=(*f)(xs,x,z);
		double r1=(y1m+y1p)/2-ys;
		double r2=(y2m+y2p)/2-ys;
			lcalg5+=5;
		if( fabs((r2*a*a-r1*b*b)*10*h)<dokl)     
			return h*(f0*ys+fa*(y1m+y1p)+fb*(y2m+y2p));  
		else
		{
			dokl/=2;h/=2;
			return calg5x_3int(f,x,z,x1,x1+h,dokl)+calg5x_3int(f,x,z,x1+h,xs,dokl)
				+calg5x_3int(f,x,z,xs,xs+h,dokl)+calg5x_3int(f,x,z,xs+h,x2,dokl);    
		}  
	}
	else
	{
		double d=(x2-x1)/ile;
		dokl/=sqrt(ile);
		double suma=0;
		while(ile-->0)
		{
			suma+=calg5x_3int(f,x,z,x1+ile*d,x1+ile*d+d,dokl,1);
		}
		return suma;
	}
}

double calg5x_4int(double (*f)(double,double,int,int,int,int),double x,
               int int1,int int2,int int3,int int4,double x1,double x2,double dokl=defdokl,int ile=1)
{//double static const sqrt2=sqrt(2);
	if(ile==1)
	{ 
		static double const a=0.53846931010568309105; //sqrt(245-14*sqrt(70))/21
		static double const b=0.90617984593866399282; //sqrt(245+14*sqrt(70))/21
		static double const f0=128.0/225;
		static double const fb=(1.0/3 -a*a*(1-f0/2))/(b*b-a*a); // 0.23692689
		static double const fa=1-f0/2-fb; // 0.47862867
		double h=(x2-x1)/2;
		if(dokl==defdokl)
		lcalg5=0;
		double xs=x1+h;
		double   y1m=(*f)(xs-h*a,x,int1, int2,int3, int4);
		double   y1p=(*f)(xs+h*a,x,int1, int2,int3, int4);
		double   y2m=(*f)(xs-h*b,x,int1, int2,int3, int4);
		double   y2p=(*f)(xs+h*b,x,int1, int2,int3, int4);
		double   ys=(*f)(xs,x,int1, int2,int3, int4);
		double r1=(y1m+y1p)/2-ys;
		double r2=(y2m+y2p)/2-ys;
		lcalg5+=5;
		if( fabs((r2*a*a-r1*b*b)*10*h)<dokl)     
			return h*(f0*ys+fa*(y1m+y1p)+fb*(y2m+y2p));  
		else
		{	dokl/=2;h/=2;
			return 
			calg5x_4int(f,x,int1, int2,int3, int4,x1,x1+h,dokl)+calg5x_4int(f,x,int1, int2,int3, int4, x1+h,xs,dokl)
			+   calg5x_4int(f,x,int1, int2,int3, int4,xs,xs+h,dokl)+calg5x_4int(f,x,int1, int2,int3, int4, xs+h,x2,dokl);    
		}  
	}
	else
	{
		double d=(x2-x1)/ile;
		dokl/=sqrt(ile);
		double suma=0;
		while(ile-->0)
		{
			suma+=calg5x_4int(f,x,int1, int2,int3, int4,x1+ile*d,x1+ile*d+d,dokl,1);
		}
		return suma;
	}
}


double calg5xy(double (*f)(double,double,double),double x,double y,double x1,double x2,double dokl=defdokl,int ile=1)
{//double static const sqrt2=sqrt(2);
	if(ile==1)
	{ 
		static double const a=0.53846931010568309105; //sqrt(245-14*sqrt(70))/21
		static double const b=0.90617984593866399282; //sqrt(245+14*sqrt(70))/21
		static double const f0=128.0/225;
		static double const fb=(1.0/3 -a*a*(1-f0/2))/(b*b-a*a); // 0.23692689
		static double const fa=1-f0/2-fb; // 0.47862867
		double h=(x2-x1)/2;
		if(dokl==defdokl)
			lcalg5=0;
		double xs=x1+h;
		double   y1m=(*f)(xs-h*a,x,y);
		double   y1p=(*f)(xs+h*a,x,y);
		double   y2m=(*f)(xs-h*b,x,y);
		double   y2p=(*f)(xs+h*b,x,y);
		double   ys=(*f)(xs,x,y);
		double r1=(y1m+y1p)/2-ys;
		double r2=(y2m+y2p)/2-ys;
			lcalg5+=5;
		if( fabs((r2*a*a-r1*b*b)*10*h)<dokl)     
			return h*(f0*ys+fa*(y1m+y1p)+fb*(y2m+y2p));  
		else
		{
			dokl/=2;h/=2;
			return calg5xy(f,x,y,x1,x1+h,dokl)+calg5xy(f,x,y,x1+h,xs,dokl)
				  +calg5xy(f,x,y,xs,xs+h,dokl)+calg5xy(f,x,y,xs+h,x2,dokl);    
		}  
		}
		else
		{
			double d=(x2-x1)/ile;
			dokl/=sqrt(ile);
			double suma=0;
			while(ile-->0)
			{
				suma+=calg5xy(f,x,y,x1+ile*d,x1+ile*d+d,dokl,1);
			}
			return suma;
	}
}



double calg5x_int(double (*f)(double,int),int x,double x1,double x2,double dokl=defdokl,int ile=1)
{//double static const sqrt2=sqrt(2);
	if(ile==1)
	{ 
		static double const a=0.53846931010568309105; //sqrt(245-14*sqrt(70))/21
		static double const b=0.90617984593866399282; //sqrt(245+14*sqrt(70))/21
		static double const f0=128.0/225;
		static double const fb=(1.0/3 -a*a*(1-f0/2))/(b*b-a*a); // 0.23692689
		static double const fa=1-f0/2-fb; // 0.47862867
		double h=(x2-x1)/2;
		if(dokl==defdokl)
			lcalg5=0;
		double xs=x1+h;
		double   y1m=(*f)(xs-h*a,x);
		double   y1p=(*f)(xs+h*a,x);
		double   y2m=(*f)(xs-h*b,x);
		double   y2p=(*f)(xs+h*b,x);
		double   ys=(*f)(xs,x);
		double r1=(y1m+y1p)/2-ys;
		double r2=(y2m+y2p)/2-ys;
		lcalg5+=5;
		if( fabs((r2*a*a-r1*b*b)*10*h)<dokl)     
			return h*(f0*ys+fa*(y1m+y1p)+fb*(y2m+y2p));  
		else
		{
			dokl/=2;h/=2;
			return calg5x_int(f,x,x1,x1+h,dokl)+calg5x_int(f,x,x1+h,xs,dokl)
				  +calg5x_int(f,x,xs,xs+h,dokl)+calg5x_int(f,x,xs+h,x2,dokl);    
		}  
	}
	else
	{
		double d=(x2-x1)/ile;
		dokl/=sqrt(ile);
		double suma=0;
		while(ile-->0)
		{
			suma+=calg5x_int(f,x,x1+ile*d,x1+ile*d+d,dokl,1);
		}
		return suma;
	}
}


/*
double calg5x(double (*f)(double,double),double x,double x1,double x2,int depth=0)
{double h=(x2-x1)/2;
 double xs=x1+h;
 if( depth>5)
  { double   y1m=(*f)(xs-h*0.53846931,x);
    double   y1p=(*f)(xs+h*0.53846931,x);
    double   y2m=(*f)(xs-h*0.90617985,x);
    double   y2p=(*f)(xs+h*0.90617985,x);
    double   ys=(*f)(xs,x);
   if(errno) return errno=0;
   return h*(128.0/225*ys+0.47862867*(y1m+y1p)+0.23692689*(y2m+y2p)); 
  } 
 else
   return calg5x(f,x,x1,xs,depth+1)+calg5x(f,x,xs,x2,depth+1);    
}
*/
double calg3x(double (*f)(double,double),double x,double x1,double x2,double dokl=1e-3)
{
	double static const s6=sqrt(0.6);
	double static const sqrt2=sqrt(2);
	double h=(x2-x1)/2;
	double xs=x1+h;
	double y1=(*f)(xs+h*s6,x);
	double y2=(*f)(xs-h*s6,x);
	double ys=(*f)(xs,x);
	double sigma=fabs(y1+y2-2*ys);
	if(sigma*2*h<dokl)
		return h*(5*(y1+y2)+8*ys)/9;  
	else
		return calg3x(f,x,x1,xs,dokl/sqrt2)
	 		  +calg3x(f,x,xs,x2,dokl/sqrt2); 
}

double calgx(double (*f)(double,double),double x,double x1,double x2,int depth=0)
{
	double static const sqrt3=sqrt(3);
	double h=(x2-x1)/2;
	double xs=x1+h;
	double   y1=(*f)(xs+h/sqrt3,x);
	double   y2=(*f)(xs-h/sqrt3,x);
	//double ys=(*f)(xs);
	if(errno) 
		return errno=0;
	if(   depth>5)
		return h*(y1+y2); 
	else
		return calgx(f,x,x1,xs,depth+1)+calgx(f,x,xs,x2,depth+1);    
}

//calkowanie po pierwszym argumencie
double calkx(double (*f)(double,double),double x, double x1,double x2,double y1=0,double y2=0,int depth=0)
{
	if(y1==0)     
		y1=(*f)(x1,x);
	if(y2==0)     
	 	y2=(*f)(x2,x);
	double xs=(x1+x2)/2;
	double ys=(*f)(xs,x);
	if(errno) 
		return errno=0;
	if( (fabs(2*ys-y1-y2)<=0.0000001 *(fabs(y1)+fabs(y2)+fabs(ys)) || depth>10) && depth>5)
		return (4*ys+y1+y2)*(x2-x1)/6; //metoda Simpsona (dopasowanie parabol±)
	else
		return calkx(f,x,x1,xs,y1,ys,depth+1)+calkx(f,x,xs,x2,ys,y2,depth+1);    
}
}
#endif
