#include "nuwro.h"

inline double sqr(vect a)
{return a.x*a.x+a.y*a.y+a.z*a.z+a.t*a.t;
}

inline void test_boost()
{
	init_genrand(time(0));

    double e[5]={0,0,0,0,0};
    double E[5]={0,0,0,0,0};
	for(int i=0;i<100000;i++)
	{
		vect a(80,50,1,-30*frandom()),a1=a,a2=a,a3=a,a4=a;
		vec v=vec(frandom(),frandom(),frandom())/2;
		//cout <<a <<endl;
		//cout <<a*a <<endl;
		a1.boost1(v);
		a3.boost3(v);
		a4.boost4(v);
		a1.boost1(-v);
		a3.boost3(-v);
		a4.boost4(-v);
        e[1]+=sqr(a1-a);
        E[1]+=abs(a1*a1-a*a);
        e[3]+=sqr(a3-a);
        E[3]+=abs(a3*a3-a*a);
        e[4]+=sqr(a4-a);
        E[4]+=abs(a4*a4-a*a);
//		cout <<a1 <<' '<<(e[1]+=sqr(a1-a))<< (E[1]+=abs(a1*a1-a*a))<<endl;
//		cout <<a3<<' '<<(e[3]+=sqr(a3-a)) <<(E[3]+=abs(a3*a3-a*a))<<endl;
//		cout <<a4<<' '<<(e[4]+=sqr(a4-a)) <<(E[4]+=abs(a4*a4-a*a))<<endl;
	}
	for(int i=0;i<5;i++)
      cout<<"boost"<<i<<" errors = "<<e[i]<<' '<<E[i]<<endl;
}

inline void test_ort()
{ double x=0,y=0;
  frandom_init(0);
	for(int i=0;i<200000;i++)
	{  
		vec a(frandom(),frandom(),frandom());
		a*=2;
		a-=vec(1,1,1);
		vec b=rand_ort(a);
//		vec c=rand_ort2(a);
//		cout<<a*a<<' '<<b*b<<' '<<a*b<<endl;
//		cout<<a*a<<' '<<c*c<<' '<<a*c<<endl;
		x+=abs(a*b);
//		y+=abs(a*c);
    }
    cout<< x <<' '<< y<<endl;
    
    
    int n=3000;
    n*=n;
    double s0=time(0),s1,s2;
    for(int i=1;i<n;i++)
    {
    	rand_ort(vec(1,3,4));
    	
	}
	s1=time(0);
    for(int i=1;i<n;i++)
    {
//    	rand_ort2(vec(1,3,4));
    	
	}
	s2=time(0);
    cout<<s1-s0<<' '<<s2-s1<<endl;
    cout<<sin(Pi)<<endl;
    cout<<sin(2*Pi)<<endl;
}
