#include "nucleus.h"
#include "flatnucleus.h"
#include "anynucleus.h"
#include "calg5.h"
#include "util2.h"
#include <iostream>
#include <fstream>

using namespace std;

params p;

nucleus *t;


void makenucleus(int i, int j, int kind )
{
  p.nucleus_p = i;
  p.nucleus_n = j;
  if(kind ==0)  t=new flatnucleus (p);
  if(kind ==1)  t=new anynucleus (p);
}

double r2(double x)
{ 
  return x*x;
}

double one(double x)
{
  return 1;
}

double
kfp (double x)
{
  return t->localkf (pdg_proton, x);
}

double
kfn (double x)
{
  return t->localkf (pdg_neutron, x);
}

double
density (double x)
{
  return t->density (x);
}


void
wykres (int i, int j, int k)
{
  makenucleus (i,j,k);
  double R = t->radius();
  cout << "radius("<<i<<','<<j<<")="<<R/fermi<<endl;
  stringstream s;
  s << "wyk" << t->Z() << ',' << t->N() << ".txt\0" << flush;
  ofstream f (s.str ().c_str ());
  for (double r = 0; r < R; r += 0.001 * fermi)
    {
      f << r / fermi << ' ' << t->density (r) * fermi3 << endl;
    }
  delete t;
}

void
test_random_r (int i, int j, int k)
{
  makenucleus (i,j,k);
  double R = t->radius();
  cout << "radius("<<i<<','<<j<<")="<<R/fermi<<endl;
  stringstream s;
  s << "wrr" << t->Z() << ',' << t->N() << ".txt\0" << flush;
  ofstream f (s.str ().c_str ());
  for(int i=0;i<10000;i++)
    f<< t->get_random_r()/fermi<<endl;
  delete t;
}


double 
f(double r)
{
  return t->density (r) * r * r * 4 * Pi;
}

double 
norma (int p, int n)
{
  makenucleus(p,n,1);
  double dr = 0.0001;
  double suma = 0;
  for (double r = 0; r < 12; r += dr)
    suma += f (r)*dr;

  double suma1=calg5a(f,0,8,1000);
  cout << "norma(" << p << ',' << n << ")=" << suma  << " "<<suma1<<endl;
  delete t;
  return suma ;
}

double 
mean (double (*g) (double))
{
  double dr = 0.001;
  double suma = 0;
  for (double r = 0; r < 12; r += dr)
    suma += f (r) * g (r);
  return suma * dr /  t->A();
}

void
stat (int p, int n)
{
  makenucleus(p,n,0);
  cout << "stat(" << p << ',' << n << ") " <<
    " <r>=" << sqrt (mean (r2)) <<
    " <kfp>=" << mean (kfp) <<
    " kfp(0)=" << kfp (0) <<
    " <kfn>=" << mean (kfn) <<
    " kfn(0)=" << kfn (0) << " <1>=" << mean (one) << endl;
  delete t;
}

bool d[60][60];

void
set (int i, int j)
{
  d[i][j] = 1;
}

void
show ()
{
  for (int i = 0; i < 60; i++)
    {
      for (int j = 0; j < 60; j++)
	cout << d[i][j];
      cout << endl;
    }
}

void
all (int i, int j)
{
  wykres (i, j, 1);
  norma (i, j);
  stat (i, j);
}




int znane ()
{
  p.read ("kaskada.txt");
#define funkcja(a,b) wykres(a,b,1)
//#define funkcja(a,b) test_random_r(a,b, 1)
  cout << "H0 model" << endl;
  funkcja (3, 4);
  funkcja (4, 5);
  funkcja (5, 5);
  funkcja (5, 6);
  funkcja (8, 8);
  funkcja (8, 9);
  funkcja (8, 10);
  cout << "MH0 model" << endl;
  funkcja (6, 7);
  funkcja (6, 8);
  cout << "2pf model" << endl;
  funkcja (9, 10);
  funkcja (10, 10);
  funkcja (10, 12);
  funkcja (12, 14);
  funkcja (13, 14);
  funkcja (18, 18);
  funkcja (18, 22);
  funkcja (22, 26);
  funkcja (23, 28);
  funkcja (24, 26);
  funkcja (24, 28);
  funkcja (24, 29);
  funkcja (30, 25);
  funkcja (26, 28);
  funkcja (26, 30);
  funkcja (26, 32);
  funkcja (27, 32);
  funkcja (29, 34);
  funkcja (29, 36);
  funkcja (30, 34);
  funkcja (30, 36);
  funkcja (30, 38);
  funkcja (30, 40);
  funkcja (32, 38);
  funkcja (32, 40);
  funkcja (38, 50);
  funkcja (39, 50);
  funkcja (41, 52);
  cout << "3pf model" << endl;
  funkcja (11, 3);
  funkcja (11, 4);
  funkcja (12, 12);
  funkcja (12, 13);
  funkcja (14, 14);
  funkcja (14, 15);
  funkcja (14, 16);
  funkcja (15, 16);
  funkcja (17, 18);
  funkcja (17, 20);
  funkcja (19, 20);
  funkcja (20, 20);
  funkcja (20, 28);
  funkcja (28, 30);
  funkcja (28, 32);
  funkcja (28, 33);
  funkcja (28, 34);
  funkcja (28, 36);
  cout << "3pG model" << endl;
  funkcja (16, 16);
  funkcja (40, 50);
  funkcja (40, 51);
  funkcja (40, 52);
  funkcja (40, 54);
  funkcja (40, 56);
  funkcja (42, 50);
  cout << "FB model" << endl;
  funkcja (6, 6);
  funkcja (16, 18);
  funkcja (16, 20);
  funkcja (22, 28);
  funkcja (32, 42);
  funkcja (32, 44);
  funkcja (42, 52);
  funkcja (42, 54);
  funkcja (42, 56);


}

int wykresy2(int n)
{
	for(int i=1;i<n;i++)
	  for(int j=0;j<n;j++)
	    funkcja(i,j);
}

int wykresy3(int n)
{
	for(int i=1;i<n;i++)
	  for(int j=0;j<n;j++)
	     norma(i,j);
}


int main ()
{
   wykresy3(70);
   //cout<<calg5((double (*)(double))sin,0,3.141596)<<endl;
}
