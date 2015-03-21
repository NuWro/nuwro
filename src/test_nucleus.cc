#include "nucleusmaker.h"
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
	p.nucleus_model=kind;
	t=make_nucleus(p);
}

double r2(double x)
{ 
	return x*x;
}

double one(double x)
{
	return 1;
}

double kfp (double x)
{
	return t->localkf_ (pdg_proton, x);
}

double kfn (double x)
{
	return t->localkf_ (pdg_neutron, x);
}

double density (double x)
{
	return t->density (x);
}


void wykres (int i, int j, int k)
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

void test_random_r (int i, int j, int k)
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


double f(double r)
{
  return t->density (r) * r * r * 4 * Pi;
}

double norma (int p, int n,int kind=1)
{
  makenucleus(p,n,kind);

  double suma1=calg5a(f,0,t->radius(),100);
  cout << "(" << p << ',' << n << ") r="<<t->radius()/fermi<<"  norm/A=" << suma1/(p+n)<<endl;
  delete t;
  return suma1;
}

double mean (double (*g) (double))
{
  double dr = 0.001;
  double suma = 0;
  for (double r = 0; r < 12; r += dr)
    suma += f (r) * g (r);
  return suma * dr /  t->A();
}

void stat (int p, int n,int k=0)
{
  makenucleus(p,n,k);
  cout << "stat(" << p << ',' << n << ") " <<
    " <r>=" << sqrt (mean (r2))/fermi <<
    " <kfp>=" << mean (kfp) <<
    " kfp(0)=" << kfp (0) <<
    " <kfn>=" << mean (kfn) <<
    " kfn(0)=" << kfn (0) << " <1>=" << mean (one) <<
    " Eb=" <<  t->Eb() <<" MeV "<< 
    " Ef=" <<  t->Ef() <<" MeV ("<< 
    t->i->Z<<','<<t->i->N<<')'<< endl;
  delete t;
}

void all (int i, int j, int kind)
{
	wykres (i, j, kind);
	norma (i, j, kind);
	stat (i, j, kind);
}

template <class F>
void znane (F funkcja,int kind)
{
  p.read ("kaskada.txt");
  cout << "H0 model" << endl;
  funkcja (3, 4,kind);
  funkcja (4, 5,kind);
  funkcja (5, 5,kind);
  funkcja (5, 6,kind);
  funkcja (8, 8,kind);
  funkcja (8, 9,kind);
  funkcja (8, 10,kind);
  cout << "MH0 model" << endl;
  funkcja (6, 7,kind);
  funkcja (6, 8,kind);
  cout << "2pf model" << endl;
  funkcja (9, 10,kind);
  funkcja (10, 10,kind);
  funkcja (10, 12,kind);
  funkcja (12, 14,kind);
  funkcja (13, 14,kind);
  funkcja (18, 18,kind);
  funkcja (18, 22,kind);
  funkcja (22, 26,kind);
  funkcja (23, 28,kind);
  funkcja (24, 26,kind);
  funkcja (24, 28,kind);
  funkcja (24, 29,kind);
  funkcja (30, 25,kind);
  funkcja (26, 28,kind);
  funkcja (26, 30,kind);
  funkcja (26, 32,kind);
  funkcja (27, 32,kind);
  funkcja (29, 34,kind);
  funkcja (29, 36,kind);
  funkcja (30, 34,kind);
  funkcja (30, 36,kind);
  funkcja (30, 38,kind);
  funkcja (30, 40,kind);
  funkcja (32, 38,kind);
  funkcja (32, 40,kind);
  funkcja (38, 50,kind);
  funkcja (39, 50,kind);
  funkcja (41, 52,kind);
  cout << "3pf model" << endl;
  funkcja (11, 3,kind);
  funkcja (11, 4,kind);
  funkcja (12, 12,kind);
  funkcja (12, 13,kind);
  funkcja (14, 14,kind);
  funkcja (14, 15,kind);
  funkcja (14, 16,kind);
  funkcja (15, 16,kind);
  funkcja (17, 18,kind);
  funkcja (17, 20,kind);
  funkcja (19, 20,kind);
  funkcja (20, 20,kind);
  funkcja (20, 28,kind);
  funkcja (28, 30,kind);
  funkcja (28, 32,kind);
  funkcja (28, 33,kind);
  funkcja (28, 34,kind);
  funkcja (28, 36,kind);
  cout << "3pG model" << endl;
  funkcja (16, 16,kind);
  funkcja (40, 50,kind);
  funkcja (40, 51,kind);
  funkcja (40, 52,kind);
  funkcja (40, 54,kind);
  funkcja (40, 56,kind);
  funkcja (42, 50,kind);
  cout << "FB model" << endl;
  funkcja (6, 6,kind);
  funkcja (16, 18,kind);
  funkcja (16, 20,kind);
  funkcja (22, 28,kind);
  funkcja (32, 42,kind);
  funkcja (32, 44,kind);
  funkcja (42, 52,kind);
  funkcja (42, 54,kind);
  funkcja (42, 56,kind);
}

template <class F >
void wykresy2(F funkcja, int n, int kind)
{
	for(int i=1;i<n;i++)
	  for(int j=0;j<n;j++)
	    funkcja(i,j,kind);
}

template <class F >
void znane2(F funkcja, int kind)
{
	for(int i=0;dens_data[i].p()>0;i++)
	{
	    funkcja(dens_data[i].p(),dens_data[i].n(),kind);
	}
}



int main ()
{
//	wykresy2(norma,70,1);
	znane2(stat	,1);
}
