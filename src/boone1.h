#ifndef _boone_h_
#define _boone_h_

#include "event1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "util2.h"

template <class T>
double findmin( T f,double a, double b)
{ 
  double fa=f(a);
  double fb=f(b);
  double s=(a+b)/2;
  double pops=a;
  double fs= f(s);
  while(pops!=s)
  { if(fa<fb)
      { if(fs<fb)
          b=(s+b)/2;
		else   
		  //a-=(b-a);	   
		  break;
	  }
	  else  // fb<fa 
	  {if(fs<fa)
          a=(s+a)/2;
		else   
		  //b+=(b-a);	   
		  break;
	   
	  }
  fa=f(a);
  fb=f(b);
  pops=s;
  s=(a+b)/2;
//  cout<<"a="<<a<< "s="<<s<<" b=b"<<b<< " diff="<<pops-s<<endl;
  fs=f(s);   
  }
  cout<<"\t"<< s;
  return s;
}



template <class T>
T Pow2(T x)
{return x*x;
}
static bool active[20][18]=
{
    0,    0,    0,    0,    0,    0,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
    0,    0,    0,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
    0,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1};

static bool actives123123[20][18]=
{
    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1};

    
/// 2D array of the value of the double differential cross section in each bin in units of 10-41 cm2/GeV/nucleon. 
/// The muon kinetic energy increases from left to right, and the cosine of the muon scattering angle decreases from top to bottom (Table VI)
static double boonesection[20][18]=
{
  190.0,  326.5,  539.2,  901.8,   1288,   1633,   1857,   1874,   1803,   1636,   1354,   1047,  794.0,  687.9,  494.3,  372.5,  278.3,  227.4,
  401.9,  780.6,   1258,   1714,   2084,   2100,   2035,   1620,   1118,  783.6,  451.9,  239.4,  116.4,  73.07,  41.67,  36.55,      0,      0,
  553.6,  981.1,   1501,   1884,   1847,   1629,   1203,  723.8,  359.8,  156.2,  66.90,  26.87,  1.527,  19.50,      0,      0,      0,      0,
  681.9,   1222,   1546,   1738,   1365,  909.6,  526.7,  222.8,  81.65,  35.61,  11.36,  0.131,      0,      0,      0,      0,      0,      0,
  765.6,   1233,   1495,   1289,  872.2,  392.3,  157.5,  49.23,  9.241,  1.229,  4.162,      0,      0,      0,      0,      0,      0,      0,
  871.9,   1279,   1301,  989.9,  469.1,  147.4,  45.02,  12.44,  1.012,      0,      0,      0,      0,      0,      0,      0,      0,      0,
  910.2,   1157,   1054,  628.8,  231.0,  57.95,  10.69,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,
  992.3,   1148,  850.0,  394.4,  105.0,  16.96,  10.93,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,
  1007,   970.2,  547.9,  201.5,  36.51,  0.844,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,
  1003,   813.1,  404.9,  92.93,  11.63,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,
  919.3,  686.6,  272.3,  40.63,  2.176,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,
  891.8,  503.3,  134.7,  10.92,  0.071,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,
  857.5,  401.6,  79.10,  1.947,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,
  778.1,  292.1,  33.69,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,
  692.3,  202.2,  17.42,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,
  600.2,  135.2,  3.624,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,
  497.6,  85.80,  0.164,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,
  418.3,  44.84,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,
  348.7,  25.82,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,
  289.2,  15.18,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0
};

/// 2D array of the shape uncertainty of the double differential cross section in each bin in units of 10-42 cm2/GeV/nucleon. 
/// The total normalization error is 10.7% (Table VII)
static double boonesectionerror[20][18]=
{
  684.3,   1071,   1378,   1664,   1883,   2193,   2558,   3037,   3390,   3320,   3037,   3110,   2942,   2424,   2586,   2653,   3254,   3838,
  905.0,   1352,   1754,   2009,   2222,   2334,   2711,   2870,   2454,   1880,   1391,   1036,  758.7,  544.3,  505.5,  359.6,      0,      0,
   1134,   1557,   1781,   1845,   1769,   1823,   1873,   1464,  963.8,  601.6,  339.6,  184.1,  170.1,  230.6,      0,      0,      0,      0,
   1435,   1455,   1581,   1648,   1791,   1513,   1068,  598.2,  267.2,  155.1,  69.28,  89.01,      0,      0,      0,      0,      0,      0,
   1380,   1372,   1434,   1370,   1201,  870.2,  432.3,  162.2,  71.88,  49.10,  54.01,      0,      0,      0,      0,      0,      0,      0,
   1477,   1273,   1365,   1369,   1021,  475.5,  161.6,  55.58,  16.32,      0,      0,      0,      0,      0,      0,      0,      0,      0,
   1267,   1154,   1155,  965.3,  574.7,  149.2,  53.26,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,
   1293,   1105,   1041,  742.5,  250.6,  77.66,  110.3,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,
   1351,   1246,   1048,  415.1,  114.3,  41.02,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,
   1090,   1078,  695.5,  238.2,  45.96,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,
  980.4,  783.6,  515.7,  114.6,  20.92,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,
  917.7,  746.9,  337.5,  50.92,  3.422,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,
  922.7,  586.4,  215.6,  55.88,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,
  698.0,  553.3,  135.3,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,
  596.9,  482.6,  57.73,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,
  520.8,  360.7,  34.63,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,
  450.2,  236.6,  31.22,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,
  408.8,  184.4,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,
  339.7,  107.6,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,
  349.8,  63.32,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0
};

/// 2D array of the predicted CCQE-like background double differential cross section in each bin in units of 10-41 cm2/GeV/nucleon (Table VIII)
static double boonesectionbackground[20][18]=
{
  83.60,  199.8,  285.3,  364.2,  391.1,  403.7,  384.3,  349.2,  301.4,  232.7,  179.2,  136.1,  102.0,  90.73,  76.55,  52.36,  41.47,  54.50,
  111.6,  257.4,  351.0,  364.3,  353.2,  288.9,  233.8,  169.5,  106.6,  59.81,  31.21,  20.89,  10.10,  6.008,  2.376,  2.859,      0,      0,
  118.4,  270.4,  312.6,  280.3,  211.7,  135.7,  81.47,  40.97,  21.56,  9.247,  3.284,  0.875,  0.057,      0,      0,      0,      0,      0,
  118.9,  260.0,  252.8,  183.4,  101.8,  52.52,  19.75,  7.978,  2.716,  0.281,      0,      0,      0,      0,      0,      0,      0,      0,
  109.0,  215.2,  181.4,  104.6,  41.87,  16.33,  3.643,  0.492,  0.004,      0,      0,      0,      0,      0,      0,      0,      0,      0,
  109.2,  182.0,  122.4,  51.26,  19.76,  4.193,  0.183,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,
  104.0,  140.2,  73.71,  24.54,  4.613,  0.151,  0.002,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,
  93.84,  107.6,  48.56,  10.78,  0.812,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,
  76.55,  80.94,  29.02,  3.049,  0.030,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,
  67.81,  52.89,  13.71,  0.392,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,
  58.91,  37.46,  5.565,  0.011,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,
  50.47,  22.49,  1.048,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,
  39.03,  12.58,  0.118,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,
  32.41,  7.575,  0.061,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,
  25.72,  2.529,  0.080,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,
  16.78,  1.063,  0.009,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,
  9.963,  0.280,  0.002,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,
  5.005,  0.244,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,
  4.877,  0.067,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,
  3.092,  0.013,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0
};

/// 1D array of bin boundaries partitioning the reconstructed four momentum transfer, Q2
double booneq2borders[18]={0.00,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.60,0.70,0.80,1.00,1.20,1.50,2.00};

/// 1D array of the value of the single differential cross section in each bin in units of cm2/GeV2/nucleon (Table IX)
static double booneq2[17]={7.681e-39, 1.457e-38, 1.684e-38, 1.703e-38, 1.589e-38, 1.449e-38, 1.329e-38, 1.172e-38, 1.030e-38, 8.852e-39, 7.164e-39, 5.425e-39, 4.032e-39, 2.713e-39, 1.620e-39, 9.915e-40, 5.474e-40};

/// 1D array of the shape uncertainty of the single differential cross section in each bin in units of cm2/GeV2/nucleon. The total normalization error is 10.7% (Table IX)
static double booneq2error[17]={1.493e-39, 1.180e-39, 9.720e-40, 8.216e-40, 5.134e-40, 3.983e-40, 3.386e-40, 2.629e-40, 2.457e-40, 2.975e-40, 3.193e-40, 3.212e-40, 3.442e-40, 2.885e-40, 2.250e-40, 1.407e-40, 2.504e-40}; 

/// 1D array of the predicted CCQE-like background single differential cross section in each bin in units of cm2/GeV2/nucleon (Table IX)
static double booneq2background[17]={3.876e-39, 3.961e-39, 3.671e-39, 3.064e-39, 2.522e-39, 2.040e-39, 1.633e-39, 1.290e-39, 1.018e-39, 7.874e-40, 5.524e-40, 3.532e-40, 2.302e-40, 1.339e-40, 6.398e-41, 2.466e-41, 3.645e-42}; 

///////////////////////////////////////////////////////////////////////
/// class for the boone results analysis
/// 
///////////////////////////////////////////////////////////////////////
class boone
{
	double Exp[20][18];
	double ExpSum;
	long double Nuw[20][18];
	long double Nuw2[20][18];
	long double NuwSum;
	long double NuwCount;
	double Tmin,Tmax;
	double Cosmin,Cosmax;
	double binarea;
	
	double Expq2[17];	
	double Expq2Sum;	
	long double Nuwq2[17];	
	long double Nuw2q2[17];	
	long double Nuwq2Sum;	
	long double Nuwq2Count;	
	long double Nuwmcq2[17];	
	long double Nuw2mcq2[17];	
	long double Nuwmcq2Sum;	
	long double Nuwmcq2Count;	
	long double Total;
	long double Count;
	double binlen[17];    
	double sigma;
	public:
	TH1D  h1,h1mc,b1;
	TH2D  h2,h2lowq,b2;
  
  
public:
/*
double findmin1(double (boone::*f)(double&),double a, double b)
{ 
  double fa= (this->*f)(a);
  double fb=(this->*f)(b);
  double s=(a+b)/2;
  double pops=a;
  double fs=(this->*f)(s);
  while(pops!=s)
  { if(fa<fb)
      { if(fs<fb)
          b=(s+b)/2;
		else   
		  //a-=(b-a);	   
		  break;
	  }
	  else  // fb<fa 
	  {if(fs<fa)
          a=(s+a)/2;
		else   
		  //b+=(b-a);	   
		  break;
	   
	  }
  fa=(this->*f)(a);
  fb=(this->*f)(b);
  pops=s;
  s=(a+b)/2;
//  cout<<a<< ' '<<s<<' '<<b<< ' '<<pops-s<<endl;
  fs=(this->*f)(s);   
  }
  cout<<"\t"<< s;
  return s;
  	
}
*/
  
    boone():
		h1("h1","dsigma/dq2 (nuwro Q2rec)", 17,booneq2borders),
		h1mc("h1mc","dsigma/dq2 (nuwro Q2mc)", 17,booneq2borders),
		b1("b1","dsigma/dq2 (miboone)", 17,booneq2borders),
		h2("h2","dsigma/dcos theta dT (nuwro)", 20,-1,1, 18, 0.2,2),
		h2lowq("h2lowq","dsigma/dcos theta dT (nuwro |q|<350)", 20,-1,1, 18, 0.2,2),
		b2("b2","dsigma/dcos theta dT (miboone)", 20,-1,1, 18, 0.2,2),
        Cosmin(-1),
        Cosmax(1),
        Tmin(0.2),
        Tmax(2),
	    binarea((Tmax-Tmin)*(Cosmax-Cosmin)/(20*18)),
	    sigma(0.107)
        {reset();}
        
      double chi2(long double x[],long double x2[],double y[],double ys[], int n)
       {  double sumx=0,sumy=0,res=0;
	      for(int i=0;i<n;i++)
	      {sumx+=x[i];
	       sumy+=y[i];
          }
			for(int i=0;i<n;i++)
			{
				res+=Pow2(n*x[i]/sumx -n*y[i]/sumy)
			      /( (x2[i]-x[i]*x[i]/Count)/Count*n*n/sumx/sumx + ys[i]*ys[i]*n*n/sumy/sumy);
			     }
		    return res;
        }
  
	inline void reset();
	inline void addval(int i, int j, double x);
	inline void addCosT(double Cos, double T, double x);
	inline void addq2(double q2, double x, int mc=0); // 0 - reconstructed; 1 - mc 
	inline double wdiff2(double &mul);
    inline double wdiff2g(double &mul);
	inline double wdiff2q2(double &mul);
	inline double wdiff2q2g(double &mul);
	inline void add(event &e, double w);
	inline long double avg(){return Total/Count;}
	inline long double total(){return Total;}
	inline long double count(){return Count;}
	inline void printtable(ostream& o, double mul=0);
	inline void printq2table(ostream& o, double mul=0);
    inline void printmcq2table(ostream& o, double mul=0);
	inline void printdiffs(ostream& o, double mul=0);
	inline double exp(int i, int j){ return boonesection[i][j]*1e-41;}
	inline double err(int i, int j){ return boonesectionerror[i][j]*1e-42;}
	inline double nuw(int i, int j)    { return Nuw[i][j]/Count/binarea;}
	inline double expq2(int i)     { return booneq2[i];}
	inline double errq2(int i)     { return booneq2error[i];}
	inline double nuwq2(int i) 	   { return Nuwq2[i]/binlen[i]/Count;}
	inline double nuwmcq2(int i)    { return Nuwmcq2[i]/binlen[i]/Count;}
	inline double scale(){return ExpSum*1e-41*binarea/(NuwSum/Count);}
    inline double rootchi2();
	inline double scaleij()
	{ double a=0,b=0;
      for(int i=0;i<20;i++)
        for(int j=0;j<18;j++)
          {a+=exp(i,j);
           b+=nuw(i,j);
	      }
      return a/b;
    }
	
	inline double scaleq2(){return Expq2Sum/(Nuwq2Sum/Count);}
	inline double scalemcq2(){return Expq2Sum/(Nuwmcq2Sum/Count);}
	inline double scaleq2ij()
	{ double a=0,b=0;
      for(int i=0;i<17;i++)
          {a+=expq2(i)*binlen[i];
           b+=nuwq2(i)*binlen[i];
          }
      return a/b;
    }
    inline double avg9(int i ,int j);
    inline void print();
    inline void save(const char* filename);
};

void boone::reset()
  {
   ExpSum=0;
   NuwSum=NuwCount=Count=Total=0;
   Expq2Sum=0;
   Nuwq2Sum=Nuwq2Count=0;
   Nuwmcq2Sum=Nuwmcq2Count=0;
   
   for(int i=0;i<20;i++)
    for(int j=0;j<18;j++)
       {
        Nuw[i][j]=Nuw2[i][j]=0;
    ExpSum+=Exp[i][j]=boonesection[i][j];//-boonesectionbackground[i][j];
        b2.SetBinContent(i+1,j+1,boonesection[i][j]*1e-41);
        b2.SetBinError(i+1,j+1,boonesectionerror[i][j]*1e-42);
	   }
   for(int i=0;i<17;i++)
     {Nuwq2[i]=0;Nuwmcq2[i]=0;
      Nuw2q2[i]=0;Nuw2mcq2[i]=0;
      binlen[i]=booneq2borders[i+1]-booneq2borders[i];
      Expq2Sum+=expq2(i)*binlen[i];
        b1.SetBinContent(i+1,booneq2[i]);
        b1.SetBinError(i+1,booneq2error[i]);
     }
//     h1.Sumw2();
//     h1mc.Sumw2();
//     h2.Sumw2();
  }	   
  
void boone::addval(int i, int j, double x)
  { if(i>=0 and i<20 and j>=0 and j<18) 
       {Nuw[i][j]+=x;
        Nuw2[i][j]+=x*x;
        NuwSum+=x;
        NuwCount++;
	   }
  }

void boone::addCosT(double Cos, double T, double x)
   {int i= (Cosmax-Cos)/(Cosmax-Cosmin)*20;
    if(Cos>Cosmax) i=-1;
    int j =(T/GeV-Tmin)/(Tmax-Tmin)*18;
    if(T/GeV<Tmin) j=-1;
    addval(i,j,x);
   }

void boone::addq2(double q2, double x, int mc)
   {  
	 if(q2<2 and q2 >0)
     {int i=0;
      while(booneq2borders[i+1]<q2)
	      i++;
      if(mc)
      {
	      Nuwmcq2[i]+=x;
	      Nuw2mcq2[i]+=x*x;
	      Nuwmcq2Sum+=x;
	      Nuwmcq2Count++;
      }
      else
      {
	      Nuwq2[i]+=x;
	      Nuw2q2[i]+=x*x;
	      Nuwq2Sum+=x;
	      Nuwq2Count++;
      }
      if(Nuwq2Sum==Nuwq2Sum and Nuwq2[i]==Nuwq2[i])
        ;
      else
       
        {cout<<q2<<' '<<x<<endl;
         throw("nan");
        } 
        
     }
   }

double boone::wdiff2(double &mul)
  {
    int dof=-1;
    if(mul==0)
      mul=scale();
    if(mul==-1)
      mul=findmin(makefun(*this,&boone::wdiff2),0.1,10);
    double res=0;
    for(int i=0;i<20;i++)
    for(int j=0;j<18;j++)
      if(active[i][j])
      {
       	double y=exp(i,j)-nuw(i,j)*mul;
        if(err(i,j)!=0 )
          { y/=err(i,j);
            res+=y*y;
              dof++;
           } 
          else if(y!=0)
            { double a=avg9(i,j);
              if(a!=0)
              {
              y/=a; 
              res+=y*y; 
              dof++;
		      }
	       }
	      
	   }
    cout<<dof<<endl;          
	  return res;   
  }

double boone::wdiff2g(double &mul)
  {int dof=-1;
   if(mul==0)
     mul=scale();
   if(mul==-1)
     mul=findmin(makefun(*this,&boone::wdiff2g),0.1,10);
   double res=0;
   for(int i=0;i<20;i++)
    for(int j=0;j<18;j++)
       if(active[i][j])
       {
       	double y=exp(i,j)-nuw(i,j)*mul;
        if(err(i,j)!=0 )
          { y/=err(i,j);
            res+=y*y;
            dof++;
           } 
          else if(y!=0)
            { double a=avg9(i,j);
              if(a!=0)
              {
              y/=a; 
              res+=y*y; 
              dof++;
		      }
	       }
//	     if(y!=0)
//	     cout<<"y="<<y<<endl;
	   }
              dof++;
//    cout<<"dof="<<dof<<endl;          
	return res+(mul-1)*(mul-1)/(sigma*sigma);   
   }

double boone::wdiff2q2(double &mul)
  {int dof=-1;
   if(mul==0)
     mul=scaleq2();
   if(mul==-1)
     mul=findmin(makefun(*this,&boone::wdiff2q2),0.1,10);
	 double res=0;
	
	 for(int i=0;i<17;i++)
      {
       	double a=expq2(i);
        double b=nuwq2(i)*mul;
       	double y=a-b;
        if(errq2(i) >0 )
           {y/=errq2(i);
             res+=y*y;
              dof++;
             cout<<"y="<<y<<endl;
	        }
	  }
	   
//    cout<<dof<<endl;          
	 return res;   
   }

double boone::wdiff2q2g(double &mul)
  {int dof=-1;
   if(mul==0)
     mul=scaleq2();
   if(mul==-1)
     mul=findmin(makefun(*this,&boone::wdiff2q2g),0.1,10);
	 double res=0;
	
	 for(int i=0;i<17;i++)
      {
       	double a=expq2(i);
        double b=nuwq2(i)*mul;
       	double y=a-b;
        if(errq2(i) >0 )
           {y/=errq2(i);
             res+=y*y;
              dof++;
             cout<<"y="<<y<<endl;
	        }
	  }
	   
              dof++;
//    cout<<dof<<endl;          
	return res+(mul-1)*(mul-1)/(sigma*sigma);   
   }


   
void boone::add(event &e,double w)
   { Total=Total+w;
	 Count=Count+1;
	 if(e.weight>0)
      {
	   double Enu=e.in[0].E();
	   double Emu=e.out[0].E();
	   double Mmu=e.out[0].mass();
	   double Tmu=e.out[0].Ek();
	   double Cos=e.in[0].p().dir()*e.out[0].p().dir();	  
	   addCosT(Cos, Tmu, w);              /////
	   h2.Fill(-Cos, Tmu/GeV, w/binarea); /////
	   vect q=e.in[0]-e.out[0];
	   if (vec(q).length()<350*MeV)
  	     h2lowq.Fill(-Cos, Tmu/GeV, w/binarea); /////
	   {
	   double Q2mc=-q*q/GeV/GeV;
	   addq2(Q2mc,w,1);            /////
	   
	   int bin=h1mc.FindBin(Q2mc);
	   double width=h1mc.GetBinWidth(bin);
	   
       h1mc.Fill(Q2mc,w/width);   ////
	   }
	   
	   double Mn=e.in[1].mass();
	   double Mp=e.out[1].mass();
	   double B=34*MeV;
	   double Mnc=Mn-B;
	   double recEnu= 
	   (2*Mnc*Emu-(Mnc*Mnc+Mmu*Mmu-Mp*Mp))/
	   (2*(Mnc-Emu+sqrt(Emu*Emu-Mmu*Mmu  )*Cos));
	   	
	   	   
	   double q2rec=-Mmu*Mmu+2 *recEnu*
	                (Emu-sqrt(Emu*Emu-Mmu*Mmu)*Cos);
	   addq2(q2rec/GeV/GeV, w);    ////
	   
	   int bin=h1.FindBin(q2rec/GeV/GeV);
	   double width=h1.GetBinWidth(bin);            
	   
	   h1.Fill(q2rec/GeV/GeV,w/width);  /////
	  }
   }

void boone::printtable(ostream& o, double mul)
   {if(mul==0)
     mul=scale();
    for(int i=0;i<20;i++)
    {  
		for(int j=0;j<18;j++)
		   {
       	   o<<20*i+j <<" "<<setw(8)<<nuw(i,j)<<" "<<setw(8)<<nuw(i,j)*mul<<" "<<setw(8)<<exp(i,j)<<endl;
	       }
	 }  
   }


void boone::printq2table(ostream& o, double mul)
   { if(mul==0)
       mul=scaleq2();
		for(int j=0;j<17;j++)
		   {o<<setw(6)<<(booneq2borders[j]+booneq2borders[j+1])/2<<"  ";
       	    o<<setw(6)<<nuwq2(j)<<"  ";
       	    o<<setw(6)<<nuwq2(j)*mul<<"  ";
       	    o<<setw(6)<<expq2(j)<<"  ";
	        o<<endl;
	       }
       
   }

void boone::printmcq2table(ostream& o, double mul)
   { if(mul==0)
       mul=scalemcq2();
		for(int j=0;j<17;j++)
		   {o<<setw(6)<<(booneq2borders[j]+booneq2borders[j+1])/2<<"  ";
       	    o<<setw(6)<<nuwmcq2(j)<<"  ";
       	    o<<setw(6)<<nuwmcq2(j)*mul<<"  ";
       	    o<<setw(6)<<expq2(j)<<"  ";
	        o<<endl;
	       }
       
   }


void boone::printdiffs(ostream& o, double mul)
  {
    if(mul==0)
      mul=scale();
    for(int i=0;i<20;i++)
    {   
      o<<setw(2)<<i<<": ";
		  for(int j=0;j<18;j++)
		  {
        o<<setw(8)<<nuw(i,j)*mul-exp(i,j)<<"  ";
	    }
      o<<endl;   
    } 
  }



double boone::avg9(int i ,int j)
   { int n=0;
     double s=0;
   	 for(int k=max(i-1,0);k<min(i+1,20);k++)
    	for(int l=max(j-1,0);l<min(j+1,18);l++)
    	  { n++;
    	  	s+=boonesectionerror[i][j];
     	  }
     return s/n;
   }

void boone::print()
{ 
  for(int i=0;i<20;i++)
  for(int j=0;j<18;j++)
    {
    	cout<< i<<"  "<<j<< "  ========================================="<<endl;
    	cout<< h2.GetBinContent(i+1,j+1)<<endl;
    	cout<< Nuw[i][j]<<endl;
  	}
  for(int i=0;i<17;i++)
    {
    	cout<< i<< "  ========================================="<<endl;
    	cout<< h1.GetBinContent(i+1)<<endl;
    	cout<< Nuwq2[i]/binlen[i]<<endl;
    	cout<< h1mc.GetBinContent(i+1)<<endl;
    	cout<< Nuwmcq2[i]/binlen[i]<<endl;
    }
}



void boone::save(const char* filename)
{     
      TFile tf(filename,"recreate");
      double s=1.0/Count; 
      h1.Scale(s); 
      h1mc.Scale(s); 
      h2.Scale(s);
      h2lowq.Scale(s);
      cout<<endl<<"Procent zdarzeń z niskim przekazem pędu (<350MeV)"<<endl ;    
      for(int i=1;i<=20;i++) 
      {
        for(int j=1;j<=18;j++)     
        {
          double x=h2.GetBinContent(i,j);
          double y=h2lowq.GetBinContent(i,j);
          int z=(x>0)? 100*y/x : 0;
          cout<<setw(5)<<z;
	      }
	      cout<<endl;    
	    }
      h1.Write(); 
      h1mc.Write(); 
      h2.Write();
      h2lowq.Write();
      b1.Write(); 
      b2.Write();
      
      tf.Write();
      tf.Close();
      
//      cout<<"1 dim Chi2="<<h1.Chi2Test(&b1,"UU")<<endl;
//      cout<<"1 dim Chi2="<<h1.Chi2Test(&b1,"UW")<<endl;
//      cout<<"1 dim Chi2="<<h1.Chi2Test(&b1,"WW")<<endl;
//      cout<<"1 dim Chi2="<<h1.Chi2Test(&b1,"P")<<endl;
//      cout<<"1 dim Chi2="<<h1.Chi2Test(&b1,"CHI2")<<endl;
//      cout<<"1 dim Chi2="<<h1.Chi2Test(&b1,"CHI2/NDF")<<endl;
//      cout<<"2 dim Chi2="<<h2.Chi2Test(&b2,"UU")<<endl;
//      cout<<"2 dim Chi2="<<h2.Chi2Test(&b2,"UW")<<endl;
//      cout<<"2 dim Chi2="<<h2.Chi2Test(&b2,"WW")<<endl;
//      cout<<"2 dim Chi2="<<h2.Chi2Test(&b2,"P")<<endl;
//      cout<<"2 dim Chi2="<<h2.Chi2Test(&b2,"CHI2")<<endl;
//      cout<<"2 dim Chi2="<<h2.Chi2Test(&b2,"CHI2/NDF")<<endl;      
  
  //cout<<"chi2 2="<<chi2(Nuw[0],Nuw2[0],boonesection[0],boonesectionerror[0],20*18)<<endl;
  //cout<<"chi2 1="<<chi2(Nuwq2,Nuw2q2,booneq2,booneq2error,17)<<endl;
  //cout<<"chi2 1mc="<<chi2(Nuwmcq2,Nuw2mcq2,booneq2,booneq2error,17)<<endl;
}




#endif
