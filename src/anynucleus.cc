#include "anynucleus.h"
#include "particle.h"
#include <algorithm>  
#include <map>  
   
using namespace std;

static  int best_pn(int p, int n);

static inline double bessj0(double x)
{
 if(fabs(x)<0.000001) 
   return cos(x);
 else
   return sin(x)/x;
}

template < class T > static inline T
pow2 (T x)
{
  return x * x;
}


///////////////////////////////////////////////////////////////////////////////
/// Construct nucleus from params given in the input file 
///////////////////////////////////////////////////////////////////////////////
anynucleus::anynucleus(params &par):
		nucleus(par)
	{ 
		using namespace PDG;
		if(p+n==1) Eb=kf=local_kf=0;
		pr=p;
		nr=n;	 
		double A = n+p;			
		double Ap = pow(A, 1.0 / 3);
		_V = sqrt(mass_proton*mass_proton + kf*kf) - mass_proton;// + 7*MeV;// + 5*MeV;
		//potencjal, wartosc do ustalenia
		_r=calculate_radius();
	}
	
///////////////////////////////////////////////////////////////////////////////
/// calculate radius  
///////////////////////////////////////////////////////////////////////////////
double anynucleus::calculate_radius()
{ 
  static map<int,double> res; ///< Memorizing
  map<int,double>::iterator it;
  int key=p*1000+n;
  it=res.find(key);  
  if(it!=res.end())
    return it->second;///<return memorized result
  
  double r0=0, r1=20.0*fermi,rs=0;
  double d0=0.0001*density(r0);
  double oldrs;
  do
  { oldrs=rs;
    rs=(r0+r1)/2;
    double ds=density(rs);
    if(ds>d0)
      r0=rs;
    else
      r1=rs;
//    cout<<r0/fermi<<" "<< rs/fermi<< " " <<r1/fermi<<" "<< ds*fermi3<<" "<<d0*fermi3<<endl;
  }
  while(rs!=oldrs); 

//return rs;

 return res[key]=rs; ///<memorize and return
}



static double MI2(double rms2, double q,double r)
{ double qr=q*r;
double r2=r*r;
double cosqr=cos(qr);
double sinqr=sin(qr);
return (((((cosqr*qr/6-sinqr/2)*qr-cosqr)*qr+sinqr)*rms2/r2 + sinqr)/r - cosqr*q)/r2;
}

static double MI(double rms2, double q,double r)
{ double qr=q*r;
double r2=r*r;
double r3=r2*r;
double cosqr=cos(qr);
double sinqr=sin(qr);
return (-cosqr*qr+sinqr)/r3
    //   + (((cosqr*qr/6-sinqr/2)*qr -cosqr)*qr+sinqr)*rms2/r2/r3;
    +(cosqr*qr*qr-2*cosqr-2*sinqr*qr)/r2/r2/6; 
}


static double prfMI(double rms,double q1,double q2, double r)
{ static double const twoPi2=2*Pi*Pi;
if(r==0) r=0.00001; 
rms*=rms; 
return max(0.0, (MI(rms,q2,r)-MI(rms,q1,r))/twoPi2);
}


///density profile for the HO model 
static double prfH0(double a,double alfa,double rho0, double r)
{double x=r/a,x2=x*x; 
return rho0*(1+alfa*x2)*exp(-x2);
}

///density profile for the MHO model
static double prfMH0(double a,double alfa,double rho0, double r)
{double x=r/a,x2=x*x; 
return (rho0*(1+alfa*x2+alfa)*exp(-x2));
}

///density profile for the 2pF model
static double prf2pf(double a,double alfa,double rho0, double r)
{//double x=r/a,x2=x*x;
return rho0/(1+exp((r-a)/alfa));
}

///density profile for the 3pF model
static double prf3pf(double a,double alfa, double w, double rho0, double r)
{double x=r/a,x2=x*x;
return max(0.,rho0*(1+w*x2)/(1+exp(r-a)/alfa)); 
}

///density profile for the 3pG model
static double prf3pg(double a,double alfa, double w, double rho0, double r)
{double x=r/a, x2=x*x;
return rho0*(1+w*x2)/(1+exp((r*r-a*a)/(alfa*alfa)));
}

///density profile for the FB model
static double prfFB(double fb[],  double r)
{ double R=fb[0];
if(r>R)
    return 0;
  double suma=0; 
  double x=Pi*r/R;
  for(int j=1; j<=17 and fb[j] ; j++)
    {   suma+=fb[j]*bessj0(j*x);
    }	
  return max(0.,suma);	
}


///density profile for the SOG model
static double prfSOG(double sog[],double r, int A)
{ static const double twoPi32=2*Pi*sqrt(Pi);
  double rms=sog[0];
  double g=sog[1]/sqrt(1.5);
  double coef=twoPi32*g*g*g;
  double suma=0;
  for(int j=2; j<=26 ; j+=2)
    {  double Ri=sog[j];
       double Qi=sog[j+1];
       double Ai=Qi/(coef*(1+2*pow2(Ri/g)));

       suma+=Ai*(exp(-pow2((r-Ri)/g))+exp(-pow2((r+Ri)/g)));
    }	
  return max(0.,suma*A);	
}



static double density(double r,int p, int n)
{  
	// He3
	static double fb201[18]={5, 0.20020e-1, 0.41934e-1, 0.36254e-1, 0.17941e-1, 0.46608e-2,
								0.46834e-2, 0.52042e-2, 0.38280e-2, 0.25661e-2, 0.14182e-2,
								0.61390e-3, 0.22929e-3};

	//N15
	static double fb708[18]={7, 0.25491e-1, 0.50618e-1, 0.29822e-1, -0.55196e-2, -0.15913e-1, 
							   -0.76184e-2,-0.23992e-2,-0.47940e-3, 0, 0, 0, 0, 0, 0, 0, 0, 
				  0};

	//C12
	static double fb606[18]={8, 0.15721e-1, 0.38732e-1, 0.36808e-1, 0.14671e-1,-0.43277e-2, 
							   -0.97752e-2,-0.68908e-2,-0.27631e-2,-0.63568e-3, 0.71809e-5,
								0.18441e-3, 0.75066e-4, 0.51069e-4, 0.14308e-4, 0.23170e-5,
								0.68465e-6, 0};
	//
	static double fb1618[18]={8,0.37036e-1, 0.58506e-1, 0.12082e-1,-0.19022e-1,-0.83421e-2,
								0.45434e-2, 0.28346e-2,-0.52304e-3, 0.27754e-4, 0.59403e-4,
							   -0.42794e-4, 0.20407e-4,-0.79934e-5, 0.27354e-5,-0.83914e-6,
								0, 0};

	static double fb1620[18]={8,0.37032e-1, 0.57939e-1, 0.10049e-1,-0.19852e-1,-0.67176e-2,
								0.61882e-2, 0.37795e-2,-0.55272e-3,-0.12904e-3, 0.15845e-3,
							   -0.84063e-4, 0.34010e-4,-0.11663e-4, 0.35204e-5, 0.95135e-6,
								0, 0};

	static double fb2228[18]={9.5, 0.31818e-1, 0.58556e-1, 0.19637e-1,-0.24309e-1,-0.18748e-1,
								   0.33741e-2, 0.89961e-2, 0.37954e-2,-0.41238e-3, 0.12540e-3,
								   0, 0, 0, 0, 0, 
								   0, 0};
	//Cr54
	static double fb2430[18]={9.0, 0.39002e-1, 0.60305e-1, 0.45845e-2,-0.30723e-1,-0.91355e-2,
					   0.93251e-2, 0.60583e-2,-0.15602e-2,-0.76809e-3, 0.76809e-3, 
					   -0.34804e-3, 0, 0,0,0,0,0};


	static double fb3242[18]={10,  0.37989e-1, 0.58298e-1, 0.27406e-2,-0.30666e-1,-0.81505e-2,
								   0.10231e-1, 0.49382e-2,-0.16270e-2,-0.13937e-2, 0.15476e-3,
								   0.14396e-3,-0.73075e-4, 0.31998e-4,-0.12822e-4, 0.48406e-5,
								   0, 0};

	static double fb3244[18]={10,  0.37951e-1, 0.57876e-1, 0.15303e-2,-0.31822e-1,-0.76875e-2,
								   0.11237e-1, 0.50780e-2,-0.17293e-2,-0.15523e-2, 0.72439e-4,
								   0.16560e-3,-0.86631e-4, 0.39159e-4,-0.16259e-4, 0.63681e-5,
								   0, 0};

	static double fb4252[18]={12,  0.30661e-1, 0.58828e-1, 0.20396e-1,-0.28830e-1,-0.25077e-1, 
								   0.44768e-2, 0.13127e-1, 0.19548e-2,-0.61403e-2,-0.35825e-2, 
								   0.73790e-3, 0.61882e-3,-0.40556e-3,-0.55748e-5,-0.12453e-3, 
								  -0.57812e-4,-0.21657e-4};

	static double fb4254[18]={12,  0.30564e-1, 0.58013e-1, 0.19255e-1,-0.28372e-1,-0.23304e-1, 
								   0.49894e-2, 0.12126e-1, 0.10496e-2,-0.62592e-2,-0.32814e-2, 
								   0.89668e-3, 0.50636e-3,-0.43412e-3, 0.71531e-4, 0.76745e-4,
								  -0.54316e-4, 0.23386e-6};

	static double fb4256[18]={12,  0.30483e-1, 0.57207e-1, 0.17888e-1,-0.28388e-1,-0.21778e-1, 
								   0.56780e-2, 0.11236e-1, 0.82176e-3,-0.50390e-2,-0.23877e-2, 
								   0.71492e-3, 0.29839e-3,-0.31408e-3, 0.80177e-3, 0.43682e-4,
								  -0.51394e-4, 0.22293e-4};
	//Mo100
	static double fb4258[18]={12,  0.30353e-1, 0.56087e-1, 0.16057e-1,-0.28767e-1,-0.20683e-1,
					   0.62429e-2, 0.11058e-1, 0.11502e-2,-0.39395e-2,-0.14978e-2, 
					   0.76350e-3, 0.10554e-3,-0.25658e-3, 0.10964e-3, 0.10015e-4,
					  -0.40341e-4, 0.25744e-4};
	//Pd104
	static double fb4658[18]={11,  0.41210e-1, 0.62846e-1,-0.21202e-2,-0.38359e-1,-0.44693e-2,
					   0.16656e-1, 0.36873e-2,-0.57534e-2,-0.32499e-2, 0.69844e-3,
					   0.16304e-2, 0.59882e-3, 0, 0, 0,
					   0, 0};
	//Pd106
	static double fb4660[18]={11,  0.41056e-1, 0.61757e-1,-0.29891e-2,-0.37356e-1,-0.35348e-2,
					   0.16085e-1, 0.28502e-2,-0.55764e-2,-0.15433e-2, 0.22281e-2,
					   0.13160e-2, 0.16508e-4, 0, 0, 0,
					   0,0};
	//Pd108
	static double fb4662[18]={11,  0.40754e-1, 0.59460e-1,-0.54077e-2,-0.36305e-1,-0.21987e-2,
					   0.15418e-1, 0.25927e-2,-0.52781e-2,-0.19757e-2, 0.10339e-2,
					   0.22891e-3,-0.33464e-3,0, 0, 0,
					   0, 0};

	// Sum Of Gausians  { rms, RP, R1, Q1,R2, Q2, ...,R12 ,Q12}
	// H3
	static double sog102[26]={1.764,0.8,
							  0.0,0.035952,
							  0.2,0.027778,
							  0.5,0.131291,
							  0.8,0.221551,
							  1.2,0.253691,
							  1.6,0.072905,
							  2.0,0.152243,
							  2.5,0.051564,
							  3.0,0.053023,
							  0,0,
							  0,0,
							  0,0};


	static double sog202[26]={1.6768,1.00,
							  0.2, 0.034724,
							  0.6, 0.430761,
							  0.9, 0.203166,
							  1.4, 0.192986,
							  1.9, 0.083866,
							  
							  2.3, 0.033007,
							  2.6, 0.014201,
							  3.1, 0.000000,
							  3.5, 0.006860,
							  4.2, 0.000000,
							  
							  4.9, 0.000438,
							  5.2, 0.000000
							  };
	  double A=p+n; // the size of the nucleus will not change (only the density)
      if(p>46) 
      {    // for big nuclei use rescaled  Palladium   
		  return ::density(r*pow((46+62)/A,1.0/3.0),46,62);
	  }


	 r/=fermi;
     double s=1/fermi3;
	 double r2=r*r; 
	 double r5=r2*r2*r;
	 
	 double m=1;
	 int kod=p*1000+n;
	 

	while(true) 
	{int a=kod/1000;
	//cout<<"DENSITY "<<kod<<endl;
	switch(kod)
	{
	case 0: return 0;
	case 1: return 0;
	case 1000: return 0;
	case 3004: return m*s*prfH0(1.772, 0.327, 0.151583,r);
	case 4005: return m*s*prfH0(1.776, 0.631, 0.148229,r);
	case 5005: return m*s*prfH0(1.718, 0.837, 0.157023,r);
	case 5006: return m*s*prfH0(1.698, 0.811, 0.182048,r);
	case 8008: return m*s*prfH0(1.83314, 1.544, 0.140667,r);
	case 8009: return m*s*prfH0(1.79818, 1.498, 0.161712,r);
	case 8010: return m*s*prfH0(1.84114, 1.513, 0.158419,r);

	//Modified harmonic-oscillator model (MHO):liczba protonow, liczba neutronow, a, alfa, rho0
	case 6007: return m*s*prfMH0(1.6359, 1.403, 0.118308, r);
	case 6008: return m*s*prfMH0(1.734, 1.3812, 0.108294, r);

	//Two - parameter Fermi model (2pF):liczba protonow, liczba neutronow, c=a, z=alfa, rho0
	case  9010: return m*s*prf2pf(2.58,   0.567,  0.178781,r);
	case 10010: return m*s*prf2pf(2.740,  0.572,  0.162248,r);
	case 10012: return m*s*prf2pf(2.782,  0.549,  0.176167,r);
	case 12014: return m*s*prf2pf(3.05,   0.523,  0.16955, r);
	case 13014: return m*s*prf2pf(2.84,   0.569,  0.201502,r);
	case 18018: return m*s*prf2pf(3.54,   0.507,  0.161114,r);
	case 18022: return m*s*prf2pf(3.53,   0.542,  0.176112,r);
	case 22026: return m*s*prf2pf(3.84315,0.5885, 0.163935,r);
	case 23028: return m*s*prf2pf(3.94,   0.505,  0.17129,r);
	case 24026: return m*s*prf2pf(3.941,  0.5665, 0.161978,r);
	case 24028: return m*s*prf2pf(4.01015,0.5785, 0.159698,r);
	case 24029: return m*s*prf2pf(4.000,  0.557,  0.165941,r);
	case 25030: return m*s*prf2pf(3.8912, 0.567,  0.184036,r);//dopisane
	case 30025: return m*s*prf2pf(3.8912, 0.567,  0.184243,r);
	case 26028: return m*s*prf2pf(4.074,  0.536,  0.162833,r);
	case 26030: return m*s*prf2pf(4.111,  0.558,  0.162816,r);
	case 26032: return m*s*prf2pf(4.027,  0.576,  0.176406,r);
	case 27032: return m*s*prf2pf(4.158,  0.575,  0.164824,r);
	case 29034: return m*s*prf2pf(4.218,  0.596,  0.167424,r);
	case 29036: return m*s*prf2pf(4.252,  0.589,  0.169715,r);
	case 30034: return m*s*prf2pf(4.285,  0.584,  0.164109,r);
	case 30036: return m*s*prf2pf(4.340,  0.559,  0.165627,r);
	case 30038: return m*s*prf2pf(4.393,  0.544,  0.166314,r);
	case 30040: return m*s*prf2pf(4.426,  0.551,  0.167171,r);
	case 32038: return m*s*prf2pf(4.442,  0.5857, 0.162741,r);
	case 32040: return m*s*prf2pf(4.452,  0.5737, 0.167365,r);
	case 38050: return m*s*prf2pf(4.83,   0.49611,0.168863,r);
	case 39050: return m*s*prf2pf(4.76,   0.57129,0.172486,r);
	case 41052: return m*s*prf2pf(4.875,  0.57329,0.16862,r);

	 //Three - parameter Fermi model (3pF): liczba protonow, liczba neutronow, c=a, z=alfa, w, rho0
	case 11003: return m*s*prf3pf(2.57000, 0.5052, -0.18070, 0.179227*14/11.2289,r);
	case 11004: return m*s*prf3pf(2.33430, 0.4985,  0.13930, 0.165689*15/20.8773,r);
	case 12012: return m*s*prf3pf(3.19230, 0.6046, -0.24920, 0.178889*24/19.4862,r);  //           norma(12,12)=19.4862
	case 12013: return m*s*prf3pf(3.22500, 0.5840, -0.23634, 0.178949*25/20.5491,r);  //           norma(12,13)=20.5491
	case 14014: return m*s*prf3pf(3.34090, 0.5803, -0.23390, 0.181216*28/23.298,r);   //           norma(14,14)=23.298
	case 14015: return m*s*prf3pf(3.33812, 0.5472, -0.20312, 0.183256*29/24.5853,r);  //           norma(14,15)=24.5853
	case 14016: return m*s*prf3pf(3.25221, 0.5532, -0.07822, 0.175686*30/29.0314,r);  //           norma(14,16)=29.0314
	case 15016: return m*s*prf3pf(3.36925, 0.5826, -0.17324, 0.181196*31/27.7566,r);  //           norma(15,16)=27.7566
	case 17018: return m*s*prf3pf(3.47632, 0.5995, -0.102,   0.171426*35/33.5399,r);  //           norma(17,18)=33.5399
	case 17020: return m*s*prf3pf(3.55427, 0.5885, -0.132,   0.177794*37/34.4144,r);  //           norma(17,20)=34.4144
	case 19020: return m*s*prf3pf(3.74325, 0.5856, -0.2012,  0.176237*39/34.2749,r);  //           norma(19,20)=34.2749
	case 20020: return m*s*prf3pf(3.76623, 0.5865, -0.16123,  0.16983*40/36.259,r);   //           norma(20,20)=36.259
	case 20028: return m*s*prf3pf(3.7369, 0.5245, -0.030,    0.188785*48/45.6756,r);  //           norma(20,28)=45.6756
	case 28030: return m*s*prf3pf(4.3092, 0.5169, -0.1308,   0.169202*58/50.8799,r);  //           norma(28,30)=50.8799
	case 28032: return m*s*prf3pf(4.4891, 0.5369, -0.2668,   0.176259*60/49.7148,r);  //           norma(28,32)=49.7148
	case 28033: return m*s*prf3pf(4.4024, 0.5401, -0.1983,   0.176952*61/52.5212,r);  //           norma(28,33)=52.5212
	case 28034: return m*s*prf3pf(4.4425, 0.5386, -0.2090,   0.177167*62/53.0451,r);  //           norma(28,34)=53.0451
	case 28036: return m*s*prf3pf(4.5211, 0.5278, -0.2284,   0.177731*64/53.855,r);   //           norma(28,36)=53.855

	 //Three - parameter Gaussian model (3pG):liczba protonow, liczba neutronow, c=a, z=alfa, w, rho0
	case 16016: return m*s*prf3pg(2.549,2.19110, 0.16112,0.224317*32/32,r);       //     norma(16,16)=32
	case 40050: return m*s*prf3pg(4.434,2.5283,  0.35025,0.189844*90/102.817,r);  //     norma(40,50)=102.817
	case 40051: return m*s*prf3pg(4.32520,2.5813,0.43325,0.200597*91/110.438,r);  //     norma(40,51)=110.438
	case 40052: return m*s*prf3pg(4.45520,2.5503, 0.33425,0.19084*92/104.035,r);  //     norma(40,52)=104.035
	case 40054: return m*s*prf3pg(4.49420,2.5853,0.29625,0.189365*94/103.64,r);   //    norma(40,54)=103.64
	case 40056: return m*s*prf3pg(4.50320,2.6023,0.34125,0.191673*96/109.176,r);  //     norma(40,56)=109.176
	case 42050: return m*s*prf3pg(4.6110,2.527,0.1911,   0.176603*92/94.0081,r);  //     norma(42,50)=94.0081

	  //Fourier-Bessel Coefficients
	 case  2001: return m*s*prfFB(fb201,   r)*A/a;
	 case  6006: return m*s*prfFB(fb606,   r)*A/a; 
	 case  7008: return m*s*prfFB(fb708,   r)*A/a;
	 case 16018: return m*s*prfFB(fb1618,  r)*A/a;
	 case 16020: return m*s*prfFB(fb1620,  r)*A/a;
	 case 22028: return m*s*prfFB(fb2228,  r)*A/a;
	 case 24030: return m*s*prfFB(fb2430,  r)*A/a;
	 case 32042: return m*s*prfFB(fb3242,  r)*A/a;
	 case 32044: return m*s*prfFB(fb3244,  r)*A/a;
	 case 42052: return m*s*prfFB(fb4252,  r)*A/a;
	 case 42054: return m*s*prfFB(fb4254,  r)*A/a;
	 case 42056: return m*s*prfFB(fb4256,  r)*A/a;
	 case 42058: return m*s*prfFB(fb4258,  r)*A/a;
	 case 46058: return m*s*prfFB(fb4658,  r)*A/a;
	 case 46060: return m*s*prfFB(fb4660,  r)*A/a;
	 case 46062: return m*s*prfFB(fb4662,  r)*A/a;

	//Sum of gaussians model
	 case 1002: return m*s*prfSOG(sog102, r, A);
	 case 2002: return m*s*prfSOG(sog202, r, A);

	//charge distribution analysis
	case 1001: return m*s*prfMI(2.1166,0.21,0.77,r); 
	//case 1001: return s*prfMI(2.116,0.21,0.77,r); 
	//case 2002: return s*prfMI(1.67114,0.45,1.92,r); 

	default:  kod=best_pn(p, n); // szukamy najbliższego zaimplementowanego modelu
			  m=(kod/1000+kod%1000)/A; 
		  break; 
	}}
	return 0; //tu nigdy nie dojdzie 
}

///////////////////////////////////////////////////////////////////////////////
/// calculates nuclear density of proton. neutrons at radius r from the center
///////////////////////////////////////////////////////////////////////////////
double anynucleus::density(double r)
{
	//double dens=::density(r,pr,nr); //size of nucleus changes with density
	//return dens;
	
	double dens=::density(r,p,n);
	return dens*(pr+nr)/(p+n); // is there are missing nucleons assume 
	                           // same size but lower density
}




///////////////////////////////////////////////////////////////////////////////
/// for given p, n find the nearest known  p, n   
///////////////////////////////////////////////////////////////////////////////
static  int best_pn(int p, int n)
    { 
      const int N=85;
      /// known density profiles
      int znane[N]={0,1,1000,1001,1002,2001,2002,3004,4005,5005,5006,6006,6007,6008,7008,8008,8009,8010,
		 9010,10010,10012,11003,11004,12012,12013,12014,13014,14014,14015,14016,15016,
		16016,16018,16020,17018,17020,18018,18022,19020,20020,20028,22026,22028,23028,
		24026,24028,24029,24030,25030,26028,26030,26032,27032,28030,28032,28033,28034,
		28036,29034,29036,30025,30034,30036,30038,30040,32038,32040,32042,32044,38050,
		39050,40050,40051,40052,40054,40056,41052,42050,42052,42054,42056,42058,46058,
		46060,46062};
      int besti=0;
      int best=pow2(znane[0]/1000-p)+pow2(znane[0]%1000-n);
      for(int i=1; i<N; i++)
        {int x=pow2(znane[i]/1000-p)+pow2(znane[i]%1000-n);
		 if(x==0) return znane[i];
		 if(x<best) {best=x;besti=i;
		 }
	}
    return znane[besti]; 	
   }
  
///////////////////////////////////////////////////////////////////////////////
/// random distance from the center weighted with nuclear density
///////////////////////////////////////////////////////////////////////////////
  double anynucleus::get_random_r()
  { double max_density=0.3/fermi3;//correction needed!!!
    double r0=radius();
    double r1;
    do
      {
      r1 = r0*frandom_sqr(); // quadratic distribution
      }
    while ( density(r1) < frandom ()*max_density*(pr+nr+1)/A());
//    cout<< "r1="<<r1<<endl;
    return r1;
  } 
  
