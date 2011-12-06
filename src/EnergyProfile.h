#ifndef _EnergyProfile_h_
#define _EnergyProfile_h_
#include <string>
#include <sstream>
#include <iostream>
#include <cstdlib>
#include "generatormt.h"
#include "jednostki.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
/// EnergyProfile shoots a double precision number whose probability density is
/// given by the parameters to the constructor
////////////////////////////////////////////////////////////////////////////////

class EnergyProfile
{ private:
     int n;
     double Emin;
     double Emax;
     double *prob;
     double *prob2; // for dis
  public:
     ~EnergyProfile(){delete [] prob;delete [] prob2;}

     /// Beam of Zero Energy
     EnergyProfile():n(0),Emin(0),Emax(0),prob(NULL),prob2(NULL){}

     /// Monoenergetic beam
     EnergyProfile(double E):n(0),Emin(E),Emax(E),prob(NULL),prob2(NULL){}

     /// uniform probability density betweeen E1 and E2 
     EnergyProfile(double E1,double E2):n(1),Emin(E1),Emax(E2),prob(NULL),prob2(NULL){}

     /// energy profile given by a histogram encoded in a string
     /// the string contains space separated values:
     /// Emin Emax b1 b2 ... bn 
     /// where b1 b2 ... bn are the heightsof the bars of the histogram
     /// the widths are assumed equal to (Emax-Emin)/n 
     EnergyProfile(string s):n(0),Emin(0),Emax(0),prob(NULL),prob2(NULL) {read(s);}

     /// the auxiliary helper function for parsing the string paramater to 
     /// the EnergyProfile constructor  
     void read(string s);

     /// a random number conforming to the probability profile
     /// specified by construcor parameters
     double shoot(bool dis);

     /// diagnostic function for verifying if the beam 
     /// has been correctly constructed
     double print();
     double minE();
     /// how many times weighting by E enhances this beam
     double disratio(){if(n==0) return 1;else return prob2[n-1]/prob[n-1];}
     
};  // Energy Profile class

////////////////////////////////////////////////////////////////////////////////
///              I M P L E M E N T A I O N
////////////////////////////////////////////////////////////////////////////////

/// parse the string parameter to the EnergyProfile constructor  
inline void EnergyProfile::read(string s)
     {
      if(prob) 
        delete[] prob;
      if(prob2) 
        delete[] prob2;
      prob2=prob=NULL;
      Emin=Emax=n=0; 
      
      stringstream in(s);
      in>>Emin>>Emax;
      Emin*=MeV;
      Emax*=MeV;
      if(!in) 
        return;
      n=1;	
      double * bufor=new double[5000];	
      int i=0;
      while(in>>bufor[i])
         i++;
      if(i>1)
		{
			n=i;
			prob=new double[n];
			prob2=new double[n];
			double prev1=0,prev2=0;
			for(int i=0;i<n;i++) 
			  {prob[i]=prev1+=bufor[i];   
			   prob2[i]=prev2+=bufor[i]*(Emin+(0.5+i)*(Emax-Emin)/n);
		      }
		}   
		delete []bufor;
     }

/// yield a random number with the probability profile
/// specified by the constructor parameters
inline double EnergyProfile::shoot(bool dis) 		 
     { if(n==0) return Emin;
       if(n==1) return Emin+frandom()*(Emax-Emin);	  
       double *pro=dis?prob2:prob;
	      double x=frandom()*pro[n-1];
	      int i=0,j=n-1;
	      while(i<j)
	      {
			  int s=(i+j)/2;
			  if(x<pro[s])
			    j=s;
			  else
			    i=s+1;  
		  }
		  double z=frandom();
		  if(dis)
		    {double a=Emin+i*(Emax-Emin)/n;
		     double b=a+(Emax-Emin)/n;
		     return sqrt(a*a+z*(b*b-a*a));
			}
		  return Emin+(i+z)*(Emax-Emin)/n;      
     }

/// diagnostic function for verifying if the beam 
/// has been correctly constructed
inline double EnergyProfile::print() 		 
     { 
       cout<<"Emin="<<Emin<<endl;
       cout<<"Emax="<<Emax<<endl;
       cout<<"n="<<n<<endl;
       if(n>1)
       for(int i=0;i<n;i++)
         cout<<"prob["<<i<<"]="<<prob[i]-(i?prob[i-1]:0)<<endl;
     }
     
inline double EnergyProfile::minE()
{
	return Emin;
}
////////////////////////////////////////////////////////////////////////////////
///        END OF  I M P L E M E N T A I O N
////////////////////////////////////////////////////////////////////////////////

#endif
