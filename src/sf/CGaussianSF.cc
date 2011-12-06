#include "CGaussianSF.h"
#include <iostream>
using namespace std;
//Gaussian approximation of the spectral function

CGaussianSF::CGaussianSF(const MomDistribs i_momDistrib, 
						 const IsospinOfSF i_isospin,
						 TargetData atom,
//						 int ZZ, int NN, double E2, double Pauli,
						 double pmin, double pmax, double emin, double emax,
  	                 	 const double *levP,const double * levN)
						 :
			m_isospin(i_isospin),atomZ(atom.Z),atomN(atom.N),
			m_targetBeta( atom.Beta ),
			m_targetE2(atom.E2)
{
			m_pBlock=(atom.PBlock);
			m_pStart=(pmin);
			m_pStop=(pmax);
			m_eStart=(emin);
			m_eStop=(emax);
	m_momDistrib = createDistrib(i_momDistrib, i_isospin);
    
	set_targetAlpha(m_targetBeta);

	m_coeff = M*sqrt(m_targetAlpha/pi);
    
    if(m_isospin==proton)
	  levels=levP;
	else  
	  levels=levN;
	
}


double CGaussianSF::evalSF(const double p, const double removalE) const
{

	double mf=0;
	for(const double *l=levels;l[0];l+=3)
	   mf+=l[0]*Gauss(l[1],l[2],removalE);

	if (m_isospin == proton)
	  mf*=m_momDistrib->MF(p)/atomZ;
	else
	  mf*=m_momDistrib->MF(p)/atomN;

	if ( p == 0.0 )
		return mf;

	double squareRoot = 2.0*M*m_targetBeta*(removalE - m_targetE2);

	if(squareRoot<0.0)
		return mf;

	squareRoot = sqrt(squareRoot);

	const double pMinimum(m_targetBeta*p - squareRoot);
	const double pMaximum(m_targetBeta*p + squareRoot);

	const double corr = m_coeff/p*m_momDistrib->Corr(p)*( exp(-m_targetAlpha*pMinimum*pMinimum) 
	                                                    - exp(-m_targetAlpha*pMaximum*pMaximum) 
														);
    return mf + corr;

}

double CGaussianSF::norm() const
{
    double mf=0;
	for(const double *l=levels;l[0];l+=3)
	   {mf+=l[0]*(GaussInt(l[1],l[2],m_eStop)-GaussInt(l[1],l[2],m_eStart));
//	   std::cout<<m_coeff<<std::endl;
       }
    return mf;
}



double CGaussianSF::generateE(const double p) const
{
	double a=m_eStart;
	double b=m_eStop;
	double x;
	double mf;
	
	double y1=0,y2=0;
	
	for(const double *l=levels;l[0];l+=3)
	   {y1+=l[0]*GaussInt(l[1],l[2],a);
	    y2+=l[0]*GaussInt(l[1],l[2],b);
	   } 
//	cout<<"y1="<<y1<<endl;   
//	cout<<"y2="<<y2<<endl;   
//	cout<<"y2-y1="<<y2-y1<<endl;   
//	cout<<"norm="<<this->norm()<<endl;   
	
  	   
	double des=y1+(y2-y1)*frandom();  
//	cout<<"des="<<des<<endl; 
	
	while(true)
	{ 
		x=(a+b)/2;
		mf=0;
		for(const double *l=levels;l[0];l+=3)
		   mf+=l[0]*GaussInt(l[1],l[2],x);
//		cout<<"mf="<<mf<<endl;   
		if(mf<des)
		  {if(a==x) break; 
		   a=x;
		   }
		else
		   {if(b==x) break;
		   b=x;
	       }
     };
	 return x;
}
