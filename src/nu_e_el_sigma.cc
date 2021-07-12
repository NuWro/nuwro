#include "nu_e_el_sigma.h"
#include <cstdlib>  // rand, abs, ...
#include <iostream>
#include <cmath>
#include "pdg.h"
//#include "ff.h"
//#include "scatter.h"
#include "jednostki.h" // G, sin2thetaW, cm2

using namespace std;

// elastic neutrino - electron scattering cross section 

double nu_e_el_sigma ( double E0nu, // neutrino energy in the target frame
                    //double q2, // squared momentum transfer (Mandelstam t)
					double t_ratio, // final-to-initial neutrino energy ratio
                    int kind,  // process type: 12 - nu-e elastic scattering
                    bool anty, // true for antineutrinos
                    double me,  //, (target) electron mass
                    int switch_sigma, // process type
                    double m_prime, // final charged lepton mass
                    double delta_m2
                    //int what  // what to return
                  )
{    
    double xsec = 0; // cross section is 0 until the reaction occurs


	double sin2W = sin2thetaW; // Weinberg angle

	double deltall;
	if (kind == 12 || (kind == -12 && switch_sigma==1))
		deltall = 1.0;
	else	deltall = 0.0;
	double c0 = 2.0*sqrt(2.0)*G; // G is the Fermi constant
	double cL = c0*(sin2W-0.5+deltall);
	double cR = c0*sin2W;

	double EFactor = E0nu*me/4.0/Pi;
	double AA,BB,CC;
    if (anty==false)
    {
	    AA = cL*cL;
	    BB = cR*cR;
    } 
    else 
    {
    	AA = cR*cR;
	    BB = cL*cL;
    }
	    CC = cL*cR;

    if (switch_sigma==1)
        xsec = EFactor * ( AA + BB*t_ratio*t_ratio + CC*(t_ratio-1.0)*me/E0nu ) / cm2;
    else if (switch_sigma==2 && kind != -12)
    {
        double IL  = 1.0 - (m_prime*m_prime - me*me)/2.0/me/E0nu;
        xsec = EFactor *  c0*c0 * IL / cm2;
    }
    else if (switch_sigma==3 || (switch_sigma==2 && kind == -12))
    {
        double IR  = t_ratio*t_ratio + t_ratio*(m_prime*m_prime - me*me)/2.0/me/E0nu;
        xsec = EFactor *  c0*c0 * IR / cm2;
    }

    xsec = xsec * (2.0*E0nu/(me+2.0*E0nu) - delta_m2/me/(me+2.0*E0nu)); // averaged <sigma>*(b-a)


    return xsec;
};
