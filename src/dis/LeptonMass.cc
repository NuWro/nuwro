#include "pdg_name.h"
#include <cmath>
#include <cassert>
#include "vec.h"
#include "jednostki.h"
#include "vect.h"
#include "generatormt.h"

double
lepton_mass (int lepton_in, bool cur)
{ 
  int lepton = lepton_in;
  if (lepton< 0)
    lepton= -lepton;
  double m = 0;

  if (cur == true)
    {
      switch (lepton)
	{
	case nu_e:
	  m = mass_e;
	  break;
	case nu_mu:
	  m = mass_mu;
	  break;
	case nu_tau:
	  m = mass_tau;
	  break;
	}
    }
  return m;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////This function produces random unit vector d with a known cosine with respect to a vector  b
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
kinfinder (vec b, vec & d, double kos)
{
  if (fabs (kos) > 1)
    cout << "Error !!! cos must be between -1 and 1 !!!" << endl;

  vec dr, trz, czw;
  if (b.x != 0 || b.y != 0)
    dr = vec (0, 0, 1);
  else
    dr = vec (1, 0, 0);		//dr is an arbitrary vector not parallel to b

  trz = vecprod (dr, b);	//veeeery stupid way of defining vector product but the rules of the game were unclear for JTS
  czw = vecprod (trz, b);	// vectors b, trz and czw form orthogonal system

  double phi = 2 * Pi * frandom ();
  double sphi = sin (phi);
  double cphi = cos (phi);
  double ssin = sqrt (1 - kos * kos);
  b = b / b.length ();
  trz = trz / trz.length ();
  czw = czw / czw.length ();	//now they are orthonormal

  d = kos * b + ssin * sphi * trz + ssin * cphi * czw;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////Decay of hadronic mass W into 2 particles; no correlations
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
kin2part (double hama, int nukleon2, int meson, vect & finnuk, vect & finpion)
{
  double nukmass, pionmass;
  double W2 = hama * hama;
  double W4 = W2 * W2;

  if (meson == 211 || meson == -211)	//this part must be done better !!!
    pionmass = 139.57;
  if (meson == 111)
    pionmass = 134.98;
  if (nukleon2 == 2212)
    nukmass = 938.27;
  if (nukleon2 == 2112)
    nukmass = 939.565;

  double nukmass2 = nukmass * nukmass;
  double pionmass2 = pionmass * pionmass;

  if (hama < 1080)
    cout << "error!!! Delta decay is kinematically impossible !!!" << endl;

  double momentum =
    sqrt (W4 - 2 * W2 * (pionmass2 + nukmass2) +
	  (nukmass2 - pionmass2) * (nukmass2 - pionmass2)) / 2.0 / hama;
  double nukenergy = (W2 + nukmass2 - pionmass2) / 2.0 / hama;
  double pionenergy = (W2 - nukmass2 + pionmass2) / 2.0 / hama;

  vec kierunek = rand_dir ();

  finnuk =
    vect (nukenergy, momentum * kierunek.x, momentum * kierunek.y,
	  momentum * kierunek.z);
  finpion =
    vect (pionenergy, -momentum * kierunek.x, -momentum * kierunek.y,
	  -momentum * kierunek.z);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////// Rotation of the particle produced by PYTHIA acoording to the direction of the momentum transfer
//////////////////////// ( in PYTHIA it is assumed that this direction is the Z axis)
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
rotation (vect & cztero, vec trzy)
{
  vec wybrany = vec (1, 0, 0);	//

  vec nowy1 = vecprod (trzy, wybrany);	//the vector orthogonal to trzy and odcztero
  vec nowy2 = vecprod (nowy1, trzy);	//vectors trzy, nowy1 and nowy2 are orthogonal

  trzy = trzy / trzy.length ();
  nowy1 = nowy1 / nowy1.length ();
  nowy2 = nowy2 / nowy2.length ();	//they are now orthonormal

  vec piaty = cztero.z * trzy + cztero.x * nowy1 + cztero.y * nowy2;

  cztero.x = piaty.x;
  cztero.y = piaty.y;
  cztero.z = piaty.z;

}

////////////////////////////////////////////////////////////////////////
/////////// Binding energy in the effective approach
///////////////////////////////////////////////////////////////////////
double
binen (vec mom, int p, int n)
{
  double x = mom.length ();

if (p==8 && n==8)
{
  if (x >= 0 && x < 330)
    return 4.748678e+01 - 1.722751e-01*x + 1.399298e-04*x*x + 5.699710e-06*x*x*x - 2.908039e-08*x*x*x*x + 1.266971e-10*x* x*x*x*x - 1.634620e-13*x*x*x*x*x*x;


  if (x >= 330 && x <= 800)
    return -1.144537e+02 + 9.459352e-01*x -4.987450e-04*x*x +3.857279e-07*x*x*x;

}

if (p==6 && n==6)
{
if (x >= 0 && x < 330)
    return 4.192962e+01  -2.513311e-01*x + 5.997493e-03*x*x - 8.607908e-05*x*x*x +6.315994e-07*x*x*x*x - 2.080878e-09*x*x*x*x*x + 2.539746e-12*x*x*x*x*x*x;


 if (x >= 330 && x <= 800)
    return 1.201904e+01 + 5.318283e-02*x + 1.279047e-03*x*x -5.507682e-07*x*x*x;

}
return 0;
}

double
binen2 (double x, int p, int n)
{

if (p==8 && n==8)//oxygen
{
  if (x >= 0 && x < 330)
    return 4.748678e+01 - 1.722751e-01*x + 1.399298e-04*x*x + 5.699710e-06*x*x*x - 2.908039e-08*x*x*x*x + 1.266971e-10*x* x*x*x*x - 1.634620e-13*x*x*x*x*x*x;


  if (x >= 330 && x <= 800)
    return -1.144537e+02 + 9.459352e-01*x -4.987450e-04*x*x +3.857279e-07*x*x*x;

}

if (p==6 && n==6)//carbon
{
if (x >= 0 && x < 330)
    return 4.192962e+01  -2.513311e-01*x + 5.997493e-03*x*x - 8.607908e-05*x*x*x +6.315994e-07*x*x*x*x - 2.080878e-09*x*x*x*x*x + 2.539746e-12*x*x*x*x*x*x;


 if (x >= 330 && x <= 800)
    return 1.201904e+01 + 5.318283e-02*x + 1.279047e-03*x*x -5.507682e-07*x*x*x;

}
return 0;
}

double
deuter_binen (vec mom)
{
double pp = mom.length ();
return 2.224 + 2* (sqrt(pp*pp + 938.5*938.5) - 938.5);
}

double
deuter_binen2 (double pp)
{
return 2.224 + 2* ( sqrt(pp*pp + 938.5*938.5) - 938.5 );
}
