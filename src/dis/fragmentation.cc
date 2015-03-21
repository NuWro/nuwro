#include<TMCParticle.h>
#include<TPythia6.h>
//#include "sobek.h"
#include "generatormt.h"
#include "fragmentation_cc.h"
#include "fragmentation_nc.h"
#include "pdg_name.h"

//routines called from pythia6


void
hadronization (double E, double W, double nu, double m, int lepton_in,
	       int nukleon_in, bool current)
{

  if (current == true && lepton_in > 0)
    {
      switch (nukleon_in)
	{
	case proton:
	  hadronization_cc_nu_p (E, W, nu, m);
	  break;
	case neutron:
	  hadronization_cc_nu_n (E, W, nu, m);
	  break;
	}
    }

  if (current == true && lepton_in < 0)
    {
      switch (nukleon_in)
	{
	case proton:
	  hadronization_cc_anu_p (E, W, nu, m);
	  break;
	case neutron:
	  hadronization_cc_anu_n (E, W, nu, m);
	  break;
	}
    }

  if (current == false && lepton_in > 0)
    {
      m = 0;
      switch (nukleon_in)
	{
	case proton:
	  hadronization_nc_nu_p (E, W, nu, m);
	  break;
	case neutron:
	  hadronization_nc_nu_n (E, W, nu, m);
	  break;
	}
    }

  if (current == false && lepton_in < 0)
    {
      m = 0;
      switch (nukleon_in)
	{
	case proton:
	  hadronization_nc_anu_p (E, W, nu, m);
	  break;
	case neutron:
	  hadronization_nc_anu_n (E, W, nu, m);
	  break;
	}
    }


}
