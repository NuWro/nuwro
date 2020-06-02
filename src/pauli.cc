#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include "params.h"
#include "jednostki.h"
#include "dis/pdg_name.h"
#include "dis/parameters.h"
#include "jednostki.h"
#include "vect.h"
#include "nucleus.h"
#include "event1.h"

////////////////////////////////////////////////////////////////////////////////////
//      Pauli blocking
/////////////////////////////////////////////////////////////////////////////////
using namespace NSNWRO;

void
mypauli_spp (event & e, nucleus & t)
{
  if (e.out.size () < 4)
    {
      for (int i = 1; i < e.out.size (); i++)
	{
	  if ( t.pauli_blocking_old (e.out[i], e.in[1].length() ) ) 
	  {
	    e.weight = 0;
	    return;
	  }
	}
    }
}

//      End Pauli blocking
/////////////////////////////////////////////////////////////////////////////////////
