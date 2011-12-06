#ifndef _beam_h_
#define _beam_h_

#include "params.h"
#include "particle.h"

class beam
{
public:
	virtual particle shoot() = 0;
};

/// The beam shoots with identical particles (whose energy is defined by some 
/// energy profile but the direction of momentum is identical)
/// beam: this class provides the input of the initial beam particle 
/// into the generator.The beam particle is defined by providing 
/// its energy, direction and its PDG code to the constructor
/// of the beam class 

#include "beam_uniform.h"
#include "beamRF.h"
#include "beamHist.h"
#include "beam_mixed.h"

inline beam * create_beam(params &p)
{   
	int pdg = p.beam_particle;
	if (!(pdg == 12 or pdg == 14 or pdg == 16 or pdg == -12 or pdg == -14 or pdg == -16)) throw "the PDG code of the beam particle is not neutrino code";
	cout<<"BEAM"<<endl;
	switch(p.beam_type)
     { 
	   case 0:	  return new beam_uniform(p);
	   case 1:	  return new beam_mixed(p);
	   case 2:	  return new BeamRF(p);
	   case 3:	  {
					int dimsizes[6] = { 5, 5, 5, 5, 5, 5 };
					auto_ptr< NArray< int > >  histogram( new NArray<int>(dimsizes, 6) );
					for( int i = 0; i < histogram->FlatCount(); ++i )
						histogram->Flat( i ) = 0;

					ifstream file( "statsout.txt" );
					ND5Statistics * stats = new ND5Statistics;
					if( stats->Load( file ) )
					{			
						histogram->Load( "histout.txt"  ); 
						histogram->DoSum();
						vector< NArrayItem > * hist_sum = histogram->DoUniqueElemsArray();
						return new BeamHist< int >( hist_sum, stats, histogram->FlatCount(), true );
					}
					else delete stats;
				  }

		          break;
	   default:  {throw "no beam defined";return NULL;}
     }
}

#endif
