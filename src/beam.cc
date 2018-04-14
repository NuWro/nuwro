#include "beam_uniform.h"
#include "beam_singleroothist.h"
#include "beam_mixedroothist.h"
#include "beamRF.h"
#include "beamHist.h"
#include "beam_mixed.h"
#include "nozero_array.h"
#include "geomy.h"


void CreateNewHistogram( int dimsizes[5], string hist_out );


beam * create_beam(params &p, geomy *detector)
{   
	int pdg = p.beam_particle;
	if (!(pdg == 12 or pdg == 14 or pdg == 16 or pdg == -12 or pdg == -14 or pdg == -16 or pdg==11))
		throw "the PDG code of the beam particle is not neutrino or electron code";

	switch( p.beam_type )
	{
	case 0:   return new beam_uniform(p);
	case 1:   return new beam_mixed(p);
	case 2:   return new BeamRF(p,detector);
	case 3:
	{
		int dimsizes[5] = { 17, 17, 17, 17, 17 };
		string hist_file( "histout.txt" );

		return new BeamHist( dimsizes, hist_file );
	}
	case 4:
	{
		int dimsizes[5] = { 17, 17, 17, 17, 17 };
		string histout( "histout.txt" );

		/// create histogram and save it into a file
		CreateNewHistogram( dimsizes, histout );
		
		/// return NULL, rest of the simulation will stop here
		return NULL;
	}
	case 5:
	{
	  cout << "Using single flux from root file." << endl;
	  return new beam_singleroothist( p );
	}
	case 6:
	{
	  cout << "Using combined flux from root file." << endl;
	  return new beam_mixedroothist( p );
	}
	default:
	{
		throw "no beam defined";
		return NULL;
	}
	} /// switch end
}

void CreateNewHistogram( int dimsizes[5], string hist_out )
{
	/// create 5-dimensional matrix
	NArray< int >  histogram( dimsizes, 5 );
	for( int i = 0; i < histogram.FlatCount(); ++i )
		histogram.Flat( i ) = 0;

	/// create new histogram using root files with neutrins
	HistMaker< Nd280Element, Nd280Statistics > histmaker( string("../flux/"), string("h3002") );
	histmaker.Create( histogram );
	
	/// convert to speed-optimized array
	NonZeroArray< int > opt_histogram( histogram );

	/// save histogram and detail description (statistics) into files
	ofstream histfile( hist_out.c_str() );
	histfile << histmaker.Stats()->String();
	opt_histogram.Save( histfile );
}
