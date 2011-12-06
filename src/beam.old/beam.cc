
#include "beam_uniform.h"
#include "beamRF.h"
#include "beamHist.h"
#include "beam_mixed.h"


void CreateNewHistogram( int dimsizes[5], string stats_out, string hist_out );


beam * create_beam(params &p)
{   
	int pdg = p.beam_particle;
	if (!(pdg == 12 or pdg == 14 or pdg == 16 or pdg == -12 or pdg == -14 or pdg == -16))
		throw "the PDG code of the beam particle is not neutrino code";

	switch( p.beam_type )
	{
	case 0:   return new beam_uniform(p);
	case 1:   return new beam_mixed(p);
	case 2:   return new BeamRF("../flux");
	case 3:
	{
		int dimsizes[5] = { 17, 17, 17, 17, 17 };
		string stats("statsout.txt");
		string histout("histout.txt");

		return new BeamHist(dimsizes, stats, histout);
	}
	case 4:
	{
		int dimsizes[5] = { 17, 17, 17, 17, 17 };
		string statsout("statsout.txt");
		string histout("histout.txt");

		/// create histogram and save it into a file
		CreateNewHistogram( dimsizes, statsout, histout );
		
		/// return NULL, rest of the simulation will stop here
		return NULL;
	}
	default:
	{
		throw "no beam defined";
		return NULL;
	}
	} /// switch end
}


void CreateNewHistogram( int dimsizes[5], string stats_out, string hist_out )
{
	/// create 5-dimensional matrix
	NArray< int >  histogram( dimsizes, 5 );
	for( int i = 0; i < histogram.FlatCount(); ++i )
		histogram.Flat( i ) = 0;

	/// create new histogram using root files with neutrins
	HistMaker< Nd280Element, Nd280Statistics > histmaker( string("../flux/"), string("h3002") );
	histmaker.Create( histogram );

	/// save histogram and detail description (statistics) into files 
	histogram.Save( hist_out );
	ofstream statsfile( stats_out.c_str() );
	statsfile << histmaker.Stats()->SummaryStr();
}
