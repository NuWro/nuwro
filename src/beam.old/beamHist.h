

#ifndef _BEAMHIST_H_
#define _BEAMHIST_H_

#include "rootf.h"
#include "beam.h"
#include "generatormt.h"
#include "particle.h"
#include "narray.h"
#include "makeHist.h"
#include "nd280stats.h"




class BeamHist : public beam
{
	/// number of dimesnions (xnu, ynu, nnu[3])
	enum
	{ 
		EDims = 5
	};
	

	NArray< int >     *_hist;  /// histogram
	Nd280Statistics   _stats;  /// detail description for histogram, like min & max
	vector<double>    _inf;    /// array of infinitizimal values,
	                           /// width of each element
	vector<int>       _f;      /// this factor says about n-dim array resolution,
	                           /// each dimension has own number of elements
	
	//vector< NArrayItem > 	*	_hist;
	//RootFSaver< Nd280Element > *  _rootout;
	
public:
	BeamHist( int dim_res[EDims], string stats_fname, string histout_fname );
	BeamHist( const BeamHist & copy ) { throw "err: BeamHist<T> copying is not allowed"; }
	virtual ~BeamHist();
	
	/// create new particle for simulation.
	/// each one is defined using the histogram
	virtual particle shoot();
};


#endif // _BEAMHIST_H_


