

#ifndef _ND280STATISTICS_H_
#define _ND280STATISTICS_H_

#include "TFile.h"
#include "TTree.h"
#include "narray.h"



class Nd280Element
{
public:
	Float_t Enu;       /// energia
	Float_t nnu[3];    /// x, y, z kierunek pedu
	Float_t xnu, ynu;  /// wsp. gdzie neutrino przecina plaszczyzne Z=0
	
	Nd280Element( ) : Enu(0), xnu(0), ynu(0)
	{
		nnu[0] = 0;
		nnu[1] = 0;
		nnu[2] = 0;
	}

	Nd280Element( const Nd280Element & e ) : Enu(e.Enu), xnu(e.xnu), ynu(e.ynu)
	{
		nnu[0] = e.nnu[0];
		nnu[1] = e.nnu[1];
		nnu[2] = e.nnu[2];
	}
	
	void CreateBranches( TTree * tree )
	{
		tree->Branch( "Enu", &Enu, "Enu/F" );
		tree->Branch( "nnu", &nnu, "nnu[3]/F" );
		tree->Branch( "xnu", &xnu, "xnu/F" );
		tree->Branch( "ynu", &ynu, "ynu/F" );
	}
	
	void OpenBranches( TTree * tree )
	{
		tree->SetBranchAddress( "Enu", &Enu );
		tree->SetBranchAddress( "nnu", &nnu );
		tree->SetBranchAddress( "xnu", &xnu );
		tree->SetBranchAddress( "ynu", &ynu );
	}
};




class Nd280Statistics
{
	typedef const double cdouble;
	typedef const unsigned int cuint;

	Nd280Element      _min;
	Nd280Element      _max;
	
	/// Number of elements under porocess.
	/// It's -1 if object is initialised from a file
	int               _count;

public:
	Nd280Statistics() : _count(0) {}
	
	Nd280Statistics( string fname );
	
	Nd280Statistics(Nd280Element minimum, Nd280Element maximum)
	:_min(minimum),
	_max(maximum),
	_count(0)
	{}
	
	/// get min and maximum value in histogram
	Nd280Element GetMin() const { return _min; }
	Nd280Element GetMax() const { return _max; }
	
	/// get number of elements in histogram.
	/// this is equal -1 if object was init from a file
	int Count() const { return _count; }

	/// an element might be a min or maximum,
	/// use the function to check it
	void FindExtremes( const Nd280Element & entry )
	{
		if( _count == 0 )
		{
			Init(entry);
		}
		else
		{
			CheckMin( entry );
			CheckMax( entry );
		}
		++_count;
	}

	/// increment one pick in the histogram, another words
	/// icrease number of element in one histogram cell
	void Increment( const  Nd280Element & x,  NArray< int > & hist )
	{
            double E = x.Enu * 1000;

		/// Increase one cell in n-dimensional array (histogram).
		++hist( 
			Idx(x.xnu, _min.xnu, _max.xnu, hist.Count(0)),
			Idx(x.ynu, _min.ynu, _max.ynu, hist.Count(1)),
			Idx(x.nnu[0] * E, _min.nnu[0], _max.nnu[0], hist.Count(2)),
			Idx(x.nnu[1] * E, _min.nnu[1], _max.nnu[1], hist.Count(3)),
			Idx(x.nnu[2] * E, _min.nnu[2], _max.nnu[2], hist.Count(4)) );
	}

	/// return the string ready for saving in a file
	string SummaryStr() const;
	
	/// fill the object using the file (input file stream)
	bool Load( ifstream & file );
	
	
private:
	/// check if x element is smaller then current minimum,
	/// if yes then overwrite old value
	void CheckMin( const Nd280Element & x );

	/// check if x element is larger then current maximum,
	/// if yes then overwrite old value	
	void CheckMax( const Nd280Element & x );
	
	/// use it at the begining, run with first element
	void Init( const Nd280Element & entry );
	
	/// count one of indexes in N-dim array.
	/// bunch of these indexes discribe one cell which value is incremented
	inline unsigned int Idx( cdouble & val, cdouble & min, cdouble & max, cuint & count );
};





unsigned int Nd280Statistics::Idx( cdouble & val, cdouble & min, cdouble & max, cuint & count )
{
	double fraction = (val-min) * count;
	double scope = max - min;
	unsigned int result = static_cast<int>(fraction/scope);
	if( result >= count )
		result = count - 1;

	return result;
}
	
	
#endif // _ND280STATISTICS_H_
