

#ifndef _MAKEHIST_H_
#define _MAKEHIST_H_


#include "rootf.h"
#include "narray.h"
#include "iostream"
#include <sstream>
#include "nd280stats.h"
using namespace std;



template< typename TEvent, typename TStats > class HistMaker
{
	/// folder with root files which handle specific type of event, e.g. Nd280Statistics
	RootFolder<  RootFReader< TEvent >  >     _folder;
	/// object knows what kind of root data we are looking for and how to 
	/// fetch min, max values and also how to count index of histogram, 
	/// e.g. Nd280Statistics
	TStats                                    _stats;
	
	/// looks for min and max value of each parameters inside event
	void GetStats();
	
public:
	HistMaker( string foldername, string treename )
      :_folder( foldername, treename )
      {}
		
	/// Create histogram
	void Create( NArray< int > & hist );
	
	/// Get access to statistics
	TStats * Stats();
};


/// looks for min and max value of each parameters inside event
template< typename TEvent, typename TStats >
void HistMaker<TEvent, TStats>::GetStats()
{
	cerr << "HistMaker preparing statistics.";
	if( _folder.Count() <= 0 )
	{
		cerr << endl;
		return;
	}
	else
		cerr << endl << "Current file: ";
	
	int i = 0, j = 0;
	/// initialize parameters in first step
	RootFReader< TEvent > * file = _folder.File( i );

	++j;
	
	/// go through all root files in folder
	while(  i < _folder.Count()  )
	{
		cerr << i << " ";
		
		/// and now... for all events in one currently opened root file
		for( ; j < file->Count(); ++j )
		{
			const TEvent * entry = file->GetEntry( j );
			_stats.FindExtremes( *entry );
		}
		j = 0;
		
		/// take next file
		file = _folder.File( ++i );
	}
	cerr << endl;
}
	
	

/// Create histogram
template< typename TEvent, typename TStats >
void HistMaker<TEvent, TStats>::Create( NArray< int > & hist )
{
	if( _stats.Count() <= 0 )
		GetStats();
	
	if( _stats.Count() <= 0 )
	{
		cerr << "error, number of elements in statistic's object is 0 or less" << endl;
		throw 0;
	}
	
	cerr << "HistMaker creates histogram, goes again through all files." << endl << "Current file: ";
	for( int i = 0; i < _folder.Count(); ++i )
	{
		RootFReader< TEvent > * file = _folder.File( i );
		cerr << i << " ";
		
		for( int j = 0; j < file->Count(); ++j )
		{
			const TEvent * entry = file->GetEntry( j );
		    
			/// let TEvent decide which histogram cell should be increased
			_stats.Increment( *entry, hist );
		}
	}
	cerr << endl;
}


/// Get access to statistics
template< typename TEvent, typename TStats >
TStats * HistMaker<TEvent, TStats>::Stats()
{
	if( _stats.Count() > 0 )
	{
		return &_stats;
	}
	else
	{
		cerr << "error, statistics are not prepared yet, run GetStats() method first." << endl;
	}
	return 0;
}

#endif // _MAKEHIST_H_


