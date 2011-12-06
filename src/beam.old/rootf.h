#ifndef _ROOTF_H_
#define _ROOTF_H_

#include <string>
#include <dirent.h>
#include <iostream>
#include "TFile.h"
#include "TTree.h"


using namespace std;


/***********************************************************************************
 * Save data to the root file.
 */
template< typename T> class RootFSaver
{
	TFile	* 	fH;
	TTree	*	tree;
	T			event;
	
public:
	RootFSaver( string fname, string tname, string treeDescript ) : fH( 0 ), tree( 0 ) 
		{
			fH =  new TFile( fname.c_str(), "recreate" ) ;
			tree = new TTree( tname.c_str(), treeDescript.c_str() );
			event.CreateBranches( tree );
		}
	
	RootFSaver( string fname, string tname ) : fH( 0 ), tree( 0 ) 
		{
			fH =  new TFile( fname.c_str(), "recreate" ) ;
			tree = new TTree( tname.c_str(), fname.c_str() );
			event.CreateBranches( tree );
		}
		
	~RootFSaver()
	{
		Close();
	}
	
	void Close()
	{
		if( fH != 0 )
		{
			fH->Write();
			fH->Close(); 
			delete fH;
			fH = 0;
		}
	}
	
	void Append( const T & entry )
	{
		event = entry;
		tree->Fill();
	}
	
};


/***********************************************************************************
* root file reader.
*/
template< typename T> class RootFReader
{
	TFile	* 	fH;
	TTree	*	tree;
	T			event;
	int			count;
	
public:	

	RootFReader( string fname, string treename ) : fH( 0 ), tree( 0 ), count( 0 )	
	{	
		fH = new TFile( fname.c_str() );
		tree = static_cast<TTree*>(  fH->Get(treename.c_str())  );
		event.OpenBranches( tree );
		count = tree->GetEntries();
	}
	
	~RootFReader()
	{
		Close();
	}
	
	void Close()
	{
		if( fH != 0 )
		{
			fH->Close(); 
			delete fH;
			fH = 0;
		}		
	}
	
	const int Count() { return count; }
	
	const T * GetEntry( const int & idx )
	{
		tree->GetEntry( idx );
		return &event;
	}
	
};


/***********************************************************************************
* Object of this class looks for root files in the selected path.
* The second input variable, called treename, describes name of the tree inside
* a root file.
*/
template<typename T> class RootFolder
{
	T 		*	_file;
	vector<string> 	_fnames;
	DIR 	*	_dp;
	string	_treename;

public:
	RootFolder( string path, string treename );
	
	~RootFolder() 
	{
		delete _file;
		if( _dp )
			closedir( _dp );
	}

	const int Count() const { return _fnames.size(); }
	
	T * File( int i )
	{
		delete _file;
		_file = 0;
		if( i >= 0 && i < _fnames.size() )
			_file = new T(  _fnames[i], _treename );
		return _file;
	}
};


template<typename T>
RootFolder<T>::RootFolder( string path, string treename ) : _file( 0 ), _dp( 0 ), _treename( treename )
{
	if( path[path.size() - 1] != '/' )
		path += '/';
		
	_dp = opendir( path.c_str() );
	
	while( true )
	{
		dirent * dirp = readdir( _dp );
		if( dirp == 0 )
			break;
		string name( dirp->d_name );
		 if(   name.find( string(".root") )  !=  string::npos   )
			_fnames.push_back( path + name );
	}
	if( _fnames.size() > 0 )
		cout << "List of root files ready to open:\n";
	for( int i = 0; i < _fnames.size(); ++i )
		cout<< i+1 << '\t' << _fnames[i] << endl;
}


#endif // _ROOTF_H_
