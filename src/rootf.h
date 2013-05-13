#ifndef _ROOTF_H_
#define _ROOTF_H_

#include <string>
#include <dirent.h>
#include <iostream>
#include <cstdlib>
#include "TFile.h"
#include "TTree.h"


using namespace std;

struct ND5Event
{
	
	Float_t Enu; 		/// energia
	Float_t nnu[3]; 	/// x, y, z kierunek pedu
	Float_t xnu, ynu; 	/// wsp. gdzie neutrino przecina plaszczyzne Z=0
	Int_t mode;	        /// reaction mode
//	UChar_t gipart; 	/// pdg neutrino mionowa 14
//	Int_t idfd; 		/// sprawdzic czy jest to wartosc 5, identyfikacje detektora
	Float_t norm; 	    /// norma neutrina
	
	ND5Event(  ) : Enu(0), xnu(0), ynu(0), mode(0), /*idfd(0),*/ norm(1)
	{
		nnu[0] = 0;
		nnu[1] = 0;
		nnu[2] = 0;
	}

	void CreateBranches( TTree * tree )
	{
		tree->Branch( "Enu", &Enu, "Enu/F" );
		tree->Branch( "nnu", &nnu, "nnu[3]/F" );
		tree->Branch( "xnu", &xnu, "xnu/F" );
		tree->Branch( "ynu", &ynu, "ynu/F" );
		tree->Branch( "mode", &mode, "mode/I" );
//		tree->Branch( "gipart", &gipart, "gipart/I" );
//		tree->Branch( "idfd", &idfd, "idfd/I" );
		tree->Branch( "norm", &ynu, "norm/F" );
	}
	
	


	
	void OpenBranches( TTree * tree )
	{
		tree->SetBranchAddress( "Enu", &Enu );
		tree->SetBranchAddress( "nnu", &nnu );
		tree->SetBranchAddress( "xnu", &xnu );
		tree->SetBranchAddress( "ynu", &ynu );
		tree->SetBranchAddress( "mode", &mode );
//		tree->SetBranchAddress( "gipart", &gipart );
//		tree->SetBranchAddress( "idfd", &idfd );
		tree->SetBranchAddress( "norm", &norm );
	}
	

};


/***********************************************************************************
 * Save short version of root file.
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
* Ordinary root file reader.
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
		if(tree==NULL) 
		  {cerr<< "tree \""<<treename<<"\" not found in file \""<<fname<<"\""<<endl;
		   exit(1);
	      }

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
	RootFolder( string directory, string treename ) : _file( 0 ), _dp( 0 ), _treename( treename )
	{
		if( directory[directory.size() - 1] != '/' )
			directory += '/';
			
		_dp = opendir( directory.c_str() );
		if(_dp==NULL)
			{cerr << "Directory \""<<directory<<"\" not found."<<endl;
			 exit (1);
		    }
		while( true )
		{
			dirent * dirp = readdir( _dp );
			if( dirp == 0 )
				break;
			string name( dirp->d_name );
			 if(   name.find( string(".root") )  !=  string::npos   )
				_fnames.push_back( directory + name );
		}
		if( _fnames.size() > 0 )
			cout << _fnames.size()<<" root files ready to open.\n";
		else
			{cerr << "No root files fount in directory \""<<directory<<"\""<<endl;
			 exit (1);
		    }
//		for( int i = 0; i < _fnames.size(); ++i )
//			cout<< _fnames[i] << endl;
	}
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
#endif // #ifndef _ROOTF_H_
