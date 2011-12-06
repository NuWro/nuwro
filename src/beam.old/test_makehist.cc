

#include <memory>
#include "makeHist.h"
#include "fstream"
#include "beam.h"
#include "particle.h"

using namespace std;

#define HISTOGRAM_FILE_PATH "histout.txt"
#define SQUEEZED_HISTOGRAM_FILE_PATH "squeezed_histout.txt"
#define STATISTICS_FILE_PATH "statsout.txt"


// it takes folder with *.root files
void In( const string & path )
{
	int dimsizes[6] = { 5, 5, 5, 5, 5, 5 };
	NArray< int >  histogram( dimsizes, 6 );
	for( int i = 0; i < histogram.FlatCount(); ++i )
		histogram.Flat( i ) = 0;
	
	HistMaker< Nd280Element, Nd280Statistics >  histmaker( path, string("h3002") );
	histmaker.Create( histogram );
	histogram.Save( HISTOGRAM_FILE_PATH );
	
	ofstream statsfile( STATISTICS_FILE_PATH );
	statsfile << histmaker.Stats()->SummaryStr();
}



void Out()
{
/*
	int dimsizes[6] = { 5, 5, 5, 5, 5, 5 };
	NArray< int >  histogram( dimsizes, 6 );
	for( int i = 0; i < histogram.FlatCount(); ++i )
		histogram.Flat( i ) = 0;
		
	ifstream file( STATISTICS_FILE_PATH );
	ND5Statistics stats;
	if( stats.Load( file ) )
	{			
		histogram.Load( HISTOGRAM_FILE_PATH );
		histogram.DoSum();
		BeamHist< int >  beam( &histogram, &stats );
		
		const int max = 100;
		cout << "\ngenerate banch of " << max << " particles\n";
		
		for( int ii = 0; ii < max; ++ii )
		{
			particle p = beam.shoot();
			cout << p << endl;
		}
		cout << "end of particle generating\n\n";
	}	
*/
}


//void MatrixSqueezeTest( )
//{
	//int dimsizes_out[2] = { 2, 2 };
	//int dimsizes_in[3] = { 2, 2, 2 };
	//NArray< int > out( dimsizes_out, 2 );
	//NArray< int > in( dimsizes_in, 3 );
	
	//cout << endl << endl;
	//cout << "squeeze test - init\n";
	//for( int i = 0; i < dimsizes_in[0]; ++i )
		//for( int j = 0; j < dimsizes_in[1]; ++j )
			//for( int k = 0; k < dimsizes_in[2]; ++k )
			//{
				//in( i, j, k ) = 1;
			//}
	
	//cout << "squeeze test - start\n";
	
	//in.Squeeze( out, 2 );
	
	//cout << "squeeze test - results\n";
	
	//for( int i = 0; i < dimsizes_out[0]; ++i )
		//for( int j = 0; j < dimsizes_out[1]; ++j )
		//{
			//cout << "(" << i << "," << j << ")=" << out(i,j) << "\n";
		//}
//}


//void DecDimensions()
//{
	//enum { maxdim = 6 };
	//int dimsizes[maxdim] = { 5, 5, 5, 5, 5, 5 };
	//auto_ptr< NArray< int > > 		old_hist(  new	NArray< int >(dimsizes, maxdim)  );
	//for( int i = 0; i < old_hist->FlatCount(); ++i )
		//old_hist->Flat( i ) = 0;
		
	//cout << "Loading histogram file, " << HISTOGRAM_FILE_PATH << '\n';
		
	//old_hist->Load( HISTOGRAM_FILE_PATH );
	
	///// take a look into <<particle BeamHist<T>::shoot()>>
	///// where first dimension is energy profile. Squeeze all other
	///// dimensions and keep first one.
	
	//cout << "Squeezing histogram from " << maxdim << " dimmensions to one (Energy)\n";
	
	//for(  int i = maxdim - 1;  i >= 1;  --i  )
	//{
		//auto_ptr< NArray< int > >		new_hist(  new NArray< int >(dimsizes, i)  );
		//old_hist->Squeeze( *new_hist, i );
		//old_hist = new_hist;
	//}	
	
	//cout << "Save one dimmensional matrix to " << SQUEEZED_HISTOGRAM_FILE_PATH << '\n';
	//old_hist->Save( SQUEEZED_HISTOGRAM_FILE_PATH );
//}


int main( int argc, char * * argv )
{
	cout << "\n\n[BEGIN]\n";
	string path;
	
	/// NArray dimensions' squeeze test
	//MatrixSqueezeTest( );
	
	if( argc > 1 )
		path = argv[1];
	else
		path = ".";		
	
	cerr << "looking for path: \"" << path << "\"" << endl;
	
	/// on empty 6D matrix,
	/// using data from root files for h3002
	/// 1. prepare histogram
	/// 2. save it as HISTOGRAM_FILE_PATH
	/// 3. prepare statistics, like min and max of energy
	/// 4. save it as STATISTICS_FILE_PATH
	In( path );
	
	/// statistics and histogram load from files (see above)
	/// create beam object
	/// shoot some particles
	Out();
	
	/// squeeze histogram matrix from 6D to 1D (Energy)
	/// save it in SQUEEZED_HISTOGRAM_FILE_PATH
	///DecDimensions();
	
	cout << "\n[END]\n\n\n";
	return 0;
}
