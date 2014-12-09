

#include <iostream>
#include <string>
#include "beam.h"
#include "beamRF.h"


/// funkcjonalność tej klasy sprowadza sie do ustawienia zmiennej bool 
/// na false.
/// Dzieje sie to gdy BeamRF wlasnie zwraca ostatni element z ostatniego
/// znalezionego pliku
/// Nic nie stoi na przeszkodzie aby nie przerywać pętli, wtedy BeamRF skoczy
/// do pliku pierwszego itd



int main( int argc, char * * argv )
{
	bool run = true;
	/// folder z plikami root
	/// testowalem na plikach 
	/// nu.nd5_horn250ka.1.root
	/// nu.nd5_horn250ka.10.root
	/// nu.nd5_horn250ka.100.root
	string path;
	
	if( argc > 1 )
		path = argv[1];
	else
		path = ".";		
	
	cerr << "looking for path: " << path << endl;
	
	/// drugi string to nazwa drzewa w pliku root
	params p;
	BeamRF beam( p);
	
	/// et voila, generujemy elementy na podstawie wszystkich plikow
	int ii=0;
	while( run )
	{ 
	    particle p = beam.shoot(false);
	    //if(ii++ %1000000==0 )cout<< ii<< ' '<<flush;
	    //	cout<<	p <<' '<<p.momentum()<<endl;
	}
	
	return 0;
}
