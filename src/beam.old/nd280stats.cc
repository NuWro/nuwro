

#include <sstream>
#include "nd280stats.h"


/// object is created and initialised from the file
Nd280Statistics::Nd280Statistics( string fname )
{
	ifstream in( fname.c_str() );
	if( Load( in ) == false )
		throw "err: Nd280Statistics can't load the file\n";
}


/// check if x element is smaller then current minimum,
/// if yes then overwrite old value
void Nd280Statistics::CheckMin( const Nd280Element & x )
{
	/// searching of minimum values
	if( _min.xnu > x.xnu )
		_min.xnu = x.xnu;
		
	if( _min.ynu > x.ynu )
		_min.ynu = x.ynu;
	
	for( int i = 0; i < 3; ++i )
	{
		double p = x.nnu[i] * x.Enu * 1000;
		if( _min.nnu[i] > p ) 
			_min.nnu[i] = p;
	}	
}


/// check if x element is larger then current maximum,
/// if yes then overwrite old value	
void Nd280Statistics::CheckMax( const Nd280Element & x )
{
	/// searching of maximum values
	if( _max.xnu < x.xnu )
		_max.xnu = x.xnu;
		
	if( _max.ynu < x.ynu )
		_max.ynu = x.ynu;

	for( int i = 0; i < 3; ++i )
	{
		double p = x.nnu[i] * x.Enu * 1000;
		if( _max.nnu[i] < p ) 
			_max.nnu[i] = p;
	}
}


/// use it at the begining, run with first element
void Nd280Statistics::Init( const Nd280Element & entry )
{
	_min = entry;
	_max = entry;
	
	double E = _min.Enu*1000;
	_min.nnu[0] *= E;
	_min.nnu[1] *= E;
	_min.nnu[2] *= E;
	_max.nnu[0] *= E;
	_max.nnu[1] *= E;
	_max.nnu[2] *= E;		
}


/// return the string ready for saving in a file
string Nd280Statistics::SummaryStr() const
{
	stringstream out;
	out << "extremes\n";
	out << "min.xnu ";						out << _min.xnu;		out << '\n';
	out << "min.ynu ";						out << _min.ynu;		out << '\n';
	out << "min.nnu0 ";					out << _min.nnu[0];	out << '\n';
	out << "min.nnu1 ";					out << _min.nnu[1];	out << '\n';
	out << "min.nnu2 ";					out << _min.nnu[2]; 	out << '\n';
	out << "max.xnu ";					out << _max.xnu;		out << '\n';
	out << "max.ynu ";					out << _max.ynu;		out << '\n';
	out << "max.nnu0 ";					out << _max.nnu[0];	out << '\n';
	out << "max.nnu1 ";					out << _max.nnu[1];	out << '\n';
	out << "max.nnu2 ";					out << _max.nnu[2]; 

	return out.str();
}


/// fill the object using the file (input file stream)
bool Nd280Statistics::Load( ifstream & file )
{
	string name;
	if( !file.eof() )
		file >> name;
	bool result = false;
	double value = 0;
	
	if( name == "extremes" )
	{
		_count = -1;
		file >> name >> _min.xnu;
		file >> name >> _min.ynu;
		file >> name >> _min.nnu[0];
		file >> name >> _min.nnu[1];
		file >> name >> _min.nnu[2];
		file >> name >> _max.xnu;	
		file >> name >> _max.ynu;	
		file >> name >> _max.nnu[0];
		file >> name >> _max.nnu[1];
		file >> name >> _max.nnu[2]; 
		if( name == "max.nnu2" )
			result = true;
	}

	return result;
}
