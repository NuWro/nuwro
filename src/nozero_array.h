

#ifndef _NOZERO_ARRAY_H_
#define _NOZERO_ARRAY_H_

#include <vector>
#include "narray.h"

using namespace std;
typedef const int & cir;

template< typename T >  struct NArrayItem
{
	int idx;
	T value;
	NArrayItem( cir i, const T & v ) : idx(i), value(v) {}
	NArrayItem() : idx(-1) {}
};


template< typename T > ofstream & operator<<( ofstream & out, const NArrayItem< T > & item )
{
	out << item.idx << ' ' << item.value;
	return out;
}

template< typename T > ifstream & operator>>( ifstream & in, NArrayItem< T > & item )
{
	in >> item.idx >> item.value;
	return in;
}





template< typename T> class NonZeroArray
{
	vector< NArrayItem< T > >   _arr;
	vector< int >               _D;

	/// Make sum array with unique elements (1,3,3,5,7)->((0,1),(1,3),(3,5),(4,7))
	/// so the example neglects double value 3.
	/// Notice that walking throught second array is faster then first one.
	/// This method can be run after DoSum() only.
	void Convert( const NArray<T> & narray );

public:
	NonZeroArray( NArray<T> & narray );
	
	/// default constructor
	NonZeroArray() {}
	
	/// Save matrix to file
	void Save( ofstream & file_out ) const;
	
	/// Load matrix from file,
	/// notice that it destroys previous data and array properties
	void Load( ifstream & file );
	
	/// matrix dimesion
	int Dim() const { return _D.size(); }
	
	/// number of cells for particular dimension
	int Count( cir dim ) const { return _D[dim]; } 
	
	/// access to elements of the array
	inline  NArrayItem< T > & operator()( cir i ) { return _arr[i]; }
	inline  const NArrayItem< T > & operator()( cir i ) const { return _arr[i]; }

	/// array max dimesnion
	int MaxDim() const { return EMaxDim; }
	
	/// Get quasi random element from the array.
	/// It's not strict random value because of dispersion inside the histogram
	void RandomIdx( int idx[] );
	
	/// jumps on flat array, is searching index of int x value
	int Bisection( cir x );
};



template< typename T> NonZeroArray<T>::NonZeroArray( NArray<T> & narray )
{
	if( narray.SumDone() == false )
	{
		cout << "Notice that NonZeroArray did conversion on NArray object using his DoSum() method" << endl;
		cout << "You can revert it with NArray< T >::UndoSum()" << endl;
		narray.DoSum();
	}
	Convert( narray );
	_D.resize( narray.Dim() );
	for( int i = 0; i < _D.size(); ++i )
		_D[i] = narray.Count(i);
}



/// Make sum array with unique elements (1,3,3,5,7)->((0,1),(1,3),(3,5),(4,7))
/// so the example neglects double value 3.
/// Notice that walking throught second array is faster then first one.
/// This method can be run after DoSum() only.
template<typename T> void NonZeroArray< T >::Convert( const NArray<T> & narray )
{
      if( narray.FlatCount() < 1 )
            throw string("Zero array size");

      /// add first element here because the loop condition omits him
      NArrayItem< T > first( 0, narray.Flat(0) );
      _arr.push_back( first );

	for( int i = 1; i < narray.FlatCount(); ++i )
	{
		if( narray.Flat(i-1) != narray.Flat(i) )
		{
			NArrayItem< T > item( i, narray.Flat(i) );
			_arr.push_back( item );
		}
	}
}



/// Save matrix to file
template<typename T> void NonZeroArray< T >::Save( ofstream & file_out ) const
{
	file_out << endl << string("nozero_array ");
	file_out << Dim() << ' ';
	for( int i = 0; i < Dim(); ++i )
		file_out << _D[i] << ' ';
	
	file_out << _arr.size() << ' ';
	
	for( int i = 0; i < _arr.size(); ++i )
		file_out << _arr[i] << ' ';
}


	
/// Load matrix from file, notice that it destroys previous data and array properties
template<typename T> void NonZeroArray< T >::Load( ifstream & file )
{
	string nozero_array;
	file >> nozero_array;
	if( nozero_array == string("nozero_array") )
	{
		int dim = 0;
		file >> dim;
		_D.resize( dim );
		
		for( int i = 0; i < Dim(); ++i )
		{
			file >> _D[i];
		}
		
		int count;
		file >> count;
		_arr.resize( count );
		
		for( int i = 0; i < _arr.size(); ++i )
			file >> _arr[i];
	}
	else
	{
		_arr.resize(0);
		_D.resize(0);
		throw "err: NonZeroArray<T> doesn't support load for this file";
	}
}



/// jumps on flat array, is searching index of int x value
template<typename T> int NonZeroArray< T >::Bisection( cir x )
{
	register int a = 0;
	register int b = _arr.size()-1;
	register int i;

	while( b > a )
	{
		i = (a+b)/2;

		if( x<_arr[i].value )
			b = i;
		else
			a = i + 1;
	}

	return _arr[a].idx;
}


/// Get quasi random element from the array.
/// It's not strict random value because of dispersion inside the histogram 
template<typename T> void NonZeroArray< T >::RandomIdx( int idx[] )
{
	int x = frandom() * _arr[_arr.size() - 1].value;
	int i = Bisection( x );
	
	for( int a = 0; a < _D.size(); ++a )
	{
		idx[a] = i%_D[a];
		i/= _D[a];
	}
}

#endif	// _NOZERO_ARRAY_H_
