

#ifndef _NARRAY_H_
#define _NARRAY_H_

#include <iostream>
#include <fstream>
#include <memory>
#include <vector>

#include "generatormt.h"

using namespace std;

typedef const int &   cir;

enum
{
	EMaxDim = 7
};

template< typename T > class NArray
{
	int	* _D;          /// size of dimensions
	T	* _arr;        /// array with elements
	int   _count;        /// number of all elems
	int   _dim;          /// number of dimensions
	bool  _accumulated;  /// element is the sum of his & previous value

public:

	/// takes number of cells for particular dimension
	NArray( const int * ds, const int n );
	
	/// copy constructor is not allowed
	NArray( const NArray & /*x*/ );
	
	/// destroy 1d array, it handles all elems
	~NArray() { delete[] _arr; delete[] _D; }
	
	/// array max dimesnion
	int MaxDim() const { return EMaxDim; }
	
	/// array dimesion
	int Dim() const { return _dim; }
	
	/// number of cells for particular dimension
	int Count( cir dim ) const { return _D[dim]; } 
	
	/// direct access to array
	T & Flat( cir i ) { return _arr[i]; }
	const T & Flat( cir i ) const { return _arr[i]; }
	cir FlatCount() const { return _count; }

	/// access to elements of the array
	inline  T & operator()( cir i, cir j = 0, cir k = 0, cir l = 0, cir m = 0, cir n = 0, cir o = 0 );
	inline  const T operator()( cir i, cir j = 0, cir k = 0, cir l = 0, cir m = 0, cir n = 0, cir o = 0  ) const;
	
	/// Save matrix to file
	void Save( ofstream & file_out );
	
	/// Load matrix from file, notice that it destroys previous data and array properties
	void Load( ifstream & file );
	
	/// Conversion to sum array, next element is (prevoius value + this value)
	void DoSum();
	
	/// Conversion array to state before make it sum of elements
	void UndoSum();

	/// Get quasi random element from the array.
	/// It's not strict random value because of dispersion inside the histogram 
	void RandomIdx( int idx[] );
	
	/// Check if object is after accumulation or not
	/// If DoSum method was run
	bool SumDone() const { return _accumulated; }
	
private:

	/// jumps on flat array, is searching index of int x value
	int Bisection( cir x );
};




/// takes number of cells for particular dimension
template<typename T> NArray< T >::NArray( const int * ds, const int n )
: _arr( 0 ), _count( 1 ), _dim( n ), _accumulated(false)
{
	_D = new int[Dim()]; 
 
	for( int i = 0; i < Dim(); ++i )
		_count *= (_D[i] = ds[i]);

	_arr = new T[_count];
}



/// copy constructor is not allowed
template<typename T> NArray< T >::NArray( const NArray & /*x*/ )
{
	cerr << "copy constructor is not allowed" << endl;
	throw 0;
}



/// direct access to the element
template<typename T> 
T & NArray< T >::operator()( cir i, cir j, cir k, cir l, cir m, cir n, cir o )
{
	/// notice that UndoSum() has N complicity, it may take time
	UndoSum();
	
	return _arr[ i + _D[0] * ( j + _D[1] * ( k + _D[2] * ( l + _D[3] * (m + _D[4] * (n + _D[5] * o))))) ];
}



/// direct access to the element
template<typename T> 
const T NArray< T >::operator()( cir i, cir j, cir k, cir l, cir m, cir n, cir o ) const
{   
	int idx =  i + _D[0] * ( j + _D[1] * ( k + _D[2] * ( l + _D[3] * (m + _D[4] * (n + _D[5] * o))))) ;  
	
	if( idx && _accumulated )
		return _arr[idx]-_arr[idx-1];
	else
		return _arr[idx];
}



/// Save matrix to file
template<typename T> void NArray< T >::Save( ofstream & file_out )
{
	UndoSum();
	file_out << Dim() << ' ';
	for( int i = 0; i < Dim(); ++i )
		file_out << _D[i] << ' ';
	
	file_out << _count << ' ';
	
	for( int i = 0; i < _count; ++i )
		file_out << _arr[i] << ' ';
}


	
/// Load matrix from file, notice that it destroys previous data and array properties
template<typename T> void NArray< T >::Load( ifstream & file )
{
	_accumulated = false;

	delete[] _arr;		_arr = 0;
	delete[] _D;		_D = 0; 
	
	file >> _dim;

	_D = new int[ Dim() ];
	
	for( int i = 0; i < Dim(); ++i )
		file >> _D[i];
	
	file >> _count;
	_arr = new T[_count];
	
	for( int i = 0; i < _count; ++i )
		file >> _arr[i];
}



/// Conversion to sum array, element is sum of his and previous value
template<typename T> void NArray< T >::DoSum()
{
	if( !_accumulated )
	{
		for( int i = 1; i < _count; ++i )
			_arr[i] += _arr[i-1];
		_accumulated = true;
	}
}



/// Conversion array to state before make it sum of elements
template<typename T> void NArray< T >::UndoSum()
{
	if( _accumulated )
	{
		for( int i = _count - 1; i > 0; --i )
			_arr[i] -= _arr[i-1]; 
		_accumulated = false;
	}
}


/// jumps on flat array, is searching index of int x value
template<typename T> int NArray< T >::Bisection( cir x )
{
	register int a = 0;
	register int b = _count-1;
	register int i;

	while( b > a )
	{
		i = (a+b)/2;

		if( x<_arr[i] )
			b = i;
		else
			a = i + 1;
	}

	return a;
}


/// Get quasi random element from the array.
/// It's not strict random value because of dispersion inside the histogram 
template<typename T> void NArray< T >::RandomIdx( int idx[] )
{
	DoSum();
	int x = frandom() * _arr[_count - 1];
	int i = Bisection( x );
	
	for( int a = 0; a < _dim; ++a )
	{
		idx[a] = i%_D[a];
		i/= _D[a];
	}
}

#endif		/// _NARRAY_H_
