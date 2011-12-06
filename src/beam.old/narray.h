

#ifndef _NARRAY_H_
#define _NARRAY_H_

#include <iostream>
#include <fstream>
#include <memory>
#include <vector>

#include "generatormt.h"

using namespace std;


struct NArrayItem
{
	int i;
	int sumIdx;
};


template< typename T > class NArray
{
	typedef const int &	cir;
	enum
	{
		EMaxDim = 7
	};
	
	int	* _D;          /// size of dimensions
	T	* _arr;        /// array with elements
	int _count;          /// number of all elems
	int _dim;            /// number of dimensions
	bool _accumulated;// TODO: co dla przypadku DoUniqueElemsArray?
	
public:

	/// takes number of cells for particular dimension
	NArray( const int * ds, const int n );
	
	/// copy constructor is not allowed
	NArray( const NArray & /*x*/ );
	
	/// destroy 1d array, it handles all elems
	~NArray() { delete[] _arr; delete[] _D; }
	
	/// array max dimesnion
	const int MaxDim() const { return EMaxDim; }
	
	/// array dimesion
	const int Dim() const { return _dim; }
	
	/// number of cells for particular dimension
	const int Count( cir dim ) const { return _D[dim]; } 
	
	/// direct access to array
	T & Flat( cir i ) { return _arr[i]; }
	const T & Flat( cir i ) const { return _arr[i]; }
	cir FlatCount() const { return _count; }

	/// access to elements of the array
	inline  T & operator()( cir i, cir j = 0, cir k = 0, cir l = 0, cir m = 0, cir n = 0, cir o = 0 );
	inline  const T operator()( cir i, cir j = 0, cir k = 0, cir l = 0, cir m = 0, cir n = 0, cir o = 0  ) const;
	
	
	/// Save matrix to file
	void Save( string fname ) const;
	
	/// Load matrix from file, notice that it destroys previous data and array properties
	void Load( string fname );
	
	/// Conversion to sum array, next element is (prevoius value + this value)
	void DoSum();
	
	/// Conversion array to state before make it sum of elements
	void RedoSum();

	/// Get quasi random element from the array.
	/// It's not strict random value because of dispersion inside the histogram 
	void RandomIdx( int idx[] );

	/// Make sum array with unique elements (1,3,3,5,7)->((0,1),(1,3),(3,5),(4,7))
	/// Notice that walking throught second array is faster then first one.
	/// This method can be run after DoSum() only.
	/// It allocates resource which must be fried later
	vector< NArrayItem > * DoUniqueElemsArray();
	
private:

	/// jumps on flat array, is searching index of int x value
	const int Bisection( const int & x );
};




/// takes number of cells for particular dimension
template<typename T> 
NArray< T >::NArray( const int * ds, const int n )  : _arr( 0 ), _count( 1 ), _dim( n ), _accumulated(false)
{
	_D = new int[Dim()]; 
 
	for( int i = 0; i < Dim(); ++i )
		_count *= (_D[i] = ds[i]);

	_arr = new T[_count];
}


/// copy constructor is not allowed
template<typename T> 
NArray< T >::NArray( const NArray & /*x*/ )
{
	cerr << "copy constructor is not allowed" << endl;
	throw 0;
}
	
	

template<typename T> 
T & NArray< T >::operator()( cir i, cir j, cir k, cir l, cir m, cir n, cir o )
{
	RedoSum();
	return _arr[ i + _D[0] * ( j + _D[1] * ( k + _D[2] * ( l + _D[3] * (m + _D[4] * (n + _D[5] * o))))) ];
}



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
template<typename T> 
void NArray< T >::Save( string fname ) const
{
	// TODO: rozwiązac problem z ta czynnoscia RedoSum();
	ofstream file( fname.c_str() );
	file << Dim() << ' ';
	for( int i = 0; i < Dim(); ++i )
		file << _D[i] << ' ';
	
	file << _count << ' ';
	
	for( int i = 0; i < _count; ++i )
		file << _arr[i] << ' ';
}


	
/// Load matrix from file, notice that it destroys previous data and array properties
template<typename T> 
void NArray< T >::Load( string fname )
{
	_accumulated = false;
	ifstream file( fname.c_str() );

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



/// Conversion to sum array, next element is (prevoius value + this value)
template<typename T> 
void NArray< T >::DoSum()
{
	if( !_accumulated )
	{
		for( int i = 1; i < _count; ++i )
			_arr[i] += _arr[i-1];
		_accumulated = true;
	}
}
	


/// Conversion array to state before make it sum of elements
template<typename T> 
void NArray< T >::RedoSum()
{
	if( _accumulated )
	{
		for( int i = _count - 1; i > 0; -- i )
			_arr[i] -= _arr[i-1]; 
		_accumulated = false;
	}
}



/// Make sum array with unique elements (1,3,3,5,7)->((0,1),(1,3),(3,5),(4,7))
/// Notice that walking throught second array is faster then first one.
/// This method can be run after DoSum() only.
/// It allocates resource which must be fried later
template<typename T> 
vector< NArrayItem > * NArray< T >::DoUniqueElemsArray()
{
      if( _count < 1 )
            throw string("Zero array size");
            
	auto_ptr<  vector< NArrayItem >  >	vec(new vector< NArrayItem >());
      
      /// add first element here because the loop condition omits him
      NArrayItem first = { 0, _arr[0] };
      vec->push_back( first );
	
	for( int i = 1; i < _count; ++i )
	{
		if( _arr[i-1] != _arr[i] )
		{
			NArrayItem item = { i, _arr[i] };
			vec->push_back( item );
		}
	}

	return vec.release();
}


/// jumps on flat array, is searching index of int x value
template<typename T> 
const int NArray< T >::Bisection( const int & x )
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
template<typename T> 
void NArray< T >::RandomIdx( int idx[] )
{
	DoSum();
	int x = frandom() * _arr[_count - 1];
	int i = Bisection( x );
	
	for( int a = 0; a < _dim; ++a )
	{
		idx[a] = i%_D[a];
		i/= _D[a];
	}
	
	/*
	idx[4] = i/_f[3];
	idx[3] = (i%_f[3]) / _f[2];
	idx[2] = ((i%_f[3]) % _f[2]) / _f[1];
	idx[1] = (((i%_f[3]) % _f[2]) % _f[1]) / _f[0];
	idx[0] = (((i%_f[3]) % _f[2]) % _f[1]) % _f[0];
	*/
}


///// Get element from the array.
///// idx handles 
//template<typename T> 
//void NArray< T >::Idx( int idx[],int i )
//{
	//for( int a = 0; a < _dim; ++a )
	//{
		//idx[a] = i%_D[a];
		//i/= _D[a];
	//}
	
	///*
	//idx[4] = i/_f[3];
	//idx[3] = (i%_f[3]) / _f[2];
	//idx[2] = ((i%_f[3]) % _f[2]) / _f[1];
	//idx[1] = (((i%_f[3]) % _f[2]) % _f[1]) / _f[0];
	//idx[0] = (((i%_f[3]) % _f[2]) % _f[1]) % _f[0];
	//*/
//}

#endif		/// _NARRAY_H_
