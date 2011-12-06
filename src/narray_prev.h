

#ifndef _NARRAY_H_
#define _NARRAY_H_

#include <iostream>
#include <fstream>
#include <memory>
#include "indexer.h"
using namespace std;


struct NArrayItem
{
	int i;
	int sumIdx;
};


template< typename T > class NArray
{
	enum { EMaxDim = 7 };
	
	int * _D;							/// size of dimensions
	T 	* _arr;							/// array with elements
	int _count;							/// number of all elems
	NArrayIndexerBase  * _idx;			/// class encapsulates array indexes
	int _dim;							/// number of dimensions
	
public:

	/// takes number of cells for particular dimension
	NArray( int * ds, int n );
	
	/// copy constructor is not allowed
	NArray( const NArray & /*x*/ );
	
	/// destroy 1d array, it handles all elems
	~NArray() { delete[] _arr; delete[] _D; delete _idx; }
	
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
	inline  const T & operator()( cir i, cir j = 0, cir k = 0, cir l = 0, cir m = 0, cir n = 0, cir o = 0  ) const;
	
	/// Save matrix to file
	void Save( string fname ) const;
	
	/// Load matrix from file, notice that it destroys previous data and array properties
	void Load( string fname );
	
	/// Conversion to sum array, next element is (prevoius value + this value)
	void DoSum();
	
	/// Conversion array to state before make it sum of elements
	void RedoSum();

	/// Make sum array with unique elements (1,3,3,5,7)->((0,1),(1,3),(3,5),(4,7))
	/// Notice that walking throught second array is faster then first one.
	/// This method can be run after DoSum() only.
	/// It allocates resource which must be fried later
	vector< NArrayItem > * DoUniqueElemsArray();
	
	template < typename TT > struct TWrapper
	{
		TWrapper( TT & toZero ) { toZero = static_cast<TT>(0) ; }
	};
	
	/// Removes one of the dimesnions. Notice that it takes
	/// already created array so all properties have to be set properly before.
	/// New array has decreased dimesion.
	/// d - dim index to remove
	void Squeeze( NArray<T> & result, int dim ) const;

private:

	void Recursive( NArray<T> & result, int dim, T out[EMaxDim], int trx[EMaxDim] ) const;
};




/// takes number of cells for particular dimension
template<typename T> 
NArray< T >::NArray( int * ds, int n )  : _arr( 0 ), _count( 1 ), _idx( 0 ), _dim( n )
{
	_D = new int[Dim()]; 
	IndexFactory factory;
	_idx = factory.Get( Dim(), _D );
 
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
    return _arr[ _idx->Idx(i,j,k,l,m,n,o) ];
}



template<typename T> 
const T & NArray< T >::operator()( cir i, cir j, cir k, cir l, cir m, cir n, cir o ) const
{
    return _arr[ _idx->Idx(i,j,k,l,m,n,o) ];
}



/// Save matrix to file
template<typename T> 
void NArray< T >::Save( string fname ) const
{
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
	ifstream file( fname.c_str() );

	delete[] _arr;		_arr = 0;
	delete[] _D;		_D = 0; 
	delete _idx;		_idx = 0;
	
	file >> _dim;
	IndexFactory factory;
	_D = new int[ Dim() ];
	_idx = factory.Get( Dim(), _D );
	
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
	for( int i = 1; i < _count; ++i )
		_arr[i] += _arr[i-1];
}
	


/// Conversion array to state before make it sum of elements
template<typename T> 
void NArray< T >::RedoSum()
{
	for( int i = _count - 1; i > 0; -- i )
		_arr[i] -= _arr[i-1]; 
}



/// Make sum array with unique elements (1,3,3,5,7)->((0,1),(1,3),(3,5),(4,7))
/// Notice that walking throught second array is faster then first one.
/// This method can be run after DoSum() only.
/// It allocates resource which must be fried later
template<typename T> 
vector< NArrayItem > * NArray< T >::DoUniqueElemsArray()
{
	auto_ptr<  vector< NArrayItem >  >	vec(new vector< NArrayItem >());
	
	for( int i = 1; i < _count; ++i )
	{
		if( _arr[i-1] != _arr[i] )
		{
			NArrayItem item = { i-1, _arr[i-1] };
			vec->push_back( item );
		}
	}

	return vec.release();
}



/// Removes one of the dimesnions. Notice that it takes
/// already created array so all properties have to be set up properly before.
/// New array has decreased dimesion.
/// d - dim index to remove
template<typename T> 
void NArray<T>::Squeeze( NArray<T> & result, int dim ) const
{
	T out[MaxDim()] ;
	int trx[MaxDim()];
	
	for( int i = 0;   i < MaxDim();   ++i )
	{
		TWrapper<T> zero_it( out[i] );
		trx[i] = 0;
	}

	for( int i = 0,  j = 0;   i < MaxDim();   ++i,  ++j )
	{
		if( i == dim )
			++i;
		trx[j] =  ( j >= dim  ?  1 : 0 );
	}
	
	Recursive( result, dim, out, trx );
}




template<typename T> 
void NArray<T>::Recursive( NArray<T> & result, int dim, T out[EMaxDim], int trx[EMaxDim] ) const
{
	T sum;
	TWrapper<T> zero_it( sum );

	 for( out[dim] = 0; out[dim] < Count(dim); ++out[dim] )
	 {
	 	 sum += (*this)( out[0], out[1], out[2], out[3], out[4], out[5], out[6] ) ;
	 }
	result( out[0+trx[0]], out[1+trx[1]], out[2+trx[2]], out[3+trx[3]], out[4+trx[4]], out[5+trx[5]], out[6+trx[6]] ) = sum;
	
	bool next = true;
	for( int i = 0; i < Dim(); ++i )
	{
		if( i == dim )
			continue;
		if( next )
		{
			if( out[i] < Count(i) -1 )
			{
				++out[i];
				next = false;
			}
			else
			{
				out[i] = 0;
				next = true;
			}
		}
	}
	if( next == false )
		Recursive( result, dim, out, trx );
}

#endif // _NARRAY_H_

