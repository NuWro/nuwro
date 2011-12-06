
#ifndef _INDEXER_H_
#define _INDEXER_H_

#include <iostream>
using namespace std;

typedef int const & cir;
typedef int * const & ptrRef;


class NArrayIndexerBase
{
public:
	virtual const int Idx( cir i ) const = 0;
	virtual const int Idx( cir i, cir j ) const = 0;
	virtual const int Idx( cir i, cir j, cir k ) const = 0;
	virtual const int Idx( cir i, cir j, cir k, cir l ) const = 0;
	virtual const int Idx( cir i, cir j, cir k, cir l, cir m ) const = 0;
	virtual const int Idx( cir i, cir j, cir k, cir l, cir m, cir n ) const = 0;
	virtual const int Idx( cir i, cir j, cir k, cir l, cir m, cir n, cir o ) const = 0;
};


template< int N > class NArrayIndexer : public NArrayIndexerBase
{
	ptrRef _D;
	
public:
	NArrayIndexer( ptrRef d ) : _D(d) {}
	
	const int Idx( cir i ) const { return i; }
	const int Idx( cir i, cir j ) const  { return Idx(i) * _D[0] + j; }
	const int Idx( cir i, cir j, cir k ) const  { return Idx(i,j) * _D[1] + k; }
	const int Idx( cir i, cir j, cir k, cir l ) const  { return Idx(i,j,k) * _D[2] + l; }
	const int Idx( cir i, cir j, cir k, cir l, cir m ) const  { return Idx(i,j,k,l) * _D[3] + m; }
	const int Idx( cir i, cir j, cir k, cir l, cir m, cir n ) const  { return Idx(i,j,k,l,m) * _D[4] + m; }
	const int Idx( cir i, cir j, cir k, cir l, cir m, cir n, cir o ) const  { return Idx(i,j,k,l,m,n) * _D[5] + o; }
};

template <> class NArrayIndexer< 6 > : public NArrayIndexerBase
{
	ptrRef _D;
	
public:
	NArrayIndexer( ptrRef d ) : _D(d) {}
	
	const int Idx( cir i ) const { return i; }
	const int Idx( cir i, cir j ) const  { return Idx(i) * _D[0] + j; }
	const int Idx( cir i, cir j, cir k ) const  { return Idx(i,j) * _D[1] + k; }
	const int Idx( cir i, cir j, cir k, cir l ) const  { return Idx(i,j,k) * _D[2] + l; }
	const int Idx( cir i, cir j, cir k, cir l, cir m ) const  { return Idx(i,j,k,l) * _D[3] + m; }
	const int Idx( cir i, cir j, cir k, cir l, cir m, cir n ) const  { return Idx(i,j,k,l,m) * _D[4] + m; }
	const int Idx( cir i, cir j, cir k, cir l, cir m, cir n, cir /*o*/ ) const  { return Idx(i,j,k,l,m,n); }
};

template <> class NArrayIndexer< 5 > : public NArrayIndexerBase
{
	ptrRef _D;
	
public:
	NArrayIndexer( ptrRef d ) : _D(d) {}
	
	const int Idx( cir i ) const { return i; }
	const int Idx( cir i, cir j ) const  { return Idx(i) * _D[0] + j; }
	const int Idx( cir i, cir j, cir k ) const  { return Idx(i,j) * _D[1] + k; }
	const int Idx( cir i, cir j, cir k, cir l ) const  { return Idx(i,j,k) * _D[2] + l; }
	const int Idx( cir i, cir j, cir k, cir l, cir m ) const  { return Idx(i,j,k,l) * _D[3] + m; }
	const int Idx( cir i, cir j, cir k, cir l, cir m, cir /*n*/ ) const  { return Idx(i,j,k,l,m); }
	const int Idx( cir i, cir j, cir k, cir l, cir m, cir /*n*/, cir /*o*/ ) const  { return Idx(i,j,k,l,m); }
};

template <> class NArrayIndexer< 4 > : public NArrayIndexerBase
{
	ptrRef _D;
	
public:
	NArrayIndexer( ptrRef d ) : _D(d) {}
	
	const int Idx( cir i ) const { return i; }
	const int Idx( cir i, cir j ) const  { return Idx(i) * _D[0] + j; }
	const int Idx( cir i, cir j, cir k ) const  { return Idx(i,j) * _D[1] + k; }
	const int Idx( cir i, cir j, cir k, cir l ) const  { return Idx(i,j,k) * _D[2] + l; }
	const int Idx( cir i, cir j, cir k, cir l, cir /*n*/ ) const  { return Idx(i,j,k,l); }
	const int Idx( cir i, cir j, cir k, cir l, cir /*n*/, cir /*n*/ ) const  { return Idx(i,j,k,l); }
	const int Idx( cir i, cir j, cir k, cir l, cir /*n*/, cir /*n*/, cir /*o*/ ) const  { return Idx(i,j,k,l); }
};

template <> class NArrayIndexer< 3 > : public NArrayIndexerBase
{
	ptrRef _D;
	
public:
	NArrayIndexer( ptrRef d ) : _D(d) {}
	
	const int Idx( cir i ) const { return i; }
	const int Idx( cir i, cir j ) const  { return Idx(i) * _D[0] + j; }
	const int Idx( cir i, cir j, cir k ) const  { return Idx(i,j) * _D[1] + k; }
	const int Idx( cir i, cir j, cir k, cir /*l*/ ) const  { return Idx(i,j,k); }
	const int Idx( cir i, cir j, cir k, cir /*l*/, cir /*n*/ ) const  { return Idx(i,j,k); }
	const int Idx( cir i, cir j, cir k, cir /*l*/, cir /*n*/, cir /*n*/ ) const  { return Idx(i,j,k); }
	const int Idx( cir i, cir j, cir k, cir /*l*/, cir /*n*/, cir /*n*/, cir /*o*/ ) const  { return Idx(i,j,k); }
};

template <> class NArrayIndexer< 2 > : public NArrayIndexerBase
{
	ptrRef _D;
	
public:
	NArrayIndexer( ptrRef d ) : _D(d) {}
	
	const int Idx( cir i ) const { return i; }
	const int Idx( cir i, cir j ) const  { return Idx(i) * _D[0] + j; }
	const int Idx( cir i, cir j, cir /*k*/ ) const  { return Idx(i,j); }
	const int Idx( cir i, cir j, cir /*k*/, cir /*l*/ ) const  { return Idx(i,j); }
	const int Idx( cir i, cir j, cir /*k*/, cir /*l*/, cir /*n*/ ) const  { return Idx(i,j); }
	const int Idx( cir i, cir j, cir /*k*/, cir /*l*/, cir /*n*/, cir /*n*/ ) const  { return Idx(i,j); }
	const int Idx( cir i, cir j, cir /*k*/, cir /*l*/, cir /*n*/, cir /*n*/, cir /*o*/ ) const  { return Idx(i,j); }
};

template <> class NArrayIndexer< 1 > : public NArrayIndexerBase
{
	ptrRef _D;
	
public:
	NArrayIndexer( ptrRef d ) : _D(d) {}
	
	const int Idx( cir i ) const { return i; }
	const int Idx( cir i, cir /*j*/ ) const  { return Idx(i); }
	const int Idx( cir i, cir /*j*/, cir /*k*/ ) const  { return Idx(i); }
	const int Idx( cir i, cir /*j*/, cir /*k*/, cir /*l*/ ) const  { return Idx(i); }
	const int Idx( cir i, cir /*j*/, cir /*k*/, cir /*l*/, cir /*n*/ ) const  { return Idx(i); }
	const int Idx( cir i, cir /*j*/, cir /*k*/, cir /*l*/, cir /*n*/, cir /*n*/ ) const  { return Idx(i); }
	const int Idx( cir i, cir /*j*/, cir /*k*/, cir /*l*/, cir /*n*/, cir /*n*/, cir /*o*/ ) const  { return Idx(i); }
};


class IndexFactory
{
public:
	NArrayIndexerBase * Get( int n, ptrRef d )
	{
		if( n > 7  ||  n < 1 )
		{
			cerr << n << " is wrong array dim" << endl;
			throw 0;
		}
		
		NArrayIndexerBase * ptr = 0;
		switch( n )
		{
			case 1 : ptr = new NArrayIndexer< 1 >( d ); break;
			case 2 : ptr = new NArrayIndexer< 2 >( d ); break;
			case 3 : ptr = new NArrayIndexer< 3 >( d ); break;
			case 4 : ptr = new NArrayIndexer< 4 >( d ); break;
			case 5 : ptr = new NArrayIndexer< 5 >( d ); break;
			case 6 : ptr = new NArrayIndexer< 6 >( d ); break;
			default:
				ptr = new NArrayIndexer< 7 >( d );
		}
		return ptr;
	} 
		
};


#endif // _INDEXER_H_

