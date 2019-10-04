#ifndef MATRIX_H     // use preprocessor wrapper (Deitel p. 747)
#define MATRIX_H

#include <complex>
using std::complex;
#include <iostream>
using std::ostream;
/********************************************************************

The Matrix class stores (4, 4) complex matrices

PRIVATE
the viewpoint we take is that of a matrix M represented by a pointer
to pointers (see web-link). a "first level" pointer is allocated; it
points to an array of #1 pointers, one for each row. next, we allocate
each row-pointer to point to an array of #2 elements. in general, #1
stands for nrows and #2 for ncolumns.

PUBLIC
a copy-constructor is defined, next to a constructor with and w/o 
parameters.      
                                                                    
*********************************************************************/

class Spr;  // use forward declaration 
class TSpr;

class Matrix{
  friend Matrix ConjMatrix( Matrix );
  friend Matrix GammaMatrixGamma( Matrix );
  friend complex<double> Trace( Matrix );

  friend void Trace_product(const Matrix &M1, const Matrix &M2, complex<double> &result);
  friend Matrix operator*( double, Matrix );
  friend Matrix operator*( complex<double>, Matrix );
  friend Matrix operator*( Matrix, Matrix );
  friend Matrix operator+( Matrix, Matrix );
  friend Matrix operator-( Matrix, Matrix );
  friend ostream &operator<<( ostream &, Matrix );

  friend Spr operator*( Matrix, Spr ); 
  friend Matrix operator*( Spr, TSpr );
  friend TSpr operator*( TSpr, Matrix );

 public:
  Matrix( const Matrix & );  // copy constructor
  Matrix( complex<double> A[4][4] );     
  Matrix();
  ~Matrix();

  const Matrix &operator=( const Matrix & );

//  private:
  complex<double> M[4][4]={};
};



void Construct_GammaMatrices( Matrix &, Matrix [] );

// Matrix Dagger( double [], Matrix [] );


// // // // Nuevo

class Spr{
  complex<double> N[4]={};
public:
  Spr(const Spr &n);
  Spr(complex<double> B[4]);
  Spr();
  ~Spr();
  friend Spr operator*(complex<double> op1, Spr op2);
  friend Spr operator+(Spr op1, Spr op2);
  friend ostream &operator<<(ostream &stream, Spr n);

  friend TSpr ConjSpr1(Spr n);
  friend Spr ConjSpr2(TSpr n);
  friend complex<double> operator*(TSpr op1, Spr op2);

  friend Spr operator*(Matrix op1, Spr op2);
  friend Matrix operator*(Spr op1, TSpr op2);
};

class TSpr{
  complex<double> Nt[4]={};
public:
  TSpr(const TSpr &n);
  TSpr(complex<double> B[4]);
  TSpr();
  ~TSpr();
  friend TSpr operator*(complex<double> op1, TSpr op2);
  friend TSpr operator+(TSpr op1, TSpr op2);
  friend ostream &operator<<(ostream &stream, TSpr n);

  friend TSpr ConjSpr1(Spr n);
  friend Spr ConjSpr2(TSpr n);
  friend complex<double> operator*(TSpr op1, Spr op2);

  friend Matrix operator*(Spr op1, TSpr op2);
  friend TSpr operator*(TSpr op1, Matrix op2);

  const TSpr &operator=( const TSpr & );
};

#endif
