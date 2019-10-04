#include <iostream>
#include <complex>

using namespace std;

#include "Matrix.h"


// // // 
// Elements of a Matrix M are called as M.[i][j] with i->row and j->column
// By default a Matrix is filled with zeroes everywhere
// // // 


Matrix::Matrix(const Matrix &m)
{

  for(int i=0; i<4; i++){
    for(int j=0; j<4; j++){
      M[i][j] = m.M[i][j];
    } 
  }  
  
//   cout << "Copy constructor Matrix" << endl;

}


Matrix::Matrix()
{

}


Matrix::Matrix(complex<double> A[4][4])
{
  
  
  for(int i=0; i<4; i++){
    for(int j=0; j<4; j++){
      M[i][j] = A[i][j];
    }
  }
  
  //cout << "Constructor with argts Matrix" << endl;

}


Matrix::~Matrix()
{

  //cout << "Kaputt Matrix!" << endl;

}


/* 
In order to assign one matrix (right) to another ('this'), one must overload the
assignment operator for the Matrix class.  For explanation about the used syntax, we 
refer to p.848 of Deitel, "C, how to program. 5th edition" 
*/

const Matrix &Matrix::operator=( const Matrix &right)
{

  if ( &right != this)
    {

      for(int i=0; i<4; i++){
	for(int j=0; j<4; j++){
	  M[i][j] = right.M[i][j];
    	}
      }

    }

  return *this;

}


Matrix ConjMatrix(Matrix m)
{
  
  Matrix temp;

  for(int i=0; i<4; i++){
    for(int j=0; j<4; j++){ 
      temp.M[i][j] = conj(m.M[j][i]); 
    }
  }

  return temp; 

}

Matrix GammaMatrixGamma(Matrix m)
{

  Matrix temp;

  temp.M[0][0] = m.M[0][0];
  temp.M[0][1] = m.M[0][1];
  temp.M[0][2] = -m.M[0][2];
  temp.M[0][3] = -m.M[0][3];

  temp.M[1][0] = m.M[1][0];
  temp.M[1][1] = m.M[1][1];
  temp.M[1][2] = -m.M[1][2];
  temp.M[1][3] = -m.M[1][3];

  temp.M[2][0] = -m.M[2][0];
  temp.M[2][1] = -m.M[2][1];
  temp.M[2][2] = m.M[2][2];
  temp.M[2][3] = m.M[2][3];

  temp.M[3][0] = -m.M[3][0];
  temp.M[3][1] = -m.M[3][1];
  temp.M[3][2] = m.M[3][2];
  temp.M[3][3] = m.M[3][3];

  return temp;

}


complex<double> Trace(Matrix A){

  complex<double> temp;

  for(int i=0; i<4; i++){
    
    temp += A.M[i][i]; 
    
  }

  return temp;

}

void Trace_product(const Matrix &M1, const Matrix &M2, complex<double> &result)
{
	//We use this function to remove some matrix multiplication, where 4*4 loops would be done to do the product, and then an additional 4 to compute the trace of the result
 	for (int i = 0 ; i < 4 ; i++)
	{
		for (int j=0 ; j < 4 ; j++)
		{
			result += M1.M[i][j]*M2.M[j][i];
		}
	} 
}


Matrix operator*(double op1, Matrix op2)
{

  Matrix temp;
  
  for(int i=0; i<4; i++){
    for(int j=0; j<4; j++){ 
      temp.M[i][j] = op1*op2.M[i][j]; 
    }
  }

  return temp;

}


Matrix operator*(complex<double> op1, Matrix op2)
{

  Matrix temp;
  
  for(int i=0; i<4; i++){
    for(int j=0; j<4; j++){ 
      temp.M[i][j] = op1*op2.M[i][j]; 
    }
  }

  return temp;

}


Matrix operator*(Matrix op1, Matrix op2)
{

  Matrix temp;
  
  for(int i=0; i<4; i++){
    for(int j=0; j<4; j++){ 
      for(int k=0; k<4; k++){
      temp.M[i][j] += op1.M[i][k]*op2.M[k][j];
      } 
    }
  }

  return temp;

}


Matrix operator+(Matrix op1, Matrix op2)
{

  Matrix temp;
  
  for(int i=0; i<4; i++){
    for(int j=0; j<4; j++){
      temp.M[i][j] = op1.M[i][j] + op2.M[i][j]; 
    }
  }

  return temp;

}


Matrix operator-(Matrix op1, Matrix op2)
{

  Matrix temp;
  
  for(int i=0; i<4; i++){
    for(int j=0; j<4; j++){
      temp.M[i][j] = op1.M[i][j] - op2.M[i][j]; 
    }
  }

  return temp;

}


ostream &operator<<(ostream &stream, Matrix m)
{

  for(int i=0; i<4; i++){
    stream << "|" << " " ;
    for(int j=0; j<4; j++){
      stream << m.M[i][j] << " " ;
    }
    stream << "|" << endl;
  }  

  return stream;

}


// /*
// The procedure below defines the identity matrix and the gamma-Dirac matrices in the contravariant notation ('upper indices').
// */
// 
// void Construct_GammaMatrices(Matrix &Id, Matrix Gamma[])
// {
// 
//   complex<double> id[4][4] = { {complex<double> (1, 0), complex<double> (0, 0), complex<double> (0, 0), complex<double> (0, 0)},
// 			     {complex<double> (0, 0), complex<double> (1, 0), complex<double> (0, 0), complex<double> (0, 0)},
// 			     {complex<double> (0, 0), complex<double> (0, 0), complex<double> (1, 0), complex<double> (0, 0)},
// 			     {complex<double> (0, 0), complex<double> (0, 0), complex<double> (0, 0), complex<double> (1, 0)} };
// 
//   complex<double> g0[4][4] = { {complex<double> (1, 0), complex<double> (0, 0), complex<double> (0, 0), complex<double> (0, 0)},
// 			     {complex<double> (0, 0), complex<double> (1, 0), complex<double> (0, 0), complex<double> (0, 0)},
// 			     {complex<double> (0, 0), complex<double> (0, 0), complex<double> (-1, 0), complex<double> (0, 0)},
// 			     {complex<double> (0, 0), complex<double> (0, 0), complex<double> (0, 0), complex<double> (-1, 0)} };
// 
//   complex<double> g1[4][4] = { {complex<double> (0, 0), complex<double> (0, 0), complex<double> (0, 0), complex<double> (1, 0)},
// 			     {complex<double> (0, 0), complex<double> (0, 0), complex<double> (1, 0), complex<double> (0, 0)},
// 			     {complex<double> (0, 0), complex<double> (-1, 0), complex<double> (0, 0), complex<double> (0, 0)},
// 			     {complex<double> (-1, 0), complex<double> (0, 0), complex<double> (0, 0), complex<double> (0, 0)} };
// 
//   complex<double> g2[4][4] = { {complex<double> (0, 0), complex<double> (0, 0), complex<double> (0, 0), complex<double> (0, -1)},
// 			     {complex<double> (0, 0), complex<double> (0, 0), complex<double> (0, 1), complex<double> (0, 0)},
// 			     {complex<double> (0, 0), complex<double> (0, 1), complex<double> (0, 0), complex<double> (0, 0)},
// 			     {complex<double> (0, -1), complex<double> (0, 0), complex<double> (0, 0), complex<double> (0, 0)} };
// 
//   complex<double> g3[4][4] = { {complex<double> (0, 0), complex<double> (0, 0), complex<double> (1, 0), complex<double> (0, 0)},
// 			     {complex<double> (0, 0), complex<double> (0, 0), complex<double> (0, 0), complex<double> (-1, 0)},
// 			     {complex<double> (-1, 0), complex<double> (0, 0), complex<double> (0, 0), complex<double> (0, 0)},
// 			     {complex<double> (0, 0), complex<double> (1, 0), complex<double> (0, 0), complex<double> (0, 0)} };
// 
//   Id = Matrix(id);
//   
//   Gamma[0] = Matrix(g0);
//   Gamma[1] = Matrix(g1);
//   Gamma[2] = Matrix(g2);
//   Gamma[3] = Matrix(g3);
// 
// }

// Matrix Dagger(double k[], Matrix Gamma[]){
// 
//   Matrix temp;
// 
//   temp = k[0]*Gamma[0] - k[1]*Gamma[1] - k[2]*Gamma[2] - k[3]*Gamma[3];
// 
//   return temp;
//   
// }
  
  
  
// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //   
// // // // // // // // // // // // // // //   Nuevo  // // // // // // // // // // // // //
// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //   
  
// class Spr{
//   complex<double>* N;
// public:
//   Spr(const Spr &n);
//   Spr(complex<double> B[4]);
//   Spr();
//   ~Spr();
//   friend Spr operator*(complex<double> op1, Spr op2);
//   friend Spr operator+(Spr op1, Spr op2);
//   friend ostream &operator<<(ostream &stream, Spr n);
// 
//   friend TSpr ConjSpr1(Spr n);
//   friend Spr ConjSpr2(TSpr n);
//   friend complex<double> operator*(TSpr op1, Spr op2);
// 
//   friend Spr operator*(Matrix op1, Spr op2);
//   friend Matrix operator*(Spr op1, TSpr op2);
// };


Spr::Spr(const Spr &n)
{


  for(int i=0; i<4; i++){
    N[i] = n.N[i]; 
  }  
  
  //cout << "Copy constructor Spr" << endl;

}


Spr::Spr()
{

  
  //cout << "Constructor w/o argts Spr" << endl;

}


Spr::Spr(complex<double> B[4])
{
  
  
  for(int i=0; i<4; i++){
    N[i] = B[i];
  }
  
  //cout << "Constructor with argts Spr" << endl;

}


Spr::~Spr()
{


  //cout << "Kaputt Spr!" << endl;

}


Spr operator*(complex<double> op1, Spr op2)
{

  Spr temp;

  for(int i=0; i<4; i++){
    temp.N[i] = op1*op2.N[i]; 
  }

  return temp;

}


Spr operator+(Spr op1, Spr op2)
{

  Spr temp;

  for(int i=0; i<4; i++){
    temp.N[i] = op1.N[i] + op2.N[i]; 
  }

  return temp;

}


ostream &operator<<(ostream &stream, Spr n)
{

  for(int i=0; i<4; i++){
    stream << "|" << " " << n.N[i] << " " << "|" << endl;
  }    

  return stream;

}


Spr operator*(Matrix op1, Spr op2)
{

  Spr temp;

  for(int i=0; i<4; i++){
    for(int k=0; k<4; k++){
      temp.N[i] += op1.M[i][k]*op2.N[k];
    }
  }

  return temp;

}


// class TSpr{
//   complex<double>* Nt;
// public:
//   TSpr(const TSpr &n);
//   TSpr(complex<double> B[4]);
//   TSpr();
//   ~TSpr();
//   friend TSpr operator*(complex<double> op1, TSpr op2);
//   friend TSpr operator+(TSpr op1, TSpr op2);
//   friend ostream &operator<<(ostream &stream, TSpr n);
// 
//   friend TSpr ConjSpr1(Spr n);
//   friend Spr ConjSpr2(TSpr n);
//   friend complex<double> operator*(TSpr op1, Spr op2);
// 
//   friend Matrix operator*(Spr op1, TSpr op2);
//   friend TSpr operator*(TSpr op1, Matrix op2);
// };


TSpr::TSpr(const TSpr &n)
{

  for(int i=0; i<4; i++){
    Nt[i] = n.Nt[i]; 
  }  
  
  //cout << "Copy constructor TSpr" << endl;

}


TSpr::TSpr()
{

  //cout << "Constructor w/o argts TSpr" << endl;

}


TSpr::TSpr(complex<double> B[4])
{
  
  
  for(int i=0; i<4; i++){
    Nt[i] = B[i];
  }
  
  //cout << "Constructor with argts TSpr" << endl;

}


TSpr::~TSpr()
{


  //cout << "Kaputt TSpr!" << endl;

}

// // Raul: I added this "assignment operator"
const TSpr &TSpr::operator=(const TSpr &right)
{

  if ( &right != this )
  {
    for(int i=0; i<4; i++){
      Nt[i] = right.Nt[i];
    }
  }
  
  return *this;
  
}
// // // // // // // // 


TSpr operator*(complex<double> op1, TSpr op2)
{

  TSpr temp;

  for(int i=0; i<4; i++){
    temp.Nt[i] = op1*op2.Nt[i]; 
  }

  return temp;

}


TSpr operator+(TSpr op1, TSpr op2) //Raul: this works
{

  TSpr temp;

  for(int i=0; i<4; i++){
    temp.Nt[i] = op1.Nt[i] + op2.Nt[i]; 
  }

  return temp;

}


ostream &operator<<(ostream &stream, TSpr n)
{

  stream << "|" << " ";  
  for(int i=0; i<4; i++){
    stream << n.Nt[i] << " ";
  }    
  stream << "|" << endl;

  return stream;

}


TSpr ConjSpr1(Spr n)
{

  TSpr temp;

  for(int i=0; i<4; i++){
    temp.Nt[i] = conj(n.N[i]);
  }

  return temp;

}


Spr ConjSpr2(TSpr n)
{
  
  Spr temp;

  for(int i=0; i<4; i++){
    temp.N[i] = conj(n.Nt[i]);
  }

  return temp;

}


complex<double> operator*(TSpr op1, Spr op2)
{

  complex<double> temp;

  for(int i=0; i<4; i++){
    temp += op1.Nt[i]*op2.N[i];
  }

  return temp;

}


Matrix operator*(Spr op1, TSpr op2)
{

  Matrix temp;

  for(int i=0; i<4; i++){
    for(int j=0; j<4; j++){
      temp.M[i][j] = op1.N[i]*op2.Nt[j];
    }
  }

  return temp;

}


TSpr operator*(TSpr op1, Matrix op2) //Raul: now... works
{

  TSpr temp;

  for(int i=0; i<4; i++){
    for(int k=0; k<4; k++){
      temp.Nt[i] += op1.Nt[k]*op2.M[k][i];
//   cout << "Hello" << endl;
    }
  }

  return temp;

}  
 




