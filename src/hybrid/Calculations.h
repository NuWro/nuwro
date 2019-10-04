#ifndef CALCULATIONS_H
#define CALCULATIONS_H
#include "Constants.h"

double CS_interpolate(double QQ,double W,CS_tabulated* CS_table, Lepton_kin* Lep, Reaction_parameters* Reac);

double Do_Fivefold_Calc( Reaction_parameters *Reac, Lepton_kin* Lep, Hadron_prime * Had);

// void Leptonic_Tensor( const int &, const double &, const double &, const double &, const double &, const double &, complex<double> [][4], complex<double> [][4] );
// 
// Matrix Free_Nucleon_Projection_Operator( double [], Matrix [], Matrix );
// 
// // // // // // // // // // // // // // // // 
// // // // // DELTA // // // // // // // // // 
// 
// void Gamma_WNDelta( const int &, double, double, double [], double [], double [], Matrix [], Matrix, Matrix, Matrix [][4] );
// 
// void S_Delta( int, double, double [], Matrix [], Matrix, Matrix [][4], Matrix [][4] );	
// 
// void DP_current( int , int , int , int , double , double , double , double , double [], double [], double [], double [], Matrix , Matrix [], Matrix , Matrix [][4], Matrix [] );
// 
// void CDP_current( int , int , int , int , double , double , double , double , double [], double [], double [], double [], Matrix , Matrix [], Matrix , Matrix [][4], Matrix [] );
// // // // // // // D13 // // // // // // // // 
// 
// void Gamma_WND13( int, const int &, int, int, double, double, double [], double [], double [], Matrix [], Matrix, Matrix, Matrix [][4] );
// 
// void S_D13( double, double [], Matrix [], Matrix, Matrix [][4], Matrix [][4] );
// 
// void D13P_current( int , int , int , int , double , double , double , double , double [], double [], double [], double [], Matrix , Matrix [], Matrix , Matrix [][4], Matrix [] );
// 
// void CD13P_current( int , int , int , int , double , double , double , double , double [], double [], double [], double [], Matrix , Matrix [], Matrix , Matrix [][4], Matrix [] );

// // // // // // // S11 // // // // // // // // 
// 
// void Gamma_WNS11( int, int, double, double, double [], double [], Matrix [], Matrix, Matrix, Matrix [], Matrix [], Matrix [][4] );
// 
// void S_S11( double, double [], Matrix [], Matrix , Matrix & );
// 
// void Gamma_S11Npi( double, Matrix, Matrix, Matrix & );
// 
// // // // // // // // P11 resonance // // // // // // // // 
// 
// void Gamma_WNP11( int, int, double, double, double [], double [], Matrix [], Matrix, Matrix, Matrix [], Matrix [], Matrix [][4] );
// 
// void S_P11( double, double [], Matrix [], Matrix , Matrix & );
// 
// void Gamma_P11Npi( double, Matrix, Matrix, Matrix & );

// // // // // // // NP and CNP // // // // // // // // 
// 
// void Gamma_WNN( int, int, int, int, double, double, double [], Matrix [], Matrix, Matrix, Matrix [], Matrix [], Matrix [][4] );
// 
// void S_Nprop( double, double [], Matrix [], Matrix , Matrix & );
// 
// void Gamma_NNpi( double, double [], Matrix, Matrix [], Matrix, Matrix & );
// 
// void NP_current( int , int , int , int , double , double , double , double , double [], double [], Matrix , Matrix [], Matrix [], Matrix , Matrix [][4], Matrix [] );
// 
// void CNP_current( int , int , int , int , double , double , double , double , double [], double [], Matrix , Matrix [], Matrix [], Matrix , Matrix [][4], Matrix [] );
// 
// // // // // // CT PF PP  // // // // //
// 
// void formfactor(double , double , double [], double [], double , double &, double &, double &);
// 
// void CT_PF_PP( int, int, int, double , double , double [], double [], double , double , double , double , Matrix , Matrix [], Matrix [], Matrix , Matrix [][4], Matrix, Matrix [], Matrix [], Matrix [] );
// 
// // // // // // // // // // // // // 

#endif



