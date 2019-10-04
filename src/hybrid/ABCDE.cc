#include <complex>
using namespace std;
#include "Constants.h"

void ABCDE(const std::complex<double> Hadron[4][4], const std::complex<double> Lepton[4][4], double *ABC_factors, double *SL_factors, double* e_responses, Lepton_kin* Lep, Lepton_kin* Lep_b)
{
//We calculate the factors A B C D and E from lepton and hadron tensor as described in eq. 31 from Sobscyk paper.
//We store them in the array ABC factors
	double Hadron_s[4][4];
	double Hadron_a[4][4]; //Symmetric and antisymmetric hadron tensor elements determined by eq. 15, not taking into account the complex unit. 
	for (int mu = 0 ; mu < 4 ; mu++)
	{
		for (int nu =0 ; nu < 4 ; nu++)
		{
			Hadron_s[mu][nu] = real(0.5*(Hadron[mu][nu]+Hadron[nu][mu]));
			Hadron_a[mu][nu] = real(-I*0.5*(Hadron[mu][nu]-Hadron[nu][mu]));
		}
	}


	//A:
	ABC_factors[0] = real(Lepton[0][0]*Hadron_s[0][0] + 2.*Lepton[3][0]*Hadron_s[3][0] + Lepton[3][3]*Hadron_s[3][3] + 0.5*(Lepton[1][1]+Lepton[2][2])*(Hadron_s[1][1]+Hadron_s[2][2]) + 2.*I*Lepton[1][2]*Hadron_a[1][2]);
//	cout << Lepton[0][0]*Hadron_s[0][0] << "  " << 2.*Lepton[0][3]*Hadron_s[0][3] << "  " << Lepton[3][3]*Hadron_s[3][3] << "  " << 0.5*(Lepton[1][1]+Lepton[2][2])*(Hadron_s[1][1] + Hadron_s[2][2]) << "  " << 2.*I*Lepton[1][2]*Hadron_a[1][2] << endl;


//	cout << "RT_symmetrix :  Sato Lee : " << pow(Lep->kl_inc[0]*Lep->kl[1]/Lep->q, 2) + 0.5*(Lep->Qsq + muonmass*muonmass) << "  tensor : " << 0.5*(Lepton[1][1]+Lepton[2][2]) << endl; 
//        cout << "RT_anti : Sato Lee : " << ((Lep->kl_inc[0]+Lep->kl[0])*Lep->Qsq - Lep->Q[0]*muonmass*muonmass)/Lep->q << "  " << x << "  " << " tensor   : " << 2.*I*Lepton[1][2] << endl;
//
//	cout << "L00, Sato : " << A*0.5/pow(qc,2) + wc*B/pow(qc,2) + 0.5*pow(wc,2)*C/pow(qc,2) << " Lepton = " << Lepton[0][0] << endl;
//
//	cout << "L30 Sato : " << -B/qc - wc*C/qc << " Lepton*2 " << Lepton[3][0]*2. <<  "  " << Lepton[0][3] << endl;
//
//	cout << "L33 Sato : " << 0.5*C << " Lepton : " << Lepton[3][3] << endl;
	
	//B:
	ABC_factors[1] = real( 2.*(Lepton[0][1]*Hadron_s[0][1]+Lepton[1][3]*Hadron_s[1][3]  + I*Lepton[0][2]*Hadron_a[0][2] + I*Lepton[2][3]*Hadron_a[2][3]));
//	ABC_factors[1] = real(2.*I*Lepton[1][2]*Hadron_a[1][2] ); //AS term

	//C:
	ABC_factors[2] = real(0.5*(Lepton[1][1]-Lepton[2][2])*(Hadron_s[1][1]-Hadron_s[2][2]));
//	ABC_factors[2] = real(Lepton[0][0]*Hadron_s[0][0] +(-B/qc - wc*C/qc)*Hadron_s[3][0] + Lepton[3][3]*Hadron_s[3][3] + 0.5*(Lepton[1][1]+Lepton[2][2])*(Hadron_s[1][1]+Hadron_s[2][2])); //Symm term

	//D:
	ABC_factors[3] = real(2.*(-Lepton[0][1]*Hadron_s[0][2] -Lepton[1][3]*Hadron_s[2][3] + I*Lepton[0][2]*Hadron_a[0][1]+ I*Lepton[2][3]*Hadron_a[1][3]));

	//E:
	ABC_factors[4] = real((Lepton[2][2] - Lepton[1][1])*Hadron_s[1][2]);

	//Terms for calculation of the inclusive cross section on a grid of Q and W
	SL_factors[0] = Hadron_s[0][0];
	SL_factors[1] = Hadron_s[3][0];
	SL_factors[2] = Hadron_s[3][3];
	SL_factors[3] = Hadron_s[1][1]+Hadron_s[2][2];
	SL_factors[4] = Hadron_a[1][2];

	//eq (D8) without any prefactors
	e_responses[0] = 0.5*real(Hadron[1][1]+Hadron[2][2]) ;
	e_responses[1] = real(Hadron[3][3]);
	e_responses[2] = 0.5*real(Hadron[1][1]-Hadron[2][2]);
	e_responses[3] = real(Hadron[1][3]); 
	e_responses[4] = imag(Hadron[1][3]);


}


