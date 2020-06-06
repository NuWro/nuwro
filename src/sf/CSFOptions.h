#ifndef C_SF_OPTIONS_H
#define C_SF_OPTIONS_H

#include <iostream>
#include "pdg.h"
#include "generatormt.h"
#include "GConstants.h"
#include "FormFactors.h"
#include "CSpectralFunc.h"
#include "params.h"

typedef double (*DFormFactorPtr)(const double);

class CSFOptions
{
public:
	CSFOptions(NSNWRO::params &p, bool cc,bool is_proton,bool anty);
	~CSFOptions()
	{//delete f;
	}
/*	
	CSFOptions(const TargetElement i_target,
		const IsospinOfSF i_isospin,
		const MomDistribs i_momDistrib,
		const int i_neutrino_pdg,
		const FormFactors i_formFactors,
//		const bool i_switchToEM,
		bool cc=1);
*/		
    void set_formFactors(int choice);


    CSpectralFunc *get_SF()const { return f;}


	double evalLH(const double q4til2, 
				  const double p4k4, 
				  const double p4kPrime4, 
				  const double p4q4til,
				  const double k4q4til,
				  const double kPrime4q4til,
				  const double k4kPrime4) const;

	double evalLHnc(const double q4til2, 
				  const double p4k4, 
				  const double p4kPrime4, 
				  const double p4q4til,
				  const double k4q4til,
				  const double kPrime4q4til,
				  const double k4kPrime4) 
				  const;


private:
    bool m_cc;
    bool m_proton;

	mutable DFormFactorPtr FFGE;
	mutable DFormFactorPtr FFGM;
	mutable DFormFactorPtr FFGEp;
	mutable DFormFactorPtr FFGMp;
	mutable DFormFactorPtr FFGEn;
	mutable DFormFactorPtr FFGMn;
	
	CSpectralFunc *f;
	
	bool m_switchAntineut;
	int m_qel_new;
	
};

#endif
