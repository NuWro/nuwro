#include "CSFOptions.h"
#include "form_factory_nnff.h"
#include "ff.h"

extern params p;


void CSFOptions::set_formFactors(int choice)
{	switch(choice)
	{
		case dipoleFF:
			FFGE = &DipoleGE;
			FFGM = &DipoleGM;
			FFGEp = &ProtonDipoleGE;
			FFGMp = &ProtonDipoleGM;
			FFGEn = &NeutronDipoleGE;
			FFGMn = &NeutronDipoleGM;
			break;
			
		case BBBA05FF:
			FFGE = &BBBA05_GE;
			FFGM = &BBBA05_GM;
			FFGEp = &ProtonBBBA05_GE;
			FFGMp = &ProtonBBBA05_GM;
			FFGEn = &NeutronBBBA05_GE;
			FFGMn = &NeutronBBBA05_GM;
			break;
			
		case BBA03FF:
			FFGE = &BBA03_GE;
			FFGM = &BBA03_GM;
			FFGEp = &ProtonBBA03_GE;
			FFGMp = &ProtonBBA03_GM;
			FFGEn = &NeutronBBA03_GE;
			FFGMn = &NeutronBBA03_GM;
			break;
			
		case JLabFF:
			FFGE = &JLabGE;
			FFGM = &JLabGM;
			FFGEp = &ProtonJLabGE;
			FFGMp = &ProtonJLabGM;
			FFGEn = &NeutronJLabGE;
			FFGMn = &NeutronJLabGM;
			break;

		default:
			std::cout<<"There are no such parameterization of the form factors"<<std::endl;
			return;

	};
}


const  inline int fun(int a,int b,int c, int d)
{return ((a*1000+b)*10+c)*10+d;
}

	
CSFOptions::CSFOptions(params &p, bool cc,bool is_proton,bool anty)
                      :m_proton(is_proton),
                       m_switchAntineut(anty)
{  
	m_cc=cc;
	m_qel_new=1;//p.qel_new;
    TargetElement te;
	MomDistribs md;
	
	int method=p.sf_method;      // 0 - none , 1 - grid, 2 - GSF
	//int mdkind=p.sf_mom_distrib; // 1 - Ben,  2 - CdA, 3 - CGo?? 
	int mdkind=method;          //force match


	switch(fun(p.nucleus_p,p.nucleus_n,method,mdkind))
	{ case  600611: te= targC_Ben;  md= md_C12_Ben;  break;
	  case  800811: te= targO_Ben;  md= md_O16_Ben;  break;	
	  case  800822: te= targO_GSF;  md= md_O16_CdA;  break;
	  case 2002022: te= targCa_GSF; md= md_Ca40_CdA; break;
	  case 1802211: te= targAr_Ben; md= md_Ca40_GCo; break;	
 	  case 1802222: te= targAr_GSF; md= md_Ca40_CdA; break;	
 	  case 2603011: te= targFe_Ben; md= md_Fe56_Ben; break;	
	  default: 
	    cerr <<"Illegal SF configuration: "<<fun(p.nucleus_p,p.nucleus_n,method,mdkind)<<endl;
		exit(24);
	}
	
	f=createSF(te,md,(is_proton? proton : neutron));
	if(!m_qel_new)
	  set_formFactors(p.qel_vector_ff_set);
	
//    m_switchToEM=p.sf_switchToEM;
	m_switchAntineut=anty;
		
}

/*
CSFOptions::CSFOptions(const TargetElement i_target,
		const IsospinOfSF i_isospin,
		const MomDistribs i_momDistrib,
		const int i_neutrino_pdg,
		const FormFactors i_formFactors,
//		const bool i_switchToEM,
		bool cc)
		 :
//		 m_switchToEM(i_switchToEM),
		 m_cc(cc)
{	m_switchAntineut=i_neutrino_pdg<0;
	set_formFactors(i_formFactors);
	m_proton=(i_isospin==proton);
	f=createSF(i_target,i_momDistrib,i_isospin);
}

*/


double CSFOptions::evalLH(const double q4til2, 
						   const double p4k4, 
						   const double p4kPrime4, 
						   const double p4q4til,
						   const double k4q4til,
						   const double kPrime4q4til,
						   const double k4kPrime4) const
{

	const double qM=q4til2/M2;
	const double tau=-qM/4.0;
	
    double a1= 2.0*M2*k4kPrime4;
    double a2= 2.0*p4k4*p4kPrime4-M2*k4kPrime4+p4k4*kPrime4q4til+p4kPrime4*k4q4til-k4kPrime4*p4q4til;
    double a3= p4kPrime4*k4q4til-p4k4*kPrime4q4til;
    double a4= k4kPrime4*q4til2-2.0*k4q4til*kPrime4q4til;

	double f1,f2,fa,fp;
	double f1x,f2x,fax,fpx;
	if(m_qel_new!=1)
	{
	const double ge=FFGE(q4til2);
	const double gm=FFGM(q4til2);
	f1x=f1= (ge + tau*gm)/(1 + tau) ;
	f2x=f2= (gm - ge)/(1 + tau) ;
	fax=fa= FA(q4til2) ;
	fpx=fp= 2.0*M2*fa/(piMass2 - q4til2) ;
	}

  if(m_qel_new)
  {
    list(f1,f2)=f12(q4til2,0);
    list(fa,fp)=fap(q4til2,0);
    if(m_qel_new==2)
    {
      cout<<f1<< '\t'<<f2<<'\t'<<fa<< '\t'<<fp<<'\t'<<endl;
       cout<<f1x<< '\t'<<f2x<<'\t'<<fax<< '\t'<<fpx<<'\t'<<endl<<endl;
    }
  }

	const double f11= f1*f1;
	const double f22= f2*f2;
	const double fa2= fa*fa;
	const double fp2= fp*fp;

	const double ff= f1 + f2;

	const double h1= fa2*(1.0+tau) + ff*ff*tau ;
	const double h2= fa2 + f11 + f22*tau ;
		    double h3= 2.0*fa*ff ;

	if (m_switchAntineut) 
		h3= -h3;
	const double h4= 0.25*f22*(1.0-tau) + 0.5*f1*f2 + fa*fp - fp2*tau;

	const double lh= 2.0*(a1*h1 + a2*h2 + a3*h3 + a4*h4) ;

	return  lh;

/*		const double a = M2*(m_leptMass2 - q4til2)/4.0*
 *                       (fa2*(4.0 - qM) - f11*(4.0 + qM) - qM*f22*(1.0 - tau) - 4.0*qM*f1*f2 
 *                        - m_leptMass2/M2*(ff*ff + fa2 + 4*fa*fp + qM*fp2));
		const double b = -q4til2*ff*fa;//M2*B(q2)
		const double c = 0.25*(f11 + tau*f22 + fa2);
	
		const double su = 4.0*p4k4 + q4til2 - m_leptMass2;		
		const int nuclNumb( m_N );
		return  nuclNumb*(a + su*( c*su -(!m_switchAntineut)*b +(m_switchAntineut)*b ));*/
}

double CSFOptions::evalLHnc(const double q4til2, 
						   const double p4k4, 
						   const double p4kPrime4, 
						   const double p4q4til,
						   const double k4q4til,
						   const double kPrime4q4til,
						   const double k4kPrime4
						   ) const
{

	const double qM=q4til2/M2;
	const double tau=-qM/4.0;
    double a1= 2.0*M2*k4kPrime4;
    double a2= 2.0*p4k4*p4kPrime4-M2*k4kPrime4+p4k4*kPrime4q4til+p4kPrime4*k4q4til-k4kPrime4*p4q4til;
    double a3= p4kPrime4*k4q4til-p4k4*kPrime4q4til;
    double a4= k4kPrime4*q4til2-2.0*k4q4til*kPrime4q4til;
	
	int sign=(m_proton?1:-1);
	//double ge,gm;
    //list(ge,gm)=Gem(q4til2,(proton?1:2));
    //cout<<ge<< '\t'<<gm<<endl;
    //cout<<gex<< '\t'<<gmx<<endl<<endl;
    double f1,f2,fa,fp;
    double f1x,f2x,fax,fpx;
	if(m_qel_new!=1)
	{
		const double ge=0.5*sign*FFGE(q4til2)-2*sin2thetaW*(m_proton?FFGEp(q4til2):FFGEn(q4til2));/// - 0.5*GEs(q4til2);
		const double gm=0.5*sign*FFGM(q4til2)-2*sin2thetaW*(m_proton?FFGMp(q4til2):FFGMn(q4til2));/// - 0.5*GMs(q4til2);
		f1x=f1= (ge + tau*gm)/(1 + tau);
		f2x=f2= (gm - ge)/(1 + tau);
		fax=fa= sign*0.5*FA(q4til2);//- 0.5*GAs(q4til2);
		fpx=fp= 2.0*M2*fa/(piMass2 - q4til2);
    }
	if(m_qel_new)
	{
    list(f1,f2)=f12(q4til2,(m_proton?1:2));
    list(fa,fp)=fap(q4til2,(m_proton?1:2));
		if(m_qel_new==2)
		{
			cout<<endl;
			cout<<f1<< '\t'<<f2<<'\t'<<fa<< '\t'<<fp<<'\t'<<endl;
			cout<<f1x<< '\t'<<f2x<<'\t'<<fax<< '\t'<<fpx<<'\t'<<endl<<endl;
		}
    }

	const double f11= f1*f1;
	const double f22= f2*f2;
	const double fa2= fa*fa;
	const double fp2= fp*fp;

	const double ff= f1 + f2;

	const double h1= fa2*(1.0+tau) + ff*ff*tau;
	const double h2= fa2 + f11 + f22*tau;
	double h3= 2.0*fa*ff ;
	if (m_switchAntineut)  
	   h3= -h3;
	
	const double h4= 0.25*f22*(1.0-tau) + 0.5*f1*f2 + fa*fp - fp2*tau;

	const double lh= 2.0*(a1*h1 + a2*h2 + a3*h3 + a4*h4);

	return  lh;
}

double CSFOptions::evalLHel(const double q4til2, 
						   const double p4k4, 
						   const double p4kPrime4, 
						   const double p4q4til,
						   const double k4q4til,
						   const double kPrime4q4til,
						   const double k4kPrime4) const
{
	const double qM=q4til2/M2;
	const double tau=-qM/4.0;

  double a1= 2.0*M2*k4kPrime4;
  double a2= 2.0*p4k4*p4kPrime4-M2*k4kPrime4+p4k4*kPrime4q4til+p4kPrime4*k4q4til-k4kPrime4*p4q4til;
  double a3= p4kPrime4*k4q4til-p4k4*kPrime4q4til;
  double a4= k4kPrime4*q4til2-2.0*k4q4til*kPrime4q4til;

  double ge,gm,lh=0;

  //for(int proton=0;proton<2;proton++)
  // {
  // if(proton)
  // {
  // ge= FFGEp(q4til2) ;
	// gm= FFGMp(q4til2) ;
  //  }
  //  else
  //  {
  // ge = FFGEn(q4til2);
 	//  gm = FFGMn(q4til2);
  //  }

	// double f1= (ge + tau*gm)/(1 + tau) ;
	// double f2= (gm - ge)/(1 + tau) ;
//    list(F1,F2)=f12(q2,cc,proton,anty);

  double f1, f2;
  list(f1,f2)=f12(q4til2,(m_proton?10:11));

	double f11= f1*f1 ;
	double f22= f2*f2 ;

	double ff= f1 + f2 ;

	double h1= ff*ff*tau ;
	double h2= f11 + f22*tau ;
	double h4= 0.25*f22*(1.0-tau) + 0.5*f1*f2 ;

	      lh+= 2.0*(a1*h1 +  a2*h2 + a4*h4 );
	// }
	return lh;
}
