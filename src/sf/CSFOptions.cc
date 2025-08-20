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


CSFOptions::CSFOptions(params &par, bool cc, bool is_proton, bool anty)
    : m_proton(is_proton),
      m_switchAntineut(anty)
{
    m_cc       = cc;
    m_qel_new  = 1;   // forced ; could be: par.qel_new

    const int method  = par.sf_method;
    const int mdkind  = method; // for now: force match

    if (method == 1) return; // For method 1 we do NOT need createSF(te,md,...) /
    
    TargetElement te;
    MomDistribs   md;

    switch (fun(par.nucleus_p, par.nucleus_n, method, mdkind)) {
        case  600622: te = targC_Ben;   md = md_C12_Ben;   break; // C, Benhar
        case  800833: te = targO_GSF;   md = md_O16_CdA;   break; // O, GSF
        case 2002033: te = targCa_GSF;  md = md_Ca40_CdA;  break; // Ca, GSF
        case 1802222: te = targAr_Ben;  md = md_Ca40_GCo;  break; // Ar, Benhar
        case 1802233: te = targAr_GSF;  md = md_Ca40_CdA;  break; // Ar, GSF
        default:
            std::cerr << "Illegal SF configuration: "
                      << fun(par.nucleus_p, par.nucleus_n, method, mdkind)
                      << std::endl;
            std::exit(24);
    }

    f = createSF(te, md, (is_proton ? proton : neutron));

    if (!m_qel_new)
        set_formFactors(par.qel_vector_ff_set);

    // m_switchToEM=p.sf_switchToEM;
    m_switchAntineut=anty;

}

//new constant name change : M2->Mass2
double CSFOptions::evalLH(const double q4til2,
			  const double p4k4,
			  const double p4kPrime4,
	                  const double p4q4til,
		          const double k4q4til,
			  const double kPrime4q4til,
                          const double k4kPrime4) const
{ // Calculating kinematic variables
  const double qM = q4til2 / Mass2;
  const double tau = -qM / 4.0;
    double a1 = 2.0 * Mass2 * k4kPrime4;
    double a2 = 2.0 * p4k4 * p4kPrime4 - Mass2 * k4kPrime4 + p4k4 * kPrime4q4til + p4kPrime4 * k4q4til - k4kPrime4 * p4q4til;
    double a3 = p4kPrime4 * k4q4til - p4k4 * kPrime4q4til;
    double a4 = k4kPrime4 * q4til2 - 2.0 * k4q4til * kPrime4q4til;

  double f1(0.0),f2(0.0),fa(0.0),fp(0.0),f1x(0.0),f2x(0.0),fax(0.0),fpx(0.0);

  if(m_qel_new!=1)
	{
	const double ge = FFGE(q4til2);
	const double gm = FFGM(q4til2);
	f1x = f1 = (ge + tau * gm) / (1 + tau) ;
	f2x = f2 = (gm - ge) / (1 + tau) ;
	fax = fa = FA(q4til2) ;
	fpx = fp = 2.0 * Mass2 * fa / (piMass2 - q4til2) ;
	std::cout<<ge<<gm<<std::endl;
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
}
double CSFOptions::evalLHnc(const double q4til2,
			    const double p4k4,
			    const double p4kPrime4,
			    const double p4q4til,
			    const double k4q4til,
			    const double kPrime4q4til,
			    const double k4kPrime4 ) const
{
	const double qM=q4til2/Mass2;
	const double tau=-qM/4.0;
    double a1= 2.0*Mass2*k4kPrime4;
    double a2= 2.0*p4k4*p4kPrime4-Mass2*k4kPrime4+p4k4*kPrime4q4til+p4kPrime4*k4q4til-k4kPrime4*p4q4til;
    double a3= p4kPrime4*k4q4til-p4k4*kPrime4q4til;
    double a4= k4kPrime4*q4til2-2.0*k4q4til*kPrime4q4til;

    int sign=(m_proton?1:-1);
    double f1(0.0),f2(0.0),fa(0.0),fp(0.0);
    double f1x(0.0),f2x(0.0),fax(0.0),fpx(0.0);
	if(m_qel_new!=1)
	{
		const double ge=0.5*sign*FFGE(q4til2)-2*sin2thetaW*(m_proton?FFGEp(q4til2):FFGEn(q4til2));/// - 0.5*GEs(q4til2);
		const double gm=0.5*sign*FFGM(q4til2)-2*sin2thetaW*(m_proton?FFGMp(q4til2):FFGMn(q4til2));/// - 0.5*GMs(q4til2);
		f1x=f1= (ge + tau*gm)/(1 + tau);
		f2x=f2= (gm - ge)/(1 + tau);
		fax=fa= sign*0.5*FA(q4til2);//- 0.5*GAs(q4til2);
		fpx=fp= 2.0*Mass2*fa/(piMass2 - q4til2);
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
	const double qM=q4til2/Mass2;
	const double tau=-qM/4.0;
  double a1= 2.0*Mass2*k4kPrime4;
  double a2= 2.0*p4k4*p4kPrime4-Mass2*k4kPrime4+p4k4*kPrime4q4til+p4kPrime4*k4q4til-k4kPrime4*p4q4til;
  double a3= p4kPrime4*k4q4til-p4k4*kPrime4q4til;
  double a4= k4kPrime4*q4til2-2.0*k4q4til*kPrime4q4til;
  double ge,gm,lh=0;
  double f1(0.0), f2(0.0);
  list(f1,f2)=f12(q4til2,(m_proton?10:11));
	double f11= f1*f1 ;
	double f22= f2*f2 ;
	double ff= f1 + f2 ;
	double h1= ff*ff*tau ;
	double h2= f11 + f22*tau ;
	double h4= 0.25*f22*(1.0-tau) + 0.5*f1*f2 ;
	      lh+= 2.0*(a1*h1 +  a2*h2 + a4*h4 );
	return lh;
}
