#include "crossection.h"
//channel you want to do
int chan=1;
int PV=1;
int FP=0;
//how many integration subregions
int no1=1;
int ffset=4;
D4V<DM> proton_em_v_mu_(0,0,0,0);
D4V<DM> neutron_em_v_mu_(0,0,0,0);
const D4V<double> p_(M,0,0,0);//initial nucleon at rest
double Q2_=0;
double q2_=0;
double q0_=0;
double q_=0;
double W_=0;
double W2_=0;
double F1V=0;
double DW=0;
double qp=0;
DM qslash=1;
double qpd=0;
DM pdslash=1;
DM pqslash=1;
DM gfivpqslash=1;
D4V<double> q4_(0,0,0,0);
D4V <double> pd_(0,0,0,0);
D4T<double> Lept_(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
double F1P=0;
double F2P=0;
double F1N=0;
double F2N=0;

bool selfenergy=true;

double dwidth_(double s)
{
	double result=0;
	if(s>((M+mpi)*(M+mpi)))
	{
		double lamb=sqrt(s*s+mpi2*mpi2+M2*M2-2*(s*mpi2+mpi2*M2+s*M2));
		double W=sqrt(s);
		double qcm=0.5*lamb/W;
		double rho=qcm*qcm*qcm/(0.57*0.57*GeV2+qcm*qcm)/W;
		result=0.13*GeV*rho/rhodelpin0;
	}
	//if ((result<0)||(result!=result)) std::cerr<<"Warning in dwidth_!\n";
	return result;
}
//Oset model coefficients, linear interpolation
double Cgammax(unsigned int i, double omega)
{
	int j=int(omega/100/MeV);
	double wynik=0;
	if((j<2)&&((i==0)||(i==2)))
	{
		double x1;
		if(i==0) x1=160.0;
		else x1=115.0;
		double y2=dataset[2][i+1];
		wynik=(omega-x1)*y2/(200.0-x1);
		if (wynik<0) wynik =0;
	}
	else if(j<5)
	{
		double x1=j*100.0;
		double y1=dataset[j][i+1];
		double y2=dataset[j+1][i+1];
		wynik=y1+(omega-x1)*(y2-y1)/100.0;
		if (wynik<0) wynik=0;
	}
	else if(j>4)
	{
		double x1=400.0;
		double y1=dataset[4][i+1];
		double y2=dataset[5][i+1];
		wynik=y1+(omega-x1)*(y2-y1)/100.0;
		if (wynik<0) wynik=0;
	}
	return wynik;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
///////                                 CHANNELS:                          ///////////////////////////////
///////     1: e + p -> e + p + \pi^\0                                     ///////////////////////////////
///////     2: e + n -> e + p + \pi^-                                      ///////////////////////////////
///////     3: e + p -> e + n + \pi^+                                      ///////////////////////////////
///////     4: e + n -> e + n + \pi^0                                      ///////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////

double deltaffv[3]={0,0,0};


double dsigma_domega_dEprime_depi2_(double Epi)
{
	double result=0;
	double kpi=sqrt(Epi*Epi-mpi2);
	double cospi=0.5*(2*Epi*(M+q0_)+M2-W_*W_-mpi2)/kpi/q_;
	//Lept_=D4T<double>(q4_,q4_);
	if(fabs(cospi)<=1)
	{
		
		double sinpi=sqrt(-cospi*cospi+1.0);
		D4V<double> k(Epi,kpi*sinpi,0,kpi*cospi);
		D4V<double> pp=p_+q4_;
		pp-=k;
		DM pprimeslash=pp*dirac;
		pprimeslash+=M;
		
		D4V<DM> smu(0,0,0,0);
		//4-vector describing the vertex choice:
		
		
		//Delta Pole
		std::complex<double> denom(0,0);
		D4T<DM> tmp=gmunu;
		tmp.Add((pd_/(-1.5*Md2)),pd_);
		tmp.Add((dirac/-3),dirac);
		tmp.Add((pd_/(3*Md)),dirac);
		tmp.Add((dirac/(-3*Md)),pd_);	
		tmp&=pdslash;
		smu=k*tmp;
		denom=W2_-Md2+delusive*Md*DW;
		denom=1.0/denom;
		
		tmp=gmunu;
		tmp*=(qslash*deltaffv[0] + qpd*deltaffv[1] + qp*deltaffv[2]);
		tmp.Add((q4_*(-deltaffv[0])),dirac);
		tmp.Add((q4_*(-deltaffv[1])),pd_);
		tmp.Add((q4_*(-deltaffv[2])),p_);
		tmp*=gamma5;
		
		
		smu=smu*tmp;
		smu*=denom*fstar/mpi*CGDP[chan];
		D4V<DM> tmpv(0,0,0,0);
		D4V<double> pcr=p_;
		pcr-=k;
		D4V<double> pcrq=pcr;
		pcrq+=q4_;
		double W2cr=pcr*pcr;
		
		//crossed Delta pole
		
		denom=W2cr-Md2;
		denom=1.0/denom;
		DM pcrs=((dirac*pcr)+Md);
		pcrs*=-1;
		tmp=gmunu;
		tmp.Add((pcr/(-1.5*Md2)),pcr);
		tmp.Add((dirac/(-3)),dirac);
		tmp.Add((pcr/(3*Md)),dirac);
		tmp.Add((dirac/(-3*Md)),pcr);
		tmp&=pcrs;
		tmpv=tmp*k;
		double qpcrq=-1.0;
		qpcrq*=q4_*pcrq;
		
		double qpcr=q4_*pcr;
		qpcr*=-1;
		tmp=gmunu;
		tmp*=(qslash*(-deltaffv[0]) + qpcr*deltaffv[1] + qpcrq*deltaffv[2]);
		tmp.Add(dirac,(q4_*(deltaffv[0])));
		tmp.Add(pcr,(q4_*(deltaffv[1])));
		tmp.Add(pcrq,(q4_*(deltaffv[2])));
		tmp*=gamma5;
		
		tmpv=tmp*tmpv;
		tmpv*=denom*fstar/mpi*CGDPC[chan];
		smu+=tmpv;
			
		DM kslash=k*dirac;
		
		double t=(q4_-k)*(q4_-k);
		
		//optional pion form factor a'la HNV PRD76
		double F_pion=1.0;
		if(FP)
		{
			F_pion=(lambdapi2-mpi2)/(lambdapi2-t);
		}
		//nucleon pole
		if(chan%2==1) tmpv=proton_em_v_mu_;
		else tmpv=neutron_em_v_mu_;
		switch(PV)
		{
			case 0:
			{
				tmpv&=gfivpqslash;
				tmpv*=(M*gA/fpi*sqrt2*CGNP[chan]*F_pion/(W2_-M2) );
				break;
			}
			case 1:
			{
				tmpv&=(kslash*gfivpqslash);
				tmpv*=(-gA/sqrt2*CGNP[chan]*F_pion/( fpi * (W2_-M2) ));
				break;
			}
			default:
			{
				tmpv&=(kslash*gfivpqslash);
				tmpv*=(-gA/sqrt2*CGNP[chan]*F_pion/( fpi * (W2_-M2) ) );
				break;
			}
			
		}
		
		smu+=tmpv;
		
		//nucleon pole crossed
		
		if(chan<3) tmpv=proton_em_v_mu_;
		else tmpv=neutron_em_v_mu_;
		
		DM pcrslash=pcr*dirac;
		pcrslash+=M;
		
		switch(PV)
		{
			case 0:
			{
				pcrslash*=gamma5;
				tmpv*=pcrslash;
				tmpv*=(M*gA/fpi*sqrt2*CGNPC[chan]*F_pion/(W2cr-M2) );
				break;
			}
			case 1:
			{
				pcrslash*=(kslash*gamma5);
				tmpv*=pcrslash;
				tmpv*=(-gA/sqrt2*CGNPC[chan]*F_pion/( fpi * (W2cr-M2) ) );
				break;
			}
			default:
			{
				pcrslash*=(kslash*gamma5);
				tmpv*=pcrslash;
				tmpv*=(-gA/sqrt2*CGNPC[chan]*F_pion /( fpi * (W2cr-M2)) );
				break;
			}
		}
	
		smu+=tmpv;
		//contact term exists only in 2 and 3 (charged pi channels)
		if((chan!=1)&&(chan!=4))
		{
			
			if(PV)
			{
				tmpv= dirac;
				tmpv*=((-1/sqrt2*CGOTH[chan]*F_pion/fpi)*(gamma5*gA*F1V));
				smu+=tmpv;
			}
			//pion-in-flight
			D4V<double> k2q=k;
			k2q*=2;
			k2q-=q4_;
			tmpv= k2q*gamma5;
			
			tmpv*=(-(1/sqrt2*F1V*2*M*gA/fpi*CGOTH[chan]*F_pion)/( (t-mpi2) ) );
			
			smu+=tmpv;
		}
		D4V<DM> snu=smu;
	
	//////Here you do (\slash{p'}+M)s^\mu(\slash{p}+M)
		smu &=pprimeslash;	
		const DM pslash=((p_*dirac)+M);
		smu*=pslash;
		
	///////here you do \gamma^0(s^{\nu})^\dag\gamma^0
		snu.hermit();
		snu &=gamma0;
		snu *=gamma0;
		
		//you construct hadronic tensor
		D4T<DM> Amunu_(smu,snu);
				
		//contraction of the two tensors: little trick to avoid \int d\phi_\pi
		DM contr=(Amunu_(0,0)*Lept_(0,0)+Amunu_(3,3)*Lept_(3,3)-Amunu_(3,0)*Lept_(3,0)-Amunu_(0,3)*Lept_(0,3)+
				0.5*(Amunu_(1,1)+Amunu_(2,2))*(Lept_(1,1)+Lept_(2,2))+Amunu_(1,2)*Lept_(1,2)+Amunu_(2,1)*Lept_(2,1));
		result+=real(contr.CTrace());
		result*=2*Pi;
		
	}
	return result;
}

double dsigma_domega_dEprime_( double E, double Eprime, double coslep, double *CV)
{
	double result=0;
	//basic kinematics
	
	//energy transfer
	q0_=E-Eprime;
	//momentum transfer^2
	q2_=E*E+Eprime*Eprime-2*E*Eprime*coslep;
	//Q2
	Q2_=q2_-q0_*q0_;
	
	W2_=M2+2*M*q0_-Q2_;
	W_=sqrt(W2_);
	q_=sqrt(q2_);
	//std::cerr<<"Q2= "<<Q2_/GeV2<<" GeV^2 x= "<<x<<" q0= "<<q0_/GeV<<" GeV W= "<<W_/GeV<<" GeV q_= "<<q_/GeV<<" GeV\n";
	if((W2_>w2min)&&(q2_>0))
	{
		//Delta ff set uncomment to set from CV
		/*deltaffv[0]=CV[0]/M;
		deltaffv[1]=CV[1]/M2;
		deltaffv[2]=CV[2]/M2;
		*/
		//comment this out if you want to use CV
		switch(ffset)
		{
				case 0://model 0A
					deltaffv[0]=2.03/(1.0 + 3.69*Q2_/GeV2 + 1.82*Q2_/GeV2*Q2_/GeV2 + 1.13*Q2_/GeV2*Q2_/GeV2*Q2_/GeV2)/M;
					deltaffv[1]=-Mp*deltaffv[0]/(M*W_);
					deltaffv[2]=0;
					break;
				case 1://model 0B
					deltaffv[0]=2.60*(1+0.89*Q2_/GeV2)/(1.0 + 8.08*Q2_/GeV2 + 0.72*Q2_/GeV2*Q2_/GeV2 + 9.03*Q2_/GeV2*Q2_/GeV2*Q2_/GeV2)/M;
					deltaffv[1]=-2.60*Mp*(1+1.45*Q2_/GeV2)/(1.0 + 8.08*Q2_/GeV2 + 0.72*Q2_/GeV2*Q2_/GeV2 + 9.03*Q2_/GeV2*Q2_/GeV2*Q2_/GeV2)/(M2*W_);
					deltaffv[2]=0;
					break;
				case 2://model II
					deltaffv[0]=2.10*(1+0.13*Q2_/GeV2)/(1.0 + 4.73*Q2_/GeV2 -0.39*Q2_/GeV2*Q2_/GeV2 + 5.59*Q2_/GeV2*Q2_/GeV2*Q2_/GeV2)/M;
					deltaffv[1]=-2.10*Mp*(1+1.68*Q2_/GeV2)/(1.0 + 4.73*Q2_/GeV2 -0.39*Q2_/GeV2*Q2_/GeV2 + 5.59*Q2_/GeV2*Q2_/GeV2*Q2_/GeV2)/(M2*W_);
					deltaffv[2]=0.62/(pow(1.0+Q2_/Mv2,2)*M2);
					break;
				case 3://model I (Our fit of Olga)
					deltaffv[0]=2.00/(1.0+0.68*Q2_/Mv2)/pow(1.0+1.15*Q2_/Mv2,2)/M;
					deltaffv[1]=-6.77/2.00*deltaffv[0]/M;
					deltaffv[2]=5.95/((pow(1.0+1.15*Q2_/Mv2,2)*M2)*(1.0+1.40*Q2_/Mv2));
					break;
				default://Olga
					deltaffv[0]=2.13/(1.0+0.25*Q2_/Mv2)/pow(1.0+Q2_/Mv2,2)/M;
					deltaffv[1]=-1.51/2.13*deltaffv[0]/M;
					deltaffv[2]=0.48/((pow(1.0+Q2_/Mv2,2)*M2)*(1.0+1.29*Q2_/Mv2));
					break;
		}
		//set global 4-mometum transfer	
		q4_=D4V<double>(q0_,0,0,q_);
		//set global Delta 4-momentum and some invariants
		pd_=p_+q4_;
		qp=q4_*p_;
		qslash=dirac*q4_;
		qpd=q4_*pd_;
		pdslash=dirac*pd_;
		pdslash+=Md;
		pdslash*=-1;
		pqslash=dirac*pd_;
		pqslash+=M;
		gfivpqslash=gamma5*pqslash;
		//set global Delta width (we know Q2, q0 and W already...)
		DW=0;
		if(W2_>w2min)
		{
			double lamb=sqrt(W2_*W2_+M4+mpi2*mpi2-2*(W2_*mpi2+W2_*M2+mpi2*M2));
			
			double qcm=0.5*lamb/W_;
			
			double qcm2=qcm*qcm;
			double qcm3=qcm2*qcm;
			double nominator=fstar*fstar*qcm3*(sqrt(M2+qcm2)+M);
			
			double denominator = 12*Pi*mpi2*W_;
			DW=nominator/denominator;			
		}
		
		//set global nucleon ff
		double tau=Q2_/Mn24;
		double dipole=1.0/(1+Q2_/Mv2)/(1+Q2_/Mv2);
		F1P=(1+mup*tau)/(1+tau)*dipole;
		F2P=(mup-1.0)/(1+tau)*dipole;
		F1N=lambdan*mun*tau*tau/((1+tau)*(1+lambdan*tau))*dipole;
		F2N=mun*(tau+1+lambdan*tau)/((1+tau)*(1+lambdan*tau))*dipole;
		F1V=F1P-F1N;
		//set global nucleon em currents (we know Q2 and q0 already...)
		proton_em_v_mu_=q4_*sigma;
		neutron_em_v_mu_=proton_em_v_mu_;
		proton_em_v_mu_*=(-0.5*delusive*F2P/M);
		neutron_em_v_mu_*=(-0.5*delusive*F2N/M);
		proton_em_v_mu_+=(dirac*F1P);
		neutron_em_v_mu_+=(dirac*F1N);
		//leptonic tensor construction (global)
		double cosl1=0.5*(E*E+q2_-Eprime*Eprime)/(E*sqrt(q2_));
		double cosl2=-0.5*(q2_+Eprime*Eprime-E*E)/(Eprime*sqrt(q2_));
		if((fabs(cosl1)>1)||(fabs(cosl2)>1)) std::cerr<<"Warning in lepton kinematics!\n";
		//incoming electron variables
		double kl11=sqrt(1.0-cosl1*cosl1)*E;
		double kl13=cosl1*E;
		//outgoing electron variables
		double kl21=sqrt(1.0-cosl2*cosl2)*Eprime;
		double kl23=cosl2*Eprime;
		//assumed scattering in x-z plane
		D4V<double> l1(E,kl11,0,kl13);
		D4V<double> l2(Eprime,kl21,0,kl23);
		//leptonic tensor construction
		double l1l2=l1*l2;
		//first comes g^{\mu\nu} of the type double:
		Lept_=gmunud;
		//now it's -g^{\mu\nu} l1_\mu l2^\nu
		Lept_ *= -l1l2;
		//you add l1^\mu l2^\nu part
		Lept_.Add(l1,l2);
		//you add l2^\nu l2^\mu part:
		Lept_.Add(l2,l1);
		//electron leptonic tensor is ready
		
		//integration limits in pion energy
		double del=sqrt(W2_*W2_+M4+mpi2*mpi2-2*(W2_*M2+W2_*mpi2+M2*mpi2));
		double Emin=0.5*((q0_+M)*(W2_+mpi2-M2)-q_*del)/W2_;
		double Emax=0.5*((q0_+M)*(W2_+mpi2-M2)+q_*del)/W2_;
		
		if(Emax>Emin)
		{
			result=calg5(dsigma_domega_dEprime_depi2_,Emin,Emax,no1);
		}
		result*=(Eprime)/(E*sqrt(q2_));
		result*=alpha*alpha/(M*Q2_*Q2_*Pi*Pi*Pi*64*cm2); //so you get cm^2/GeV Sr (MeV=0.001) or cm^2/MeV Sr (MeV=1)
	}
	return result;
}

double dsigma_dQ2_dW_( double E, double Q2, double W, double *CV2)
{
	double result=0;
	//basic kinematics
	
	W_=W;
	W2_=W*W;
	//Q2
	Q2_=Q2;
	
	
	//energy transfer
	q0_=0.5*(W2_+Q2_-M2)/M;
	double Eprime=E-q0_;
	double coslep=0.5*(E*E+Eprime*Eprime-Q2_-q0_*q0_)/(E*Eprime);
	//momentum transfer^2
	q2_=E*E+Eprime*Eprime-2*E*Eprime*coslep;
	
	
	q_=sqrt(q2_);
	//std::cerr<<"Q2= "<<Q2_/GeV2<<" GeV^2 x= "<<x<<" q0= "<<q0_/GeV<<" GeV W= "<<W_/GeV<<" GeV q_= "<<q_/GeV<<" GeV\n";
	if((W2_>w2min)&&(q2_>0))
	{
		//Delta ff set uncomment to set from CV
		/*deltaffv[0]=CV[0]/M;
		deltaffv[1]=CV[1]/M2;
		deltaffv[2]=CV[2]/M2;
		*/
		//comment this out if you want to use CV
		switch(ffset)
				{
						case 0://model 0A
							deltaffv[0]=2.03/(1.0 + 3.69*Q2_/GeV2 + 1.82*Q2_/GeV2*Q2_/GeV2 + 1.13*Q2_/GeV2*Q2_/GeV2*Q2_/GeV2)/M;
							deltaffv[1]=-Mp*deltaffv[0]/(M*W_);
							deltaffv[2]=0;
							break;
						case 1://model 0B
							deltaffv[0]=2.60*(1+0.89*Q2_/GeV2)/(1.0 + 8.08*Q2_/GeV2 + 0.72*Q2_/GeV2*Q2_/GeV2 + 9.03*Q2_/GeV2*Q2_/GeV2*Q2_/GeV2)/M;
							deltaffv[1]=-2.60*Mp*(1+1.45*Q2_/GeV2)/(1.0 + 8.08*Q2_/GeV2 + 0.72*Q2_/GeV2*Q2_/GeV2 + 9.03*Q2_/GeV2*Q2_/GeV2*Q2_/GeV2)/(M2*W_);
							deltaffv[2]=0;
							break;
						case 2://model II
							deltaffv[0]=2.10*(1+0.13*Q2_/GeV2)/(1.0 + 4.73*Q2_/GeV2 -0.39*Q2_/GeV2*Q2_/GeV2 + 5.59*Q2_/GeV2*Q2_/GeV2*Q2_/GeV2)/M;
							deltaffv[1]=-2.10*Mp*(1+1.68*Q2_/GeV2)/(1.0 + 4.73*Q2_/GeV2 -0.39*Q2_/GeV2*Q2_/GeV2 + 5.59*Q2_/GeV2*Q2_/GeV2*Q2_/GeV2)/(M2*W_);
							deltaffv[2]=0.62/(pow(1.0+Q2_/Mv2,2)*M2);
							break;
						case 3://model I (Our fit of Olga)
							deltaffv[0]=2.00/(1.0+0.68*Q2_/Mv2)/pow(1.0+1.15*Q2_/Mv2,2)/M;
							deltaffv[1]=-6.77/2.00*deltaffv[0]/M;
							deltaffv[2]=5.95/((pow(1.0+1.15*Q2_/Mv2,2)*M2)*(1.0+1.40*Q2_/Mv2));
							break;
						default://Olga
							deltaffv[0]=2.13/(1.0+0.25*Q2_/Mv2)/pow(1.0+Q2_/Mv2,2)/M;
							deltaffv[1]=-1.51/2.13*deltaffv[0]/M;
							deltaffv[2]=0.48/((pow(1.0+Q2_/Mv2,2)*M2)*(1.0+1.29*Q2_/Mv2));
							break;
				}
		//set global 4-mometum transfer	
		q4_(0)=q0_;
		q4_(1)=0;
		q4_(2)=0;
		
		q4_(3)=q_;
		//D4V<double> p_(M,0,0,0);
		//set global Delta 4-momentum and some invariants
		pd_=p_+q4_;
		qp=q4_*p_;
		qslash=dirac*q4_;
		qpd=q4_*pd_;
		pdslash=dirac*pd_;
		pdslash+=Md;
		pdslash*=-1;
		pqslash=dirac*pd_;
		pqslash+=M;
		gfivpqslash=gamma5*pqslash;
		//set global Delta width (we know Q2, q0 and W already...)
		DW=0;
		if(W2_>w2min)
		{
			double lamb=sqrt(W2_*W2_+M4+mpi2*mpi2-2*(W2_*mpi2+W2_*M2+mpi2*M2));
			
			double qcm=0.5*lamb/W_;
			
			double qcm2=qcm*qcm;
			double qcm3=qcm2*qcm;
			double nominator=fstar*fstar*qcm3*(sqrt(M2+qcm2)+M);
			
			double denominator = 12*Pi*mpi2*W_;
			DW=nominator/denominator;			
		}
		
		//set global nucleon ff
		double tau=Q2_/Mn24;
		double dipole=1.0/(1+Q2_/Mv2)/(1+Q2_/Mv2);
		F1P=(1+mup*tau)/(1+tau)*dipole;
		F2P=(mup-1.0)/(1+tau)*dipole;
		F1N=lambdan*mun*tau*tau/((1+tau)*(1+lambdan*tau))*dipole;
		F2N=mun*(tau+1+lambdan*tau)/((1+tau)*(1+lambdan*tau))*dipole;
		F1V=F1P-F1N;
		//set global nucleon em currents (we know Q2 and q0 already...)
		proton_em_v_mu_=q4_*sigma;
		neutron_em_v_mu_=proton_em_v_mu_;
		proton_em_v_mu_*=(-0.5*delusive*F2P/M);
		neutron_em_v_mu_*=(-0.5*delusive*F2N/M);
		proton_em_v_mu_+=(dirac*F1P);
		neutron_em_v_mu_+=(dirac*F1N);
		//leptonic tensor construction (global)
		double cosl1=0.5*(E*E+q2_-Eprime*Eprime)/(E*sqrt(q2_));
		double cosl2=-0.5*(q2_+Eprime*Eprime-E*E)/(Eprime*sqrt(q2_));
		if((fabs(cosl1)>1)||(fabs(cosl2)>1)) std::cerr<<"Warning in lepton kinematics!\n";
		//incoming electron variables
		double kl11=sqrt(1.0-cosl1*cosl1)*E;
		double kl13=cosl1*E;
		//outgoing electron variables
		double kl21=sqrt(1.0-cosl2*cosl2)*Eprime;
		double kl23=cosl2*Eprime;
		//assumed scattering in x-z plane
		D4V<double> l1(E,kl11,0,kl13);
		D4V<double> l2(Eprime,kl21,0,kl23);
		//leptonic tensor construction
		double l1l2=l1*l2;
		//first comes g^{\mu\nu} of the type double:
		Lept_=gmunud;
		//now it's -g^{\mu\nu} l1_\mu l2^\nu
		Lept_ *= -l1l2;
		//you add l1^\mu l2^\nu part
		Lept_.Add(l1,l2);
		//you add l2^\nu l2^\mu part:
		Lept_.Add(l2,l1);
		//electron leptonic tensor is ready
		
		//integration limits in pion energy
		double del=sqrt(W2_*W2_+M4+mpi2*mpi2-2*(W2_*M2+W2_*mpi2+M2*mpi2));
		double Emin=0.5*((q0_+M)*(W2_+mpi2-M2)-q_*del)/W2_;
		double Emax=0.5*((q0_+M)*(W2_+mpi2-M2)+q_*del)/W2_;
		
		if(Emax>Emin)
		{
			result=calg5(dsigma_domega_dEprime_depi2_,Emin,Emax,no1);
		}
		result*=(Eprime)/(E*sqrt(q2_));
		result*=alpha*alpha/(M*Q2_*Q2_*Pi*Pi*Pi*64*cm2)*Pi*W/(E*Eprime*M); //so you get cm^2/GeV^3 (MeV=0.001) or cm^2/MeV^3 (MeV=1)
	}
	
	return result;
}

//laboratory electron energy, outgoing electron 4-momentum, hadronic CMS pion solid angle w.r.t hadronic cms delta direction
//Done in laboratory frame: electron coming along the z-axis, any outgoing electron angle and energy any direction of nucleon
double dsigma_dq0_dOmegal_dOmegaCMS(double E, double mass, D4V<double> lprimemu, double costhetacms, double phicms, D4V<double> pn)
{
	double result=0;
	
	const double mass2=mass*mass;
	
	//laboratory frame electron 4-momentum
	D4V <double> lmu(E,0,0,sqrt(E*E-mass2));
	//hadronic invariant mass squared
	W2_=((lmu-lprimemu+pn)*(lmu-lprimemu+pn));
	//first kinematic cut:can you produce a pion?
	if(W2_>(M+mpi)*(M+mpi))
	{
		W_=sqrt(W2_);
		//sum of electron and nucleon 4-momenta
		D4V <double> lpn=lmu+pn;
		//lepton-nucleon invariant mass squared
		double s=lpn*lpn;
		double sqrts=sqrt(s);
		
		//lepton-nucleon CMS "delta" momentum, massive electron assumption
		double pdelta=(s*s+mass2*mass2+W2_*W2_-2*s*W2_-2*mass2*s-2*mass2*W2_)/4/s;
		//second kinematic cut
		if(pdelta>0)
		{
			//lepton-nucleon CMS lepton energy
			double Ecms=(lpn*lmu+mass2)/sqrts;
			double lcms=sqrt(Ecms*Ecms-mass2);
			//outgoing lepton energy
			double Epcms=sqrt(pdelta+mass2);
			pdelta=sqrt(pdelta);
			// 4-momentum transfer squared
			Q2_=-((lmu-lprimemu)*(lmu-lprimemu));
			//neutrino laboratory energy, Q2, W, hadronic CMS pion solid angle
			//lepton scattering angle cosine in the lepton-nucleon CMS, massless electron
			double costhetal=(-Q2_+2*Ecms*Epcms-2*mass2)/2/pdelta/lcms;
			//third kinematic cut, just in case
			if(fabs(costhetal)<=1)
			{
				//since we work in CMS we have incoming electron and incoming nucleon (colinear)
				//and outgoing electron, outgoing hadronic system with invariant mass W and momentum=pdelta
				//and we want the z axis to go along outgoing hadronic system momentum and x-z lepton scattering plane always
				//and now with one boost one goes to hadronic CMS, where we define pionic angles
				double p0del=sqrt(W2_+pdelta*pdelta);
				//Incoming electron 4-momentum
				double sinthetal=sqrt(1.0-costhetal*costhetal);
				D4V <double> _lcms(((Ecms*p0del+lcms*pdelta*costhetal)/W_),lcms*sinthetal,0,-lcms*costhetal-pdelta*Ecms/W_-lcms*pdelta*pdelta*costhetal/W_/(W_+p0del));
				//Outgoing electron 4-momentum
				D4V <double> _lprimecms((Epcms*p0del+pdelta*pdelta)/W_,0,0,-pdelta-pdelta*Epcms/W_-pdelta*pdelta*pdelta/W_/(p0del+W_));
				//hadronic CMS 4-momentum transfer
				q4_=_lcms-_lprimecms;
				//result=-(q4_*q4_)/Q2_/GeV/4/Pi/1e33;
				//first comes g^{\mu\nu} of the type double:
				Lept_=gmunud;
				//now it's -g^{\mu\nu} Q^2/2
				Lept_ *= -Q2_/2;
				//you add l1^\mu l2^\nu part
				Lept_.Add(_lcms,_lprimecms);
				//you add l2^\nu l2^\mu part:
				Lept_.Add(_lprimecms,_lcms);
				//electron leptonic tensor is ready
				
				//now the rest of these bastards
				double lamb=sqrt(W2_*W2_+M4+mpi2*mpi2-2*(W2_*mpi2+W2_*M2+mpi2*M2));
			
				double qcm=0.5*lamb/W_;
				
				double qcm2=qcm*qcm;
				double qcm3=qcm2*qcm;
				double nominator=fstar*fstar*qcm3*(sqrt(M2+qcm2)+M);
				
				double denominator = 12*Pi*mpi2*W_;
				DW=nominator/denominator;
				
				//hadronic CMS Delta 4-momentum
				pd_=D4V<double>(W_,0,0,0);
				//hadronic CMS nucleon momentum
				D4V<double> p_=pd_-q4_;
				// some other 4-vectors and invariants
				qp=q4_*p_;
				qslash=dirac*q4_;
				qpd=q4_*pd_;
				pdslash=dirac*pd_;
				pdslash+=Md;
				pdslash*=-1;
				pqslash=dirac*pd_;
				pqslash+=M;
				gfivpqslash=gamma5*pqslash;
				
				switch(ffset)
				{
						case 0://model 0A
							deltaffv[0]=2.03/(1.0 + 3.69*Q2_/GeV2 + 1.82*Q2_/GeV2*Q2_/GeV2 + 1.13*Q2_/GeV2*Q2_/GeV2*Q2_/GeV2)/M;
							deltaffv[1]=-Mp*deltaffv[0]/(M*W_);
							deltaffv[2]=0;
							break;
						case 1://model 0B
							deltaffv[0]=2.60*(1+0.89*Q2_/GeV2)/(1.0 + 8.08*Q2_/GeV2 + 0.72*Q2_/GeV2*Q2_/GeV2 + 9.03*Q2_/GeV2*Q2_/GeV2*Q2_/GeV2)/M;
							deltaffv[1]=-2.60*Mp*(1+1.45*Q2_/GeV2)/(1.0 + 8.08*Q2_/GeV2 + 0.72*Q2_/GeV2*Q2_/GeV2 + 9.03*Q2_/GeV2*Q2_/GeV2*Q2_/GeV2)/(M2*W_);
							deltaffv[2]=0;
							break;
						case 2://model II
							deltaffv[0]=2.10*(1+0.13*Q2_/GeV2)/(1.0 + 4.73*Q2_/GeV2 -0.39*Q2_/GeV2*Q2_/GeV2 + 5.59*Q2_/GeV2*Q2_/GeV2*Q2_/GeV2)/M;
							deltaffv[1]=-2.10*Mp*(1+1.68*Q2_/GeV2)/(1.0 + 4.73*Q2_/GeV2 -0.39*Q2_/GeV2*Q2_/GeV2 + 5.59*Q2_/GeV2*Q2_/GeV2*Q2_/GeV2)/(M2*W_);
							deltaffv[2]=0.62/(pow(1.0+Q2_/Mv2,2)*M2);
							break;
						case 3://model I (Our fit of Olga)
							deltaffv[0]=2.00/(1.0+1.15*Q2_/Mv2)/pow(1.0+1.40*Q2_/Mv2,2)/M;
							deltaffv[1]=-6.77/(1.0+1.15*Q2_/Mv2)/pow(1.0+1.40*Q2_/Mv2,2)/M2;
							deltaffv[2]=5.95/((pow(1.0+1.40*Q2_/Mv2,2)*M2)*(1.0+0.68*Q2_/Mv2));
							break;
						default://Olga
							deltaffv[0]=2.13/(1.0+0.25*Q2_/Mv2)/pow(1.0+Q2_/Mv2,2)/M;
							deltaffv[1]=-1.51/2.13*deltaffv[0]/M;
							deltaffv[2]=0.48/((pow(1.0+Q2_/Mv2,2)*M2)*(1.0+1.29*Q2_/Mv2));
							break;
				}
				//set global nucleon ff
				double tau=Q2_/Mn24;
				double dipole=1.0/(1+Q2_/Mv2)/(1+Q2_/Mv2);
				F1P=(1+mup*tau)/(1+tau)*dipole;
				F2P=(mup-1.0)/(1+tau)*dipole;
				F1N=lambdan*mun*tau*tau/((1+tau)*(1+lambdan*tau))*dipole;
				F2N=mun*(tau+1+lambdan*tau)/((1+tau)*(1+lambdan*tau))*dipole;
				F1V=F1P-F1N;
				//set global nucleon em currents (we know Q2 and q0 already...)
				proton_em_v_mu_=q4_*sigma;
				neutron_em_v_mu_=proton_em_v_mu_;
				proton_em_v_mu_*=(-0.5*delusive*F2P/M);
				neutron_em_v_mu_*=(-0.5*delusive*F2N/M);
				proton_em_v_mu_+=(dirac*F1P);
				neutron_em_v_mu_+=(dirac*F1N);
				//now comes the pion in hadronic CMS
				double Epi=sqrt(mpi2+qcm2);
				double sinthetacms=sqrt(1.0-costhetacms*costhetacms);
				D4V<double> k(Epi,qcm*sinthetacms*cos(phicms),qcm*sinthetacms*sin(phicms),qcm*costhetacms);
				D4V<double> pp=pd_;
				pp-=k;
				DM pprimeslash=pp*dirac;
				pprimeslash+=M;
				
				D4V<DM> smu(0,0,0,0);
				//4-vector describing the vertex choice:
				
				
				//Delta Pole
				std::complex<double> denom(0,0);
				D4T<DM> tmp=gmunu;
				tmp.Add((pd_/(-1.5*Md2)),pd_);
				tmp.Add((dirac/-3),dirac);
				tmp.Add((pd_/(3*Md)),dirac);
				tmp.Add((dirac/(-3*Md)),pd_);	
				tmp&=pdslash;
				smu=k*tmp;
				denom=W2_-Md2+delusive*Md*DW;
				denom=1.0/denom;
				
				tmp=gmunu;
				tmp*=(qslash*deltaffv[0] + qpd*deltaffv[1] + qp*deltaffv[2]);
				tmp.Add((q4_*(-deltaffv[0])),dirac);
				tmp.Add((q4_*(-deltaffv[1])),pd_);
				tmp.Add((q4_*(-deltaffv[2])),p_);
				tmp*=gamma5;
				
				
				smu=smu*tmp;
				smu*=denom*fstar/mpi*CGDP[chan];
				D4V<DM> tmpv(0,0,0,0);
				D4V<double> pcr=p_;
				pcr-=k;
				D4V<double> pcrq=pcr;
				pcrq+=q4_;
				double W2cr=pcr*pcr;
				
				//crossed Delta pole
				
				denom=W2cr-Md2;
				denom=1.0/denom;
				DM pcrs=((dirac*pcr)+Md);
				pcrs*=-1;
				tmp=gmunu;
				tmp.Add((pcr/(-1.5*Md2)),pcr);
				tmp.Add((dirac/-3),dirac);
				tmp.Add((pcr/(3*Md)),dirac);
				tmp.Add((dirac/(-3*Md)),pcr);
				tmp&=pcrs;
				tmpv=tmp*k;
				double qpcrq=-1.0;
				qpcrq*=q4_*pcrq;
				
				double qpcr=q4_*pcr;
				qpcr*=-1;
				tmp=gmunu;
				tmp*=(qslash*(-deltaffv[0]) + qpcr*deltaffv[1] + qpcrq*deltaffv[2]);
				tmp.Add(dirac,(q4_*(deltaffv[0])));
				tmp.Add(pcr,(q4_*(deltaffv[1])));
				tmp.Add(pcrq,(q4_*(deltaffv[2])));
				tmp*=gamma5;
				
				tmpv=tmp*tmpv;
				tmpv*=denom*fstar/mpi*CGDPC[chan];
				smu+=tmpv;
					
				DM kslash=k*dirac;
				
				double t=(q4_-k)*(q4_-k);
				
				//optional pion form factor a'la HNV PRD76
				double F_pion=1.0;
				if(FP)
				{
					F_pion=(lambdapi2-mpi2)/(lambdapi2-t);
				}
				//nucleon pole
				if(chan%2==1) tmpv=proton_em_v_mu_;
				else tmpv=neutron_em_v_mu_;
				switch(PV)
				{
					case 0:
					{
						tmpv&=gfivpqslash;
						tmpv*=(M*gA/fpi*sqrt2*CGNP[chan]*F_pion/(W2_-M2) );
						break;
					}
					case 1:
					{
						tmpv&=(kslash*gfivpqslash);
						tmpv*=(-gA/sqrt2*CGNP[chan]*F_pion/( fpi * (W2_-M2) ));
						break;
					}
					default:
					{
						tmpv&=(kslash*gfivpqslash);
						tmpv*=(-gA/sqrt2*CGNP[chan]*F_pion/( fpi * (W2_-M2) ) );
						break;
					}
					
				}
				
				smu+=tmpv;
				
				//nucleon pole crossed
				
				tmpv= (chan<3) ? proton_em_v_mu_ 
				               : neutron_em_v_mu_;
				
				DM pcrslash=pcr*dirac;
				pcrslash+=M;
				
				switch(PV)
				{
					case 0:
					{
						pcrslash*=gamma5;
						tmpv*=pcrslash;
						tmpv*=(M*gA/fpi*sqrt2*CGNPC[chan]*F_pion/(W2cr-M2) );
						break;
					}
					case 1:
					{
						pcrslash*=(kslash*gamma5);
						tmpv*=pcrslash;
						tmpv*=(-gA/sqrt2*CGNPC[chan]*F_pion/( fpi * (W2cr-M2) ) );
						break;
					}
					default:
					{
						pcrslash*=(kslash*gamma5);
						tmpv*=pcrslash;
						tmpv*=(-gA/sqrt2*CGNPC[chan]*F_pion /( fpi * (W2cr-M2)) );
						break;
					}
				}
			
				smu+=tmpv;
				//contact term exists only in 2 and 3 (charged pi channels)
				if((chan!=1)&&(chan!=4))
				{
					
					if(PV)
					{
						tmpv= dirac;
						tmpv*=((-1/sqrt2*CGOTH[chan]*F_pion/fpi)*(gamma5*gA*F1V));
						smu+=tmpv;
					}
					//pion-in-flight
					D4V<double> k2q=k;
					k2q*=2;
					k2q-=q4_;
					tmpv= k2q*gamma5;
					
					tmpv*=(-(1/sqrt2*F1V*2*M*gA/fpi*CGOTH[chan]*F_pion)/( (t-mpi2) ) );
					
					smu+=tmpv;
				}
				D4V<DM> snu=smu;
			
			//////Here you do (\slash{p'}+M)s^\mu(\slash{p}+M)
				smu &=pprimeslash;	
				const DM pslash=((p_*dirac)+M);
				smu*=pslash;
				
			///////here you do \gamma^0(s^{\nu})^\dag\gamma^0
				snu.hermit();
				snu &=gamma0;
				snu *=gamma0;
				
				//you construct hadronic tensor
				D4T<DM> AMUNU(smu,snu);
				result=(AMUNU.contrd(AMUNU,Lept_).Trace());
				result*=qcm*alpha*alpha*lprimemu(0)/(64*Pi*Pi*Pi*W_*sqrt(pow(lmu*pn,2.0)-mass2*M2)*Q2_*Q2_*cm2);
			}
		}
	}
	return result;
}

//laboratory electron energy, outgoing electron 4-momentum, hadronic CMS pion solid angle w.r.t hadronic cms delta direction
//Done in laboratory frame: electron coming along the z-axis, any outgoing electron angle and energy any direction of nucleon
//MASSIVE lepton and Oset Delta Selfenergy effects included.
double dsigma_dq0_dOmegal_dOmegaCMS_Oset(double E, double mass, D4V<double> lprimemu, double costhetacms, double phicms, D4V<double> pn, double KF, double rhorel)
{
	double result=0;
	
	const double mass2=mass*mass;
	
	//laboratory frame electron 4-momentum
	D4V <double> lmu(E,0,0,sqrt(E*E-mass2));
	//hadronic invariant mass squared
	W2_=((lmu-lprimemu+pn)*(lmu-lprimemu+pn));
	//first kinematic cut:can you produce a pion?
	if(W2_>(M+mpi)*(M+mpi))
	{
		W_=sqrt(W2_);
		//sum of electron and nucleon 4-momenta
		D4V <double> lpn=lmu+pn;
		//lepton-nucleon invariant mass squared
		double s=lpn*lpn;
		double sqrts=sqrt(s);
		
		//lepton-nucleon CMS "delta" momentum, massive electron assumption
		double pdelta=(s*s+mass2*mass2+W2_*W2_-2*s*W2_-2*mass2*s-2*mass2*W2_)/4/s;
		//second kinematic cut
		if(pdelta>0)
		{
			//lepton-nucleon CMS lepton energy
			double Ecms=(lpn*lmu+mass2)/sqrts;
			double lcms=sqrt(Ecms*Ecms-mass2);
			//outgoing lepton energy
			double Epcms=sqrt(pdelta+mass2);
			pdelta=sqrt(pdelta);
			// 4-momentum transfer squared
			Q2_=-((lmu-lprimemu)*(lmu-lprimemu));
			//neutrino laboratory energy, Q2, W, hadronic CMS pion solid angle
			//lepton scattering angle cosine in the lepton-nucleon CMS, massless electron
			double costhetal=(-Q2_+2*Ecms*Epcms-2*mass2)/2/pdelta/lcms;
			//third kinematic cut, just in case
			if(fabs(costhetal)<=1)
			{
				//since we work in CMS we have incoming electron and incoming nucleon (colinear)
				//and outgoing electron, outgoing hadronic system with invariant mass W and momentum=pdelta
				//and we want the z axis to go along outgoing hadronic system momentum and x-z lepton scattering plane always
				//and now with one boost one goes to hadronic CMS, where we define pionic angles
				double p0del=sqrt(W2_+pdelta*pdelta);
				//Incoming electron 4-momentum
				double sinthetal=sqrt(1.0-costhetal*costhetal);
				D4V <double> _lcms(((Ecms*p0del+lcms*pdelta*costhetal)/W_),lcms*sinthetal,0,-lcms*costhetal-pdelta*Ecms/W_-lcms*pdelta*pdelta*costhetal/W_/(W_+p0del));
				//Outgoing electron 4-momentum
				D4V <double> _lprimecms((Epcms*p0del+pdelta*pdelta)/W_,0,0,-pdelta-pdelta*Epcms/W_-pdelta*pdelta*pdelta/W_/(p0del+W_));
				//hadronic CMS 4-momentum transfer
				q4_=_lcms-_lprimecms;
				//result=-(q4_*q4_)/Q2_/GeV/4/Pi/1e33;
				//first comes g^{\mu\nu} of the type double:
				Lept_=gmunud;
				//now it's -g^{\mu\nu} Q^2/2
				Lept_ *= -Q2_/2;
				//you add l1^\mu l2^\nu part
				Lept_.Add(_lcms,_lprimecms);
				//you add l2^\nu l2^\mu part:
				Lept_.Add(_lprimecms,_lcms);
				//electron leptonic tensor is ready
				
				//now the rest of these bastards
				double lamb=sqrt(W2_*W2_+M4+mpi2*mpi2-2*(W2_*mpi2+W2_*M2+mpi2*M2));
			
				double qcm=0.5*lamb/W_;
				
				double qcm2=qcm*qcm;
				double qcm3=qcm2*qcm;
				double nominator=fstar*fstar*qcm3*(sqrt(M2+qcm2)+M);
				
				double denominator = 12*Pi*mpi2*W_;
				DW=nominator/denominator;
				
				//hadronic CMS Delta 4-momentum
				pd_=D4V<double>(W_,0,0,0);
				//hadronic CMS nucleon momentum
				D4V<double> p_=pd_-q4_;
				// some other 4-vectors and invariants
				qp=q4_*p_;
				qslash=dirac*q4_;
				qpd=q4_*pd_;
				pdslash=dirac*pd_;
				pdslash+=Md;
				pdslash*=-1;
				pqslash=dirac*pd_;
				pqslash+=M;
				gfivpqslash=gamma5*pqslash;
				
				switch(ffset)
				{
						case 0://model 0A
							deltaffv[0]=2.03/(1.0 + 3.69*Q2_/GeV2 + 1.82*Q2_/GeV2*Q2_/GeV2 + 1.13*Q2_/GeV2*Q2_/GeV2*Q2_/GeV2)/M;
							deltaffv[1]=-Mp*deltaffv[0]/(M*W_);
							deltaffv[2]=0;
							break;
						case 1://model 0B
							deltaffv[0]=2.60*(1+0.89*Q2_/GeV2)/(1.0 + 8.08*Q2_/GeV2 + 0.72*Q2_/GeV2*Q2_/GeV2 + 9.03*Q2_/GeV2*Q2_/GeV2*Q2_/GeV2)/M;
							deltaffv[1]=-2.60*Mp*(1+1.45*Q2_/GeV2)/(1.0 + 8.08*Q2_/GeV2 + 0.72*Q2_/GeV2*Q2_/GeV2 + 9.03*Q2_/GeV2*Q2_/GeV2*Q2_/GeV2)/(M2*W_);
							deltaffv[2]=0;
							break;
						case 2://model II
							deltaffv[0]=2.10*(1+0.13*Q2_/GeV2)/(1.0 + 4.73*Q2_/GeV2 -0.39*Q2_/GeV2*Q2_/GeV2 + 5.59*Q2_/GeV2*Q2_/GeV2*Q2_/GeV2)/M;
							deltaffv[1]=-2.10*Mp*(1+1.68*Q2_/GeV2)/(1.0 + 4.73*Q2_/GeV2 -0.39*Q2_/GeV2*Q2_/GeV2 + 5.59*Q2_/GeV2*Q2_/GeV2*Q2_/GeV2)/(M2*W_);
							deltaffv[2]=0.62/(pow(1.0+Q2_/Mv2,2)*M2);
							break;
						case 3://model I (Our fit of Olga)
							deltaffv[0]=2.00/(1.0+0.68*Q2_/Mv2)/pow(1.0+1.15*Q2_/Mv2,2)/M;
							deltaffv[1]=-6.77/2.00*deltaffv[0]/M;
							deltaffv[2]=5.95/((pow(1.0+1.15*Q2_/Mv2,2)*M2)*(1.0+1.40*Q2_/Mv2));
							break;
						default://Olga
							deltaffv[0]=2.13/(1.0+0.25*Q2_/Mv2)/pow(1.0+Q2_/Mv2,2)/M;
							deltaffv[1]=-1.51/2.13*deltaffv[0]/M;
							deltaffv[2]=0.48/((pow(1.0+Q2_/Mv2,2)*M2)*(1.0+1.29*Q2_/Mv2));
							break;
				}
				//set global nucleon ff
				double tau=Q2_/Mn24;
				double dipole=1.0/(1+Q2_/Mv2)/(1+Q2_/Mv2);
				F1P=(1+mup*tau)/(1+tau)*dipole;
				F2P=(mup-1.0)/(1+tau)*dipole;
				F1N=lambdan*mun*tau*tau/((1+tau)*(1+lambdan*tau))*dipole;
				F2N=mun*(tau+1+lambdan*tau)/((1+tau)*(1+lambdan*tau))*dipole;
				F1V=F1P-F1N;
				//set global nucleon em currents (we know Q2 and q0 already...)
				proton_em_v_mu_=q4_*sigma;
				neutron_em_v_mu_=proton_em_v_mu_;
				proton_em_v_mu_*=(-0.5*delusive*F2P/M);
				neutron_em_v_mu_*=(-0.5*delusive*F2N/M);
				proton_em_v_mu_+=(dirac*F1P);
				neutron_em_v_mu_+=(dirac*F1N);
				//now comes the pion in hadronic CMS
				double Epi=sqrt(mpi2+qcm2);
				double sinthetacms=sqrt(1.0-costhetacms*costhetacms);
				D4V<double> k(Epi,qcm*sinthetacms*cos(phicms),qcm*sinthetacms*sin(phicms),qcm*costhetacms);
				D4V<double> pp=pd_;
				pp-=k;
				DM pprimeslash=pp*dirac;
				pprimeslash+=M;
				
				D4V<DM> smu(0,0,0,0);
				//4-vector describing the vertex choice:
				
				
				//Delta Pole
				std::complex<double> denom(0,0);
				D4T<DM> tmp=gmunu;
				tmp.Add((pd_/(-1.5*Md2)),pd_);
				tmp.Add((dirac/-3),dirac);
				tmp.Add((pd_/(3*Md)),dirac);
				tmp.Add((dirac/(-3*Md)),pd_);	
				tmp&=pdslash;
				smu=k*tmp;
				
				double SQEL = 0;
				double SA2  = 0;
				double SA3  = 0;
				double RS=0;
				double PBF=1.0;
				//Oset Selfenergy and Pauli blocking effects on Delta width
				if(selfenergy)
				{
					double _C_qel=0;
					double _C_a2=0;
					double _C_a3=0;
					double _alpha=1;
					double _beta=1;
					double _gamma=2*_beta;
					double _omegax=(W2_-M2)/(2*sqrt(M*M+0.6*KF*KF));
					if (_omegax<0) _omegax=0;
					if(_omegax>0)
					{
						_C_qel=Cgammax(0,_omegax);
						_C_a2=Cgammax(1,_omegax);
						_C_a3=Cgammax(2,_omegax);
						_alpha=Cgammax(3,_omegax);
						_beta=Cgammax(4,_omegax);
						_gamma=2*_beta;
					}
					RS=2*40*MeV*rhorel;
					SQEL = 2*_C_qel*pow(rhorel,_alpha);
					SA2  = 2*_C_a2*pow(rhorel,_beta);
					SA3  = 2*_C_a3*pow(rhorel,_gamma);
					//Phase-space factor for Pauli Blocking
					D4V<double> pdlab=lmu-lprimemu+pn;
					double vpdlab=sqrt(pdlab(1)*pdlab(1)+pdlab(2)*pdlab(2)+pdlab(3)*pdlab(3));
					double cpicmsmax=((pn(0)+E-lprimemu(0))*p_(0)-sqrt(KF*KF+M2)*W_)/(vpdlab*qcm);
					if(cpicmsmax>1) cpicmsmax=1;
					PBF=0.5*(cpicmsmax+1);
				}
				double STOT = SQEL + SA2 + SA3;
				
				
				
				denom=W2_-Md2+delusive*Md*(DW+STOT)-Md*RS;
				denom=1.0/denom;
				
				tmp=gmunu;
				tmp*=(qslash*deltaffv[0] + qpd*deltaffv[1] + qp*deltaffv[2]);
				tmp.Add((q4_*(-deltaffv[0])),dirac);
				tmp.Add((q4_*(-deltaffv[1])),pd_);
				tmp.Add((q4_*(-deltaffv[2])),p_);
				tmp*=gamma5;
				
				
				smu=smu*tmp;
				smu*=denom*fstar/mpi*CGDP[chan];
				D4V<DM> tmpv(0,0,0,0);
				D4V<double> pcr=p_;
				pcr-=k;
				D4V<double> pcrq=pcr;
				pcrq+=q4_;
				double W2cr=pcr*pcr;
				
				//crossed Delta pole
				
				denom=W2cr-Md2;
				denom=1.0/denom;
				DM pcrs=((dirac*pcr)+Md);
				pcrs*=-1;
				tmp=gmunu;
				tmp.Add((pcr/(-1.5*Md2)),pcr);
				tmp.Add((dirac/-3),dirac);
				tmp.Add((pcr/(3*Md)),dirac);
				tmp.Add((dirac/(-3*Md)),pcr);
				tmp&=pcrs;
				tmpv=tmp*k;
				double qpcrq=-1.0;
				qpcrq*=q4_*pcrq;
				
				double qpcr=q4_*pcr;
				qpcr*=-1;
				tmp=gmunu;
				tmp*=(qslash*(-deltaffv[0]) + qpcr*deltaffv[1] + qpcrq*deltaffv[2]);
				tmp.Add(dirac,(q4_*(deltaffv[0])));
				tmp.Add(pcr,(q4_*(deltaffv[1])));
				tmp.Add(pcrq,(q4_*(deltaffv[2])));
				tmp*=gamma5;
				
				tmpv=tmp*tmpv;
				tmpv*=denom*fstar/mpi*CGDPC[chan];
				smu+=tmpv;
					
				DM kslash=k*dirac;
				
				double t=(q4_-k)*(q4_-k);
				
				//optional pion form factor a'la HNV PRD76
				double F_pion=1.0;
				if(FP)
				{
					F_pion=(lambdapi2-mpi2)/(lambdapi2-t);
				}
				//nucleon pole
				if(chan%2==1) tmpv=proton_em_v_mu_;
				else tmpv=neutron_em_v_mu_;
				switch(PV)
				{
					case 0:
					{
						tmpv&=gfivpqslash;
						tmpv*=(M*gA/fpi*sqrt2*CGNP[chan]*F_pion/(W2_-M2) );
						break;
					}
					case 1:
					{
						tmpv&=(kslash*gfivpqslash);
						tmpv*=(-gA/sqrt2*CGNP[chan]*F_pion/( fpi * (W2_-M2) ));
						break;
					}
					default:
					{
						tmpv&=(kslash*gfivpqslash);
						tmpv*=(-gA/sqrt2*CGNP[chan]*F_pion/( fpi * (W2_-M2) ) );
						break;
					}
					
				}
				
				smu+=tmpv;
				
				//nucleon pole crossed
				
				tmpv= (chan<3) ? proton_em_v_mu_ 
				               : neutron_em_v_mu_;
				
				DM pcrslash=pcr*dirac;
				pcrslash+=M;
				
				switch(PV)
				{
					case 0:
					{
						pcrslash*=gamma5;
						tmpv*=pcrslash;
						tmpv*=(M*gA/fpi*sqrt2*CGNPC[chan]*F_pion/(W2cr-M2) );
						break;
					}
					case 1:
					{
						pcrslash*=(kslash*gamma5);
						tmpv*=pcrslash;
						tmpv*=(-gA/sqrt2*CGNPC[chan]*F_pion/( fpi * (W2cr-M2) ) );
						break;
					}
					default:
					{
						pcrslash*=(kslash*gamma5);
						tmpv*=pcrslash;
						tmpv*=(-gA/sqrt2*CGNPC[chan]*F_pion /( fpi * (W2cr-M2)) );
						break;
					}
				}
			
				smu+=tmpv;
				//contact term exists only in 2 and 3 (charged pi channels)
				if((chan!=1)&&(chan!=4))
				{
					
					if(PV)
					{
						tmpv= dirac;
						tmpv*=((-1/sqrt2*CGOTH[chan]*F_pion/fpi)*(gamma5*gA*F1V));
						smu+=tmpv;
					}
					//pion-in-flight
					D4V<double> k2q=k;
					k2q*=2;
					k2q-=q4_;
					tmpv= k2q*gamma5;
					
					tmpv*=(-(1/sqrt2*F1V*2*M*gA/fpi*CGOTH[chan]*F_pion)/( (t-mpi2) ) );
					
					smu+=tmpv;
				}
				D4V<DM> snu=smu;
			
			//////Here you do (\slash{p'}+M)s^\mu(\slash{p}+M)
				smu &=pprimeslash;	
				const DM pslash=((p_*dirac)+M);
				smu*=pslash;
				
			///////here you do \gamma^0(s^{\nu})^\dag\gamma^0
				snu.hermit();
				snu &=gamma0;
				snu *=gamma0;
				
				//you construct hadronic tensor
				D4T<DM> AMUNU(smu,snu);
				result=(AMUNU.contrd(AMUNU,Lept_).Trace());
				result*=qcm*alpha*alpha*lprimemu(0)/(64*Pi*Pi*Pi*W_*sqrt(pow(lmu*pn,2.0)-mass2*M2)*Q2_*Q2_*cm2);
			}
		}
	}
	return result;
}
