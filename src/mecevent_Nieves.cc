#include "mecevent_Nieves.h"
//int noev=0;
void mecevent_Nieves(params & p, event & e, nucleus & t, bool cc)
{
	e.par = p;				//save params in the event
	e.flag.cc = cc;				//set flags for the event
	e.flag.nc = !cc;
	e.flag.dis = false;
	e.flag.qel = false;
	e.flag.coh = false;
	e.flag.mec = true;
	
	
	//sadly, only \nu_\mu CC o ^{12}C available so far...
	if((t.p!=6)&&(t.n!=6)&&(e.flag.nc)&&(e.in[0].pdg!=14))
	{
		cerr<<" MEC error: Wrong Settings!\n";
		e.weight = 0;
		return;
	}
	
	double potwell  = t.Ef();		//Fermi energy
	double ebinding = t.Eb();		//Binding energy

	particle meclepton;

	meclepton.pdg = 13;	//only muon available

	meclepton.set_mass (PDG::mass (meclepton.pdg));	//set mass coresponding to pdg
	
	width_T=e.in[0].E()-mmu;
	
	if (blocked(e.in[0].E(), meclepton.mass()))	//check if neutrino energy is high enough
	{
		e.weight = 0;
		return;
	}

	particle mecnucleon[4]; //0,1 - in, 2, 3 - out
	
	
	double weight;
	
	weight=Nieves_kin_and_weight (e.in[0].E(), meclepton, mecnucleon, t);
	e.weight = weight*1e-41*width_T/12;
	if (weight>0) Nieves_do_cc ( mecnucleon,  p.mec_ratio_pp);
	e.in.push_back (mecnucleon[0]);
	e.in.push_back (mecnucleon[1]);
	e.out.push_back (meclepton);
	e.out.push_back (mecnucleon[2]);
	e.out.push_back (mecnucleon[3]);
	//noev++;
	//cerr<<"              event "<<noev<<"\r";
}

double Nieves_kin_and_weight (double E, particle &meclep, particle *nucleon, nucleus &t)
{
	//it is assumed that neutrino direction is (0,0,1); but transition to other direction in nuwro.cc!
	
	double mecm  = mmu;
	double mecm2 = mecm * mecm;
	
	//here we extrapolate the cross section from Nieves's data file
	double result=0;
	//muon kinetic energy
	double Tmu=width_T*frandom();
	double w = E-Tmu-mmu;	//energy transfer
	double muonmom = sqrt (Tmu*Tmu+2*Tmu* mmu);
	double kf2=t.kF();
	kf2*=kf2;
	double vqmax2=8*w*w+4*(sqrt(MN2+4*kf2)*w+4*kf2);
	double cosmin=0.5*(muonmom*muonmom+E*E-vqmax2)/E/muonmom;
	if(cosmin<-1) cosmin=-1;
	//double cosmin=-1;
	double cosmuon=cosmin+(1-cosmin)*frandom();
	double q = sqrt(muonmom*muonmom+E*E-2*muonmom*E*cosmuon);
	if(cosmin<1)
	{	
		//Are we in the data range etc?
		if((Tmu>=0)&&(E>154.8959*MeV)&&(E<2995.104*MeV))
		{
			//above angular binning ->move the point inside border		
			if(cosmuon<=-0.9965643) cosmuon=-0.99656429;
			if(cosmuon>=0.9965643) cosmuon=0.99656429;
			
			int pos_E=0;
			int pos_E0=0;
			//search for nearest smaller energy		
			for(int i=0; i<40; i++)
			{
				if(List_E_[i]>E)
				{
					pos_E=(i-1)*1600;
					pos_E0=i-1;
					break;
				}
			}
			int pos_cos=0;
			int pos_cos0=0;
			//search for nearest smaller angle	
			for(int i=0; i<40; i++)
			{
				if(List_cos_[i]>cosmuon)
				{
					pos_cos=(i-1)*40;
					pos_cos0=i-1;
					break;
				}
			}
			int pos_T1=0;
			int border_T1=0;
			//search for nearest smallest muon energy at lower neutrino energy	
			for(int i=0; i<40; i++)
			{
				if(List_T_[pos_E0*40+i]>Tmu)
				{
					pos_T1=i-1;
					break;
				}
			}
			//above limits?
			if(List_T_[pos_E0*40+39]<Tmu)
			{
				border_T1=1;
				pos_T1=39;
			}
			//below limits?
			if(pos_T1<0)
			{
				border_T1=-1;
				pos_T1=0;
			}
			//because of different Tmu binnings...
			int pos_T2=0;
			int border_T2=0;
			//search for nearest smallest muon energy at highe neutrino energy	
			for(int i=0; i<40; i++)
			{
				if(List_T_[(pos_E0+1)*40+i]>Tmu)
				{
					pos_T2=i-1;
					break;
				}
			}
			//above limits?
			if(List_T_[(pos_E0+1)*40+39]<Tmu)
			{
				border_T2=1;
				pos_T2=39;
			}
			//below limits?
			if(pos_T2<0)
			{
				border_T2=-1;
				pos_T2=0;
			}
			
			//positions of the cube vertices
			unsigned int p[8]={pos_E+pos_cos+pos_T1          ,pos_E+pos_cos+40+pos_T1       ,
							   pos_E+pos_cos+pos_T1+1        ,pos_E+pos_cos+40+pos_T1+1     ,
							   pos_E+1600+pos_cos+pos_T2     ,pos_E+1600+pos_cos+40+pos_T2  ,
							   pos_E+1600+pos_cos+pos_T2+1  ,pos_E+1600+pos_cos+40+pos_T2+1};
	
			//inside data grid_ case
			if((border_T1==0)&&(border_T2==0))
			{
				double downcosmin=Nieves_XS_grid_[p[0]]+(cosmuon-List_cos_[pos_cos0])/(List_cos_[pos_cos0+1]-List_cos_[pos_cos0])*(Nieves_XS_grid_[p[1]]-Nieves_XS_grid_[p[0]]);
				double downcosmax=Nieves_XS_grid_[p[2]]+(cosmuon-List_cos_[pos_cos0])/(List_cos_[pos_cos0+1]-List_cos_[pos_cos0])*(Nieves_XS_grid_[p[3]]-Nieves_XS_grid_[p[2]]);
				double downT=downcosmin+(Tmu-List_T_[pos_E0*40+pos_T1])/(List_T_[pos_E0*40+pos_T1+1]-List_T_[pos_E0*40+pos_T1])*(downcosmax-downcosmin);
				double lower_E=List_E_[pos_E0];
	
	
				double upcosmin=Nieves_XS_grid_[p[4]]+(cosmuon-List_cos_[pos_cos0])/(List_cos_[pos_cos0+1]-List_cos_[pos_cos0])*(Nieves_XS_grid_[p[5]]-Nieves_XS_grid_[p[4]]);
				double upcosmax=Nieves_XS_grid_[p[6]]+(cosmuon-List_cos_[pos_cos0])/(List_cos_[pos_cos0+1]-List_cos_[pos_cos0])*(Nieves_XS_grid_[p[7]]-Nieves_XS_grid_[p[6]]);
				double upT=upcosmin+(Tmu-List_T_[(pos_E0+1)*40+pos_T2])/(List_T_[(pos_E0+1)*40+pos_T2+1]-List_T_[(pos_E0+1)*40+pos_T2])*(upcosmax-upcosmin);
				double upper_E=List_E_[pos_E0+1];
				//linear interpolation in energy
				result=downT+(E-lower_E)/(upper_E-lower_E)*(upT-downT);
					
				
			}
			//above Tmu in lower energy, but inside for higher
			if((border_T1>0)&&(border_T2==0))
			{
				
				
				double lower_E=List_E_[pos_E0];
	
	
				double upcosmin=Nieves_XS_grid_[p[4]]+(cosmuon-List_cos_[pos_cos0])/(List_cos_[pos_cos0+1]-List_cos_[pos_cos0])*(Nieves_XS_grid_[p[5]]-Nieves_XS_grid_[p[4]]);
				double upcosmax=Nieves_XS_grid_[p[6]]+(cosmuon-List_cos_[pos_cos0])/(List_cos_[pos_cos0+1]-List_cos_[pos_cos0])*(Nieves_XS_grid_[p[7]]-Nieves_XS_grid_[p[6]]);
				double upT=upcosmin+(Tmu-List_T_[(pos_E0+1)*40+pos_T2])/(List_T_[(pos_E0+1)*40+pos_T2+1]-List_T_[(pos_E0+1)*40+pos_T2])*(upcosmax-upcosmin);
				double upper_E=List_E_[pos_E0+1];
				//linear interpolation in energy
				result=(E-lower_E)/(upper_E-lower_E)*upT;
			}
			if((border_T1==0)&&(border_T2<0))
			{
				
				
				double downcosmin=Nieves_XS_grid_[p[0]]+(cosmuon-List_cos_[pos_cos0])/(List_cos_[pos_cos0+1]-List_cos_[pos_cos0])*(Nieves_XS_grid_[p[1]]-Nieves_XS_grid_[p[0]]);
				double downcosmax=Nieves_XS_grid_[p[2]]+(cosmuon-List_cos_[pos_cos0])/(List_cos_[pos_cos0+1]-List_cos_[pos_cos0])*(Nieves_XS_grid_[p[3]]-Nieves_XS_grid_[p[2]]);
				double downT=downcosmin+(Tmu-List_T_[pos_E0*40+pos_T1])/(List_T_[pos_E0*40+pos_T1+1]-List_T_[pos_E0*40+pos_T1])*(downcosmax-downcosmin);
				double lower_E=List_E_[pos_E0];
				double upper_E=List_E_[pos_E0+1];
				//linear interpolation in energy
				result=(upper_E-E)/(upper_E-lower_E)*downT;
			}
			if(result<0) result=0;
				
		}
		//Why should we bother to do anything, if weight=0??
		if(result>0)
		{
			
			
			//double q = sqrt(muonmom*muonmom+E*E-2*muonmom*E*cosmuon);
			double phi = 2.0 * Pi * frandom();
		
			//momentum transfer
			vec qq (cos(phi) * muonmom * sqrt(1.0 - cosmuon * cosmuon),
					sin(phi) * muonmom * sqrt(1.0 - cosmuon * cosmuon),
					E - muonmom * cosmuon);
		
			//muon momentum
			vec kprim  (-cos(phi) * muonmom * sqrt(1.0 - cosmuon * cosmuon),
						-sin(phi) * muonmom * sqrt(1.0 - cosmuon * cosmuon),
						muonmom * cosmuon);
			
			//muon 4-momentum
			vect qqq (qq, w);
				
			vect suma;
			unsigned licz=0;
			do
			{
				nucleon[0] = t.get_nucleon ();
				nucleon[1] = t.get_nucleon ();
								
				suma = nucleon[0].p4() + nucleon[1].p4() + qqq;
				licz++;
				//if(licz>calls_max) break;	
			}								 //to be able to make Lorentz boost and to "decay"
			while ( (( suma.t < suma.length() ) || ( suma*suma < 4.0 * MN2)) &&(licz<calls_max) );
			if(licz<calls_max)
			{
				result*=(1-cosmin);
				vec trans = suma.v();
				nucleon[2].set_proton();	//in mec_do* we decide about isospin	
				nucleon[3].set_proton();
			
				suma.boost2 (trans);		 	//boost to the CM frame
				double Ecm = suma.t / 2.0;	//each nucleon get the same energy in CM
				
				vec dir_cm = rand_dir () * sqrt(Ecm * Ecm  - MN2);	//radnomly set direction of nucleons in CM
				nucleon[2].set_momentum (dir_cm);
				nucleon[3].set_momentum (-dir_cm);
		
				nucleon[2].p4().boost2 (-trans);
				nucleon[3].p4().boost2 (-trans);		 //nucleons in the LAB frame
			}
			else
			{
				result=0;
				vect zero (0,0,0);
				nucleon[2].set_momentum (zero);
				nucleon[3].set_momentum (zero);
				nucleon[2].set_energy (MN);
				nucleon[3].set_energy (MN);
			}
			//result*=(1-cosmin);
			meclep.set_momentum(kprim);
			
		}	
	}
	
	return result;
	
}

void Nieves_do_cc (particle *p, double ratio)
{
	
	// here is the isospin model; I assume that 3/5 times a pair is p-p and 2/5 times it is p-n
	if (frandom () < ratio) 
	{
		p[0].set_proton ();
		p[1].set_neutron ();
		p[2].set_proton ();
		p[3].set_proton ();
	}
	else
	{
		p[0].set_neutron ();
		p[1].set_neutron ();		
		p[2].set_proton ();
		p[3].set_neutron ();
	}
	
}