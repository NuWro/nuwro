#include "mecevent_tem.h"

void mecevent_tem (params & p, event & e, nucleus & t, bool cc)
{
	bool nu = e.in[0].pdg > 0;	//neutrino or antineutrino
	e.par = p;					//save params in the event
	e.flag.cc = cc;				//set flags for the event
	e.flag.nc = !cc;
	e.flag.dis = false;
	e.flag.qel = false;
	e.flag.coh = false;
	e.flag.mec = true;
	
	if(t.A() < 4)					//set cross section = 0 if nuclei has less nucleons than 4
	{
		e.weight = 0;
		return;
	}
	
	double potwell  = t.Ef();		//Fermi energy
	double ebinding = t.Eb();		//Binding energy

	particle meclepton;

	meclepton.pdg = e.in[0].pdg;	//for nc outgoing lepton = incoming neutrino

	if (cc)							//for charged current set outgoing lepton pdg
	{
		if (nu) 
			meclepton.pdg--;
		else
			meclepton.pdg++;
	}

	meclepton.set_mass (PDG::mass (meclepton.pdg));	//set mass coresponding to pdg
	
	if (blocked(e.in[0].E(), meclepton.mass()))	//check if neutrino energy is high enough
	{
		e.weight = 0;
		return;
	}

	particle mecnucleon[4]; //0,1 - in, 2, 3 - out
	
	tem_kin (e.in[0].E(), meclepton, mecnucleon, t);
	
	double weight;
	
	if (cc)
		weight = mec_do_cc (weight, e.in[0].E(), mecnucleon, meclepton.mass(), p.mec_ratio_pp, nu);
	else
		weight = mec_do_nc (weight, e.in[0].E(), mecnucleon, meclepton.mass(), p.mec_ratio_pp, nu);
	
	if (weight <= 0)
	{
		e.weight = 0;
		return;
	}
	
	e.weight = dif * weight / cm2;
		
	e.in.push_back (mecnucleon[0]);
	e.in.push_back (mecnucleon[1]);
	
	e.out.push_back (meclepton);
	e.out.push_back (mecnucleon[2]);
	e.out.push_back (mecnucleon[3]);
}

void tem_kin (double E, particle &meclep, particle *nucleon, nucleus &t)
{
	//it is assumed that neutrino direction is (0,0,1); but transition to other direction in nuwro.cc!
	
	double mecm  = meclep.mass();
	double mecm2 = mecm * mecm;
	
	Q2 = setQ2 (E, mecm2);

	double w = Q2 / 2.0 / M1;	//energy transfer
	double q = sqrt(Q2 + w * w);
	
	double muonmom = sqrt ((E - w) * (E - w) - mecm2);
	double cosmuon = (2.0 * E * (E - w) - mecm2 - q * q + w * w) / 2.0 / E / muonmom;
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

	do
	{
		nucleon[0] = t.get_nucleon ();
		vec pos (nucleon[0].x, nucleon[0].y, nucleon[0].z);
		nucleon[1] = t.get_nucleon (pos);
					
		//nucleon[1] = t.get_nucleon ();
		//nucleon[1].r=nucleon[0].r;//both should start from the same point !
						
		suma = nucleon[0].p4() + nucleon[1].p4() + qqq;
	}								 //to be able to make Lorentz boost and to "decay"
	while ( ( suma.t < suma.length() )   ||   ( suma*suma < 4.0 * M2) );
		
	nucleon[2].set_proton();	//in mec_do* we decide about isospin	
	nucleon[3].set_proton();
	
	vec trans = suma.v();
		
	suma.boost2 (trans);		 	//boost to the CM frame

	double Ecm = suma.t / 2.0;	//each nucleon get the same energy in CM
		
	vec dir_cm = rand_dir () * sqrt(Ecm * Ecm  - M2);	//radnomly set direction of nucleons in CM

	nucleon[2].set_momentum (dir_cm);
	nucleon[3].set_momentum (-dir_cm);
	
	nucleon[2].p4().boost2 (-trans);
	nucleon[3].p4().boost2 (-trans);		 //nucleons in the LAB frame
	
	meclep.set_momentum(kprim);	
}

double mec_do_cc (double &w, double E, particle *p, double m, double ratio, bool nu)
{
	double weight = qel_sigma (E, -Q2, 3, !nu, m, M1) - qel_sigma (E, -Q2, 6, !nu, m, M1);

	// here is the isospin model; I assume that 3/5 times a pair is p-p and 2/5 times it is p-n
	if (frandom () < ratio) 
	{
		if (nu)
		{
			p[0].set_proton ();
			p[1].set_neutron ();
			p[2].set_proton ();
			p[3].set_proton ();
		}
		else
		{
			p[0].set_proton ();
			p[1].set_neutron ();			
			p[2].set_neutron ();
			p[3].set_neutron ();				
		}
	}
	else
	{
		if (nu)
		{
			p[0].set_neutron ();
			p[1].set_neutron ();
		}
		else
		{
			p[0].set_proton ();
			p[1].set_proton ();
		}
		
		p[2].set_proton ();
		p[3].set_neutron ();
	}
	
	return weight / 2.0;
}

double mec_do_nc (double &w, double E, particle *p, double m, double ratio, bool nu)
{
	double weight_proton  = qel_sigma (E, -Q2, 4, !nu, m, M1) - qel_sigma (E, -Q2, 7, !nu, m, M1);
	double weight_neutron = qel_sigma (E, -Q2, 5, !nu, m, M1) - qel_sigma (E, -Q2, 8, !nu, m, M1);
				
	double frac = weight_proton / (weight_neutron + weight_proton);
	
	double ratio_nc = 1.0 / (2.0 * ratio + 1);
		
	if (frandom () < ratio_nc)
	{
		if (frandom () < frac)
		{
			p[0].set_proton ();
			p[1].set_proton ();
			p[2].set_proton ();
			p[3].set_proton ();
		}
		else
		{
			p[0].set_neutron ();
			p[1].set_neutron ();
			p[2].set_neutron ();
			p[3].set_neutron ();
		}
	}
	else
	{
		p[0].set_neutron ();
		p[1].set_proton ();
		p[2].set_neutron ();
		p[3].set_proton ();
	}	
	
	return (weight_proton + weight_neutron) / 2.0;
}
