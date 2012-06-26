#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include "generatormt.h"
#include "vect.h"
#include <queue>
#include "pdg.h"
#include "nucleus.h"
#include "flatnucleus.h"
#include "anynucleus.h"
#include "event1.h"
#include <TROOT.h>
#include <TTree.h>
#include "beam.h"
//#include "Metropolis.h"
#include "Interaction.h"
#include "proctable.h"
#include "nucleusmaker.h"
#include "fsi.h"

using namespace std;
using namespace PDG;

queue < particle > parts ;

#define LOG(x) cout<< x <<endl;
#define LOG2(x,y) cout<< x << y <<endl;
#define ERR(x) cerr<< x <<endl;
#define ERR2(x,y) cerr<< x << y <<endl;


			
/// process the event 
/// put reintreacition particles to tmp vector ???
/// put outgoing particles to out vector

bool check (particle & p1, particle & p2, particle *spect, int n, particle p[],int k)
{
	int ch1=p1.charge()+p2.charge();
	if(spect) ch1+=spect->charge();
	int i=n;
	while(i)
	  ch1-=p[--i].charge();
	if(ch1!=0)
	{   cout<<endl<<"Proc:"<<k<<endl;
  	    cout<<p1<<endl<<p2<<endl;
		if(spect) 
		  cout<<(*spect)<<endl;
		cout<<endl;
		while(i<n)
		  cout<<p[i++]<<endl;
		cout<<endl;	
	}
	return ch1==0;
}
bool check2 (particle & p1, particle & p2, particle *spect, int n, particle p[],int k)
{
	vect p4=p1+p2;
	if(spect) p4+=*spect;
	int i=n;
	while(i)
	  p4-=p[--i];
	double prec=0.001*MeV;  
	if(abs(p4.t)>prec ||abs(p4.x)>prec ||abs(p4.y)>prec ||abs(p4.z)>prec)
	{   cout<<endl<<"Proc:"<<k<<" delta p4 = "<<p4<<endl;
  	    cout<<p1<<endl<<p2<<endl;
		if(spect) 
		  cout<<(*spect)<<endl;
		cout<<endl;
		while(i<n)
		  cout<<p[i++]<<endl;
		cout<<endl;	
		return false;
	}
	return true;
}


void procinfo(particle p1, particle p2, int n, particle p[])
{

}

//double formation1 (particle &p, params &par, vect q, bool qel, vect p0, bool res);//, int nofpi);
//double formation2 (particle &p, params &par);//, vect q, bool qel, vect p0, bool res);//, int nofpi);

/*

double test_fun(vector<particle>& input,vect q)
    {
	    double t = 0;
		vec help=vec(q);
					
		for (int i = 0; i < input.size (); i++)
		{
  		  particle p = input[i];
		  if (nucleon_or_pion (p.pdg))
		  {
			  double pt = help.dir()*p.p();
			  pt = 0; //sqrt(p.momentum2() - pt*pt);
			 
			  double tau = 0.03;
			  if (p.pdg == 2212 or p.pdg == 2112)
				tau = 0.2;
			 
			  t += tau*p.mass()*p.E()/(p.mass2() + pt*pt);
		  }
		}
		return t;
	}

double test_fun(particle p[], int n)
    {
	    double t = 0;

		for (int i = 0; i < n; i++)
		{
			double tau = 0.03;
			if (p[i].pdg == 2212 or p[i].pdg == 2112)
				tau = 0.2;
				
			  t += tau*p[i].E()/p[i].mass();
		}
		return t;
	}

*/
	
vec test_fun(vector<particle>& input,vect q)
    {
	    double t = 0;
	    double e = 0;
	    vec ped = vec(0,0,0);
	    			
		for (int i = 0; i < input.size (); i++)
		{
  		  particle p = input[i];
		  if (nucleon_or_pion (p.pdg))
		  {
			  double tau = 0.01;
			  if (p.pdg == 2212 or p.pdg == 2112)
				tau = 0.35;
			 
			  t += tau;
			  e += p.E();
			  ped += p.p();
		  }
		}
		double mom = ped.length();
		double m = sqrt(e*e - mom*mom);
		
		return ped*t/m;
	}
	
vec test_fun(particle p[], int n)
    {
	    double t = 0;
	    double e = 0;
	    vec ped = vec(0,0,0);

		for (int i = 0; i < n; i++)
		{
			double tau = 0.01;
			if (p[i].pdg == 2212 or p[i].pdg == 2112)
				tau = 0.35;
				
			  t += tau;
			  e += p[i].E();
			  ped += p[i].p();
		}
		
		double mom = ped.length();
		double m = sqrt(e*e - mom*mom);

		return ped*t/m;
	}
	
	
	
void set_pt(particle &p, vect q, vect p0, bool res, bool dis)
{
	if (res)
	{
		vect delta = q + p0;
		double v = delta.v().length();
		double g = 1.0/sqrt(1-v*v);
		double t = 200.0/120.0;
		p.wfz(g*v*t*fermi);
	}
	else if (dis)
	{
		vec qdir(q.x, q.y, q.z);
								
		double pt = qdir.dir()*p.p();
		pt = sqrt(p.momentum2() - pt*pt);
		
		p.wpt(pt);
		
		double fz = p.momentum()*p.mass()/(p.mass2() + pt*pt);
		
		p.wfz(fz*fermi);
	}
	
	//p.krok(1000.0*fermi);
}

void set_nucl_hist(double *dhist, double &rhist, int &phist, int &nhist, nucleus &nucl)
{
	rhist = nucl.radius()/fermi;
	phist = nucl.Zr();
	nhist = nucl.Nr();
	
	for (int i = 0; i < 50; i++)
	{
		double r = (i + 1.0) * 0.2 * fermi;
		dhist[i] = nucl.density(r)*fermi3;
	}	
}

int
kaskadaevent (params & par, event&e, vector<particle> &input, vector<particle>& output,vector<particle> &all, int nod[], vect q, bool qel, vect pin, int place[][20], bool res, int nofpi, double &pabsen, bool &nopp, bool &abs, int &nofi, double *pen, double &absr, double &odl, double *density_hist, double &radius_hist, int &protons_hist, int &neutrons_hist)
{ 
  int result=0; 
 // double max_step=0.2*fermi; 
  double max_step = par.step * fermi;  
  if(not par.kaskada_on)
    return 0; 
  nucleus* pnucleus=make_nucleus(par);
  nucleus& nucl=*pnucleus;
  
  double Radius = nucl.radius();
  
  if (e.in[0].lepton()) // if there was lepton in the primary vertex (it is nuwro)
       nucl.remove_nucleon(e.in[1]); // remove nucleon used in the primary vertex
      
  bool test=false;
  bool writeall=true;
//  bool writeall=false;
	//double t = 0;
//    if(test)	
		vec t=test_fun(input,q)*fermi;
     
  // take all nucleons and pions from input vector
  // others copy directly to the output vector

  for (int i = 0; i < input.size (); i++)
    {
      particle p;
      
      if (i < input.size()) p = input[i];
      else p = nucl.get_nucleon(input[i].r);
      p.set_new();
      if (nucleon_or_pion (p.pdg))
	  {
		  set_pt(p, q, e.in[1], e.flag.res, e.flag.dis);
		  
		  double fz = formation1(p, par, q, qel, pin, res, input.size()-1, e.W());
		  //if (!qel) fz = formation3((q+e.in[1].p4()).v().length(), input.size()-1);
		  
		  //if (test) 
		  //fz = fermi*t;//*p.momentum()/p.E();
		  
		  /// only for tests 
		  //double pt = help.dir()*p.p();
		  //pt = sqrt(p.momentum2() - pt*pt);
		  //p.wpt(pt);
		  //p.wfz(fz); 
		  //p.endproc=escape;
		  //output.push_back (p);
		  ///must be commented
		  
		  //p.krok(t);
		  p.krok(fz);
		  if (nucleon (p.pdg))
			p.set_fermi(e.in[1].Ek());
		  		  
		  parts.push (p); 
	  }
      else
      {
		p.endproc=escape;
		output.push_back (p);
		if(writeall) 
			all.push_back(p);
      }
    }
    
    ////////////////////////////
    /// initialize statistics
	for (int i = 0; i<12; i++)
	{
		nod[i] = 0;
	}
		
	for (int i = 0; i<20; i++)
	{
		for (int j = 0; j < 11; j++)
		place[j][i] = 0;
	}
	
	nopp = true;
	abs = false;
	pabsen = 0;
	nofi = 0;
	absr = 100;
	odl = 0;
	
	for (int i = 0; i < 10; i++) pen[i] = 0;
 
////////for tests	
	//ofstream test("pion_energy_lost.txt");
    
    //int nofi = 0;
////////    
  /// MAIN LOOP IN CASCADE
  while (parts.size () > 0 and nucl.Ar()>0)
    {
      particle p1 = parts.front();
      parts.pop();
      
      
      	/*if (p1.pdg == 211 or p1.pdg == -211 or p1.pdg == 111)
		 {
			test << p1.r/fermi << endl;
		 }*/
	  
      double r = p1.r.length ();
      
      if (nucleon (p1.pdg))
      {
		  double V = nucl.calculate_V(p1) + par.kaskada_w;
		  
		  if (p1.Ek() < V)
		  {
			p1.endproc=jailed;				
			nucl.insert_nucleon (p1);
			continue;
		}
	  }

      double dens = nucl.density (r);
      assert(dens>=0);
      double d0 = dens * nucl.frac_neutron ();
      double d1 = dens * nucl.frac_proton ();
      
      static Interaction I(par.xsec);

      double s0, s1;// place for neutron and proton cross sections
      
      I.total_cross_sections (p1.pdg, p1.Ek (),  s0, s1, dens); //calculate cross sections s0,s1

      double stotal = d0*s0+d1*s1;
      double freepath=0;
      
      //if(stotal<=0) 
      //  cout<< "d0=" <<d0<< "d1=" <<d1<< "s0=" <<s0<< "s1=" <<s1<<endl;
      assert(stotal>=0);
      
      if(stotal>0)
      {
		freepath = -log (frandom ()) / stotal;
        p1.krok (min (max_step, freepath));	// propagate by no more than max_step
      }
      if (stotal==0 || nucl.radius()<=p1.r.length ()) // leaving nucleus
		{
		  if (nucleon (p1.pdg))	// reduce momentum when leaving nucleus or jailed the nucleon
		  {	
			  double V = p1.his_fermi + par.kaskada_w;
			  		
				p1.set_momentum (p1.p () /p1.momentum() * sqrt(p1.momentum()*p1.momentum() - 2.0*p1.E()*V  + V*V)); //p1.set_momentum (p1.p () * (1 - nucl.V() / p1.E ())); OLD!!!
			}
			
			p1.endproc=escape;
            output.push_back (p1);
				
		  if(writeall) 
		     all.push_back (p1);
		}
      else if (max_step < freepath)	// no interaction during max_step
	    {
	      parts.push (p1);	// put to end of queue for further processing
	      
	    }
	  else	// can_interact(p1) and freepath < max_step   (i.e. there was interaction )
	    {
	      particle p2;	// nucleon from nucleus  ///not initialized ????
	                    // initialized in particle_scattering!!!
	      particle p[5];	// place for the results of scatering
	      int n = 2;	// place for the number of particles after scattering
	      vect js (0, 0, 0, 0);	// check energy-momentum conservation 

	      int loop = 0;
	      static int call=0;
	      static int rep=0;
	      static int procid=0;
	      while(++call && 0 == I.particle_scattering (p1, p2, s1*d1 / stotal, n, nucl, p, dens))
	      {     
				if(loop==0)
					procid=I.process_id();
				++rep;++loop;
 			    cout<<I.process_id()<<'/'<<loop<<','<<endl;
				cout<<int(rep*1000./call)/10.<<" % repeated"<<endl;
				double suma=0;
				cout<<"["<<nucl.p<<","<<nucl.n<<"/"<<nucl.Zr()<<","<<nucl.Nr()<<"]"<<p1.pdg<<' '<<p2<<':'<<endl;
				for(int i=0;i<n;i++)
				  { suma+=p[i].mass();
					cout<<' '<<p[i].pdg;
				  }
				cout<<endl<<p1.mass()+p2.mass()<<' '<<suma<<endl;
//CJ	        cout<<" Interaction: "<<I.process_name()<<" ("<<I.process_id()<<") ";
               assert(procid==I.process_id());
               if(loop>100) break;
	      }
	      if(loop>100) 
			{parts.push(p1);
			 continue;
			}

	      for (int i = 0; i < n; i++)
		  if (!p[i].is_valid ())
		  {
		    cerr << I.process_name()<< "Interaction: error" << p[i] << endl;
                    delete pnucleus;
            exit(0);
		    return 0;
		  }

	      if ( (I.process_id() != 25) // Why absorption can not be Pauli blocked?
	          and 
	          nucl.pauli_blocking (p, n) and (p1.pdg == 2112 or p1.pdg == 2212)// and (par.xsec == 0)
	          )	
			{ // some particle was Pauli blocked
				parts.push (p1);	// Pauli blocked, so nothing happened - put p1 to end of queue
				continue;
			}
	      else
		{// Scattering happened  
		    assert(check(p1,p2,nucl.spectator,n,p,I.process_id()));
		    assert(check2(p1,p2,nucl.spectator,n,p,I.process_id()));
		    
		    
			result = 1;			
			int ktory=2*r/fermi;
			ktory=min(ktory,19);
			
			if (I.process_id() >= 20 and I.process_id() <= 25)
			{	
				vec help=vecprod(p1.r,input[0].r);
				odl = help.length();
			}
			
			//////for tests
			
			//if (I.process_id() == 24) test << "abs" << endl;
			
			if (I.process_id() >= 20 and I.process_id() < 25 and !abs)
			{
				nofi++;
				if (nofi < 10) pen[nofi] = p1.Ek();
			}
			switch(I.process_id())
			{case 12: case 13: case 22: case 23: case 24: nopp=false;break;
		     case 25:
						pabsen = p1.Ek();
						absr = p1.r.length();
						abs = true;
						break;
			}
			//////
			
//CJ		  cout << ": completed." << endl;				  
			  int k =kod(I.process_id());
			   nod[k]++;
			   if(k<11)
				 place[k][ktory]++;
			   
			   static int stat[11]={0,0,0,0,0,0,0,0,0,0,0};
			   stat[k]++;
	//		   cout<<p1.pdg<<": ";
	//		   for(int i=0;i<11;i++)
	//		      cout<<i<<" -> "<<stat[i]<<" ";	 
	//		   cout<<endl;
				
			  
			  p1.endproc=I.process_id();
			  if(writeall) 
			     all.push_back(p1);
			  
			  nucl.remove_nucleon (p2); // remove from the nuclear matter
			  if(nucl.spectator!=NULL)
			      nucl.remove_nucleon (*nucl.spectator);
////		  cout << "Interaction: " << n << "particle added to cascade." << endl;

				  vect qhelp;
				 double czas = 0;

				  for (int l = 0; l < n; l++) qhelp += p[l].p4();
					  qhelp -= p2.p4();
				 
			if (test)
			{
				vec h=vec(qhelp);
			
				for (int l = 0; l < n; l++)
				{
					  double pt = h.dir()*p[l].p();
					  pt = sqrt(p[l].momentum2() - pt*pt);
					  czas += 0.1*p[l].mass()*p[l].E()/(p[l].mass2() + pt*pt);
				}
			}
			
			vec time = test_fun(p,n)*fermi;

		  for (int i = 0; i < n; i++)
		    {
			  p[i].r = p1.r;
		      p[i].travelled = 0;
		      p[i].set_mother(p1);
		      
		      bool qel2 = false;
		      bool res2 = false;
		      
		      if (I.process_id() == 10 or I.process_id() == 11 or I.process_id() == 20 or I.process_id() == 21) qel2 = true;
		      if (I.process_id() == 12 or I.process_id() == 13 or I.process_id() == 22 or I.process_id() == 23 or I.process_id() == 24) res2 = true;

			  //if (I.process_id() != 25)
			  //{
				  double fz = formation2(p[i], par, n);//, qhelp, qel2, p2, res2);
				  //fz = time*fermi;//*p[i].momentum()/p[i].E();
				  //fz = 1.0*fermi;
				  
				  //p[i].krok(time);
				  p[i].krok(fz);
				  
//				  if (!test) p[i].krok(formation2(p[i], par, q, qel2, pin, res2));//, nofpi));
//				  else p[i].krok(fermi*czas*p[i].momentum()/p[i].E());				  
			  //}

		      js = js + p[i];
		      if (nucleon (p[i].pdg))
		      {
				  if (nucleon (p1.pdg))
					p[i].set_fermi((p2.Ek()+p1.his_fermi)/2.0);
				  else
					p[i].set_fermi(p2.Ek());
				}
/*				  if(p[i].Ek () < nucl.calculate_V(p[i]) + par.kaskada_w)	// to low energy to leave nucleus (old, now shouldn't get here)
				{
					p[i].endproc=jailed;
					if(writeall) 
						all.push_back(p[i]);
					nucl.insert_nucleon (p[i]);
				}
			}
		      else*/
			    //output.push_back (p[i]);	// only to test
				
				parts.push (p[i]);	// queue for further processing
		    }
		  procinfo(p1,p2,n,p);
////		  cout << p1 + p2 << "===?===" << js << "  p2=" << p2 << endl;
		}// scattering ended
	    }
    }  // end of the main loop
  // if nucleus has evaporated the part quequ may not be empty   
  while(! parts.empty())
   {   
	   particle p0=parts.front();
       p0.endproc=escape;
	   output.push_back(p0);
  	   if(writeall) 
			all.push_back(p0);
	   parts.pop();
   }
    
    ////for	 test was here

  set_nucl_hist(density_hist, radius_hist, protons_hist, neutrons_hist, nucl);
  e.pr=pnucleus->Zr();
  e.nr=pnucleus->Nr();
  delete pnucleus;
  return result;

};



		
