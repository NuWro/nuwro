#ifndef _event1_h_
#define _event1_h_

#include <iostream>
#include <vector>
#include "TObject.h"
#include "particle.h"
#include "params.h"
#include "dis/LeptonMass.h"


class flags
{
public:
/// primary vertex flags
  bool qel;
  bool res;
  bool dis;
  bool coh;
  bool mec;
  bool nc;
  bool cc;
  bool anty;			///< true if antineutrino 
};

using namespace std;

class event:public TObject
{
public:
  flags flag;
  params par;
    vector < particle > in;    ///< incoming particles
    vector < particle > temp;  ///< temporary particles (daughters of primary vertex)
    vector < particle > out;   ///< outgoing particles
    vector < particle > post;  ///< results of the kaskadaevent
    vector < particle > all;   ///< all particles

    double weight;             ///< cross section
    double norm;               ///< norm of the initial neutrino
    vec r;
    double density;
    int dyn;                   ///< dynamics code (from proctable.h)
	
	int nod[12]; ///number of dynamics: 0 - nucleon elastic, 1 - nucleon ce, 2 - nucleon spp, 3 - nucleon dpp, 4 - pion elastic, 5 - pion ce, 6 - pion spp, 7 - pion dpp, 8 - pion abs, 9 - jailed, 10 - escape
	
	int place[11][20];
	
	double pabsen; //kinetic energy of absorped pion
	bool nopp; //true if there was no pion production
	bool abs; //true if there was absorption
	int nofi; //number of pion interactions before cascade
	double pen[10]; //pion energy after ith scattering
	double absr; //place where pion was absorbed
	double odl; //distance between primary vertex and piN interaction
	double density_hist[50]; //density after cascade
	double radius_hist; //radius after cascade
	int protons_hist; //#protons after cascade
	int neutrons_hist; //#neutrons after cascade
	int pr,nr;
  void clear_fsi()
  {
	  post=vector<particle>();
	  all=vector<particle>();
	  for(int i=0;i<11;i++)
  	    for(int j=0;j<20;j++)
		  place[i][j]=0;
	  for(int i=0;i<12;i++)
         pen[i%10]=nod[i]=0;
  }
  /// initial neutrino
  particle nu ()
  {
    return in[0];
  }
  /// initial nucleon
  particle N0 ()
  {
    return in[1];
  }
  
  /// fourmomentum transfer
  vect q ()
  {
    vect q = in[0] - out[0];
    return q;
  }
  double q0 ()
  {
    return in[0].t - out[0].t;
  }
  double qv ()
  {
    vect q = in[0] - out[0];
    return vec(q).norm();
  }
  
  /// fourmomentum transfer squared
  double q2 ()
  {
    vect q = in[0] - out[0];
    return q * q;
  }
  /// s viariable
  double s ()
  {
    vect q = in[0] + in[1];
    return q * q;
  }
  /// cos theta lab
  double costheta ()
  {
    return cos(vec(in[0]),vec(out[0]));
  }
  /// neutrino energy
  double E ()
  {
    return in[0].t;
  }
  /// number of outgoing particles
  int n ()
  {
    return out.size ();
  }				

int f ()
  {
    return post.size ();
  }				
  

// number of outgoing particles with given pdg
  int nof (int pdg)		
  {
    int c = 0;
    for (int i = 0; i < out.size (); i++)
      if(out[i].pdg == pdg)
		c++;
    return c;
  }
  int nof (int pdg1,int pdg2)		
  {
    int c = 0;
    for (int i = 0; i < out.size (); i++)
      if(out[i].pdg == pdg1 || out[i].pdg == pdg2)
        c++;
    return c;
  }
  int nof (int pdg1,int pdg2,int pdg3)		
  {
    int c = 0;
    for (int i = 0; i < out.size (); i++)
      if (out[i].pdg == pdg1 || out[i].pdg == pdg2 || out[i].pdg == pdg3)
        c++;
    return c;
  }

  int fof (int pdg)		
  {
    int c = 0;
    for (int i = 0; i < post.size (); i++)
      if(post[i].pdg == pdg) 
		c++;
    return c;
  }
  int fof (int pdg1, int pdg2)		
  {
    int c = 0;
    for (int i = 0; i < post.size (); i++)
      if(post[i].pdg==pdg1 || post[i].pdg==pdg2)
        c++;
    return c;
  }

  int fof (int pdg1, int pdg2, int pdg3)		
  {
    int c = 0;
    for (int i = 0; i < post.size (); i++)
      if(post[i].pdg==pdg1 || post[i].pdg==pdg2 || post[i].pdg==pdg3)
        c++;
    return c;
  }
  int charge (int r)		
  { 
	vector<particle>  *set=&post;
    int c = 0;
	switch(r)
	{
		case 0: return par.nucleus_p+in[0].charge();break;
		case 1: set=&out;c=par.nucleus_p-(dyn<6)*in[1].charge();break;
		case 2: set=&post;c=pr;break;
	}
    for (int i = 0; i < set->size (); i++)
		c+=(*set)[i].charge();
    return c;
  }


  double W ()
  {
    vect h = out[1];
    for (int a = 2; a < out.size (); a++)
      h = h + out[a];
    return sqrt (h * h);
  }


  double przod ()
  {
    int licz = 0;
    vect h = out[1];
    for (int a = 2; a < out.size (); a++)
      h = h + out[a];

    vect tran = in[0] - out[0];
    vec ptr = vec (tran.x, tran.y, tran.z);

    for (int b = 1; b < out.size (); b++)
      {
	if ((out[b].pdg != 111) && (out[b].pdg != 22) && (out[b].pdg != 130)
	    && (out[b].pdg != 2112) && (out[b].pdg != -2112))

	  {
	    vect hh = out[b];
	    hh.boost (-h.v());

	    vec phh = vec (hh.x, hh.y, hh.z);

	    if (ptr * phh > 0)
	      licz = licz + 1;
	  }
      }
    return licz;
  }

  double tyl ()
  {
    int licz = 0;
    vect h = out[1];
    for (int a = 2; a < out.size (); a++)
      h = h + out[a];

    vect tran = in[0] - out[0];
    vec ptr = vec (tran.x, tran.y, tran.z);

    for (int b = 1; b < out.size (); b++)
      {
	if ((out[b].pdg != 111) && (out[b].pdg != 22) && (out[b].pdg != 130)
	    && (out[b].pdg != 2112) && (out[b].pdg != -2112))

	  {
	    vect hh = out[b];
	    hh.boost (-(h.v ()));

	    vec phh = vec (hh.x, hh.y, hh.z);

	    if (ptr * phh < 0)
	      licz = licz + 1;
	  }
      }
    return licz;
  }
  
  int number_of_nucleon_elastic ()
  {
  	return nod[0];
  }
  	
  int number_of_nucleon_ce ()
  {
  	return nod[1];
  }
  
  int number_of_nucleon_spp ()
  {
  	return nod[2];
  }
  
  int number_of_nucleon_dpp ()
  {
  	return nod[3];
  }
  
  int number_of_pion_elastic ()
  {
  	return nod[4];
  }
  
  int number_of_pion_ce ()
  {
  	return nod[5];
  }
  
  int number_of_pion_spp ()
  {
  	return nod[6];
  }
  
  int number_of_pion_dpp ()
  {
  	return nod[7];
  }
  
  int number_of_pion_tpp ()
  {
  	return nod[11];
  }
  
  int number_of_pion_abs ()
  {
  	return nod[8];
  }
  
  int number_of_jailed ()
  {
  	return nod[9];
  } 
  
  int number_of_escape ()
  {
  	return nod[10];
  }
  
  int number_of_interactions ()
  {
  	int noi = 0;
  	for (int i = 0; i<11; i++)
  	{
  		noi = noi + nod[i];
  	}
  	
  	return noi;
  }
  
  int number_of_particle (int pdg, bool fsi)	// fsi = 0 - before FSI, 1 - after FSI
  {
  	int number = 0;
  	
  	if(fsi)
  	{
  		for (int k = 0; k<post.size(); k++)
  		{
  			if (post[k].pdg == pdg) number++;
  		}
  	}
  	else if (!fsi)
  	{
  		for (int k = 0; k<out.size(); k++)
  		{
  			if (out[k].pdg == pdg) number++;
  		}
  	}
  	
  	return number;
  	
  }
	
public:
  event ():weight(0),norm(1)
  {
  }
  
  void check()
  {for(int i=0;i<out.size();i++)
     {particle p=out[i];
      if(!(p.x==p.x && p.y==p.y && p.z==p.z and p.t==p.t) )
        {cerr<<p<<endl;
         exit(-1);
	    }
     }
     if(!(weight==weight))
        {cerr<<"dyn="<<dyn<<endl;
         exit(-1);
	    }
      
  }
  ClassDef (event, 1);
};

////////////////////////////////////////////////////////////
/// add to list and get id
inline void registration (vector<particle>& list, particle &p)
{
	p.id=list.size();
	list.push_back(p);
}


inline void refresh (vector<particle>& list, particle &p)
{
	list[p.id]=p;
}


#endif
