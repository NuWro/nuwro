#ifndef _Interaction_h_
#define _Interaction_h_
#include <iostream>
#include <TGenPhaseSpace.h>
#include "pdg.h"
#include "pdg.h"
#include "particle.h"
#include "nucleus.h"
#include "jednostki.h"
#include "scatter.h"
//#include "proctable.h"
#include "seek.h"
//#define echo(x)
#include "piangle.h"

enum {nucleon_=10,pion_=20};
enum {elastic_=0, ce_=1, spp_=2, dpp_=3, tpp_=4, abs_=5};

class interaction_parameters
{
	public:
		double r; //distance from a center
		double dens; //density in r
		double dens_n; //neutron density in r
		double dens_p; //proton density in r
		double xsec; //total cross section
		double xsec_n; //xsec for neutron
		double xsec_p; //xsec for proton
		double frac_proton;
		int pdg; //particle pdg
		double Ek; //particle kinetic energy
		double Ekeff;
		double freepath;
		particle p2;	//nucleon from nucleus (initialized in particle_scattering)
		particle p[5];	//place for the results of scatering
		int n;			//place for the number of particles after scattering
};

inline int kod(int i)
{switch(i)
  {
	case 10: return 0;
	case 11: return 1;
	case 12: return 2;
	case 13: return 3;
	case 20: return 4;
	case 21: return 5;
	case 22: return 6;
	case 23: return 7;
	case 25: return 8;
	case 24: return 11; 
	case 99: return 9;
	case 100: return 10;
	default: return -1;
  }		  	
}
///////////////////////////////////////////)////////////////
/// Help to "reuse" random number "x" a few times in nested ifs.
/// Assume x - random in [0,1) and type:
/// if(below(x,val))
///       ....../// branch 1   //probability=val
///    else
///        ...../// branch 2   //probability=1-val
/// Then x is again random in [0,1) in each branch
////////////////////////////////////////////////////////////
inline bool below(double &x,double val)
{
  if(x<val) {x/=val;    return true; } ///<<renormalize x  
  else      {x/=(1-val);return false;} ///<<renormalize x  
}
///////////////////////////////////////////)////////////////
///one should use "if(above(x,val))" instead of "if(not below(x,val))"
///////////////////////////////////////////)////////////////
inline bool above(double &x,double val)
{
  if(x>=val) {x/=(1-val); return true; } ///<<renormalize x  
  else       {x/=val;     return false;} ///<<renormalize x  
}

///////////////////////////////////////////)////////////////
inline double pow2(double x) {return x*x;}

using namespace std;

static particle PiPlus(PDG::pdg_piP,PDG::mass_piP);
static particle PiMinus(-PDG::pdg_piP,PDG::mass_piP);
static particle PiZero(PDG::pdg_pi,PDG::mass_pi);
static particle Proton(PDG::pdg_proton,PDG::mass_proton);
static particle Neutron(PDG::pdg_neutron,PDG::mass_neutron);

///////////////////////////////////////////)////////////////
struct channel{double dist;const char* codes;};

///////////////////////////////////////////)////////////////
/// choose reaction channel from table a[] 
/// according to probability distribution a[?].dist
/// and fill p with particles according to codes string
inline void doit(int& n,const channel a[],particle p[])
{
  if(a->dist<1)
  {
	  double x=frandom();
	  while(x>=a->dist) 
		 ++a;
  }
  const char *c=a->codes;
  for(int i=0;true;i++)
	switch(c[i])  ///Fill p with particles at rest based on char codes
	{ case '-':p[i]=PiMinus;break;
	  case '.':p[i]=PiZero;break;
	  case '+':p[i]=PiPlus;break;
	  case 'n':p[i]=Neutron;break;
	  case 'p':p[i]=Proton;break;
	  case'\0':n=i;return;
	  default:cerr<<"doit: Invalid process: n="<<n<<" \""<<a->codes<<"\""<<endl;exit(27);
	 }
} 

///////////////////////////////////////////)////////////////
class NData
{
private:
     int k2;
     double Ek;
     int iE;
     double aE; 
     int ij;
	 const double *E,*s[2],*F[2],*A[2],*B[2],*Fp[2];
     double val(const double* v){return (1-aE)*v[iE]+aE*v[iE+1];}
public:
  NData(int i )
	{if(i) setMetropolis();
	 else  setOset();
	 }
  void setMetropolis();
  void setOset();
  void set_Ek(double NewEk)
	{   Ek=NewEk;
		iE=0;
		while(NewEk>=E[iE+1])
			++iE;
		aE=(NewEk-E[iE])/(E[iE+1]-E[iE]);
	}
  void set_particles(particle& p1,particle& p2){ ij= p1!=p2; }
//  void set_density(double dens){}
  inline void get_sij (double Ek,double &resii,double &resij);
  inline double sigma (double Ek);
  double finel(){return F[ij][iE];}    
  double a(){return A[ij][iE];}    
  double b(){return B[ij][iE];}    
  double fp(){return Fp[ij][iE];}
  inline bool nucleon_scattering (particle & p1, particle & p2, 
							     int &n,  particle p[]);
  inline bool nucleon_elastic    (particle & p1, particle & p2, 
                                 int &n,  particle p[]);
  inline bool nucleon_spp 		 (particle p1, particle p2,     
                                  int &n, particle p[]);// p[2] is pion
  inline bool nucleon_dpp         (particle p1, particle p2,
			                       int &n, particle p[]);// p[2], p[3] are pions
  inline int process_id()
	{
		return nucleon_+k2;
	}
  inline const char * process_name()
	{const char *name[4]=
	     {"nucleon elastic","nucleon ce","nucleon spp","nucleon dpp"};
	 if(k2<4) return name[k2];
     else     return NULL;
	}
};
///////////////////////////////////////////////////////////
/// Cross section for same (ii) and different (ij) nucleons
void NData::get_sij (double Ek_,double &resii,double &resij)
{   set_Ek(Ek_);
	if (Ek_ < 335 * MeV)
	{   Ek_=max(Ek, 30 * MeV); //makes virtually no difference to set 25 MeV here
		const double M = (mass_proton + mass_neutron) / 2;
		double v = sqrt (1 - pow2 (M / (Ek+M)));
		resij=((34.10 / v - 82.20) / v + 82.2)* millibarn;		
		resii=((10.63 / v - 29.92) / v + 42.9)* millibarn;
	}
	else
    { 
	  resii=max(0.0,val(s[0])) * millibarn;
	  resij=max(0.0,val(s[1])) * millibarn;
     }  
}
///////////////////////////////////////////////////////////
double NData::sigma (double Ek_)
{   set_Ek(Ek_);
	if (Ek_ < 335 * MeV)
	{ 
  	  Ek=max(Ek, 40 * MeV); 
	  const double M = (mass_proton + mass_neutron) / 2;
  	  double v = sqrt (1 - pow2 (M / (Ek+M)));
	  if(ij) return ((10.63 / v - 29.92) / v + 42.9)* millibarn;
	  else    return ((34.10 / v - 82.20) / v + 82.2)* millibarn;
	}
	else
	  { 
        return val(s[ij]) * millibarn;    
      }
}
///////////////////////////////////////////////////////////
bool NData::nucleon_scattering (particle & p1, particle & p2, 
			                   int &n, particle p[])
{
    set_particles(p1,p2);
    vec v = p2.v();
    assert ( v*v<1 && " nucleon  ");
    
    double Ek1 = p1.Ek_in_frame (v);
    double Ek2 = p1.Ek_in_frame (-v);
    double s2=sigma(Ek2); 
    double s1=sigma(Ek1); 
    double s = s1;

	if (frandom()*(s1 + s2) < s2)
	  {
		p2.x *= -1;
		p2.y *= -1;
		p2.z *= -1;
		set_Ek(Ek2);
	  }
	  
	if (frandom()>finel())
	    return nucleon_elastic(p1, p2, n, p);
	if (frandom() <  fp() ) 
		return nucleon_spp(p1, p2, n, p) 
		     ||nucleon_elastic(p1, p2, n, p); // in case of insufficiend energy fallback to elastic
	else 		    
	    return nucleon_dpp(p1, p2, n, p) 
	         ||nucleon_spp(p1, p2, n, p)      // in case of insufficiend energy fallback to spp
	         ||nucleon_elastic(p1, p2, n, p); // or even to to elastic
}
///////////////////////////////////////////////////////////
bool NData::nucleon_elastic(particle& p1,particle& p2, int &n, particle p[])
 { 
	k2=elastic_;
	n = 2;
	p[0] = p1;
	p[1] = p2;
	int res=scatterAB (p1, p2, p[0], p[1], 0, 0, 0, a(), b(), 0, 0, 1); //0*x^7 + 0*x^6 + 0*x^5 + a*x^4 + b8x^3 + 0*x^2 + 0*x + 1
		
	if(res==0) cerr<<"AB=0"<<endl;
	return res;
 }

///////////////////////////////////////////////////////////

/// nucleon single pion production
bool NData::nucleon_spp (particle p1, particle p2, int &n, particle p[])
  {
	k2=spp_;
    int canal=(p1==Proton)*2+(p2==Proton);
    static const double f1[]={0.11};
    static const double f2[]={0.43,0.815};
    static const channel cnls[4][3]=
      {{f1[0],"nn.",    1,"np-", 1,"ccc"},//nn
       {f2[0],"np.",f2[1],"nn+", 1,"pp-"},//np
       {f2[0],"pn.",f2[1],"pp-", 1,"nn+"},//pn
       {f1[0],"pp.",    1,"np+", 1,"ddd"} //pp
      };
    doit(n,cnls[canal],p); 	// p[2] is pion    
    return scatter_n (n, p1, p2, p);
  }
///////////////////////////////////////////////////////////  
/// nucleon double pion production
bool NData::nucleon_dpp (particle p1, particle p2, int &n, particle p[])	
  {
 	  k2=dpp_;
      int canal=(p1==Proton)*2+(p2==Proton);
      static const double f1[]={0.6,0.8};
      static const double f2[]={0.6,0.8,0.9};
      static const channel cnls[4][4]=
      {{f1[0],"nn..",f1[1],"nn+-",    1,"np.-",1,"    "},//nn
       {f2[0],"np..",f2[1],"np+-",f2[2],"pp.-",1,"nn.+"},//np
       {f2[0],"pn..",f2[1],"pn-+",f2[2],"nn.+",1,"pp.-"},//pn
       {f1[0],"pp..",f1[1],"pp-+",    1,"pn.+",1,"    "} //pp
      };
       doit(n,cnls[canal],p);// p[2],p[3] are pions
       return scatter_n (n, p1, p2, p);
  }


///////////////////////////////////////////////////////////
class PiData
{ 
private:
    int k2;
	double Ek;
	int ij;
	int nE;
	int iE;
	double aE;
    int nD;
    int iD;
    double aD;
    double maxdens;
	int xsec;
const double *E, *s[3], *F[3], *Fc[3],
             *A[3], *B[3], *Cel[3], *Cinel[3], *Fp[3], *F2p;
private:
	inline double dval(const double *V)
	{		
		if(nD==1)
			return V[iE];
		int j=iE*nD+iD;
//		int j=iE*nD+iD-1;//?? Sprawdzić
		return (1-aD)*V[j]+aD*V[j+1];
	}
	inline double dval2(const double *V)
	{ 
  	    if(nD==1)
			return (1-aE)*V[iE]+aE*V[iE+1];
		int j=iE*nD+iD; 
//		int j=iE*nD+iD-1; //?? Sprawdzić
		return ((1-aD)*V[j]+aD*V[j+1])*(1-aE) // V[iE][iD] , V[iE][iD+1]
		      +((1-aD)*V[j+nD]+aD*V[j+nD+1])*aE; // V[iE+1][iD] , V[iE+1][iD+1]
	}
	inline double angle_par(bool cex, int z)
	{
		int start = 0; //start and end points in piAngle table
		int end = 69;

		if (cex)
		{
			start = 210;
			end = 279;
		}
		else if (ij == 1)
		{
			start = 70;
			end = 139;
		}
		else if (ij == 2)
		{
			start = 140;
			end = 209;
		}
		
		int bin = 0;
		
		if (Ek <= piAngle[start][0]) return piAngle[start][z];
		if (Ek >= piAngle[end][0]) return piAngle[end][z];
		
		for (int i = start; start < end; i++)
		{
			if (Ek > piAngle[i][0] and Ek <= piAngle[i+1][0])
			{
				bin = i;				
				break;
			}
		}
		
		double E2 = piAngle[bin+1][0] - Ek;
		double E1 = Ek - piAngle[bin][0];
		
		return (E1*piAngle[bin+1][z] + E2*piAngle[bin][z])/(E1+E2);
	}
	void dump(const double *V)
	{if(nD>1) cout<<endl;
		  for (int i=0;i<nE;i++)
	     {if(nD>1) cout<<setw(8)<<E[i]<<": ";
		  for(int j=0;j<max(nD,1);j++)
	        cout<<setw(10)<<V[i*max(nD,1)+j]<<' ';
	       if (nD>1) cout<<endl; 
	      }
       cout<<endl; 
	}
	void dump()
	{
#define echo(x) cout<<#x<<":"; dump(x);
			  echo(s[0]);echo(s[1]);echo(s[2]);
//              echo(F[0]);echo(F[1]);echo(F[2]);
//              echo(Fc[0]);echo(Fc[1]);echo(Fc[2]);
//              echo(A[0]);echo(A[1]);echo(A[2]);
//              echo(B[0]);echo(B[1]);echo(B[2]);
//              echo(Cel[0]);echo(Cel[1]);echo(Cel[2]);
//              echo(Cinel[0]);echo(Cinel[1]);echo(Cinel[2]);
//              echo(Fp[0]);echo(Fp[1]);echo(Fp[2]); 
//              echo(F2p);
     }
public:
 
  PiData(int k);
  void setMetropolis();
  void set(double tab[11][16][26],double tab2[6][29]);
  double set_Ek(double NewEk)
	{   Ek=NewEk;
		iE=0;
		while(NewEk>=E[iE+1])
		  ++iE;
		aE=(NewEk-E[iE])/(E[iE+1]-E[iE]);
	}
  double set_density(double dens)
	  {   if(nD==1)
			 {aD=0; iD=0;}
		  else
		  {   
			double x=min(1.0,dens/maxdens)*(nD-1);
			iD=min(int(x),nD-2);
			aD=x-iD;
		   }
	//	  cout<<"iD="<<iD<<" nD="<<nD<< " frac="<<dens/maxdens<<endl;
	  }
  void set_particles(particle &p1,particle& p2)
	{
		if(p1==PiZero) ij=2; 
		else if(p1==PiPlus && p2==Proton || 
	            p1==PiMinus && p2==Neutron) ij=0;
		else ij= 1;
    }
    
    double cel(){return dval(Cel[ij]);}
	double cinel(){return dval(Cinel[ij]);}


	double a1(bool cex)
	{	
		//return 0;	
		return angle_par(cex, 1); //cex = false for elastic and true for CEX
	} 
	double a2(bool cex)
	{
		//return 0;
		return angle_par(cex, 2);
	} 
	double a3(bool cex)
	{
		//return 0;
		return angle_par(cex, 3);
	}
	double a4(bool cex)
	{		  
		//if (cex and Ek < 51) return 1.5;
		//if (cex) return dval(A[0]);
		//return dval(A[ij]);
		return angle_par(cex, 4);
	}
	double a5(bool cex)
	{
		//if (cex and Ek < 51) return -2.5;
		//if (cex) return dval(B[0]);
		//return dval(B[ij]);
		return angle_par(cex, 5);
	}
	double a6(bool cex)
	{
		//if (cex and Ek < 51) return 0;
		//if (cex) return cinel();
		//return cel();
		return angle_par(cex, 6);
	}
	double a7(bool cex)
	{
		//return 0;
		return angle_par(cex, 7);
	}
	double a8(bool cex)
	{
		//return 1;
		return angle_par(cex, 8);
	}
		
  double finel(){return dval(F[ij]);}    //dval2 or dval - skalowanie z energią lub bez
  double fce(){return dval(Fc[ij]);}    
  double fp(){return dval(Fp[ij]);}
  double f2p(){return dval(F2p);}
  double fabs_inel() // absorption/inelastic
	  {//cout<<"ij="<<ij<<" i="<<i<<" iD="<<iD<<"["<<s[ij][i*nD+iD]<<"]"<< " {"<<sij(0)<<','<<sij(1)<<','<<sij(2)<<")"<<endl;
		  switch(ij) 
			{case 0: return 0;
			 case 1: {double x=sij(2);return x/(x+sij(1));}
			 case 2: {double x=sij(2);return x/(x+sij(0)+sij(1));
				     }
			}
	   } 

	double sigma ()
     {switch(ij)
      {
        case 0: return sij(0)*millibarn;
        case 1: return (sij(1)+sij(2))*millibarn;
        case 2: return (sij(0)+sij(1)+sij(2))/2.0*millibarn;
       } 
      }

  int process_id() { return pion_+k2; }
  const char* process_name()
	  {const char* name[6]={"pion elastic","pion ce",
	                     "pion spp","pion dpp",
	                     "pion tpp","pion abs"};
	   if(k2<6) return name[k2];
	   else return NULL;                   
	  }
	inline double sij(int i);
	inline bool pion_scattering (particle & p1, particle & p2, nucleus &t,
						   int &n, particle p[], double dens); //p1-moving pion ,p2-target nucleon
	inline bool pion_abs (particle& p1, particle& p2, nucleus & t, int &n, particle p[]); 
	inline bool pion_elastic (particle& p1, particle& p2,  int &n, particle p[]);
	inline bool pion_ce (particle& p1, particle& p2,  int &n, particle p[]);
	inline bool pion_spp (particle& p1, particle& p2, int &n, particle p[]);    
	inline bool pion_dpp (particle& p1, particle& p2, int &n, particle p[]);	
	inline bool pion_tpp (particle& p1, particle& p2, int &n, particle p[]);    
};
///////////////////////////////////////////////////////////

  /// Cross section for pion - nucleon scattering   
  /// k = 0  // ii
  /// k = 1  // ij
  /// k = 2  // abs
  /// Ek = kinetic energy 
double PiData::sij (int k)
  { 
    if (Ek <= 49*MeV && xsec==0)
      {	// Low energy Metropolis formula
        double x = Ek / PDG::mass_pi;
        switch (k)
		{
		  case 0:  return (3.7 + 286 * x * x * x);
		  case 1:  return (6.5 + 23.9 * x);
		  case 2:  double p = sqrt(x*(x + 2));
					return (16.4 * (0.14 / p + p));
				   
		}
	}  
	else
		return max(0.0,dval2(s[k]));
}

///////////////////////////////////////////////////////////
bool PiData::pion_scattering (particle & p1, particle & p2, nucleus &t,
		       int &n, particle p[], double dens)
  { 
    set_particles(p1,p2);
    set_density(dens);
    int canal=((p1==PiPlus)-(p1==PiMinus))*2+(p2==Proton)+2;
    
    vec  v = p2.v();
    assert ( v*v<1 && " pion  ");
    double Ekm = p1.Ek_in_frame (-v);
    set_Ek(Ekm);
    double sm = sigma ();
    double Ekp = p1.Ek_in_frame ( v);
    set_Ek(Ekp);
    double sp = sigma ();
    if (frandom() * (sm + sp) < sm)
      {
		p2.x *= -1;
		p2.y *= -1;
		p2.z *= -1;
		set_Ek(Ekm);
      }
    if(frandom()<fabs_inel() and canal != 0 and canal != 5)
       return pion_abs (p1, p2, t, n, p);
    if(frandom()>=finel())
       return pion_elastic (p1, p2, n, p); 
    if(canal != 0 && canal != 5 && frandom()<fce())
       return pion_ce  (p1, p2, n, p)
            ||pion_elastic(p1, p2, n, p);
    if(frandom()<fp())
	     return pion_spp (p1, p2, n, p)
	        ||pion_ce(p1, p2, n, p)
            ||pion_elastic(p1, p2, n, p);
    if (frandom() < f2p()) 	
          return pion_dpp (p1, p2, n, p)
               ||pion_spp(p1, p2, n, p)
               ||pion_ce(p1, p2, n, p)
               ||pion_elastic(p1, p2, n, p);	
	else
		return pion_tpp (p1, p2, n, p)
		     ||pion_dpp (p1, p2, n, p)
		     ||pion_spp (p1, p2, n, p)
             ||pion_ce(p1, p2, n, p)	
             ||pion_elastic(p1, p2, n, p);	
			
  }
///////////////////////////////////////////////////////////
bool PiData::pion_elastic (particle& p1, particle& p2,  int &n, particle p[])
  {  
	n = 2;  
	k2=elastic_;
	p[0]=p1;
	p[1]=p2;
	return scatterAB (p1, p2, p[0], p[1], a1(0), a2(0), a3(0), a4(0), a5(0), a6(0), a7(0), a8(0));
   }
///////////////////////////////////////////////////////////
/// pion charge exchange
bool PiData::pion_ce (particle& p1, particle& p2, int &n, particle p[])
  {  
 	   k2=ce_;
	   
	   int canal=((p1==PiPlus)-(p1==PiMinus))*2+(p2==Proton)+2;
	   if(canal==0 || canal==5) return 0;
		
	   static const channel cnls[6][2]=
	   {{1,"ee"},//-n
		{1,".n"},//-p
		{1,"-p"},//.n
		{1,"+n"},//.p
		{1,".p"},//+n
		{1,"ff"} //+p
	   };
	   doit(n,cnls[canal],p);
			return scatterAB (p1, p2, p[0], p[1], a1(1), a2(1), a3(1), a4(1), a5(1), a6(1), a7(1), a8(1)); 
  }

///////////////////////////////////////////////////////////
bool PiData::pion_spp (particle& p1, particle& p2, int &n, particle p[])
  {
	k2=spp_;
	int canal=((p1==PiPlus)-(p1==PiMinus))*2+(p2==Proton)+2;
    static const double f1c[]={0.75};      // pi charged
    static const double f2o[]={1./3,2./3}; // pi0
    static const double f2c[]={0.65,0.90}; // pi charged

	static const channel cnls[6][3]=
	{ {f1c[0],"-n.",      1,"-p-", 1,"aaa"},//-n
      {f2c[0],"+n-", f2c[1],"-p.", 1,".n."},//-p
      {f2o[0],"+n-", f2o[1],"-p.", 1,".n."},//.n
	  {f2o[0],"-p+", f2o[1],"+n.", 1,".p."},//.p
	  {f2c[0],"-p+", f2c[1],"+n.", 1,".p."},//+n
	  {f1c[0],"+p.",      1,"+n+", 1,"bbb"} //+p
	};
	doit(n,cnls[canal],p);
    return scatter_n (n, p1, p2, p);
  }
  
///////////////////////////////////////////////////////////
/// pion double pion production
bool PiData::pion_dpp (particle& p1, particle& p2, int &n, particle p[])	// p[2], p[3] are pions
{
	k2=dpp_;
    int canal = ((p1==PiPlus)-(p1==PiMinus))*2+(p2==Proton)+2;    

//    static const double f1o[2] = {1.0/3.0, 2.0/3.0}; // pi0
    static const double f1[2] = {0.5,0.75};         // pi charged
    static const double f2[3] = {0.25, 0.50, 0.75};   

	static const channel cnls[6][4]=
	{{f1[0],"-n-+",f1[1],"-n..",    1,"-p-.",1,"    "},//-n
	 {f2[0],".n..",f2[1],".n-+",f2[2],"-p..",1,"-p-+"},//-p
	 {f2[0],".n..",f2[1],".n-+",f2[2],"-p..",1,"-p-+"},//.n
	 {f2[0],"+n..",f2[1],"+n+-",f2[2],"+p-.",1,".p.."},//.p
	 {f2[0],"+n..",f2[1],"+n+-",f2[2],"+p-.",1,".p.."},//+n
	 {f1[0],"+p+-",f1[1],"+p..",    1,"+n+.",1,"    "} //+p
	};	
    doit(n,cnls[canal],p);
    return scatter_n (n, p1, p2, p);
}
///////////////////////////////////////////////////////////
/// pion triple pion production
bool PiData::pion_tpp (particle& p1, particle& p2, int &n, particle p[])
{
	k2=tpp_;
    int canal = ((p1==PiPlus)-(p1==PiMinus))*2+(p2==Proton)+2;
        
    static const double f1[3] = {0.25, 0.5, 0.75};
    static const double f2[4] = {0.2, 0.4, 0.6, 0.8};
   
	static  const channel cnls[6][5]=
	{{f1[0],"-n...",f1[1],"-n-+.",f1[2],"-p-..",    1,"-p-+-",0,"     "},//-n
	 {f2[0],".n...",f2[1],".n+-.",f2[2],"+n-+-",f2[3],"-p...",1,"-p+-."},//-p
	 {f2[0],".n...",f2[1],".n+-.",f2[2],"+n-+-",f2[3],"-p...",1,"-p+-."},//.n
	 {f2[0],"+n...",f2[1],"+n+-.",f2[2],".p...",f2[3],".p+-.",1,"+p-+-"},//.p
	 {f2[0],"+n...",f2[1],"+n+-.",f2[2],".p...",f2[3],".p+-.",1,"+p-+-"},//+n
	 {f1[0],"+n+..",f1[1],"+n++-",f1[2],"+p...",    1,"+p+-.",0,"     "} //+p
	};
    doit(n,cnls[canal],p);
    return scatter_n (n, p1, p2, p);
 }
///////////////////////////////////////////////////////////
/// pion absorption on a pair of nucleons
bool PiData::pion_abs (particle& p1, particle& p2, nucleus & t, int & n, particle p[])
  { 
	static particle p2a, p3a;
	if(t.Ar()<2) return 0;
    n=2;  
	k2=abs_;
    
    //int canal = ((p1==PiPlus)-(p1==PiMinus))*2+(p2==Proton)+2;
    //switch(canal)  
    //{ case 0: case 5: return 0; // no absorbtion for (pi- n) and (pi+ p) 
    //}
    /* old version
	if(!t.remove_nucleon(p2))
		return false;       // pretend that p2 is already removed
    p2a = t.get_nucleon (p1.r);	// get nucleon from this place in Nucleus
    t.spectator=&p2a;            // remember to remove p2a later (if PB permits)
    t.insert_nucleon(p2);       // stop pretending :)
    p[0] = p2;
    p[1] = p2a;
    
    if (p1.pdg == pdg_piP && p2.pdg == pdg_neutron)
      p[0].set_proton ();
    else if (p1.pdg == -pdg_piP && p2.pdg == pdg_proton)
      p[0].set_neutron ();
    end old version
    */ 
    p2a = t.get_nucleon (p1.r);	// get nucleon from this place in Nucleus
    p2a.pdg = pdg_proton;
    p3a.pdg = pdg_neutron;
    p[0] = p2;
    p[1] = p2a;
    
    if (p1.pdg == pdg_piP)
    {
    if ( frandom ()<7.0/8.0 ) 
    { p[0].set_proton (); 
      p[1].set_proton (); 
      if (p2.pdg==pdg_proton)
      t.remove_nucleon(p3a);
      if (p2.pdg==pdg_neutron)
      t.remove_nucleon(p2a);
    }
      else 
      { p[0].set_neutron (); 
        p[1].set_proton (); 
	if (p2.pdg==pdg_proton)
	{t.remove_nucleon(p3a);
	t.remove_nucleon(p3a);
	t.insert_nucleon(p2a);
	}
        if (p2.pdg==pdg_neutron)
        t.remove_nucleon(p3a);
       }
      }
      
      if (p1.pdg == -pdg_piP)
    {
    if ( frandom ()<7.0/8.0 ) 
    { p[0].set_neutron (); 
      p[1].set_neutron (); 
      if (p2.pdg==pdg_neutron)
      t.remove_nucleon(p2a);
      if (p2.pdg==pdg_proton)
      t.remove_nucleon(p3a);
    }
      else 
      { p[0].set_neutron (); 
        p[1].set_proton (); 
	if (p2.pdg==pdg_neutron)
	{t.remove_nucleon(p2a);
	t.remove_nucleon(p2a);
	t.insert_nucleon(p3a);
	}
        if (p2.pdg==pdg_proton)
        t.remove_nucleon(p2a);
       }
      }
      
      if (p1.pdg == 111)//p2a=proton p3a=neutron
    {
    if ( frandom ()<4.0/5.0 ) 
    { p[0].set_neutron (); 
      p[1].set_proton (); 
      if (p2.pdg==pdg_neutron)
      t.remove_nucleon(p2a);
      if (p2.pdg==pdg_proton)
      t.remove_nucleon(p3a);
    }
      else 
      { 
	if ( frandom ()<0.5 ) 
	{
	  p[0].set_neutron (); 
          p[1].set_neutron ();
	  if (p2.pdg==pdg_neutron)
      t.remove_nucleon(p3a);
      if (p2.pdg==pdg_proton)
      {t.remove_nucleon(p2a);
      t.remove_nucleon(p2a);
      t.insert_nucleon(p3a);
      }
	}
	else
	{
	  p[0].set_proton (); 
          p[1].set_proton ();
	  if (p2.pdg==pdg_proton)
      t.remove_nucleon(p2a);
      if (p2.pdg==pdg_neutron)
      {t.remove_nucleon(p3a);
      t.remove_nucleon(p3a);
      t.insert_nucleon(p2a);
      }
	}
       }
      }
      /*
    if (p1.pdg == - pdg_piP)
    {
    if ( frandom ()<7.0/8.0 ) { p[0].set_neutron (); p[1].set_neutron (); }
      else { p[0].set_proton (); p[1].set_neutron (); }
    }
    
    if (p1.pdg == 111)
    {
    if ( frandom ()<0.8 ) { p[0].set_neutron (); p[1].set_proton (); }
      else { 
	
	p[0].set_proton (); p[1].set_neutron (); }
    }
    */
    
    
    //cout<<"abs"<<"  "<<endl;
    return ::decay (p1 + p2 + p2a, p[0], p[1]);
  }


///////////////////////////////////////////////////////////
/// Interaction
///////////////////////////////////////////////////////////
class Interaction
{
	int k1;
	PiData PD;
	NData ND;
  public: 
    Interaction(int xsec):PD(xsec),ND(xsec){}
  	inline void total_cross_sections(particle&, nucleus&, interaction_parameters &X);
	inline bool particle_scattering (particle & p1, nucleus &t, interaction_parameters &X);
    inline int test ();
	int process_id()
		{ 
			return k1==nucleon_ ? ND.process_id() : PD.process_id(); 
		}
	const char* process_name()
		{ 
			return k1==nucleon_ ? ND.process_name() : PD.process_name(); 
		}
};
///////////////////////////////////////////////////////////
/// scattering of pion or nucleon on nucleus
/// returns: s0 - cross setion on neutron
///          s1 - cross setion on proton
void Interaction::total_cross_sections(particle &p1, nucleus &t, interaction_parameters &X)
  {
    
     X.p2 = t.get_nucleon (p1.r);// added - JTS
	vec vvv = X.p2.v();
	p1.p4().boost2 (vvv);
	X.Ekeff = p1.Ek();
	
	p1.p4().boost2 (-vvv);
	
	double dens00 = 0.16/fermi3;
        double Masssa = (938.272+939.56533)/2.0;
	
	double beta = -116.0*X.dens/dens00;
	double lambda = (3.29 - 0.373*X.dens/dens00)/fermi;
	double effmass1 = Masssa/
	(1 - 2.0*Masssa*beta/lambda/lambda/( 1+p1.momentum2()/lambda/lambda )/( 1+p1.momentum2()/lambda/lambda ) );
	double effmass2 = Masssa/
	(1 - 2.0*Masssa*beta/lambda/lambda/( 1+X.p2.momentum2()/lambda/lambda )/( 1+X.p2.momentum2()/lambda/lambda ) );
	double effmass3 = Masssa/
	(1 - 2.0*Masssa*beta/lambda/lambda
	/( 1+ (p1.momentum2()+ X.p2.momentum2())/2.0/lambda/lambda )
	/( 1+ (p1.momentum2()+ X.p2.momentum2())/2.0/lambda/lambda ) );
		
	double k1minusk2 = ( p1.p()-X.p2.p() ).length()/Masssa;
	double k1minusk2star = ( 1.0/effmass1*p1.p() - 1.0/effmass2*X.p2.p() ).length(); 
	
	double mod = k1minusk2/k1minusk2star*effmass3/Masssa;//in medium modification of the cross section - JTS
	
	//cout<<"mod_proton  "<<X.dens<<"  "<<p1.momentum()<<"  "<<effmass1<<"  "<<mod<<endl;
	
	double resc=1;
	
    switch (X.pdg)
      {
      case pdg_neutron: ND.get_sij (X.Ek,X.xsec_n,X.xsec_p);
      
      //cout<<"neutron  "<<X.r<<"  "<<X.xsec_n<<"  "<<X.xsec_p<<"  "<<mod<<endl;
      X.xsec_n*=mod;
      X.xsec_p*=mod;
      X.xsec_n*=resc;
      X.xsec_p*=resc;
      if (X.Ek<40)
	X.xsec_p*=0.9;//effective Pauli blocking 
      //cout<<"neutron2  "<<X.r<<"  "<<X.xsec_n<<"  "<<X.xsec_p<<"  "<<mod<<endl;
      break;
	  case pdg_proton:	ND.get_sij (X.Ek,X.xsec_p,X.xsec_n);
	  
	  //cout<<"proton  "<<X.r<<"  "<<X.xsec_n<<"  "<<X.xsec_p<<"  "<<mod<<endl;
	  X.xsec_n*=mod;
         X.xsec_p*=mod;
	 X.xsec_n*=resc;
         X.xsec_p*=resc;
	 if (X.Ek<40)
	X.xsec_n*=0.9;//effective Pauli blocking 
	  //cout<<"proton2  "<<X.r<<"  "<<X.xsec_n<<"  "<<X.xsec_p<<"  "<<mod<<endl;
	  break;
	  default:
	     { 
	      PD.set_density(X.dens);
	      PD.set_Ek(X.Ek);
	      double rii=PD.sij(0)*millibarn;
	      double rij=PD.sij(1)*millibarn;
	      double rabs=PD.sij(2)*millibarn;
	      switch(X.pdg)
	      {
			  case pdg_pi:  
				X.xsec_n =(rii+rij+rabs)/2; //.n 
				X.xsec_p = X.xsec_n;              //.p
				break;
			  case pdg_piP: 
				X.xsec_n =  rij+rabs;//+n
				X.xsec_p =  rii;     //+p
				break;
			  case -pdg_piP:
				X.xsec_n =  rii;     //-n
				X.xsec_p =  rij+rabs;//-p
				break;
			  default:
				X.xsec_n=X.xsec_p=0;
		   };
          };
      }
  }

///////////////////////////////////////////////////////////
/// pion or nucleon scattering on nucleon
/// fracproton - probability that the scattering will be on proton 
/// scattering results are writen to table p 
/// and n is set to number of outgoing particles
///////////////////////////////////////////////////////////
bool Interaction::particle_scattering (particle & p1, nucleus &t, interaction_parameters &X)
  {
	t.spectator=NULL;
    X.p2 = t.get_nucleon (p1.r);
	//X.p2.set_momentum(vec(0,0,0));

    if(X.p2.v().length()>=1) 
         cout<<"t.n="<<t.n<<"  t.p="<<t.p<<"   "<<X.p2<<endl;
	assert(X.p2.v().length()<1 && "particle scattering");

    if (frandom () < X.frac_proton) X.p2.set_proton ();
    else                    	  X.p2.set_neutron ();

    switch (p1.pdg)
    {
      case pdg_proton: 
      case pdg_neutron:
						k1=nucleon_;	
						return ND.nucleon_scattering (p1, X.p2 , X.n,  X.p);
      case pdg_pi: 
      case pdg_piP: 
      case -pdg_piP:
					    k1=pion_;
					    return PD.pion_scattering (p1, X.p2, t, X.n, X.p, X.dens);
      default:	        
						return 0;
    }
  }
///////////////////////////////////////////////////////////
/// test functions present in this module
int Interaction::test ()
  {
//    echo (Interaction::test);
    init_genrand (time (NULL));
    particle proton (pdg_proton, mass_proton);
    particle neutron (pdg_pi, mass_pi);
    particle neutron2 (neutron);
    particle proton2 (proton);
    proton.set_momentum (vec (1000, 0, 0));
    neutron.set_momentum (vec (-1000, 0, 0));
    int t[1000];
    for (int i = 0; i < 1000; i++)
      {
	t[i] = 0;
      }
    for (int i = 0; i < 1000000; i++)
      {
	scatterAB (proton, neutron, proton2, neutron2, 100, 0, 0, 0, 0, 0, 0, 0);
	t[int (proton2.y + 1000) / 2]++;
      }
    for (int i = 0; i < 1000; i++)
      {
	cout << i << '\t' << t[i] << endl;
      }
  }

#undef echo

#endif
