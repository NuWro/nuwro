#include "kinsolver.h"
#include <cstdlib>

double V(double,double);
extern double kF_P;

///////////////////////////////////////////////////////////////////////////////
kinsolver::kinsolver(double E_lab, vec cms_p_dir, vec cmsspeed,
 	             particle nu, particle N0, particle &lepton, particle &N1):	      
///////////////////////////////////////////////////////////////////////////////
	      _E_lab(E_lab), _cms_p_dir(cms_p_dir), _cmsspeed(cmsspeed),
	      _nu(nu), _N0(N0), _lepton(lepton), _N1(N1)
	    {
	    }
	      
///////////////////////////////////////////////////////////////////////////////
double kinsolver::dE(double pnew)
///////////////////////////////////////////////////////////////////////////////
 {  vec p=_cms_p_dir*pnew;    
    _lepton.set_momentum(p);
    _lepton.boost(_cmsspeed);
    _N1.set_momentum(_N0.p()+_nu.p()
                     -_lepton.p());
#ifdef DEBUG_dE
    cout<<_nu<<_N0<<_lepton<<_N1<<endl;
    cout<<_nu+_N0-_N1-_lepton<<"dE = "<<_lepton.energy()+_N1.energy()+V(_N1.momentum())-_E_lab <<endl;
#endif    
    return _lepton.energy()+_N1.energy()+V(_N1.momentum(),kF_P)-_E_lab;
  }
  


/////////////////////////////////////////////////////////////////////////////
double kinsolver::findmomentum()
/////////////////////////////////////////////////////////////////////////////  
  {
   double pa=1e-12*MeV,pb=_E_lab+1*GeV,pc,pold=-1e35;
   double fa,fb,fc=0;
   pc=pa;  
   fa=dE(pa);
   fb=dE(pb);
   if(fa>0) 
      {//cerr<< "Energy below threshold:("<<fa<<","<<fb<<")"<<endl;
//!!!!!!!!!!1       cerr<<'*';    
/*       for(double p =pa*1.01; p<pb;p+=(pb-pa)/2000)
         {double f=dE(p);
	  if(f<fa) 
           cout<<p<<'\t'<<f<<"<"<<fa<<endl;
	  if(f<0)
           {cout<<p<<'\t'<<f<<endl;
	   
            exit(1);
	    }
	  } 
*/	    
       return -1;
      } 
   if(fb<0) 
      {cerr<< "Upper limit incorrect:("<<fa<<","<<fb<<")"<<endl;
       exit(1);
       return -1;
      }
//   cout<<"Solution exists"<<endl;      
   while(pc!=pold)
   { assert(fa*fb<0);
    // cout<<"fa="<<fa<<" fb="<<fb<<endl;
     pold=pc;
     pc=(pa+pb)/2;
     fc=dE(pc);
//     cout<<"dE="<<fc/MeV<<" MeV. dE("<<pa<<")="<<dE(pa)<<"dE("<<pb<<")="<<dE(pb)<<endl;
     if (fc==0) 
        return pc;
     if((fa>0)==(fc>0))       
	  { pa=pc; fa=fc;}
       else     
	  { pb=pc; fb=fc;}
   }
   if(fc/MeV<0.00000001)
     {// cout<<"stety"<<endl;
      return pc;
      }
   else 
      {cout<< "Niestety fc="<<fc/MeV<<"MeV" ;
       throw "Ola";
       return -1;      
      } 
  }
