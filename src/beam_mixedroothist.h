#ifndef _beam_mixedroothist_h_
#define _beam_mixedroothist_h_
#include "beam_uniform.h"

/// The beam shoots with identical particles (whose energy is defined by some 
/// energy profile but the direction of momentum is identical)
/// beam: this class provides the input of the initial beam particle 
/// into the generator.The beam particle is defined by providing 
/// its energy, direction and its PDG code to the constructor
/// of the beam class 

class beam_mixedroothist : public beam
{
protected:

  int throws;
   int n;
   vec dir;
   beam_singleroothist* beams[10];
   double w1[10];
   double w2[10];
   
public:   
   beam_mixedroothist(NSNWRO::params& p):n(0),dir(p.beam_direction.dir())
   {

     throws = 0;
     string fname = p.beam_inputroot;

     int fpdg[6] = {12, -12, 14, -14, 16, -16};
     
     string fneut[6];
     fneut[0] = p.beam_inputroot_nue;
     fneut[1] = p.beam_inputroot_nueb;
     fneut[2] = p.beam_inputroot_numu;
     fneut[3] = p.beam_inputroot_numub;
     fneut[4] = p.beam_inputroot_nutau;
     fneut[5] = p.beam_inputroot_nutaub;

     // Get Flux Inputs
     for (int i = 0; i < 6; i++){
       if (fneut[i].empty()) continue;
       
       // Get PDG
       int code = fpdg[i];
       
       // Setup new beam
       double m = mass(code);
       p.beam_inputroot_flux = fneut[i];
       p.beam_particle = code;
       p.beam_type = 5;
       beams[n] = new beam_singleroothist(p);
       n++;

     }

     double total = 0.0;
     for (int i = 0; i < n; i++){
       total += beams[i]->EProf().GetHist().Integral("width");
     }

     for (int i = 0; i < n; i++){
       w1[i] = beams[i]->EProf().GetHist().Integral("width")/ total;
       w2[i]=w1[i]*beams[i]->disratio();
       if(i)
	 {
	   w1[i]+=w1[i-1];
	   w2[i]+=w2[i-1];
	 }
     }

     p.beam_type = 6;
     p.beam_inputroot_flux = "";
     p.beam_particle = 0;
   }

   ~beam_mixedroothist()
     {
       while(n>0)
	 delete beams[--n];
     }
   
   /// get next particle form the beam
   virtual particle shoot(bool dis)
   { double *w=(dis?w2:w1);
     double x=frandom()*w[n-1];
     int i=0;
     throws++;
     
     while(x>w[i]){
       i++;
     }     
     return beams[i]->shoot(dis);
   }

   // Get N Beams
   int GetN() const {return n;};

   // Get Beam i
   beam_singleroothist* GetBeam(int i) const {return beams[i];};

   
};
#endif // _beam_mixedroothist_h_
