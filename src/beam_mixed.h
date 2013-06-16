#ifndef _beam_mixed_h_
#define _beam_mixed_h_
#include "beam_uniform.h"


/// The beam shoots with identical particles (whose energy is defined by some 
/// energy profile but the direction of momentum is identical)
/// beam: this class provides the input of the initial beam particle 
/// into the generator.The beam particle is defined by providing 
/// its energy, direction and its PDG code to the constructor
/// of the beam class 

class beam_mixed : public beam
{
protected:
   int n;
   vec dir;
   beam_uniform* beams[10];
   double w1[10];
   double w2[10];
public:   
   beam_mixed(params& p):n(0),dir(p.beam_direction.dir())
   { stringstream s(p.beam_content);
     while (s)
     {
	   int code;
       double ratio;
       char perc;
       string line;
       getline(s,line);
       stringstream s1(line);
       s1>> code >> ratio>> perc;
       if(!s1) break;
//       cout<< setw(4) << code <<' '<<ratio<<' '<<perc<<endl;
       if(perc!='%') 
         {cerr<< "beam_content= "<<code<<' '<<perc<< " % sign expected"<<endl;
         exit(30);}
       w1[n]=ratio;
       double m=mass(code);
       string rest;
       getline(s1,rest);
       p.beam_particle=code;
       p.beam_energy=rest;
       beams[n]=new beam_uniform(p);
       w2[n]=w1[n]*beams[n]->disratio();
       if(n)
         {
			 w1[n]+=w1[n-1];
             w2[n]+=w2[n-1];
         }
 //      cout<< w1[n]<<' '<<w2[n]<<' '<<beams[n]->disratio()<<endl;
       n++;
     } 
       p.beam_energy="";
       p.beam_particle=0;
   //    cout<<"n="<<n<<endl;
   }
   ~beam_mixed()
   {while(n>0)
      delete beams[--n];
   }
   
   /// get next particle form the beam
   virtual particle shoot(bool dis)
   { double *w=(dis?w2:w1);
     double x=frandom()*w[n-1];
     int i=0;
     while(x>w[i]) 
       i++;
     return beams[i]->shoot(dis);
   }
   
};
#endif // _beam_mixed_h_
