#include "mecevent_common.h"

//int noev=0;
double ml=0;
double ml2=0;
double width_q0=0;
double Bmec=0;
int ap=0;
int nucl=0;
bool PB=0;

////////////////////////////////////////

void mec_do_cc (particle *p, double ratio)
{
  // here is the isospin model; I assume that 3/5 times a pair is p-p and 2/5 times it is p-n
  // KM: I don't understand the comment above
  if (frandom () < ratio)
  {
    p[0].set_proton ();
    p[1].set_neutron();
    if(ap)
    {
      p[2].set_neutron();
      p[3].set_neutron();
    }
    else
    {
      p[2].set_proton();
      p[3].set_proton();
    }
  }
  else
  {
    if(ap)
    {
      p[0].set_proton();
      p[1].set_proton();
      p[2].set_neutron();
      p[3].set_proton ();
    }
    else
    {
      p[0].set_neutron();
      p[1].set_neutron();
      p[2].set_proton ();
      p[3].set_neutron();
    }
  }
}