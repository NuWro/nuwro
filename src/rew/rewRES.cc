#include "../params.h"
#include "rewparams.h"
#include "../event1.h"
#include "../nucleus.h"
#include "ff.h"


#include "../dis/alfa.h"
#include "../dis/charge.h"
#include "../dis/delta.h"
#include "../dis/dis_cr_sec.h"
#include "../dis/resevent2.h"
#include "../dis/singlepion.h"
#include "../dis/LeptonMass.h"


extern double SPP[2][2][2][3][40];

extern "C" {
void shhpythiaitokay_(void);
void youcanspeaknowpythia_(void);
}

void SetupSPP(params &param) {
  if (true) {  //! CheckSPPSetup()){ -- No way to know in general
    std::cout << "[INFO]: Setting up singlepion tables..." << std::endl;
    shhpythiaitokay_();
    singlepion(param);
    youcanspeaknowpythia_();
    std::cout << "[INFO]: Set up singlepion tables!" << std::endl;
  }
}
void SetupSPP() {
  std::cout << "[WARN]: Setting up singlepion tables with default parameter "
               "set."
            << std::endl;
  params default_p;
  SetupSPP(default_p);
}

double GetEBind(event &nuwro_event, params const &rwparams) {
  switch (rwparams.nucleus_target) {
    case 0: return 0;  // free
    case 1: return rwparams.nucleus_E_b;  // FG
    case 2: return 0;  // local FG
    case 3: return 0;  // Bodek
    case 4: return binen(nuwro_event.in[1].p(), rwparams.nucleus_p, rwparams.nucleus_n);// effective SF
    case 5: return deuter_binen(nuwro_event.in[1].p()); // deuterium
    case 6: return rwparams.nucleus_E_b;  // deuterium like Fermi gas
    default:return 0; 
  }
}


double calcRES(event & e, params &p, nucleus &t)
{
	if(!e.flag.res)
		return 1;
	return 1;
}
