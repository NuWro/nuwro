#ifndef _proctable_h_
#define _proctable_h_
#include <string>

using namespace std;


enum {
     none=-1,

     qel_cc=0,
     qel_nc=1,
     res_cc=2,
     res_nc=3,
     dis_cc=4,
     dis_nc=5,
     coh_cc=6, 
     coh_nc=7,
     mec_cc=8,
     mec_nc=9,

     nucleon_el=10, 
     nucleon_ce=11, 
     nucleon_spp=12, 
     nucleon_dpp=13, 
      
     pion_el=20, 
     pion_ce=21, 
     pion_spp=22, 
     pion_dpp=23, 
     pion_tpp=24,
     pion_abs=25, 

     jailed=99,
     escape=100
     };



#endif
