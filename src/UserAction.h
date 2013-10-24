#ifndef _UserAction_h_
#define _UserAction_h_
#include "event1.h"
#include "nucleus.h"


// typical useraction function
void NuWro::UserAction(params& p)
{params p1=p;
 event *e;
 bool active[]={p.dyn_qel_cc,p.dyn_qel_nc,
                 p.dyn_res_cc,p.dyn_res_nc,
	             p.dyn_dis_cc,p.dyn_dis_nc,
	             p.dyn_coh_cc,p.dyn_coh_nc};

/// 
///_________________________________________  
///some kind of loop with changing parameters
//for(;;) 
///_________________________________________  
  { _procesy.reset(active);
  /// prepare stuctures for this values of parameters
  
  ///________________________________________________
  for (int i = 0; i < p.number_of_test_events; i++)    // OK
    { 
       e = new event ();	                            // OK
       int k=0;                                         // OK
       e->dyn = _procesy.choose ();	// choose dynamics  // OK
       ///_________________________________________  
       /// modify event before 
          
       ///_____________________________________
       if(_mixer)
       	   _mixer->prepare(p);                            // OK
       makeevent(e,p); // OK
       ///_________________________________________  
       /// analize event fill histograms etc.
   
       ///____________________________________
       _procesy.add (e->dyn, e->weight,e->in[0].t);                   // OK
       delete e;                                          // OK
       raport(i+1,p.number_of_test_events," % of User events ready..."); // OK
     } // end of typical nuwro User Action loop
  ///_________________________________________  
  /// report before changing parameters (save histograms)
  
  
  ///_____________________________________ 
	_procesy.report();                                      //OK
	
   } // end of  loop
///_________________________________________  
///   Final report 

///_________________________________________  
   p=p1;                                                //OK
}


#endif
