#ifndef _UserAction_h_
#define _UserAction_h_
#include "event1.h"
#include "nucleus.h"


// typical useraction function
void NuWro::UserAction()
{params p1=p;
 event *e;

/// 
///_________________________________________  
///some kind of loop with changing parameters
//for(;;) 
///_________________________________________  
  { _procesy.reset(p);
  /// prepare stuctures for this values of parameters
  
  ///________________________________________________
  for (int i = 0; i < p.number_of_test_events; i++)     // OK
    { 
       e = new event ();	                            // OK
       int k=_procesy.choose ();                        // OK
       e->dyn = _procesy.dyn (k);	// choose dynamics  // OK
       ///_________________________________________  
       /// modify event before 
          
       ///_____________________________________
       if(_mixer)
       	   _mixer->prepare(p);                           // OK
       makeevent(e); // OK
       ///_________________________________________  
       /// analize event fill histograms etc.
   
       ///____________________________________
       _procesy.add (k, e->weight,e->in[0].t);            // OK
       delete e;                                          // OK
       raport(i+1,p.number_of_test_events," % of User events ready..."); // OK
     } // end of typical nuwro User Action loop
  ///_________________________________________  
  /// report before changing parameters (save histograms)
  
  
  ///_____________________________________ 
	_procesy.report();                                     //OK
	
   } // end of  loop
///_________________________________________  
///   Final report 

///_________________________________________  
   p=p1;                                                   //OK
}


#endif
