#ifndef _analyser_h_
#define _analyser_h_

#include "event1.h"
///class Analyser
/*
///class Analyser
	usage:
 	Analyser rap(init); //init 0 or 1
	
	while()(analyser.loop())
	 {
		 nuwro loop();
		 * {	....
	   analyser.prepare_event(event &e);
	   .....
	   
	   analyser.process_event(event &e);
	     }
	     * //end nuwro loop
	   analyse.partial_raport(event &e);
	  
	  }
	analysr.final_raport(event &e);
	
}

*/


class Analyser
{
    boone b;	
	
	public: 
	
	Analyser(){}
	
	init();
	
	bool loop(); 
	
	prepare_event(event &e);
	
	process_event(event &e);
	
	partial_raport(event &e);
	
	end_loop();
	
	
	final_raport(event &e);
	
	~Analyser();
	
}




#endif
