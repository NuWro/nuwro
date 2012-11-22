#ifndef _Analyser_h_
#define _Analyser_h_

#include "event1.h"

/**
	@usage example:
	
 	Analyser A(p);
	
	for(a.start(); !a.end(); a.step() ) 
	{
	   for..... //nuwro loop 
	   {	....
	   
	   A.prepare_event();   
	   .....
	    makeevent();
	   .....
	   A.process_event();
	   
	   }
	    
	   A.partial_raport();
	  
	}
	A.final_raport();
	
}

*/

class Analyser
{
	public:
		params &p;
	
	public: 
	
		Analyser(params &p0):p(p0){}
		
		virtual void start()=0;
		virtual bool end()=0;
		virtual void step()=0; 
		
		virtual void prepare_event(event&e)=0;
		virtual void process_event(event&e)=0;
		
		virtual void partial_report()=0;
		virtual void final_report()=0;
		
        virtual ~Analyser(){}	
};


#endif
