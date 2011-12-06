#ifndef _Analyser_h_
#define _Analyser_h_
#include "TMath.h"

#include "event1.h"
#include "boone1.h"
///class Analyser
/*
	@usage example:
	
-> 	Analyser a;
	
->	while()(a.loop(start,end,step)) // loop in axial mass for instance 
	 {
	   for..... //nuwro loop 
	   {	....
	   
->	   a.prepare_event(event &e);   //set axial mass for event
	   .....
	    makeevent();;
	   .....
->	   a.process_event(event &e);
	   
	   }
	    
->	   a.partial_raport(event &e);
	  
	  }
->	  a.final_raport(event &e);
	
}

*/


class Analyser
{
    ofstream plots[4][3];
    double mult[4][3];
    int scales[3];
	ofstream mafit,matables;
	double ma;
	boone *boon;
	public: 
	
	Analyser():mafit("ma_fit.txt"),matables("ma_tables.txt"),ma(-1),boon(NULL)
	{
     string names[]={"b2-","b2g-",
                     "b1-","b1g-"};
     scales[0]=1;
     scales[1]=0;
     scales[2]=-1;
     string scalenames[]={"fixed","samenorm","errmin"}; 
     stringstream s;
     
      for(int i=0;i<4;i++)
        for(int j=0;j<3;j++)
         { stringstream s;
            s<<names[i]<<"."<<scalenames[j]<<".txt\0"<<flush; 
	 	    plots[i][j].open(s.str().c_str());
	 	 }
	}
	
	bool loop(double start,double end, double step) 
	{ if(ma==-1) 
	     {ma=start;
	      delete boon;
	      boon = new boone();
	     }
	  else 
	     {ma+=step;
	      delete boon;
	      boon = new boone();
	     } 
	  return ma<=end;   
	}
	
	void prepare_event(event &e,params &p)
	{ 
		//p.sf_axial_mass=ma;
		p.qel_cc_axial_mass=ma;
    }
	void process_event(event &e,params &p)
	{
        boon->add(e,e.weight*12/6);
	}
	
	void partial_raport()
	{ cout<<boon->count()<< "events generated."<<endl;
	  cout<<"sigma="<<boon->avg()<< " cm2" <<endl;
	     { 
        double res[4][3];
		  for(int j=0;j<3;j++)
		    {
		     res[0][j]=boon->wdiff2(mult[0][j]=scales[j]);
		     res[1][j]=boon->wdiff2g(mult[1][j]=scales[j]);
		     
		     res[2][j]=boon->wdiff2q2(mult[2][j]=scales[j]);
		     res[3][j]=boon->wdiff2q2g(mult[3][j]=scales[j]);
		    }
		 
        for(int i=0;i<4;i++)
		  {for(int j=0;j<3;j++)
            {plots[i][j]<< ma <<setw(12)<<res[i][j]<<endl;
             }
		 }
	  cout<<endl<<"Ma="<<ma<<endl;
	  for(int i=0;i<4;i++)
	   { for(int j=0;j<3;j++)
	      {  int ndf=i<2?127:17;
	         ndf-=(j>0);
			 cout.precision(5);
			 cout<<"   "<<setw(12)<<res[i][j];
			 cout<<"("<<TMath::Prob(res[i][j],ndf)<<'/'<<ndf<<')';
   		     cout<<'*'<<setw(8)<<mult[i][j];
   		     
		   }  
		     cout<<endl;
	   }
		 
     }     
//     boon->printtable(matables);
     {stringstream s;
	  s<<"wyk2-"<<ma<<"\0"<<flush;
     ofstream fil(s.str().c_str());
     boon->printtable(fil,mult[3][2]);
     }

     {stringstream s;
	  s<<"wyk1-"<<ma<<"\0"<<flush;
     ofstream fil(s.str().c_str());
     boon->printq2table(fil,mult[7][2]);
     }
     {stringstream s;
	  s<<"wyk1mc-"<<ma<<"\0"<<flush;
     ofstream fil(s.str().c_str());
     boon->printmcq2table(fil,mult[7][2]);
     }
//     boon->printq2table(matables);
//     boon->printmcq2table(matables);
     
     stringstream s15;
     s15<<"ma="<<ma<<".root\0"<<flush;
     boon->save(s15.str().c_str());  

	
		
	}
	
	void final_report()
	{
    for(int i=0;i<4;i++)
       for(int j=0;j<3;j++)
        { 
	 	  plots[i][j].close();
	 	 }
	}
	
	~Analyser(){ delete boon;}
	
};




#endif
