#include <iomanip>
#include <sstream>
#include <vector>
#include "event1.h"
#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "qelevent.h"
#include "pdg.h"
#include "chooser.h"
#include "beam.h"
#include "beamHist.h"
#include "target_mixer.h"
#include "stdlib.h"
#include "pauli.h"
#include "dis/singlepion.h"
#include "dis/disevent.h"
#include "dis/resevent2.h"
#include "cohevent2.h"
#include "coh.h"
#include "args.h"
#include "kaskada7.h"
#include "sfevent.h"
#include "Analyser.h"
#include "geomy.h"
#include "ff.h"
#include "hist.h"
#include "nucleusmaker.h"

extern double SPP[2][2][2][3][40];
//extern double sppweight;

params *p1=NULL;
string data_dir;
#include "nuwro.h"

NuWro::~NuWro()
{
	delete mixer;
	delete detector;
	delete neutrino_beam;
	delete nucleuss;
}


geomy* NuWro::make_detector(params &p)
{
	if(p.geo_file.length())
	   {
		try{
			   if(p.geo_d.norm2()==0)
				 return new geomy(p.geo_file,p.geo_name);
				else 
				 return new geomy(p.geo_file,p.geo_name,p.geo_volume,p.geo_d,p.geo_o);
		    }
		catch(...)
			{
				cerr<<"Failed to make detector."<<endl;
				exit(0);
			}
        }
	else
	   return NULL;
	
}

int NuWro::init (int argc, char **argv)
{ 
//  dismode=false;
  dismode=true;
  set_dirs(argv[0]);
  a.read (argc, argv);
  p.read (a.input);
  p.read (a.params, "command line");
  p.list ();
  p1=&p;
  progress.open(a.progress);
  frandom_init(p.random_seed);
  if(p.beam_test_only==0 && p.kaskada_redo==0)
  if(p.dyn_dis_nc or p.dyn_res_nc  or p.dyn_dis_cc or p.dyn_res_cc )
  {  cout<<" Calculating the one pion functions ..."<<endl;
	  singlepion (p);
	/*  for (int c = 0; c < 40; c++)
		cout << 1210 + c * 20 << "  MeV -> " 
		 << setw(10)<< SPP[0][0][0][0][c] << "  " 
		 << setw(10)<< SPP[0][0][0][1][c] << "  "
		 << setw(10)<< SPP[0][0][0][2][c] << "  " 
		 << setw(10)<< SPP[0][0][1][0][c] << "  " 
		 << setw(10)<< SPP[0][0][1][1][c] << "  " 
			 << setw(10)<< SPP[0][0][1][2][c] << endl;
*/  }
   if(p.kaskada_redo==0)
   {
	   cout<<"Creating the beam ..."<<endl;
	   neutrino_beam=create_beam(p);  
	   if(neutrino_beam==NULL) 
		{cerr<<"No beam defined."<<endl;exit(1);}

	  nucleuss = make_nucleus(p);
	  if(p.target_type==1)
		mixer=new target_mixer(p);
	  detector=make_detector(p);
	  
	  if(p.geo_file!="" and detector==NULL) 
		 {cerr<<"Detector geometry not created."<<endl;exit(1);}
  }
  ff_configure(p);
  
  bool active[]={p.dyn_qel_cc,p.dyn_qel_nc,
                 p.dyn_res_cc,p.dyn_res_nc,
	             p.dyn_dis_cc,p.dyn_dis_nc,
	             p.dyn_coh_cc,p.dyn_coh_nc};
  //const int NPROC = sizeof(active)/sizeof(bool);
  procesy.reset(active);

}


void NuWro::makeevent(event* e, params &p)
{  static double maxweight=0,maxnorm=0;
   particle nu;
   material mat;
	int dyn = e->dyn;
	if(detector)
	{ 
       bool akceptacja=false;
       while(not akceptacja)
         { nu=neutrino_beam->shoot(dyn>1 && dyn<6 && dismode);
           if(nu.travelled>0 && p.beam_weighted==0)
              {if(nu.travelled<frandom()*maxnorm) 
                 continue;
               if(nu.travelled>maxnorm)
                 maxnorm=nu.travelled;
               nu.travelled=0;
              }
           nu.r=vec(nu.r)+p.beam_offset;
           if(nu.r.x==0 && nu.r.y==0 && nu.r.z==0)
              mat=detector->getpoint();
           else
              mat=detector->getpoint(nu.p(),nu.r);// jednostki ???
           double weight=mat.w_density;//*mat.w_length;
           if(weight>maxweight)
              maxweight=weight;
           akceptacja= (weight>=frandom()*maxweight) && (mat.Z+mat.N>0);
          }
       ///reset jądra
	   p.nucleus_p=mat.Z;
	   p.nucleus_n=mat.N;
	   p.nucleus_E_b=mat.e_bin;
	   p.nucleus_kf=mat.p_fermi;
	   if(mat.Z==0&&mat.N==0) 
	      throw "Empty isotope 00";
   }   
   else
     {
	  nu=neutrino_beam->shoot(dyn>1 && dyn<6 && dismode);
      nu.r=vec(nu.r)+p.beam_offset;
     }
   if(detector or mixer)
      {delete nucleuss;
       nucleuss= make_nucleus(p);
       //cout<<"make_nucleus "<<nucleuss->p<<" "<<nucleuss->n<<endl;
      }
   e->in.push_back (nu);	// insert neutrino
   if(dyn<6)
   {
	e->in.push_back (nucleuss->shoot ());	// insert target nucleon
	e->in[0].r=e->in[1].r;
    assert(e->in[1]*e->in[1]>0);
   }
   e->r=mat.r;
  // cout<< "mat.r="<< mat.r<<endl;

   e->weight=0;
   if(nu.travelled>0)
     e->norm=nu.travelled;
   // else e->norm remains 1;
     
	e->flag.cc  = dyn == 0 || dyn == 2 || dyn == 4 || dyn == 6;
	e->flag.nc  = dyn == 1 || dyn == 3 || dyn == 5 || dyn == 7;
  
	e->flag.qel = dyn == 0 || dyn == 1;
	e->flag.res = dyn == 2 || dyn == 3;
	e->flag.dis = dyn == 4 || dyn == 5;
	e->flag.coh = dyn == 6 || dyn == 7;
	
    if(p.beam_test_only)
      { e->weight=1;
        e->out.push_back(e->in[0]);
        return;
      }
    double factor=1.0;
     if(p.cc_smoothing and e->flag.cc and !e->flag.coh) // coherent makes own thing
      {double A=nucleuss->n+nucleuss->p;
         if(e->in[0].pdg>0)
           {
			 factor=nucleuss->n/A;
             e->in[1].set_neutron();
		   }
         else
           {
			 factor=nucleuss->p/A;
             e->in[1].set_proton();
		   }
      }
    e->par =p;
    switch (dyn)
	{
	case 0:
	  if (p.dyn_qel_cc)
	    {
		if(p.sf_method>0 and has_sf(*nucleuss))
	        sfevent2 (p, *e, *nucleuss);
		  else
		  {
	      qelevent1 (p, *e, *nucleuss, false);
	         if (p.pauli_blocking)
		        mypauli_qel (*e, *nucleuss);
	      }
	    }
	  break;		// qel cc
	case 1:
	  if (p.dyn_qel_nc)
		if(p.sf_method>0 and has_sf(*nucleuss))
	        sfevent2nc (p, *e, *nucleuss);
		  else
	    {
	      qelevent1 (p, *e, *nucleuss, true);
	      if (p.pauli_blocking)
		  mypauli_qel (*e, *nucleuss);
	    }
	  break;		// qel nc
	case 2:
	  if (p.dyn_res_cc)
	    {
	      resevent2 (p, *e, true);
	      if (p.pauli_blocking)
		  mypauli_spp (*e, *nucleuss);
	    }
	  break;		//true for dis cc
	case 3:
	  if (p.dyn_res_nc)
	    {
	      resevent2 (p, *e, false);
	      if (p.pauli_blocking)
		  mypauli_spp (*e, *nucleuss);
	    }
	  break;
	case 4:
	  if (p.dyn_dis_cc)
	    {
	      disevent (p, *e, true);
	      if (p.pauli_blocking)
		  mypauli_spp (*e, *nucleuss);
	    }
	  break;		//true for dis cc
	case 5:
	  if (p.dyn_dis_nc)
	    {
	      disevent (p, *e, false);
	      if (p.pauli_blocking)
		  mypauli_spp (*e, *nucleuss);
	    }
	  break;
	case 6:
	  if (p.dyn_coh_cc)
	    {if(p.coh_new) cohevent_cj (p, *e, *nucleuss, true);
	     else          cohevent2   (p, *e, *nucleuss, true);
	    }
	  break;
	case 7:
	  if (p.dyn_coh_nc)
	    {if(p.coh_new) cohevent_cj (p, *e, *nucleuss, false);
	     else          cohevent2   (p, *e, *nucleuss, false);
	    }
	  break;
	}
    e->weight*=factor;
//    e->r=mat.r;
//    cout <<mat.r<<endl;
//  e->density=density
        if (e->weight == 0)
          {
            e->out.clear ();
            e->out.push_back (e->in[0]);	//po co to wpisywanie ???
            e->out.push_back (e->in[1]);
          }
//      e->check();
}// end of makeevent

void NuWro::finishevent(event* e, params &p)
{
	for(int i=0;i<1/* e->in.size()*/;i++)
	    {  e->in[i].endproc=e->dyn;
		  registration(e->all,e->in[i]);
	    }
	for(int i=0;i<e->out.size();i++)
	   {e->out[i].mother=0;
	   registration(e->all,e->out[i]);
       } 
	if(p.beam_test_only)
	   return;
	    for(int j=0;j<e->out.size();j++)
			e->out[j].r=e->in[1].r;
			
        for(int j=0;j<e->in.size();j++)
		{
			particle p=e->in[j];
			p.endproc=e->dyn;
			registration(e->all,p);
		}
		
		//e->pr=nucleuss->Zr(); 	// 1. po co to?
		//e->nr=nucleuss->Nr(); 	// 2. powoduje break, segmentation fault
				
		if (e->flag.coh && p.kaskada_on)                        // copy particle from out to post if coherent interaction
		{
			for (int j = 0; j<e->out.size(); j++)
			{
				particle p = e->out[j];
				p.endproc = e->dyn;
				registration(e->all,p);
				e->post.push_back(p);
			}
		}
		else 
		  kaskadaevent(p,*e);  // runs only if p.kaskada_on is true
}//end of finishevent


void NuWro::raport(double i,double n,const char* text,int precision, int k, bool toFile)
{static int prev=-1;
 int proc=precision*i/n;
 if(proc!=prev)
   {   
	   prev=proc;
	   cerr.precision(3);
	   if(k>=0) 
	       cerr<<"Dyn["<<k<<"] ";
	   cerr.precision();
	   cerr<<showpoint<<proc*100.0/precision<<text<<'\r'<<flush;
	   if(toFile){
		progress.seekp(ios_base::beg);
		progress << proc*100.0/precision;
		progress.flush();
	   }
//	   printf("%f3.1 %s\r", proc/10.0,text);
   }
}//end of report


//////////////////////////////////////////////////////////////
//              Test events
//////////////////////////////////////////////////////////////    
void NuWro::test_events(params & p)
{
  
  if(p.number_of_test_events>0  && p.beam_test_only==0)
  {hist hq2((char*)"q2",0,2*GeV2,100);
   hist hT((char*)"T",0,2*GeV,100);
     bool active[]={p.dyn_qel_cc,p.dyn_qel_nc,
                 p.dyn_res_cc,p.dyn_res_nc,
	             p.dyn_dis_cc,p.dyn_dis_nc,
	             p.dyn_coh_cc,p.dyn_coh_nc};

  procesy.reset(active);
  for (int i = 0; i < p.number_of_test_events; i++)
    { 
       event *e = new event ();
       
       e->dyn = procesy.choose ();	// choose dynamics 
       if(mixer)
         mixer->prepare(p);
       makeevent(e,p);
       double bias=1;
       if(dismode && e->dyn>1 && e->dyn<6)
          bias=e->in[0].t;
       procesy.add (e->dyn, e->weight, bias);
       if(e->weight>0)
         {hq2.insert_value(-e->q2(),e->weight*cm2);
          hT.insert_value(e->in[0].E(),e->weight*cm2);
         }
       delete e;
       raport(i+1,p.number_of_test_events," % of test events ready...");
     } // end of nuwro loop
     cout<<endl;
    procesy.report();
    ofstream totals ("totals.txt",ios::app);
    totals<<p.beam_energy;
    for(int i=0;i<8;i++)
      totals << ' '<<procesy.avg(i);
    totals<<endl;
    procesy.set_weights_to_avg ();
    ofstream fhq2 ("q2.txt");
    hq2.wykres(fhq2,GeV2,1e-38*cm2/GeV2);
    ofstream fhT ("T.txt");
    hT.wykres(fhT,GeV,1e-38*cm2/GeV);
  } 
}


void NuWro::analizer_events(params &p)
  {  bool active[]={p.dyn_qel_cc,p.dyn_qel_nc,
                 p.dyn_res_cc,p.dyn_res_nc,
	             p.dyn_dis_cc,p.dyn_dis_nc,
	             p.dyn_coh_cc,p.dyn_coh_nc};
  
  params p1=p;
  Analyser analyser;
//   UserAction analyser;
// Useraction(procesy,mixer.e.p,neutrino_beam,detector)
  while(analyser.loop(1000,2000,50))
  { procesy.reset(active);
  for (int i = 0; i < p.number_of_test_events; i++)
    { 
       event *e = new event ();
       int k=0;
       e->dyn = procesy.choose ();	// choose dynamics 
/////////////
       analyser.prepare_event(*e,p);
       if(mixer)  
          mixer->prepare(p);
       makeevent(e,p);
   
       analyser.process_event(*e,p);   

       double bias=1;
       if(dismode && e->dyn>1 && e->dyn<6)
          bias=e->in[0].t;
   
       procesy.add (e->dyn, e->weight, bias);
       
       delete e;
       raport(i+1,p.number_of_test_events," % of analyser events ready...");
     } // end of nuwro loop
     
   analyser.partial_raport();
   procesy.report();
	
   } // end of analyser loop
   analyser.final_report(); 
   p=p1;
  }

void NuWro::real_events(params& p)
{
  if(p.number_of_events<1) 
    return;

	/// calculate desired number of events for each dynamics
	procesy.calculate_counts(p.number_of_events);
	{ /// Write cross sections and counts to screen and file
	  ofstream f((string(a.output)+".txt").c_str());
 	  procesy.short_report(cout);
	  procesy.short_report(f);
	}
	
	event *e = new event;
	
	string output=a.output;
	TFile *ff = new TFile (output.c_str(), "recreate");
	TTree *tf = new TTree ("treeout", "Tree of events");
	tf->Branch ("e", "event", &e);
	delete e;
    TH1 * xsections= new TH1D("xsections","xsections",8,0,7);
    for(int i=0;i<8;i++)
       {xsections->SetBinContent(i+1,procesy.avg(i));
        xsections->SetBinError(i+1,procesy.sigma(i)); 
	   }
    


/////////////////////////////////////////////////////////////
// The main loop in NPROC -- generate file with unweighted events
//////////////////////////////////////////////////////////////

  char filename[230];
  if(p.kaskada_redo==0)
  for (int k = 0; k < procesy.size(); k++)
    if(procesy.desired(k)>0)
    { 
      sprintf(filename,"%s.%d.part",a.output,k);
      TFile *f1 = new TFile (filename, "recreate");
      TTree *t1 = new TTree ("treeout", "Tree of events");

      e = new event ();
      t1->Branch ("e", "event", &e);
      delete e;

      while(procesy.ready(k)<procesy.desired(k))
		{
		  e = new event ();
		  e->dyn = k;

		  if(mixer)
			 mixer->prepare(p);
		  makeevent(e,p);
		  double bias=1;
          if(!p.beam_test_only && dismode & k>1 && k<6)
             bias=e->in[0].t;
          if (procesy.accept(e->dyn,e->weight,bias))
			{
				finishevent(e, p);
				t1->Fill ();
			}
		  delete e;  	  

		  raport(procesy.ready(k),procesy.desired(k)," % of events ready...",1000,k,bool(a.progress));
		}
      f1->Write ();

// elimination of spurious events for dynamics k 
// by copying only last desired[k] events to outfile 
      if(p.mixed_order==0)
      {
		  cout<<endl;
		  int nn = t1->GetEntries ();
		  int start = nn-procesy.desired(k);
		  for (int jj = start; jj < nn; jj++)
			{
			  e = new event();
			  t1->GetEntry (jj);
			  tf->Fill ();
			  delete e;
			  raport(jj-start+1,nn-start," % events copied...",100,k);
			}
		  cout<<endl;
	  }
	  else
	        cout<<endl;

      f1->Close ();
      delete f1;
      if(p.mixed_order==0)
        unlink(filename);
    };
//////////////////////////////////////////////////////////////////////////////////////
//                    end of the main loop in NPROC
//////////////////////////////////////////////////////////////////////////////////////
  if(p.mixed_order)
  { 
	  TFile *f[8];
	  TTree *t[8];
	  int n[8],u[8];
	  int ile=0;
	  e=new event();
	  for (int k = 0; k < procesy.size(); k++)
		if((u[k]=procesy.desired(k))>0)
		{ 
		  sprintf(filename,"%s.%d.part",a.output,k);
		  f[k] = new TFile (filename);
		  t[k] = (TTree*) f[k]->Get("treeout");
		  if(t[k]==NULL)
			{
				cerr<< "tree \"treeout\" not found in file \""<<filename<<"\""<<endl;
				exit(1);
			}
		     
		  t[k]->SetBranchAddress("e", &e);
		  n[k]=t[k]->GetEntries();
		  ile+=u[k];
        }
      int nn=ile;  
      while(ile>0)
      {
		int i=0;
		int x=frandom()*ile;
		while(x>=u[i]) 
		  x-=u[i++];
		t[i]->GetEntry(n[i]-u[i]);
		tf->Fill();
		ile--;
		u[i]--;	  
 	   raport(nn-ile,nn," % events copied...",100,i,false);
	  }  
      delete e;
	  for (int k = 0; k < procesy.size(); k++)
		if(procesy.desired(k))
		  {
			  f[k]->Close();
			  delete f[k];
  		      sprintf(filename,"%s.%d.part",a.output,k);
			  unlink (filename);
		  }
  }
  ff->Write ();
  ff->Close ();
  progress.close();
  remove(a.progress);
  delete ff;
  procesy.report();

  cout << "Output file: \"" << output << "\"" << endl;

}



void NuWro::kaskada_redo(string input,string output)
{ 
	event *e=new event;

	TFile *fi = new TFile (input.c_str());
	TTree* ti = (TTree*) fi->Get("treeout");
	if(ti==NULL) 
	{
		cerr<< "tree \"treeout\" not found in file \""<<input<<"\""<<endl;
		exit(1);
	}
	ti->SetBranchAddress("e", &e);

	TFile *ff= new TFile(output.c_str(),"recreate");    
	TTree *tf = new TTree("treeout","Tree of events");    
	tf->Branch("e","event",&e);


	int nn = ti->GetEntries ();
	for (int i = 0; i < nn; i++)
	{
//		e = new event();
		ti->GetEntry (i);
		e->clear_fsi();
		finishevent(e,p);
		tf->Fill ();
//		delete e;
		//if(i%1000==0)
		//cout<<i/1000<<"/"<<nn/1000<<"\r"<<endl;
		raport(i+1,nn," % events processed...",100,e->dyn,false);
	}
	cout<<endl;
	fi->Close ();
	delete fi;


	ff->Write();
	ff->Close();
	delete e;  				
	delete ff;
	cout<<"Output: \""<<output<<"\""<<endl;
}

#include "UserAction.h"

int NuWro::main (int argc, char **argv)
{
	try{
		init(argc,argv);
		if(p.kaskada_redo==1)
		   kaskada_redo(a.output,string(a.output)+".fsi.root");
	    else
	    {
			if(not p.beam_test_only)
				{
				 if(p.user_events>0)  	 
			        analizer_events(p);
			     else
			     	test_events(p);
			    }
			real_events(p);
		}
		genrand_write_state();
	}
	catch(string s){
		cout<<s<<endl;
	}
	catch(char const* s){
		cout<<s<<endl;
	}
	catch(...){
		cout<<"Nuwro failed"<<endl;
	}
}


NuWro nuwro;

