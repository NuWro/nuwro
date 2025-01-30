#ifndef _BEAMRF_H_
#define _BEAMRF_H_

#include "rootf.h"
#include "beam.h"
#include "particle.h"
#include "params.h"
#include <dirent.h>
#include <vector>
#include "geomy.h"

double beam_length_scaling=1; 

int nu_pdg_from_mode(int mode)
{
//~ #define MODE_NUMU_PI      11  /* numu from pi+  */
//~ #define MODE_NUMU_K       12  /* numu from K+(2)*/
//~ #define MODE_NUMU_MU      13  /* numu from mu-  */
//~ #define MODE_NUMU_KPLUS   14  /* numu from K+(3)*/
//~ #define MODE_NUMU_K0      15  /* numu from K0(3)*/
//~ 
//~ #define MODE_NUMUB_PI     21  /* numu_bar from pi-  */
//~ #define MODE_NUMUB_K      22  /* numu_bar from K-(2)*/
//~ #define MODE_NUMUB_MU     23  /* numu_bar from mu+  */
//~ #define MODE_NUMUB_KMINUS 24  /* numu_bar from K-(3)*/
//~ #define MODE_NUMUB_K0     25  /* numu_bar from K0(3)*/
//~ 
//~ #define MODE_NUE_KPLUS    31  /* nue from K+ (Ke3) */
//~ #define MODE_NUE_K0       32  /* nue from K0L(Ke3) */
//~ #define MODE_NUE_MU       33  /* nue from Mu+      */
//~ #define MODE_NUE_PI       34  /* nue from pi+      */
//~ 
//~ #define MODE_NUEB_KMINUS  41  /* nue_bar from K- (Ke3) */
//~ #define MODE_NUEB_K0      42  /* nue_bar from K0L(Ke3) */
//~ #define MODE_NUEB_MU      43  /* nue_bar from Mu-      */
//~ #define MODE_NUEB_PI      44  /* nue_bar from pi-      */   
    switch(mode/10)
    {
        case 1: return  14; //nu_mu
        case 2: return -14; //nu_mu_bar
        case 3: return  12; //nu_e
        case 4: return -12; //nu_e
        default: cerr<<"Unknown reaction mode "<<mode<<" reading flux files"<<endl;
                 exit(16);
    }
    
}

particle nu_from_event(ND5Event e)
{
  int pdg=nu_pdg_from_mode(e.mode);
  particle p( pdg, 0.0 );
  double E=e.Enu*1000; // from GeV to MeV

  // Length scaling compares beam and geo units, sets beam units to match those in geo 
  p.r.x = e.xnu*beam_length_scaling;
  p.r.y = e.ynu*beam_length_scaling;

  // Read z component of neutrino position from file (in case beam window is not x-y plane)
  // if not specified by flux file this should default to 0
  p.r.z = e.znu*beam_length_scaling;

  p.r.t = 0;
  p.t=E;
  p.x=e.nnu[0]*E;
  p.y=e.nnu[1]*E;
  p.z=e.nnu[2]*E;
  p.travelled=1;

  return p;
}

vector<string>  root_files(string directory)
{
  vector<string> names;

  DIR * dp= opendir( directory.c_str() );
  if(dp==NULL)
  {
    cerr << "Directory \""<<directory<<"\" not found."<<endl;
    exit(3);
  }
  dirent * dirp;
  if(dp)
    while((dirp = readdir( dp )))
    {
      string name=dirp->d_name;
      if(name.find( string(".root") )  !=  string::npos)
        names.push_back( directory + "/" + name );
    }
  if( names.size() > 0 )
    cout << names.size()<<" root files ready to open.\n";
  else
  {
    cerr << "No root files found in directory \""<<directory<<"\""<<endl;
    exit(34);
  }
  return names;
}


class BeamRF : public beam
{
  ND5Event *events[100000];
    double P_region;
  int N=0;
  double * acum;
  double * acum2;
  int curevent;
  string folder;
  int first;
  int limit;
  double minx,miny,maxx,maxy;
  double POT=0;
  int nnu=0; //total num of events  
  double sum_weights=0; //sum of norms used
  double POT_per_file=0;

public:
  
  BeamRF(params &p,geomy *detector=NULL):
    folder(p.beam_folder),
    first(p.beam_file_first),
    limit(p.beam_file_limit),
    POT_per_file(p.beam_pot_per_file)
  { 
        P_region=1;
  
  // Set length units for positioning window to match those used by the geometry
  
  if(p.beam_length_units == p.geom_length_units) beam_length_scaling = 1; //if beam and geo units are same, no scaling needed
  else if(p.beam_length_units == "mm" && p.geom_length_units == "cm") beam_length_scaling = 0.1;    
  else if(p.beam_length_units == "cm" && p.geom_length_units == "mm") beam_length_scaling = 10;
  else if(p.beam_length_units == "m" && p.geom_length_units == "cm") beam_length_scaling = 10;
  else if(p.beam_length_units == "m" && p.geom_length_units == "mm") beam_length_scaling = 100;
  else if(p.beam_length_units == "cm" && p.geom_length_units == "m") beam_length_scaling = 0.1;
  else if(p.beam_length_units == "mm" && p.geom_length_units == "m") beam_length_scaling = 0.01;
  else {
     std::cout << "Unrecognised length units: " << p.beam_length_units << " " << p.geom_length_units << std::endl
               << "Use either mm, cm or m" << std::endl;
     exit(1);
  }  
     
    read_events(detector);
    cout<<endl;
    acum=new double[N];
    acum2=new double[N];
    double prev=0;
    for(int i=0;i<N;i++)
      acum[i]=prev+=events[i/100000][i%100000].norm;


    prev=0;
    for(int i=0;i<N;i++)
    { 
      ND5Event& ev=events[i/100000][i%100000];
      acum2[i]=prev+=ev.norm*ev.Enu;
    }
      
    // Total POT used is the number of pot per file x number of files.
    // Assumes each file has equal POT - might want to consider 
    // storing the POT info in the root flux files themselves in future
  
    POT = limit*POT_per_file;
    
    cout<<" nu/POT="<<nu_per_POT()<<endl;
    cout<<" POT/nu="<<1/nu_per_POT()<<endl;
    cout<<" nfiles="<<limit<<endl;

  // This asssumed flux is produced in x-y plane - not always true
  // Replace with more suitable calculation later

  //double surf=(maxx-minx)*(maxy-miny);
  //cout<<" Beam Surface="<<maxx-minx<<" cm x "<<maxy-miny<<" cm = "<<surf<<" cm2"<<endl;
  }
  
  ~BeamRF()
  {  
    for(int i=0;i<(N+100000-1)/100000;i++)
      delete []events[i];
    delete []acum;
    delete []acum2;
  }
  void read_events(geomy *detector=NULL)
  {
        int total=0;
        int unsaved=0;
        double normsall=0;
        double normsbad=0;

    vector<string> names=root_files(folder);
    int n=names.size();
    if(first>n)
      first=0;
    if(limit>n || limit==0) // 0 means no limit
      limit=n;
    for(int j=0; j < limit; j++)
    {   
      RootFReader<ND5Event> reader( names[j], "h3002");
      for(int k=0;k<reader.Count();k++)
            {   
                total++;
                const ND5Event *e=reader.GetEntry(k);
                normsall+=e->norm;
                if(detector)
                {
                   particle nu=nu_from_event(*e);
                   if(detector->is_hit_by(nu.p(),nu.r)){
                        store(*e);
			//adding weight calc
			sum_weights += e->norm;
			}
                   else
                   {
                        unsaved++;
                        normsbad+=e->norm;
                    }
                }
                else
                    store(*e);
            if(total%1000==0)
                cout<<total<<" events read from "<<j<<" files.\r"<<flush;
            }
    }
        if(detector)    
            P_region=1-normsbad/normsall;
        cout<<total<<" events read from "<<limit<<" files.\r"<<flush;
        cout<<endl<<total-unsaved<<" neutrinos can hit the target (will be used).\r"<<endl;
        cout<<endl<<int(P_region*1e6)/1e4<<"% POTS will be used.\r"<<endl;
  nnu = total-unsaved;
  }
  
////////////////////////////////////////////////////////////////////////
  void store( const ND5Event& e)
  {     

    // C Thorpe: This does not appear to affect the POT counting calculation
    // but the assumtion made here that the flux window is the x-y plane
    // is no longer correct - this needs updating

    if(N%100000==0) 
    {
//       cerr<< "File:"<<_nextFile<<" "<< N<<" beam events read...\r"<<flush;
       events[N/100000]=new ND5Event[100000];
    }
    events[N/100000][N%100000]=e;
    if(N==0)
    {
      minx=maxx=e.xnu;//in cm
      miny=maxy=e.ynu;//in cm
    }
    else
    {
      double x=e.xnu; //in cm
      double y=e.ynu; //in cm
      maxx=max(x,maxx);
      maxy=max(y,maxy);
      minx=min(x,minx);
      miny=min(y,miny);     
    }
    N++;  
  } 
  
  double nu_per_POT()
  {
     //if POT per file has been specified
     if(POT_per_file != -1) return sum_weights/POT;  

     //otherwise assume ND280 setup
     else return (acum[N-1]/limit)/1e21;//*P_region;
  }

  ////////////////////////////////////////////////////////////////////////
  virtual particle shoot(bool dis)
  {  int n=N;//events.size();
    double *acc=(dis?acum2:acum);
    double x=frandom()*acc[n-1];
    int i=0,j=n-1;
    while(i<j)   //binary search
    {
      int s=(i+j)/2;
      if(x<acc[s])
        j=s;
      else
        i=s+1;  
    }
    return nu_from_event(events[i/100000][i%100000]);
  } 
};
#endif // _BEAMRF_H_
