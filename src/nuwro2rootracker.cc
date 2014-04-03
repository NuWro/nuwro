//nuwro to rootracker converter, written by Pawel Przewlocki vel Pafcio, 2013.

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <sys/stat.h> 
#include "event1.h"
#include "event1dict.h"
#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TBits.h"
#include "pdg.h"
#include "generatormt.h"

bool FileExists(string strFilename) {
  struct stat stFileInfo;
  bool blnReturn;
  int intStat;
  // Attempt to get the file attributes
  intStat = stat(strFilename.c_str(),&stFileInfo);
  if(intStat == 0) {
    // We were able to get the file attributes
    // so the file obviously exists.
    blnReturn = true;
  } else {
    // We were not able to get the file attributes.
    // This may mean that we don't have permission to
    // access the folder which contains this file. If you
    // need to do that level of checking, lookup the
    // return values of stat which will give you
    // more details on why stat failed.
    blnReturn = false;
  }
  
  return(blnReturn);
}

int GetNeutChannel(event *e){
  //simplified neut reaction codes
  //TODO: coh correct, what about antinus?
  if (e->flag.qel){
    if (e->flag.cc){
      return 1;
    }else{
      return 51;
    }
  }
  if (e->flag.res){//this should be improved
    if (e->flag.cc){
      return 11;
    }else{
      return 31;
    }
  }
  if (e->flag.dis){
    if (e->flag.cc){
      return 26;
    }else{
      return 46;
    }
  }
  if (e->flag.coh){//CORRECTED. to be corrected should be 16
    if (e->flag.cc){
      return 16;
    }else{
      return 36;
    }
  }
  if (e->flag.mec){
    if (e->flag.cc){
      return 70;
    }else{
      return 90;
    }
  }
  return 100;
}


int main (int argc, char *argv[]){
  if (argc < 3){    
      std::cout << "Usage: nuwro2rootracker_series nuwrofile.root rootrackerfile_template.root [#evts_in_a_file]" << std::endl;
      std::cout << "Converts nuwro output file to a series of rootracker tree format files." << std::endl;
      std::cout << "The @ character in the template will be substituted with file number." << std::endl;
      std::cout << "If #evts not specified, one output file will be created." << std::endl;

  }else{
    if (!FileExists(argv[1])){
	  std::cout << argv[1] << ": File does not exist or is inaccessible." << std::endl;
	  return 1;
    }
    event *e = new event;
    TFile *fin = new TFile(argv[1]);
    TTree *tt1 = (TTree*)fin->Get("treeout");
    tt1->SetBranchAddress ("e", &e);
    int n = tt1->GetEntries();
    int ncopy=n;//by default one big rootracker with n events in it.
    bool onefile=true;
    if (argc>=4){//we have number of events specified
      sscanf(argv[3],"%d",&ncopy);
      onefile=false;
      printf("Each file will have %d evts.\n",ncopy);
    }

      
      //vars
    int channel;

   //rooTracker
    /// The generator-specific event flags.
    TBits* fEvtFlags;

    /// The generator-specific string with the 'event code'
    TObjString* fEvtCode;

    /// The sequence number of the event (the event number).
    int fEvtNum;

    /// The cross section for the event (1E-38 cm2)
    double fEvtXSec;

    /// The differential cross section for the event kinematics 
    /// (1E-38 cm2/{K^n})
    double fEvtDXSec;

    /// The weight for the event
    double fEvtWght;

    /// The probability for the event (given the cross section, path lengths,
    double fEvtProb;

    /// The event vertex position in detector coord syst (in meters and seconds).
    double fEvtVtx[4];

    /// The number of particles in the particle arrays to track
    int fStdHepN;

    /// The maximum number of particles that can be in the particle arrays.
    static const int kNPmax = 4000;

    /// The PDG codes for the particles to track. This may include generator
    /// specific codes for pseudo particles.
    int fStdHepPdg[kNPmax]; //[fStdHepN]

    /// The a generator specific status for each particle. Particles with a
    /// status equal to 1 will be tracked.
    int fStdHepStatus[kNPmax]; //[fStdHepN]

    /// The position (x, y, z, t) (fm, second) of the particle in the nuclear frame
    double fStdHepX4[kNPmax][4]; //[fStdHepN]

    /// The 4-momentum (px, py, pz, E) of the particle in the LAB frame (GeV)
    double fStdHepP4[kNPmax][4]; //[fStdHepN]

    /// The particle polarization vector.
    double fStdHepPolz[kNPmax][3]; //[fStdHepN]

    /// The index of the first daughter of the particle in the arrays.
    int fStdHepFd[kNPmax]; //[fStdHepN]

    /// The index last daughter of the particle in the arrays.
    int fStdHepLd[kNPmax]; //[fStdHepN]

    /// The index of the first mother of the particle in there arrays.
    int fStdHepFm[kNPmax]; //[fStdHepN]

    /// The index of the last mother of the particle in the arrays.
    int fStdHepLm[kNPmax]; //[fStdHepN]

    /// The PDG code of the particle which created the parent neutrino.
    int fNuParentPdg;

    /// The interaction mode at the vertex which created the parent neutrino.
    /// This is normally the decay mode of the parent particle.
    int fNuParentDecMode;

    /// The 4 momentum of the particle at the vertex which created the parent
    /// neutrino. This is normally the momentum of the parent particle at the
    /// decay point.
    double fNuParentDecP4[4];

    /// The position (meters, seconds) of the vertex at which the neutrino was
    /// created. This uses the target as the origin.
    double fNuParentDecX4[4];

    /// The momentum (GeV) of the parent particle at it's production point.
    double fNuParentProP4[4];

    /// The position (meters, seconds) of the parent particle at it's
    /// production point. This uses the target as the origin.
    double fNuParentProX4[4];

    /// The vertex ID of the parent particle vertex.
    int fNuParentProNVtx;



      // Create a TTree  
      
    const int fEmpty=-999999;  
      

    
    printf("Number of entries in the input file: %d\n", n);
    double coef=1e38/n;
    int i0=0;

    char CodeStr[20];
    int ncopied=0;
    int i=0,fnum=0;
    string outtemplate=argv[2];
    int index=outtemplate.find("@");
    if (!onefile && index==std::string::npos){
      printf("No @ character in the output template, terminating.\n");
      exit(1);
    }
      //event loops
    while (i<n){//until end-of-file
      string fname=outtemplate;
      char buf[5];
      sprintf(buf,"%d",fnum);
      if (!onefile) fname.replace(index,1,buf);
      TFile *fout=TFile::Open(fname.c_str(),"RECREATE");
      fnum++;
      TTree *fOutputTree = new TTree("nRooTracker","RooTracker");
      fEvtFlags = NULL;
      fOutputTree->Branch("EvtFlags", "TBits",     &fEvtFlags, 32000, 1);
      fEvtCode = NULL;
      fOutputTree->Branch("EvtCode",  "TObjString", &fEvtCode, 32000, 1);
      fOutputTree->Branch("EvtNum",       &fEvtNum     ,     "EvtNum/I"                );
      fOutputTree->Branch("EvtXSec",      &fEvtXSec    ,     "EvtXSec/D"               );
      fOutputTree->Branch("EvtDXSec",     &fEvtDXSec   ,     "EvtDXSec/D"              );
      fOutputTree->Branch("EvtWght",      &fEvtWght    ,     "EvtWght/D"               );
      fOutputTree->Branch("EvtProb",      &fEvtProb    ,     "EvtProb/D"               );
      // vertex in det coord. [m],[s]
      fOutputTree->Branch("EvtVtx",        fEvtVtx     ,     "EvtVtx[4]/D"             );
      fOutputTree->Branch("StdHepN",      &fStdHepN    ,     "StdHepN/I"               );
      fOutputTree->Branch("StdHepPdg",     fStdHepPdg  ,     "StdHepPdg[StdHepN]/I"    );
      fOutputTree->Branch("StdHepStatus",  fStdHepStatus,    "StdHepStatus[StdHepN]/I" );
      fOutputTree->Branch("StdHepX4",      fStdHepX4,        "StdHepX4[StdHepN][4]/D"  );
      // px,py,pz,E in LAB, GeV
      fOutputTree->Branch("StdHepP4",      fStdHepP4,        "StdHepP4[StdHepN][4]/D"  );
      fOutputTree->Branch("StdHepPolz",    fStdHepPolz,      "StdHepPolz[StdHepN][3]/D");
      fOutputTree->Branch("StdHepFd",      fStdHepFd,        "StdHepFd[StdHepN]/I"     );
      fOutputTree->Branch("StdHepLd",      fStdHepLd,        "StdHepLd[StdHepN]/I"     );
      fOutputTree->Branch("StdHepFm",      fStdHepFm,        "StdHepFm[StdHepN]/I"     );
      fOutputTree->Branch("StdHepLm",      fStdHepLm,        "StdHepLm[StdHepN]/I"     );
      
      fOutputTree->Branch("NuParentPdg",    &fNuParentPdg,     "NuParentPdg/I"         );
      fOutputTree->Branch("NuParentDecMode",&fNuParentDecMode, "NuParentDecMode/I"     );
      fOutputTree->Branch("NuParentDecP4",   fNuParentDecP4,   "NuParentDecP4[4]/D"    );
      fOutputTree->Branch("NuParentDecX4",   fNuParentDecX4,   "NuParentDecX4[4]/D"    );
      fOutputTree->Branch("NuParentProP4",   fNuParentProP4,   "NuParentProP4[4]/D"    );
      fOutputTree->Branch("NuParentProX4",   fNuParentProX4,   "NuParentProX4[4]/D"    );
      fOutputTree->Branch("NuParentProNVtx",&fNuParentProNVtx, "NuParentProNVtx/I"     );  //obsolete
      
      ncopied=0;//reset to start copying again
      while (ncopied<ncopy){
	  if (i>=n) break; //end of file, wrapping up
	  tt1->GetEntry(i);
	  i++;
	  if (!(e->weight>0.)) continue; //we don't copy events with 0 weight, that's for weighted only
                                        //regular events don't have 0 weights
	  ncopied++;
	  fEvtFlags->Set(8,"00000000");
	  sprintf(CodeStr,"%d",GetNeutChannel(e));
	  fEvtCode->SetString(CodeStr);
	  fEvtNum=i;
	  fEvtXSec=e->weight*coef; //unimportant for regular events
	  fEvtDXSec=fEmpty;
	  fEvtWght=e->weight*1e38; //unimportant for regular events
	  fEvtProb=fEmpty;
	  fEvtVtx[0]=e->r.x/1000;
	  fEvtVtx[1]=e->r.y/1000;
	  fEvtVtx[2]=e->r.z/1000;
	  fEvtVtx[3]=0.;
	  fStdHepN=0;
	  

	  for (int nin=0; nin<e->in.size();nin++){//incoming particles
		  fStdHepPdg[fStdHepN]=e->in[nin].pdg; 
		  fStdHepStatus[fStdHepN]=0;//incoming
		  for (int k=0;k<4;k++){
			  fStdHepX4[fStdHepN][k]=fEmpty;
		  }
		  for (int k=0;k<3;k++){
			  fStdHepPolz[fStdHepN][k]=fEmpty;
		  }
		  fStdHepP4[fStdHepN][3]=e->in[nin].E()/1000;
		  fStdHepP4[fStdHepN][0]=e->in[nin].x/1000;
		  fStdHepP4[fStdHepN][1]=e->in[nin].y/1000;
		  fStdHepP4[fStdHepN][2]=e->in[nin].z/1000;
		  fStdHepFd[fStdHepN]=fEmpty;
		  fStdHepLd[fStdHepN]=fEmpty;
		  fStdHepFm[fStdHepN]=fEmpty;
		  fStdHepLm[fStdHepN]=fEmpty;
		  fStdHepN++;
	  }
	  vector <particle> & out =(e->post.size()>0 ? e->post : e->out);
	  for (int nout=0; nout< out.size();nout++){//outgoing
		  fStdHepPdg[fStdHepN]=out[nout].pdg; 
		  fStdHepStatus[fStdHepN]=1;//outgoing
		  for (int k=0;k<4;k++){
			  fStdHepX4[fStdHepN][k]=fEmpty;
		  }
		  for (int k=0;k<3;k++){
			  fStdHepPolz[fStdHepN][k]=fEmpty;
		  }
		  fStdHepP4[fStdHepN][3]=out[nout].E()/1000;
		  fStdHepP4[fStdHepN][0]=out[nout].x/1000;
		  fStdHepP4[fStdHepN][1]=out[nout].y/1000;
		  fStdHepP4[fStdHepN][2]=out[nout].z/1000;
		  fStdHepFd[fStdHepN]=fEmpty;
		  fStdHepLd[fStdHepN]=fEmpty;
		  fStdHepFm[fStdHepN]=fEmpty;
		  fStdHepLm[fStdHepN]=fEmpty;
		  fStdHepN++;
	  }
	  //random values for parent info, disregard!
	  fNuParentPdg=211;
	  fNuParentDecMode=11;
	  //for (int k=0;k<4;k++){
	  fNuParentDecP4[0]=0.;fNuParentDecP4[1]=-0.410;fNuParentDecP4[2]=6.630;fNuParentDecP4[3]=6.644;
	  fNuParentProP4[0]=0.202;fNuParentProP4[1]=-0.227;fNuParentProP4[2]=6.738;fNuParentProP4[3]=6.747;
	  fNuParentDecX4[0]=3.361;fNuParentDecX4[1]=-212.127;fNuParentDecX4[2]=3519.273;fNuParentDecX4[3]=0.;
	  fNuParentProX4[0]=0.935;fNuParentProX4[1]=31.355;fNuParentProX4[2]=-459.387;fNuParentProX4[3]=0;
	  //}
	  fNuParentProNVtx=1;

	  fOutputTree->Fill();
	  printf("\r%d%% done.",(int)((i+1)*100/n));
      }
      fout->Write();
      fout->Close();
      printf("\nFile %s written, %d events copied, at %dth event in file.\n",fname.c_str(),ncopied,i);
      //delete fOutputTree;
      //delete fout;
    }

    fin->Close();
    delete fin;
    printf("\n%d events copied.\n",i);
  }
}
