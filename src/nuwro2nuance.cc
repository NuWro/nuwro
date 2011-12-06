//nuwro to nuance converter, written by Pawel Przewlocki vel Pafcio, 2011.

#include <iostream>
#include <stdio.h>
#include <vector>
#include <sys/stat.h> 
#include "event1.h"
#include "event1dict.h"
#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
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

int GetChannel(event *e){
  //simplified nuance reaction codes
  if (e->flag.qel){
    if (e->flag.cc){
      return 1;
    }else{
      return 2;
    }
  }
  if (e->flag.res){//this should be improved
    if (e->flag.cc){
      return 16;
    }else{
      return 17;
    }
  }
  if (e->flag.dis){
    if (e->flag.cc){
      return 91;
    }else{
      return 92;
    }
  }
  if (e->flag.coh){
    if (e->flag.cc){
      return 96;
    }else{
      return 97;
    }
  }
  return 100;
}

int main (int argc, char *argv[]){
  if (argc == 1 || (argc == 2 && strcmp(argv[1],"--help")==0)){    
      std::cout << "Usage: nuwro2nuancekin nuwrofile.root [nuancefile.kin]" << std::endl;
      std::cout << "Converts nuwro output file to nuance kin text format." << std::endl;
      std::cout << "If second arg not given, prints everything onto the screen," << std::endl;
      std::cout << "so you can type nuwro2nuancekin in.root > out.kin." << std::endl;
  }else{
      if (!FileExists(argv[1])){
	std::cout << argv[1] << ": File does not exist or is inaccessible." << std::endl;
	return 1;
      }
      event *e = new event;
      TFile *fin = new TFile(argv[1]);
      FILE *fout;
      TTree *tt1 = (TTree*)fin->Get("treeout");
      tt1->SetBranchAddress ("e", &e);
      int n = tt1->GetEntries ();
      bool filemode=false;
      if (argc==3){
	fout=fopen(argv[2],"w");
	printf("Number of entries in input file: %d\n", n);
	filemode=true;
      }else{
	fout=stdout;
      }
      float mom=0;
      for (int i = 0; i < n; i++){
	  tt1->GetEntry(i);
	  fprintf(fout," $ begin\n");
	  fprintf(fout," $ nuance %d\n",GetChannel(e));
	  fprintf(fout," $ vertex %.1f %.1f %.1f %.4f\n",e->r.x/10,e->r.y/10,e->r.z/10,0.); 
	      //coordinates in cm      
	      //last param not important, set to 0
	  for (int nin=0; nin<e->in.size();nin++){
	    mom=vec(e->in[nin].x,e->in[nin].y,e->in[nin].z).length(); //momentum
	    if (mom==0.){
	      fprintf(fout," $ track %d %.4f %.5f %.5f %.5f %d\n",
		    e->in[nin].pdg,e->in[nin].E(),0.,0.,0.,-1); 
		    //last param -1 - incoming
	      
	    }else{
	      fprintf(fout," $ track %d %.4f %.5f %.5f %.5f %d\n",
		    e->in[nin].pdg,e->in[nin].E(),e->in[nin].x/mom,e->in[nin].y/mom,e->in[nin].z/mom,-1); 
		    //last param -1 - incoming
//	      fprintf(fout,"$ debug\t%d\t%f\t%f\t%f\t%f\t%d\n",
//		    e->in[nin].pdg,mom,e->in[nin].x,e->in[nin].y,e->in[nin].z,-1); 
	    }
	  }
	  fprintf(fout," $ info %d 0 0.0\n",i);//these numbers are irrelevant - we'll use the first as evt number
	  vector <particle> & out =(e->post.size()>0 ? e->post : e->out);
	  for (int nout=0; nout< out.size();nout++){
	    mom=out[nout].momentum(); //momentum
	    if (mom==0.){
	      fprintf(fout," $ track %d %.4f %.5f %.5f %.5f %d\n",
		    out[nout].pdg,out[nout].E(),0.,0.,0.,0); 
		    //last param 0 - outgoing final state
	    }else{
	      fprintf(fout," $ track %d %.4f %.5f %.5f %.5f %d\n",
		    out[nout].pdg,out[nout].E(),out[nout].x/mom,out[nout].y/mom,out[nout].z/mom,0); 
		    //last param 0 - outgoing final state
//	      fprintf(fout,"$ debug\t%d\t%f\t%f\t%f\t%f\t%d\n",
//		    out[nout].pdg,mom,out[nout].x,out[nout].y,out[nout].z,-1); 	      
	    }
	  }
	  fprintf(fout," $ end\n");
	  if (filemode) printf("\r%d%% done.",(int)((i+1)*100/n));
      }
      fprintf(fout," $ stop\n");
      printf("\n");
      fin->Close ();
      if (filemode) fclose(fout);
      delete fin;  
  }
}
