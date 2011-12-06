#ifndef _nuwro_h_
#define _nuwro_h_
#include "event1.h"
#include "params.h"
#include "args.h"
#include "geomy.h"
#include "target_mixer.h"
#include "beam.h"
#include <ostream>
#include "nucleus.h"
#include "chooser.h"


class NuWro
{ 
 public:
	geomy* make_detector(params &p);
    void makeevent(event* e, params &p);
    void finishevent(event* e, params &p);
    void raport(double i,double n,const char* text,int precision=1000, int k=-1,bool toFile=false);
    int init  (int argc, char **argv);
    void test_events(params &p);
    void analizer_events(params &p);
    void UserAction(params& p);
    void real_events(params &p);
    void kaskada_redo(string input,string output);
    int main (int argc, char **argv);
    ~NuWro();

 private:
        params p;  
	args a;
	chooser < 8 > procesy;
	ofstream progress;
	geomy* detector;
	beam* neutrino_beam;
	nucleus * nucleuss;  
	target_mixer *mixer;
	bool dismode;
};


extern NuWro nuwro;

#endif
