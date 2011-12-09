#ifndef _params_h_
#define _params_h_
#include <string>
#include <string.h>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>

#include "vec.h"
#include "dirs.h"

// all the parameters (type, name ,default_value) are defined here
#define PARAMS_ALL()\
        PARAM(int, random_seed, 1)\
        PARAM(int, number_of_events, 10000)\
        PARAM(int, number_of_test_events,10000)\
        PARAM(int, user_events,0)\
        PARAM(int, beam_type, 0)\
        PARAM(line, beam_energy,"1000")\
        PARAM(int, beam_particle, 14)\
        PARAM(vec, beam_direction, vec(0,0,1))\
        /* beam type=1 */ \
        PARAM(line, beam_content,"")\
        /* beam type=2 */ \
        PARAM(string, beam_folder, "flux")\
        PARAM(int, beam_file_first, 1)\
        PARAM(int, beam_file_limit, 0)\
        PARAM(bool, beam_weighted, 0)\
        /* beam type=3 */ \
        PARAM(string, beam_file, "beam/ND280hist.txt")\
        PARAM(vec, beam_offset, vec(0,0,0))\
        PARAM(int, beam_placement, 0)\
        PARAM(int, beam_test_only, 0)\
        PARAM(int, nucleus_p, 0)\
        PARAM(int, nucleus_n, 0)\
        PARAM(string,nucleus_density,"0") \
        PARAM(double,nucleus_E_b,0)\
        PARAM(double,nucleus_kf,0)\
		PARAM(int, nucleus_target,0)\
        PARAM(int, nucleus_model,0)\
		PARAM(int, target_type,0)\
		PARAM(line, target_content,"")\
        PARAM(string, geo_file,"")\
        PARAM(string, geo_name,"ND280Geometry")\
        PARAM(string, geo_volume,"")\
        PARAM(vec, geo_o,vec(0,0,0))\
        PARAM(vec, geo_d,vec(0,0,0))\
        PARAM(bool, dyn_qel_cc,0)\
        PARAM(bool, dyn_qel_nc,0)\
        PARAM(bool, dyn_res_cc,0)\
        PARAM(bool, dyn_res_nc,0)\
        PARAM(bool, dyn_dis_cc,0)\
        PARAM(bool, dyn_dis_nc,0)\
        PARAM(bool, dyn_coh_cc,0)\
        PARAM(bool, dyn_coh_nc,0)\
/*       PARAM(int, qel_kinematics,0)*/ \
        PARAM(int, qel_vector_ff_set, 2)\
        PARAM(int, qel_axial_ff_set, 1)\
		PARAM(double, qel_cc_vector_mass, 840)\
		PARAM(double, qel_cc_axial_mass,1030)\
        PARAM(double, qel_nc_axial_mass,1030)\
        PARAM(double, qel_s_axial_mass,1030)\
		PARAM(int, qel_strange, 0)\
		PARAM(int, qel_strangeEM, 0)\
		PARAM(double, delta_s, 0)\
		PARAM(bool, flux_correction, 1)\
        PARAM(int, delta_FF_set,1)\
		PARAM(double, pion_axial_mass, 1)\
		PARAM(double, pion_C5A, 1.2)\
        PARAM(double, res_dis_cut,0)\
        PARAM(int, spp_precision,0)\
/*      PARAM(bool, qel_new,0)\
        PARAM(bool, qel_cc_japan,0)\
        PARAM(bool, qel_relati,0)\
*/      PARAM(bool, coh_mass_correction,0)\
        PARAM(bool, coh_new,1)\
        PARAM(bool, kaskada_on,1)\
        PARAM(bool, kaskada_newangle,1)\
        PARAM(bool, kaskada_redo,0)\
        PARAM(bool, pauli_blocking,1)\
		PARAM(string, formation_zone, "nofz")\
		PARAM(bool, first_step, 1)\
		PARAM(double, step, 0.2)\
		PARAM(int, xsec, 0)\
        PARAM(int, sf_method,0)\
        PARAM(bool, cc_smoothing,1)\
        PARAM(bool, mixed_order,1)\
        PARAM(double,rmin,0)\
        PARAM(double,rmax,0)


using namespace std;

struct line :public string
{
	line(){}
	line(const string&s):string(s)
	{
	};
	line& operator=(const string&s)
	{
		*(string*)this=s;		
		return *this;
	}
};

template <class T>
inline bool read(T &x, stringstream & s,char z)
{
  s >>x;
  return true;
}

inline bool read(vec &a, stringstream & s,char z)
{
  s >>a.x>>a.y>>a.z;   
  return true;
}

inline bool read( line &x , stringstream & s,char z)
{ 
  if(z=='=')
    getline(s,x) ;
  else
    {string y;
     getline(s,y) ;
     x+="\n"+y;
	}
  return true;
}


/// the params class is used to input the run time parameters
/// to the generator from a file and to provide them to the  
/// interested objects in the generator
/// it also helps to store them in the output file for the record
/// the type of the paremeter must be string, vector or  
/// any type with friend operator>>(istream&,type&) operator defined 
/// all the parameters (type, name  are defined here
class params
{//public fields
public:
/// parameter definitions
#define PARAM(type,name,default_value) type name;
PARAMS_ALL() 
#undef PARAM

public:
string path_to_data;
/// constructor
params():path_to_data("")
{
#define PARAM(type,name,default_value) name=default_value;
PARAMS_ALL() 
#undef PARAM
}

/// read all parameter values from config file filename
inline bool read (const char *filename)
{
    ifstream dane;
    open_data_file(dane,filename);
    if (!dane)
      cerr << "Could not open: '" << filename << "' using defaults " << endl;
    return  read(dane,filename);
}

inline bool read (istream& dane,const char *filename)
{
    int line_number = 0;
    string varname, line;
    char eq;
    double value;
    while (dane.good ()) // while not end of file
      {	line_number++;   // increase line number 
	getline (dane, line);  //read next line to buffer "line" 
	stringstream ins (line); // create istream for reading from "line" 
	if (line[0] == '#')      // skip to next line if this line is commented
	  continue;
	if(line[0] == '@')
	  { char malpa;
	    string incfname;
		ins>>malpa>>incfname;
		ifstream incfile;
		if(open_data_file(incfile,incfname))
		  read(incfile,incfname.c_str());  
		else
		  {
	       cerr << filename << ':' << line_number <<". Could not open file '"<<
	       incfname <<"' for reading."<<endl;
		  throw("error");
	      }
		continue;
	  }
	ins >> ws;
	if (ins.peek () == '#')  // comments must start at the beginning of the line  
	  {
	    cerr << filename << ':' << line_number <<
	      ": '#' in the middle of the line " << endl;
	    continue;
	  }
	if (ins.eof ()) // the line was empty
	  continue;
	varname = "";  // erase the name of the variable from the previous line
	eq = ' ';
	while (!ins.eof () && ins.peek () != '=' && ins.peek () != ' '
	       && ins.peek () != '\t')
	  varname += ins.get (); // read up to the next ' ', '\t', '=' or end of line
	ins >> eq;  // read the '=' sign
	if(eq!='=')
	    if(eq=='+' && ins.peek()=='=')
	    {
	    	char z=ins.get();
	    }
		else   
	    {cerr << filename << ':' << line_number <<
	      ": '=' expected after  \"" <<varname<<"\""<< endl;
	    continue;
	    }
// check if the varname is in the fieldlist and call appiopriate read function
#define PARAM(type,name,default_value) ( varname == #name ? ::read(name,ins,eq) : false) ||
    PARAMS_ALL() 
	(cerr << filename << ':' << line_number <<": unknown parameter \""<<varname<<"\""<<endl
	, exit(-1),1)
	; 
#undef PARAM
      }
  }

void list()
{
#define PARAM(type,name,default_value) cout<< #name" = "<< name<<endl; 
    PARAMS_ALL();
#undef PARAM
}
void makeFormParams()
{
		ofstream f;
			f.open("formParam.ini");
#define PARAM(type,name,default_value) \
				f<< "form.name." #name "=" << name <<endl; \
				f<< "form.type." #name "=" << #type <<endl;
			    PARAMS_ALL();
#undef PARAM
			    	f.close();
}

};


#endif
