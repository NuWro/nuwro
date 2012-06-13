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

//#include "params_all_orig.h" // original params (hand written)
#include "params_all.h"  // auto_generated (from params.xml)

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
inline bool read(T &x, istream & s,char z)
{
  s >>x;
  return true;
}

template <class T>
inline void write(T x, ostream & s)
{
  s <<x;
}

inline bool read(vec &a, istream & s,char z)
{
  s >>a.x>>a.y>>a.z;   
  return true;
}

inline void write(vec a, ostream & s)
{
  s <<a.x<<' '<<a.y<<' '<<a.z;   
}

inline bool read( line &x , istream & s,char z)
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
	#define PARAM(type,name,default_value) ,name(default_value)
	PARAMS_ALL() 
	#undef PARAM
	{
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
	#define PARAM(type,name,default_value) cout<< #name" = "; write(name,cout); cout<<endl; 
		PARAMS_ALL();
	#undef PARAM
	}

};


#endif
