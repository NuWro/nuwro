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
#include "params_all.h"		// auto_generated (from params.xml)

using namespace std;

struct line:public string
{
  line ()
  {
  }
  line (const string & s):string (s)
  {
  };
  line & operator= (const string & s)
  {
    *(string *) this = s;
    return *this;
  }
};

template < class T > inline bool read (T & x, istream & s, char z)
{
  s >> x;
  return z == '=' && s;		// '+=' allowed only for lines
}


template < class T > inline void
write (T x, ostream & s)
{
  s << ' ' << x;		// space looks nice ofter in "a = 1"
}


inline bool
read (vec & a, istream & s, char z)
{
  s >> a.x >> a.y >> a.z;
  return z == '=' && s;		// '+=' allowed only for lines
}


inline void
write (vec a, ostream & s)
{
  s << ' ' << a.x << ' ' << a.y << ' ' << a.z;
}

inline bool
read (string & a, istream & s, char z)
{
  s >> a;
  if (!s)
    a = "";
  return true;
}

inline bool
read (line & x, istream & s, char z)
{
  switch (z)
    {
    case '=':
      getline (s, x);
      return true;
    case '+':
      string y;
      getline (s, y);
      x += "\n" + y;
      return true;
    }
  return false;
}

inline void
write (line a, ostream & s)
{
  s << a;
}


/// the params class is used to input the run time parameters
/// to the generator from a file and to provide them to the
/// interested objects in the generator
/// it also helps to store them in the output file for the record
/// the type of the paremeter must be string, vector or
/// any type with friend operator>>(istream&,type&) operator defined
/// all the parameters (type, name  are defined here
class params
{				//public fields
public:

  /// parameter definitions
#define PARAM(type,name,default_value) type name;
  PARAMS_ALL ()
#undef PARAM
public:
  string path_to_data;

  /// constructor
  params ():
#define PARAM(type,name,default_value) name(default_value),
  PARAMS_ALL ()
#undef PARAM
	path_to_data ("")
  {
  }

  /// read all parameter values from config file filename
  inline bool read (const char *filename)
  {
    ifstream dane;
    open_data_file (dane, filename);
    if (!dane)
      cerr << "Could not open: '" << filename << "' using defaults " << endl;
    return read (dane, filename);
  }

  inline bool read (istream & dane, const char *filename)
  {
    int line_number = 0;
    double value;
    while (dane.good ())	// while not end of file
      {				// increase line number
	line_number++;
	//read next line to buffer "line"
	string line;
	getline (dane, line);
	// create istream for reading from "line"
	stringstream ins (line);
	// skip to next line if this line is commented
	if (line[0] == '#')
	  continue;
	if (line[0] == '@')
	  {
	    char malpa;
	    string incfname;
	    ins >> malpa >> incfname;
	    ifstream incfile;
	    if (open_data_file (incfile, incfname))
	      read (incfile, incfname.c_str ());
	    else
	      {
		cerr << filename << ':' << line_number <<
		  ". Could not open file '" << incfname << "' for reading." <<
		  endl;
		exit (11);
		//throw("error");

	      }
	    continue;
	  }
	ins >> ws;
	// comments must start at the beginning of the line
	if (ins.peek () == '#')
	  {
	    cerr << filename << ':' << line_number <<
	      ": '#' in the middle of the line " << endl;
	    continue;
	  }
	if (ins.eof ())		// the line was empty
	  continue;
	string varname = "";
	while (!ins.eof ())
	  {
	    char b = ins.peek ();
	    if (b == '=' || b == '+' || b == ' ' || b == '\t')
	      break;
	    else
	      varname += ins.get ();
	  }
	char eq = ' ';
	ins >> eq;		// read the '=' sign (or '+' sign)
	switch (eq)
	  {
	  case '=':
	    break;
	  case '+':
	    if (ins.peek () == '=')
	      {
		ins.get ();
		break;
	      }
	  default:
	    cerr << filename << ':' << line_number <<
	      ": '=' expected after  \"" << varname << "\"" << endl;
	    exit (13);
	  }
	// check if the varname is in the fieldlist and call appiopriate read function
	int kod = 0;
#define PARAM(type,name,default_value) \
					if(kod==0) \
						if(varname == #name)\
							kod=1+::read(name,ins,eq);
	PARAMS_ALL ()
#undef PARAM
	  switch (kod)
	  {
	  case 0:
	    cerr << filename << ':' << line_number << ": unknown parameter \""
	      << varname << "\"" << endl;
	    exit (14);
	  case 1:
	    cerr << filename << ':' << line_number <<
	      ": bad value for parameter \"" << varname << "\"" << endl;
	    exit (15);
	  }
      }
    return true;
  }

  void list (ostream & out = cout)
  {
#define PARAM(type,name,default_value) \
				out<<"        "<< #name" ="; \
				write(name,out); \
				out<<endl;
    PARAMS_ALL ();
#undef PARAM
  }

  void list (string filename)
  {
    std::ofstream out (filename.c_str ());
    list (out);
  }

};

#endif
