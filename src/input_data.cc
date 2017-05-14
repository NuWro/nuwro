#include "input_data.h"

#include <iostream>
#include <dirent.h>

#include "dirs.h"

////////////////////////////////////////
// Public methods
////////////////////////////////////////

input_data::input_data( params _par )
{
  par = _par;
  cascade_xsec_NN = new data_container("kaskada_xsec_NN",2);
}

////////////////////////////////////////

input_data::~input_data()
{
  delete cascade_xsec_NN;
}

////////////////////////////////////////

bool input_data::initialize()
{
  if( initialize_input_path() && initialize_data_containers() )
    return 1;
  else
    return 0;
}

////////////////////////////////////////

bool input_data::load_data()
{
  if( read_data( *cascade_xsec_NN ) )
    return 1;
  else
    return 0;
}


////////////////////////////////////////
// Private methods
////////////////////////////////////////

bool input_data::initialize_input_path()
{
  // generate the input_path
  name_sstream.str(string());                  // clear the stringstream
  name_sstream << get_data_dir() << "input/"; // data_dir + relative folder
  input_path = name_sstream.str();

  // check if the directory exists
  DIR* dir = opendir(input_path.c_str());
  if ( dir )                    // the directory exists
  {
    closedir(dir);
    cerr << "test\n";
    return 1;
  }
  else if ( ENOENT == errno )   // the directory does not exist
  {
    return 0;
  }
  else                          // other problem
  {
    return 0;
  }
  return 0;
}

////////////////////////////////////////

bool input_data::initialize_data_containers()
{
  if ( par.kaskada_xsec_NN < cascade_xsec_NN->number_of_options )  // if the parameter is ok
  {
    cascade_xsec_NN->file_name = generate_file_name( cascade_xsec_NN->parameter_name, 
                                                     par.kaskada_xsec_NN );
  }
  else
  {
    cerr << "input_data: unknown kaskada_xsec_NN code " << par.kaskada_xsec_NN << "." << endl;
    exit(29);
  }
}

////////////////////////////////////////

string input_data::generate_file_name( string name, int option )
{
  name_sstream.str(string());                                    // clear the stringstream
  name_sstream << input_path << name << "_" << option << ".dat"; // path + name + extension
  return name_sstream.str();
}

////////////////////////////////////////

bool input_data::read_data( data_container &container )
{
  cout << container.file_name << "\t";
  return 1;
}