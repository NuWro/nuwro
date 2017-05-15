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
}

////////////////////////////////////////

input_data::~input_data()
{
  delete cascade_xsec_NN;
}

////////////////////////////////////////

void input_data::initialize()
{
  initialize_input_path();
  initialize_data_containers();
}

////////////////////////////////////////

void input_data::load_data()
{
  read_data( *cascade_xsec_NN );
}


////////////////////////////////////////
// Private methods
////////////////////////////////////////

void input_data::initialize_input_path()
{
  // generate the input_path
  name_sstream.str(string());                   // clear the stringstream
  name_sstream << get_data_dir() << "input/";   // data_dir + relative folder
  input_path = name_sstream.str();

  // check if the directory exists
  DIR* dir = opendir(input_path.c_str());
  if ( dir )                                    // the directory exists
  {
    closedir(dir);
  }
  else
  {
    throw "input_data error: Could not find the input folder.";
  }
}

////////////////////////////////////////

void input_data::initialize_data_containers()
{
  cascade_xsec_NN = new data_container("kaskada_xsec_NN",2);

  if ( par.kaskada_xsec_NN < cascade_xsec_NN->number_of_options )  // if the parameter is ok
  {
    cascade_xsec_NN->file_name = generate_file_name( cascade_xsec_NN->parameter_name, 
                                                     par.kaskada_xsec_NN );
  }
  else
  {
    throw "input_data error: Invalid parameter.";
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

void input_data::read_data( data_container &container )
{
  cout << container.file_name << "\n";
}