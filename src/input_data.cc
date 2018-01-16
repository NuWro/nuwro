#include "input_data.h"

#include <iostream>
#include <fstream>
#include <dirent.h>
#include <math.h>

#include "dirs.h"
#include "jednostki.h"


////////////////////////////////////////
// Public methods
////////////////////////////////////////

input_data::input_data()
{
}

////////////////////////////////////////

input_data::~input_data()
{
}

////////////////////////////////////////

void input_data::initialize( params _par )
{
  par = _par;
  initialize_input_path();
  initialize_data_containers();
}

////////////////////////////////////////

void input_data::load_data()
{
  for( int i = 0; i < containers.size(); i++ )
    containers[i].read_data_file();
}

////////////////////////////////////////

data_container* input_data::get_data_container( int i )
{
  return &containers[i];
}


////////////////////////////////////////
// Private methods
////////////////////////////////////////

void input_data::initialize_input_path()
{
  // generate the input_path
  stringstream name_sstream;
  name_sstream << get_data_dir() << "input/";   // data_dir + relative folder
  input_path = name_sstream.str();

  // check if the directory exists
  DIR* dir = opendir( input_path.c_str() );
  if ( dir )                                    // the directory exists
  {
    closedir( dir );
  }
  else
  {
    throw "input_data error: Could not find the input folder.";
  }
}

////////////////////////////////////////

void input_data::initialize_data_containers()
{
  // Provide the parameter that governs the data,
  // then the number of different fields in the file, their names and the method of interpolation:
  //  0 is no interpolation,
  //  1 is linear interpolation,
  // -1 means it is the input axis.

    string cascade_xsec_NN_file_name            = generate_file_name( "kaskada_xsec_NN", par.kaskada_xsec_NN );
    int    cascade_xsec_NN_number_of_fields     = 10;
    string cascade_xsec_NN_data_fields[]        = {"energy",
                                                   "xsec_ii", "xsec_ij",
                                                   "inel_ii", "inel_ij", "inel_1pi",
                                                   "angle_A_ii", "angle_A_ij", "angle_B_ii","angle_B_ij"};
    int    cascade_xsec_NN_interpolate_fields[] = {-1, 1, 1, 0, 0, 0, 0, 0, 0, 0};
    double cascade_xsec_NN_unit_fields[]              = {1, fermi, fermi, 1, 1, 1, 1, 1, 1, 1};
    containers.push_back( data_container( cascade_xsec_NN_file_name,   cascade_xsec_NN_number_of_fields,
                                          cascade_xsec_NN_data_fields, cascade_xsec_NN_interpolate_fields,
                                          cascade_xsec_NN_unit_fields ));
}

////////////////////////////////////////

string input_data::generate_file_name( string parameter_name, int parameter_option )
{
  stringstream name_sstream;
  name_sstream << input_path << parameter_name << "_" << parameter_option << ".dat"; // path + name + extension
  return name_sstream.str();
}
