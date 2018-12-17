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
  //  0 is taking floor,
  //  1 is linear interpolation,
  // -1 means it is the input axis.

    int    cascade_NN_xsec_number_of_fields     = 3;
    string cascade_NN_xsec_data_fields[]        = {"energy", "xsec_ii", "xsec_ij"};
    int    cascade_NN_xsec_interpolate_fields[] = {-1, 1, 1};
    double cascade_NN_xsec_unit_fields[]        = { 1, millibarn, millibarn};
    containers.push_back( data_container( input_path, "kaskada_NN_xsec", par.kaskada_NN_xsec,
                                          cascade_NN_xsec_number_of_fields, cascade_NN_xsec_data_fields,
                                          cascade_NN_xsec_interpolate_fields, cascade_NN_xsec_unit_fields ));

    int    cascade_NN_inel_number_of_fields     = 5;
    string cascade_NN_inel_data_fields[]        = {"energy", "inel_ii", "inel_ij", "inel_1pii", "inel_1pij"};
    int    cascade_NN_inel_interpolate_fields[] = {-1, 1, 1, 1, 1};
    if( par.kaskada_NN_inel <= 1 )  // for dataset 0 and 1 use no interpolation
    {cascade_NN_inel_interpolate_fields[1] = 0;cascade_NN_inel_interpolate_fields[2]=0;
     cascade_NN_inel_interpolate_fields[3] = 0;cascade_NN_inel_interpolate_fields[4]=0;}
    double cascade_NN_inel_unit_fields[]        = { 1, 1, 1, 1, 1};
    containers.push_back( data_container( input_path, "kaskada_NN_inel", par.kaskada_NN_inel,
                                          cascade_NN_inel_number_of_fields, cascade_NN_inel_data_fields,
                                          cascade_NN_inel_interpolate_fields, cascade_NN_inel_unit_fields ));

    int    cascade_NN_angle_number_of_fields    = 7;
    string cascade_NN_angle_data_fields[]       = {"energy", "angle_A_ii", "angle_A_ij", "angle_B_ii",
                                                   "angle_B_ij", "angle_C_ii", "angle_C_ij"};
    int    cascade_NN_angle_interpolate_fields[]= {-1, 1, 1, 1, 1, 1, 1};
    double cascade_NN_angle_unit_fields[]       = { 1, 1, 1, 1, 1, 1, 1};
    containers.push_back( data_container( input_path, "kaskada_NN_angle", par.kaskada_NN_angle,
                                          cascade_NN_angle_number_of_fields, cascade_NN_angle_data_fields,
                                          cascade_NN_angle_interpolate_fields, cascade_NN_angle_unit_fields ));
}
