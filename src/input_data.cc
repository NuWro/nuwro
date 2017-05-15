#include "input_data.h"

#include <iostream>
#include <fstream>
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
  name_sstream.str( string() );                 // clear the stringstream
  name_sstream << get_data_dir() << "input/";   // data_dir + relative folder
  input_path = name_sstream.str();

  // check if the directory exists
  DIR* dir = opendir( input_path.c_str() );
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
    generate_file_name( *cascade_xsec_NN, par.kaskada_xsec_NN );
  }
  else
  {
    throw "input_data error: Invalid parameter.";
  }
}

////////////////////////////////////////

void input_data::generate_file_name( data_container &container, int option )
{
  name_sstream.str( string() );                                                      // clear the stringstream
  name_sstream << input_path << container.parameter_name << "_" << option << ".dat"; // path + name + extension
  container.file_name = name_sstream.str();
}

////////////////////////////////////////

void input_data::read_data( data_container &container )
{
  file_ifstream.open( container.file_name.c_str() );        // open the file

  if( file_ifstream.is_open() )
  {
    file_ifstream.close();                                    // close the file
  }
  else
  {
    throw "input_data error: Could not open the data file.";
  }
}