#include "input_data.h"

#include <iostream>
#include <fstream>
#include <dirent.h>

#include "dirs.h"


////////////////////////////////////////
// data_container
////////////////////////////////////////

data_container::data_container( string _parameter_name, int _number_of_options ):
                                parameter_name(_parameter_name),
                                number_of_options(_number_of_options)
{}

////////////////////////////////////////

data_container::~data_container()
{
}

////////////////////////////////////////

void data_container::create_data_vector()
{
  energy.reserve(number_of_points);
}


////////////////////////////////////////
// input_data
////////////////////////////////////////

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
  cascade_xsec_NN = new data_container( "kaskada_xsec_NN", 2 );

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
    container.number_of_points = 0;                         // make sure its zero
    while( getline ( file_ifstream, file_line ) )           // first check the number of data points
    {
      if( file_line[0] == '-' )
      {
        container.number_of_points++;
      }
    }
    container.create_data_vector();                         // reserve proper amount of memory

    file_ifstream.clear();                                  // go back to the start of the file
    file_ifstream.seekg(0, ios::beg);
    while( getline ( file_ifstream, file_line ) )           // first check the number of data points
    {
      if( file_line[0] == '#' )                             // a comment starts with #
       continue;
      if( file_line[0] == '-' )                             // new point starts after -
      {
        getline ( file_ifstream, file_line );
        if( file_line.find('energy') )
        {
          char_position = file_line.find(':');              // find ":"
          file_line = file_line.substr(char_position+1);    // erase everything up to ":"
          container.energy.push_back( stod(file_line) );    // convert to double and add the data point
        }
      }
    }
    cout << container.energy[2] << "\n";
    cout << container.energy[9] << "\n";
    file_ifstream.close();                                  // close the file
  }
  else
  {
    throw "input_data error: Could not open the data file.";
  }
}