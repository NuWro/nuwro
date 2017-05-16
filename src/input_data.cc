#include "input_data.h"

#include <iostream>
#include <fstream>
#include <dirent.h>
#include <math.h>

#include "dirs.h"


////////////////////////////////////////
// data_container
////////////////////////////////////////

data_container::data_container( string _parameter_name, int _number_of_options, int _number_of_fields,
                                string &_data_fields, int &_interpolate_fields ):
                                parameter_name(_parameter_name),
                                number_of_options(_number_of_options),
                                number_of_fields(_number_of_fields)
{
  data_fields         = &_data_fields;
  interpolate_fields  = &_interpolate_fields;
}

////////////////////////////////////////

data_container::~data_container()
{
  delete data_fields;
  delete interpolate_fields;
}

////////////////////////////////////////

void data_container::create_data_vector()
{
  vector<double> empty_array(number_of_fields);                    // create a placeholder for data
  for(int i=0;i<number_of_fields;i++) empty_array[i] = NAN;        // fill it with NANs

  data.reserve(number_of_points);                                  // reserve proper amount of memory
  for(int i=0;i<number_of_points;i++) data.push_back(empty_array); // fill the vector with empty data
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
  stringstream name_sstream;
  name_sstream << get_data_dir() << "input/"; // data_dir + relative folder
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
  // Provide the name of the parameter that governs the data, the number of options,
  // then the number of different fields in the file, their names and the method of interpolation for each.
  int cascade_xsec_NN_number_of_fields     = 2;
  string cascade_xsec_NN_data_fields[]     = {"energy", "xsec_ii"};
  int cascade_xsec_NN_interpolate_fields[] = {0, 0};
  cascade_xsec_NN = new data_container( "kaskada_xsec_NN", 2, cascade_xsec_NN_number_of_fields,
                               *cascade_xsec_NN_data_fields, *cascade_xsec_NN_interpolate_fields );

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
  stringstream name_sstream;
  name_sstream << input_path << container.parameter_name << "_" << option << ".dat"; // path + name + extension
  container.file_name = name_sstream.str();
}

////////////////////////////////////////

void input_data::read_data( data_container &container )
{
  cerr << container.data_fields[0] << "\n";
  ifstream file_ifstream;
  cerr << container.data_fields[0] << "\n";
  cerr << "test\n";
  file_ifstream.open( container.file_name.c_str() );        // open the file

  if( file_ifstream.is_open() )
  {
    string file_line;

    // first check the number of data points
    container.number_of_points = 0;                         // make sure its zero

    while( getline ( file_ifstream, file_line ) )
    {
      if( file_line[0] == '-' )
      {
        container.number_of_points++;
      }
    }
    container.create_data_vector();                         // reserve proper amount of memory
                                                            // and fill the vector with empty data


    file_ifstream.clear();                                  // go back to the start of the file
    file_ifstream.seekg(0, ios::beg);


    // then iterate through points
    int point = -1;                                         // which data point
    size_t char_position;                                   // position of a given char
    string field;                                           // which field
    double value;                                           // what is the data

    while( getline ( file_ifstream, file_line ) )
    {
      if( file_line[0] == '#' )                             // a comment starts with #
      {
        continue;
      }
      if( file_line[0] == '-' )                             // new point starts after -
      {
        point++;
        continue;
      }
      char_position = file_line.find( ':' );                // ":" means data
      if( char_position != string::npos )                   // we have something to take care of
      {
        field = file_line.substr(0, char_position);         // erase everything up to ":"
        field.erase(0, field.find_first_not_of(" \n\r\t") );// trim from left

        value = stod( file_line.substr(char_position+1) );  // everything after ":", convert to double
        //cout << field << " " << value << "\n";
        for( int i=0; i<container.number_of_fields; i++ )   // determine the row and fill
        {
          //cout << container.data_fields[i] << "\n";
          //cout << field << "\n";
          if( container.data_fields[i] == field )
          {
            container.data[point][i] = value;
            cout << "siiiii\n";
            //break;
          }

          //if( i == container.number_of_fields-1 )
          //{
          //  throw "input_data error: Invalid syntax.";
          //}
        }
        //cout << field << " " << value << "\n";
      }
    }

    file_ifstream.close();                                  // close the file
  }
  else
  {
    throw "input_data error: Could not open the data file.";
  }
}