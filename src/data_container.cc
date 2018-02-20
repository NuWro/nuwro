#include "data_container.h"

#include <iostream>
#include <fstream>
#include <dirent.h>
#include <math.h>
#include <algorithm>


////////////////////////////////////////
// Public methods
////////////////////////////////////////

data_container::data_container( string _input_path,
                                string _param_name,       int _param_value,
                                int _number_of_fields,    string *_data_fields,
                                int *_interpolate_fields, double *_unit_fields ):
                                param_name(_param_name), param_value(_param_value),
                                number_of_fields(_number_of_fields)
{
  generate_file_name( _input_path );          // generate file name using the provided params and the input path
  copy_fields_information( _data_fields, _interpolate_fields, _unit_fields ); // copy the information to vectors
}

////////////////////////////////////////

data_container::~data_container()
{
}

////////////////////////////////////////

void data_container::read_data_file()
{
  ifstream file_ifstream;
  file_ifstream.open( file_name.c_str() );        // open the file

  if( file_ifstream.is_open() )                   // if the file exists
  {
    count_data_points( file_ifstream );           // count data points in the file

    create_data_vector();                         // reserve proper amount of memory
                                                  // and fill the vector with empty data

    read_data( file_ifstream );                   // read and store the actual data

    fill_nan_data();                              // fill in the points with no data

    input_point = data[0][input_axis];            // reset input points
    input_mid_point = 0;

    file_ifstream.close();                        // close the file
  }
  else
  {
    throw "input_data error: Could not open the data file. Make sure the params are correct and files for given options exist.";
  }
}

////////////////////////////////////////

void data_container::set_input_point( double input_value )
{
  if( fabs(input_value - input_point) > 1e-10 ) // do it if the input value changed, fixed epsilon 1e-10
  {
    if( input_value < data[0][input_axis] || input_value > data[number_of_points-1][input_axis] )
    {
      throw "input_data error: Cannot set the data-taking point, out of bounds.";
    }

    input_prev_bin = 0;
    while( input_value >= data[input_prev_bin+1][input_axis] )
    {
          input_prev_bin++;
    }

    input_point = input_value;
    input_mid_point = (input_value - data[input_prev_bin][input_axis])
                    / (data[input_prev_bin+1][input_axis] - data[input_prev_bin][input_axis]);
  }
}

////////////////////////////////////////

double data_container::get_value( int field )
{
  switch( interpolate_fields[field] )
  {
    case 0:                                       // taking floor
    {
      return data[input_prev_bin][field];
    }
    case 1:                                       // linear interpolation
    {
      return (1-input_mid_point)*data[input_prev_bin][field]
               +input_mid_point *data[input_prev_bin+1][field];
    }

    default:
    {
      if( interpolate_fields[field] < 0 )
      {
        throw "input_data error: Cannot return the value of the input axis.";
      }
      else
      {
        throw "input_data error: Cannot determine the interpolation type.";
      }
    }
  }
}

////////////////////////////////////////

void data_container::debug()
{
  cout << "\n# data_container: "
       << param_name << "_" << param_value << "\n";
  for(int i=0; i<number_of_points; i++)
  {
    cout << "# ";
    for(int j=0; j<number_of_fields; j++)
    {
      cout << data[i][j] << "\t";
    }
    cout << "\n";
  }
  cout << "\n";
}


////////////////////////////////////////
// Private methods
////////////////////////////////////////

void data_container::generate_file_name( string input_path )
{
  stringstream name_sstream;
  name_sstream << input_path << param_name << "_" << param_value << ".dat"; // path + name + extension
  file_name = name_sstream.str();
}

////////////////////////////////////////

void data_container::copy_fields_information( string *_data_fields, int *_interpolate_fields, double *_unit_fields)
{
  data_fields        = vector<string> (_data_fields, _data_fields + number_of_fields);
  interpolate_fields = vector<int> (_interpolate_fields, _interpolate_fields + number_of_fields);
  unit_fields        = vector<double> (_unit_fields, _unit_fields + number_of_fields);

  int check_minus = 0;                               // find the input axis and check if there is only one
  for( int i=0; i<interpolate_fields.size(); i++)
  {
    if( interpolate_fields[i] < 0 )
    {
      input_axis = i;
      check_minus++;
    }
  }
  if( check_minus != 1 )
  {
    throw "input_data error: Cannot specify a single input axis.";
  }
}

////////////////////////////////////////

void data_container::count_data_points( ifstream &file_ifstream )
{
  string file_line;
  number_of_points = 0;                           // make sure its zero

  while( getline ( file_ifstream, file_line ) )   // search the whole file
  {
    if( file_line[0] == '-' )                     // a line that starts with "-" gives a new point
    {
      number_of_points++;
    }
  }

  file_ifstream.clear();                          // go back to the start of the file
  file_ifstream.seekg(0, ios::beg);
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

void data_container::read_data( ifstream &file_ifstream )
{
  string file_line;
  int point = -1;                                 // which data point
  size_t char_position;                           // position of a given char
  string field;                                   // which field
  double value;                                   // what is the data

  while( getline ( file_ifstream, file_line ) )
  {
    if( file_line[0] == '#' )                     // a comment starts with #
    {
      continue;
    }

    if( file_line[0] == '-' )                     // new point starts after -
    {
      point++;
      continue;
    }

    char_position = file_line.find( ':' );        // ":" means data
    if( char_position != string::npos )           // we found ":"
    {
      field = file_line.substr(0, char_position);           // take everything up to ":"
      field.erase(0, field.find_first_not_of(" \n\r\t") );  // trim from left

      value = stod( file_line.substr(char_position+1) );    // everything after ":", convert to double

      for( int i=0; i<number_of_fields; i++ )               // determine the row and fill
      {
        if( data_fields[i] == field )
        {
          data[point][i] = value * unit_fields[i];          // save the value and convert to natural units
          break;
        }

        if( i == number_of_fields-1 )                       // couldn't find the field in data_fields
        {
          throw "input_data error: Invalid syntax.";
        }
      }

      if( ::isnan(data[point][input_axis]) )
      {
        throw "input_data error: Point specified without a value on the input axis.";
      }
    }
  }

  file_ifstream.clear();                              // go back to the start of the file
  file_ifstream.seekg(0, ios::beg);

  // sort the data
  std::sort(data.begin(), data.end(), data_compare(input_axis));
}

////////////////////////////////////////

void data_container::fill_nan_data()
{
  // fill the back and front if there are nans
  for( int field = 0; field < number_of_fields; field++ )
  {
    if( field != input_axis )                         // for data fields only
    {
      for( int point_up = 0; point_up < number_of_points; point_up++ )
                                                      // scan upwards
      {
        if( !::isnan(data[point_up][field]) )           // find the first one that is a number
        {
          for( int nan = point_up; nan >= 0; nan-- )  // fill the previous ones with this number
          {
            data[nan][field] = data[point_up][field];
          }
          break;
        }
        if( point_up == number_of_points-1 )          // if only nans in a field
        {
          throw "input_data error: One of the fields has no data.";
        }
      }
      for( int point_down = number_of_points-1; point_down >= 0; point_down-- )
                                                      // scan downwards
      {
        if( !::isnan(data[point_down][field]) )         // find the last one that is a number
        {
          for( int nan = point_down+1; nan < number_of_points; nan++ )
                                                      // fill the next ones with this number
          {
            data[nan][field] = data[point_down][field];
          }
          break;
        }
      }
    }
  }

  // interpolate nans in the middle
  float mid_point;                                    // relative position in the interpolation
  for( int field = 0; field < number_of_fields; field++ )
  {
    if( field != input_axis )                         // for data fields only
    {
      for( int nan_up = 0; nan_up < number_of_points; nan_up++ )
                                                      // scan upwards
      {
        if( ::isnan(data[nan_up][field]) )            // find the first nan
        {
          for( int point_up = nan_up+1; point_up < number_of_points; point_up++ )
                                                      // scan for the next number
          {
            if( !::isnan(data[point_up][field]) )     // find the next number
            {
              for ( int nan = nan_up; nan < point_up; nan++ )
                                                      // for every nan in a sequence
              {
                // now interpolate
                mid_point = ( data[nan][input_axis] - data[nan_up-1][input_axis] )
                          / ( data[point_up][input_axis] - data[nan_up-1][input_axis] );
                data[nan][field] = (1-mid_point)*data[nan_up-1][field]
                                 +    mid_point *data[point_up][field];
              }
              break;
            }
          }
        }
      }
    }
  }
}