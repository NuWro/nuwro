#include "data_container.h"

#include <iostream>
#include <fstream>
#include <dirent.h>
#include <math.h>

////////////////////////////////////////
// Public methods
////////////////////////////////////////

data_container::data_container( string _file_name, int _number_of_fields,
                                string *_data_fields, int *_interpolate_fields ):
                                file_name(_file_name),
                                number_of_fields(_number_of_fields)
{
  copy_fields_information( _data_fields, _interpolate_fields ); // copy the information to vectors
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

    file_ifstream.close();                        // close the file
  }
  else
  {
    throw "input_data error: Could not open the data file.";
  }
}


////////////////////////////////////////
// Private methods
////////////////////////////////////////

void data_container::copy_fields_information( string *_data_fields, int *_interpolate_fields)
{
  data_fields        = vector<string> (_data_fields, _data_fields + number_of_fields);
  interpolate_fields = vector<int> (_interpolate_fields, _interpolate_fields + number_of_fields);
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
      field = file_line.substr(0, char_position);           // erase everything up to ":"
      field.erase(0, field.find_first_not_of(" \n\r\t") );  // trim from left

      value = stod( file_line.substr(char_position+1) );    // everything after ":", convert to double

      for( int i=0; i<number_of_fields; i++ )               // determine the row and fill
      {
        if( data_fields[i] == field )
        {
          data[point][i] = value;
          break;
        }

        if( i == number_of_fields-1 )                       // couldn't find the field in data_fields
        {
          throw "input_data error: Invalid syntax.";
        }
      }
    }
  }

  file_ifstream.clear();                                  // go back to the start of the file
  file_ifstream.seekg(0, ios::beg);
}