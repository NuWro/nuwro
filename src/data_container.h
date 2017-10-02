#ifndef _DATA_CONTAINER_h_
#define _DATA_CONTAINER_h_

#include <string>
#include <sstream>
#include <vector>
#include <fstream>

#include "dirs.h"

//! A container for data.
/*! Reads and stores the data. */

class data_container
{
  string               file_name;                 //!< The name of the file with input data.
  int           number_of_points;                 //!< Number of data points.
  int           number_of_fields;                 //!< Number of fields within data file.
  vector<string>     data_fields;                 //!< Names of the fields.
  vector<int> interpolate_fields;                 //!< Interpolation type for the fields.
  vector< vector<double> >  data;                 //!< 2d vector with data.
  int                 input_axis;                 //!< An input axis.
  int             input_prev_bin;                 //!< An input bin previous to the data-taking point.
  float          input_mid_point;                 //!< Fraction between the data-taking bins.

  public:
    data_container( string _file_name, int _number_of_fields,
                  string *_data_fields, int *_interpolate_fields );
                                                  //!< The default constructor.
    ~data_container();                            //!< The default destructor.
    void read_data_file();                        //!< Read and store the data.
    void set_input_point( double input_value );
                                                  //!< Set the point where the data is taken.
    double get_value( int field );                //!< Get a value (interpolated) for a given field.

  private:
    void copy_fields_information( string *_data_fields, int *_interpolate_fields );
                                                  //!< Copy information about the fields.
    void count_data_points( ifstream &file_ifstream );
                                                  //!< Count number of points in the data file.
    void create_data_vector();                    //!< Create vector for data points.
    void read_data( ifstream &file_ifstream );    //!< Read and store the actual data.
};

#endif