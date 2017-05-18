#ifndef _DATA_CONTAINER_h_
#define _DATA_CONTAINER_h_

#include <string>
#include <sstream>
#include <vector>

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

  public:
    data_container( string _parameter_name, int _number_of_fields,
                  string *_data_fields, int *_interpolate_fields );
                                                  //!< The default constructor.
    ~data_container();                            //!< The default destructor.
    void read_data();                             //!< Read and store the data.

  private:
    void create_data_vector();                    //!< Create vector for data points.
};

#endif