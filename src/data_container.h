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
  bool                  nucl_dep;                 //!< Is the container nucleus dependent.
  int                    protons;                 //!< Proton number for nucleus dependent case.
  int                   neutrons;                 //!< Neutron number for nucleus dependent case.
  int           number_of_points;                 //!< Number of data points.
  int           number_of_fields;                 //!< Number of fields within data file.
  vector<string>     data_fields;                 //!< Names of the fields.
  vector<int> interpolate_fields;                 //!< Interpolation type for the fields.
  vector<double>     unit_fields;                 //!< Conversion to natural units.
  vector< vector<double> >  data;                 //!< 2d vector with data.
  int                 input_axis;                 //!< An input axis.
  double             input_point;                 //!< Currently set input point.
  float          input_mid_point;                 //!< Fraction between the data-taking bins.
  int             input_prev_bin;                 //!< An input bin previous to the data-taking point.

  public:
    string            param_name;                 //!< Name of the parameter governing the data.
    int              param_value;                 //!< Value of the parameter governing the data.

    data_container( string _input_path,
                    string _param_name,       int _param_value,
                    int _number_of_fields,    string *_data_fields,
                    int *_interpolate_fields, double *_unit_fields );
                                                  //!< The default constructor.
    data_container( string _input_path,
                    string _param_name,       int _param_value,
                    int _number_of_fields,    string *_data_fields,
                    int *_interpolate_fields, double *_unit_fields,
                    int _protons,             int _neutrons );
                                                  //!< Constructor for nucleus dependent case.
    ~data_container();                            //!< The default destructor.
    void   read_data_file();                      //!< Read and store the data.
    void   set_input_point( double input_value );
                                                  //!< Set the point where the data is taken.
    double get_value( int field );                //!< Get a value (interpolated) for a given field.
    void   debug();                               //!< Write down the contents of the container.

  private:
    void generate_file_name( string input_path ); //!< Generates file name.
    void copy_fields_information( string *_data_fields, int *_interpolate_fields,
                                  double *_unit_fields );
                                                  //!< Copy information about the fields.
    void count_data_points( ifstream &file_ifstream );
                                                  //!< Count number of points in the data file.
    void create_data_vector();                    //!< Create vector for data points.
    void read_data( ifstream &file_ifstream );    //!< Read and store the actual data.
    void fill_nan_data();                         //!< Interpolate (or not) the missing data.
};

////////////////////////////////////////

//! Contain rules used to sort the data.
/*! Allows for the generic sorting in terms of "input_axis". */

class data_compare
{
    int input_axis;
  public:
    data_compare( int _input_axis ) : input_axis(_input_axis) {}
                                                  //!< The default constructor.

    bool operator()( const vector<double>& row1, const vector<double>& row2 )
    {
      return row1[input_axis] < row2[input_axis];
    }
                                                  //!< The sorting function that uses the "input_axis".
};

#endif