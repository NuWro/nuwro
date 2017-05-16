#ifndef _INPUT_DATA_h_
#define _INPUT_DATA_h_

#include <string>
#include <sstream>
#include <vector>

#include "params.h"


////////////////////////////////////////
// data_container
////////////////////////////////////////

//! A container for data.
/*! Contains data and the information needed for reading of the input file. */

struct data_container
{
  int     number_of_options;                      //!< Number of possible options within the parameter.
  string     parameter_name;                      //!< Name of the parameter that controls that data.
  string          file_name;                      //!< The name of the file with input data.
  string       *data_fields;                      //!< Possible fields within data file.
  int      number_of_points;                      //!< Number of data points.

  vector<double>     energy;                      //!< Temporary

  data_container( string _parameter_name, int _number_of_options, string &_data_fields );
                                                  //!< The default constructor.
  ~data_container();                              //!< The default destructor.

  void create_data_vector();                      //!< Create vector for data points.
};


////////////////////////////////////////
// input_data
////////////////////////////////////////

//! Manager of the input data needed in NuWro.
/*! This class decides which data files should be read, reads them and stores the data. */

class input_data
{
  params par;                                     //!< Params of the simulation.
  data_container  *cascade_xsec_NN;               //!< Container for data.
  string                input_path;               //!< Path to the folder with input data.
  stringstream        name_sstream;               //!< Stringstream needed for generic names.
  ifstream           file_ifstream;               //!< Stream for handling files.
  string                 file_line;               //!< String for reading lines from files.
  size_t             char_position;               //!< Size_t for manipulating strings;

  public:
    input_data( params _par );                    //!< The default constructor.
                                                  /*!< Receives the params provided. */
    ~input_data();                                //!< The default destructor.
    void initialize();                            //!< Initialize objects, check essential things.
    void load_data();                             //!< Loads the data needed for given simulation.

  private:
    void initialize_input_path();                 //!< Create the input_path and check if it exists.
    void initialize_data_containers();            //!< Prepare data containers for reading files.
    void generate_file_name( data_container &container, int option ); //!< Generate file name.
    void read_data( data_container &container );                      //!< Read and store the input data.
};

#endif