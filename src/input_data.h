#ifndef _INPUT_DATA_h_
#define _INPUT_DATA_h_

#include <string>
#include <sstream>

#include "params.h"

////////////////////////////////////////
// data_container
////////////////////////////////////////

//! A container for data.
/*! Contains data and the information needed for reading of the input file. */

struct data_container
{
  int     number_of_options;                      //!< Number of possible options within the parameter.
  string  parameter_name;                         //!< Name of the parameter that controls that data.
  string  file_name;                              //!< The name of the file with input data.

  data_container( string _parameter_name, int _number_of_options ):
                  parameter_name(_parameter_name), number_of_options(_number_of_options)
                  {};                             //!< The default constructor.
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

  public:
    input_data( params _par );                    //!< The default constructor.
                                                  /*!< Receives the params provided. */
    ~input_data();                                //!< The default destructor.
    bool initialize();                            //!< Initialize objects, check essential things.
    bool load_data();                             //!< Loads the data needed for given simulation.

  private:
    bool   initialize_input_path();               //!< Create the input_path and check if it exists.
    bool   initialize_data_containers();          //!< Prepare data containers for reading files.
    string generate_file_name( string name, int option ); //!< Generate file name.
    bool   read_data( data_container &container );        //!< Read and store the input data.
};

#endif