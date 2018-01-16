#ifndef _INPUT_DATA_h_
#define _INPUT_DATA_h_

#include <string>
#include <vector>

#include "data_container.h"
#include "params.h"


//! Manager of the input data needed in NuWro.
/*! This class is resposible for the management of the input data files. It creates proper
    "data_containers" that read and store the data. It supplies them with the name of the file,
    the number of different data fields, their names and the method of interpolation for each. */

class input_data
{
  params                  par;                    //!< Params of the simulation.
  vector<data_container> containers;              //!< Containers for data.
  string                 input_path;              //!< Path to the folder with input data.

  public:
    input_data();                                 //!< The default constructor.
                                                  /*!< Receives the params provided. */
    ~input_data();                                //!< The default destructor.
    void initialize( params _par );               //!< Initializes objects, checks essential things.
    void load_data();                             //!< Loads the data needed for given simulation.
    data_container* get_data_container( int i );  //!< Provides with a specific data_container.

  private:
    void initialize_input_path();                 //!< Creates the input_path and checks if it exists.
    void initialize_data_containers();            //!< Prepares data containers for reading files.
    string generate_file_name( string parameter_name, int parameter_option ); //!< Generates file name.
};

#endif