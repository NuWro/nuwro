#ifndef _INPUT_DATA_h_
#define _INPUT_DATA_h_

#include <string>

#include "params.h"


struct data_holder
{
  string file_name;
};

//! Manager of the input data needed in NuWro.
/*! This class decides which data files should be read, reads them and stores the data */

class input_data
{
  params par;                                     //!< Params of the simulation.
  data_holder cascade_NN_xsec;                    //!< Container for data.

  public:
    input_data(params _par);                      //!< The default constructor.
                                                  /*!< Receives the params provided. */
    ~input_data();                                //!< The default destructor.
    int load_data();                              //!< Loads the data needed for given simulation.

  private:
    int init();
};

#endif