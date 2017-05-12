#ifndef _INPUT_DATA_h_
#define _INPUT_DATA_h_

#include "params.h"

//! Manager of the input data needed in NuWro.
/*! This class decides which data files should be read, reads them and stores the data */


class input_data
{
  params par;

  public:
    input_data(params _par);
    ~input_data();
    int load_data();

  private:
    
};

#endif