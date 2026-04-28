#ifndef _shell_sampler_h_
#define _shell_sampler_h_

#include <vector>
#include <string>

#include "generatormt.h"


class shell_sampler
{
  public:
    shell_sampler(const std::string& filename);  //!< Constructor that loads the shell data

    double r();                                  //!< Sample a random radius according to the distribution
    double dens(double r) const;                 //!< Density at point r

  private:
    struct point {
        double r;     //!< Radius
        double rho;   //!< Density
    };
    double norm;                //!< Normalization

    std::vector<point> data;    //!< Loaded shell data
    std::vector<double> cdf;    //!< Cumulative distribution function for sampling

    void load_data(const std::string& filename);  //!< Load the shell data from file
    void build_cdf();           //!< Build the cumulative distribution function
};

#endif
