#ifndef _E_EL_h_
#define _E_EL_h_

double e_el_sigma(double Enu, double q2, 
     int kind, 
     bool anty, 
     double m, 
     double M, 
     int what=0  /// what to return
);

/*
 *  The formula for d sigma / d Q2 as obtaine by K.M.G.
 */
 

double dsigma_dOmega_ksiazka(double En, double q2, 					
    int kind,  ///< process type: 10 - ep, 11 - en elastic scattering
    bool anty, ///< true for positrons
    double m,  ///< lepton mass
    double M   ///< nucleon (effective) mass
    );

#endif
