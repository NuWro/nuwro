#include "lepevent.h"
#include "params.h"
#include "event1.h"
#include "pdg.h"
#include "nu_e_el_sigma.h"
#include "jednostki.h" // G, sin2thetaW, cm2
//#include "generatormt.h" // frand

double lepevent(params& p, event& e) //, bool cc)
{// scattering of neutrino on electron at rest
  e.weight = 0;
  particle neut_in = e.in[0]; // e.in[0] is neutrino
  particle elec_in = e.in[1]; // e.in[1] is electron

  particle lept_out(PDG::pdg_e, PDG::mass_e); // for now final charged lepton = electron
  particle neut_out = neut_in; // for now neut_in = neut_out


  // generate kinematics
  int kind = neut_in.pdg;  // initial neutrino type: nu_e -> 12, nu_mu -> 14, nu_tau -> 16, anu_e -> -12,...
  double E0nu = neut_in.t; // incoming neutrino energy in target electron frame	
  double me = mass_e;
  double threshold_mu = (mass_mu*mass_mu-me*me)/2.0/me; // threshold to produce muon
  double threshold_tau = (mass_tau*mass_tau-me*me)/2.0/me; // threshold to produce tau lepton
  double m_prime;  // final charged lepton mass
  double delta_m2; // m_prime*m_prime - me*me;
  double t_ratio;  // final-to-initial neutrino energy ratio

  int switch_sigma = 1; //1 for now 
  // 1 - single subchannel, 2 - second option in double subchannel, 3 - third option in triple subchannel

  int switch_sigma_max = 1;
  if ((kind == 14 && E0nu>threshold_mu) || (kind == 16 && E0nu>threshold_tau) || (kind == -12 && E0nu>threshold_mu && E0nu<=threshold_tau))
    switch_sigma_max = 2;
  else if (kind == -12 && E0nu>threshold_tau)
    switch_sigma_max = 3;

  double rand1 = frandom();
  double rand2 = frandom();

double weight_tot = 0; // weight of given (multiple) channel
double weight_[switch_sigma_max];
for (int ii = 1; ii<=switch_sigma_max; ii++) // ii -> switch_sigma
{
  out_state(kind, ii, neut_in, m_prime, lept_out, neut_out);
  delta_m2 = m_prime*m_prime - me*me;

  if (ii==1)
    t_ratio = me/(me+2.0*E0nu) + 2.0*E0nu/(me+2.0*E0nu)*rand1;
  else if (ii==2 || ii==3)
    t_ratio = (me/(me+2.0*E0nu) - delta_m2/2.0/(me+2.0*E0nu)/E0nu) 
            + (2.0*E0nu/(me+2.0*E0nu) - delta_m2/me/(me+2.0*E0nu))*rand1;

  weight_[ii] = nu_e_el_sigma (E0nu, t_ratio, kind, neut_in.pdg<0, elec_in.mass(), ii, m_prime, delta_m2 );
  weight_tot = weight_tot + weight_[ii];
}

// *** Subchannel chooser *****
  double rand_chooser = frandom();
  if (switch_sigma_max == 2)
  {
    if (rand_chooser <= weight_[1]/weight_tot) //(weight_[1]+weight_[2]))
      switch_sigma = 1;
    else switch_sigma = 2;
  }
  else if (switch_sigma_max == 3)
  {
    if (rand_chooser <= weight_[1]/weight_tot) //(weight_[1]+weight_[2]+weight_[3]))
      switch_sigma = 1;
    else if (rand_chooser > (weight_[1]+weight_[2])/weight_tot)
      switch_sigma = 3;
    else switch_sigma = 2;
  }
// ****************************

  out_state(kind, switch_sigma, neut_in, m_prime, lept_out, neut_out);
  
  delta_m2 = m_prime*m_prime - me*me;

  if (switch_sigma==1)
    t_ratio = me/(me+2.0*E0nu) + 2.0*E0nu/(me+2.0*E0nu)*rand1;
  else if (switch_sigma==2 || switch_sigma==3)
    t_ratio = (me/(me+2.0*E0nu) - delta_m2/2.0/(me+2.0*E0nu)/E0nu) 
            + (2.0*E0nu/(me+2.0*E0nu) - delta_m2/me/(me+2.0*E0nu))*rand1;


	double Enu = E0nu*t_ratio;
  double phi = 2.0*Pi*rand2;
  double Ep = E0nu - Enu + me;
	double p_out = sqrt(Ep*Ep - m_prime*m_prime);

	//double cos_theta1 = (E0nu*Enu + me*(Enu-E0nu))/E0nu/Enu;
  double cos_theta1 = (2.0*E0nu*Enu + 2.0*me*(Enu-E0nu) + (m_prime*m_prime-me*me))/2.0/E0nu/Enu;
	double sin_theta1 = sqrt(1.0 - cos_theta1*cos_theta1);


// final neutrino energy in target electron frame
  neut_out.t = Enu;
// final neutrino momentum
  double xx = Enu*sin_theta1*cos(phi);
  double yy = Enu*sin_theta1*sin(phi);
//	double k_z = Enu * cos_theta1;
  vec P_3=vec(xx, yy, Enu*cos_theta1);
  neut_out.set_momentum (P_3);

// final charged lepton energy
  lept_out.t = Ep;
// final charged lepton momentum
//	double p_x = - p_out * sin_theta2 * cos(phi); // - k_x
//	double p_y = - p_out * sin_theta2 * sin(phi); // - k_y
//	double p_z =   p_out * cos_theta2;
//	double cos_theta2 = (E0nu + me)/E0nu / sqrt(1.0 + 2.0*me/(E0nu-Enu));
//	double sin_theta2 = sqrt(1.0 - cos_theta2*cos_theta2);
//  vec P_4=vec(-sin_theta2*cos(phi), -sin_theta2*sin(phi), cos_theta2)*p_out;
  vec P_4=vec(-xx, -yy, E0nu-Enu*cos_theta1);
  lept_out.set_momentum (P_4);


  e.weight = weight_tot*p.nucleus_p/(p.nucleus_p+p.nucleus_n);
  //p.target_content

    e.out.push_back(neut_out);
    e.out.push_back(lept_out);

  return e.weight;
}


// Returning multiple values of out-state, Using References
void out_state(int kind, int switch_sigma, particle neut_in, double &m_prime_, particle &lept_out_, particle &neut_out_) 
{ 
      // E0nu>=threshold_mu
  if(kind == 14 && switch_sigma == 2)
  {
    m_prime_ = mass_mu;
    lept_out_ = particle(PDG::pdg_mu, PDG::mass_mu); // final charged lepton = muon
    neut_out_ = particle(PDG::pdg_nu_e, PDG::mass_nu_e); // final neutrino = nu_e
  }
  else if(kind == -12 && switch_sigma == 2)
  {
    m_prime_ = mass_mu;
    lept_out_ = particle(PDG::pdg_mu, PDG::mass_mu); // final charged lepton = muon
    neut_out_ = particle(-PDG::pdg_nu_mu, PDG::mass_nu_mu); // final neutrino = anu_mu
  }
  // E0nu>=threshold_tau
  else if (kind == 16 && switch_sigma == 2)
  {
    m_prime_ = mass_tau;
    lept_out_ = particle(PDG::pdg_tau, PDG::mass_tau); // final charged lepton = tau lepton
    neut_out_ = particle(PDG::pdg_nu_e, PDG::mass_nu_e); // final neutrino = nu_e
  }
  else if (switch_sigma == 3)
  {
    m_prime_ = mass_tau;
    lept_out_ = particle(PDG::pdg_tau, PDG::mass_tau); // final charged lepton = tau lepton
    neut_out_ = particle(-PDG::pdg_nu_tau, PDG::mass_nu_tau); // final neutrino = anu_tau
  }
  // no energy threshold
  else 
  {
    m_prime_ = mass_e;
    lept_out_ = particle(PDG::pdg_e, PDG::mass_e); // final charged lepton = electron
    neut_out_ = neut_in; // final neutrino = initial neutrino
  } 
} 