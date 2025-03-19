#include "mecevent_Nieves.h"
#include "mecevent_2020Valencia.h"
#include "generatormt.h"
#include "mecevent_common.h"
#include "vecrand.h"
#include <thread>


bool flag_2p2h = true;
bool flag_2p2h_pn = false;
bool flag_3p3h = !flag_2p2h;


int ConvertCoordinate_2D_to_1D(int n, int m) {    
  // Convert the 2D co-ordinates of a N x N upper triangular matrix to a 1D matrix (which stores only non zero elements)
  // formula = n*(2N -n + 1 )/2  + m-n  where N is the dimension of the upper triangular matrix 
  return (    (int)( n*(240 - n + 1)/2 + m-n)  );
}

double Bilinear_Interpolation(double x, double y, double x1, double x2, double y1, double y2, double f11, double f12, double f21, double f22) {
  double spacing_x = x2 - x1;
  double spacing_y = y2 - y1;

  double r1 = (x2 - x)*f11/spacing_x + (x-x1)*f12/spacing_x;
  double r2 = (x2 - x)*f21/spacing_x + (x-x1)*f22/spacing_x;

  double value = (y2 - y)*r1/spacing_y + (y - y1)*r2/spacing_y;
  return value;
}

void mecevent_2020Valencia (params & p, event & e, nucleus & t, bool cc)
{
  e.par = p;            // save params in the event
  e.flag.cc = cc;       // set flags for the event
  e.flag.nc = !cc;
  e.flag.dis = false;
  e.flag.qel = false;
  e.flag.coh = false;
  e.flag.mec = true;
  int ile_pb = p.mec_pb_trials;
  double mec_central = p.mec_central_motion;
  double mec_smearing = p.mec_back_to_back_smearing;
  double binding = p.kaskada_w;
  int mecskalowanie = p.mec_scaling;

  double mc_sampling[3] = {0,0,0};          //sampling direction for pp, np, and pn respectively
  int mc_strength[3] = {1,1,1};      //strength for sampling direction pp, np, and pn respectively
  mc_sampling[0] = p.MEC_cm_direction_pp;
  mc_strength[0] = p.MEC_cm_strength_pp;

  mc_sampling[1] = p.MEC_cm_direction_np;
  mc_strength[1] = p.MEC_cm_strength_np;
  
  mc_sampling[2] = p.MEC_cm_direction_pn;
  mc_strength[2] = p.MEC_cm_strength_pn;
 

  // sadly, only CC events available so far...
  if(e.flag.nc)
  {
    cerr<<" MEC error: Wrong Settings!\n";
    e.weight = 0;
    return;
  }

  particle meclepton;
  particle meclepton_3p3h;
  ap=(e.in[0].pdg<0);

  if ( std::abs(e.in[0].pdg) != 14 ) {
    cerr << "Only muon neutrino available \n";
    e.weight = 0;
    return;
  }

  meclepton_3p3h.pdg = meclepton.pdg = e.in[0].pdg-1+2*ap;
  meclepton.set_mass (PDG::mass (meclepton.pdg)); //set mass coresponding to pdg
  meclepton_3p3h.set_mass (PDG::mass(meclepton_3p3h.pdg));
  ml=meclepton.mass();
  ml2=ml*ml;

 
  PB=p.MEC_pauli_blocking;

  // Binding energy / Correlation within the medium 

  switch(p.nucleus_p) {
    case 6: Bmec = ap ? E_corr[1] : E_corr[0]; break;
    case 8: Bmec = ap ? E_corr[3] : E_corr[2]; break;
    case 40: Bmec = ap ? E_corr[5] : E_corr[4]; break;
    default: Bmec = ap ? E_corr[1] : E_corr[0]; break;  // By default Carbon grid is set in src/nucleus.cc so Bmec is set as E_corr of Carbon
  }
 


  double q0max = e.in[0].energy() - ml - Bmec;
  width_q0 = q0max;

  if( q0max > qmax_Nieves)
  {
    q0max = qmax_Nieves;
  }

  
  particle mecnucleon_in[2];
  particle mecnucleon_out[2];

  particle mecnucleon_3p3h_in[3];
  particle mecnucleon_3p3h_out[3];

  double individual_weight[4] = {0,0,0,0};
  double weight_2p2h=0;
  double weight_3p3h=0;

  if(width_q0>0)
  {
    weight_2p2h=Valencia2020_kin_and_weight_2p2h (e.in[0].E(), individual_weight, meclepton, mecnucleon_in, mecnucleon_out, flag_2p2h_pn, t, mec_central, mec_smearing, binding, ile_pb, mc_sampling, mc_strength);
    
    if (p.MEC_3p3h_dynamics) weight_3p3h = Valencia2020_kin_and_weight_3p3h (e.in[0].E(), individual_weight, meclepton_3p3h, mecnucleon_3p3h_in, mecnucleon_3p3h_out, t, binding, ile_pb);
  }

  double weight = weight_2p2h + weight_3p3h;
  e.weight = weight;
  if(weight > 0)
  {
    double ratio_2p2h  = weight_2p2h/weight;
    flag_2p2h = (bool)(frandom() < ratio_2p2h);
    flag_3p3h = (!flag_2p2h);
  
    if(flag_2p2h)
    {
        e.in.push_back (mecnucleon_in[0]);
        e.in.push_back (mecnucleon_in[1]);
        e.out.push_back (meclepton);
        e.out.push_back (mecnucleon_out[0]);
        e.out.push_back (mecnucleon_out[1]);
    }
    else if (flag_3p3h) {
        e.in.push_back (mecnucleon_3p3h_in[0]);
        e.in.push_back (mecnucleon_3p3h_in[1]);
        e.in.push_back (mecnucleon_3p3h_in[2]);
        e.out.push_back (meclepton_3p3h);
        e.out.push_back (mecnucleon_3p3h_out[0]);
        e.out.push_back (mecnucleon_3p3h_out[1]);
        e.out.push_back (mecnucleon_3p3h_out[2]);
    }
    else {
      cout << "Impossible event topology !! \n Stopping .. \n";
      std::exit(1);
    }
  }
  // Reset the event topology flags
  flag_2p2h = false;
  flag_2p2h_pn = false;
  flag_3p3h = false;
  
}
////////////////////////////////////////

double Valencia2020_kin_and_weight_2p2h (double E, double *individual_dsdqdw, particle &meclep, particle *inc_nucleon_2p2h, particle *out_nucleon_2p2h, bool &flag_2p2h_pn, nucleus &t,
                              double mec_central, double mec_smearing, double binding,
                              int ile_PB, double* sampling, int* strength)
{
  // it is assumed that neutrino direction is (0,0,1); but transition to other direction in nuwro.cc!

  // here we extrapolate the cross section from 2020 Valencia's grid file stored in src/MEC_data/ pp.h ... and so on....
  double result=0;
    //Energy Outgoing Lepton Momentum As Defined In New Valencia Code
  double Eout_lepton = ml + width_q0*frandom();
  double q0 = E - Eout_lepton - Bmec;




  // squared final lepton momentum
  double lp2 = Eout_lepton*Eout_lepton - ml2;

  // final lepton momentum
  double lp=sqrt(lp2);
  // minimum cosine boundary from momentum transfer cut
  

  double cos_min=0.5*(E*E+lp2-qmax_Nieves*qmax_Nieves)/E/lp;
  if(cos_min<-1.0) cos_min=-1.0;
  if(cos_min<1)
  {
    double width_cos=1.0-cos_min;
    //we randomize scattering angle
    double ct=cos_min+width_cos*frandom();

    //squared momentum transfer
    double q2=E*E+lp2-2*E*lp*ct;
    if((q2>q0*q0))
    {
      
	//double differential cross section
     
      for(int outgoing_pair = 0; outgoing_pair < 3; outgoing_pair++)
      {
        individual_dsdqdw[outgoing_pair]=Valencia2020_dsdEdc(E, outgoing_pair, q0,Eout_lepton,ct, t)*width_q0*width_cos;
        result += individual_dsdqdw[outgoing_pair];
      }

      q0+=Bmec;

            // Why should we bother to do anything, if weight=0??
      if (result>0) {
        double phi = 2.0 * Pi * frandom();

        //momentum transfer
        vec qq (cos(phi) * lp * sqrt(1.0 - ct * ct),
                sin(phi) * lp * sqrt(1.0 - ct * ct),
                E - lp * ct);

        //lepton momentum
        vec kprim  (-cos(phi) * lp * sqrt(1.0 - ct * ct),
                    -sin(phi) * lp * sqrt(1.0 - ct * ct),
                    lp * ct);

        //4-momentum transfer
        vect qqq (qq, q0);
        unsigned licz=0;
        vect suma;

        
        double dsdqdw_pp = (individual_dsdqdw[0] > 0) ? individual_dsdqdw[0] : 0;
        double dsdqdw_np = (individual_dsdqdw[1] > 0) ? individual_dsdqdw[1] : 0;
        double dsdqdw_pn = (individual_dsdqdw[2] > 0) ? individual_dsdqdw[2] : 0;

        double ratio_pp = dsdqdw_pp/result;
        double ratio_np = ((dsdqdw_np) > 0 || (dsdqdw_pn > 0)) ? dsdqdw_np/(dsdqdw_np + dsdqdw_pn) : 0 ;
        
        Isospin_model_2p2h_2020Valencia(inc_nucleon_2p2h, out_nucleon_2p2h, ratio_pp, ratio_np, flag_2p2h_pn);
        Generate_nucleon_kinematics_2p2h(inc_nucleon_2p2h, out_nucleon_2p2h, flag_2p2h_pn, t, result, qqq, mec_central, mec_smearing, binding, ile_PB, sampling, strength);
        
        meclep.set_momentum(kprim);
        
      }
      else {
        result = 0;
      }
    }
  }
  return result;
}
////////////////////////////////////////

double Valencia2020_kin_and_weight_3p3h (double E, double *individual_dsdqdw, particle &meclep, particle *inc_nucleon_3p3h, particle *out_nucleon_3p3h, nucleus &t, double binding, int ile_PB)
{
  // it is assumed that neutrino direction is (0,0,1); but transition to other direction in nuwro.cc!

  // here we extrapolate the cross section from 2020 Valencia's data file
  double result=0;
    
  //Energy Outgoing Lepton Momentum As Defined In New Valencia Code
  double Eout_lepton = ml + (E-ml)*frandom();
  double q0 = E - Eout_lepton;




  // squared final lepton momentum

  double lp2 = Eout_lepton*Eout_lepton - ml2;

  // final lepton momentum
  double lp=sqrt(lp2);
  // minimum cosine boundary from momentum transfer cut
  

  double cos_min=0.5*(E*E+lp2-qmax_Nieves*qmax_Nieves)/E/lp;
  if(cos_min<-1.0) cos_min=-1.0;
  if(cos_min<1)
  {
    double width_cos=1.0-cos_min;
    //we randomize scattering angle
    double ct=cos_min+width_cos*frandom();

    //squared momentum transfer
    double q2=E*E+lp2-2*E*lp*ct;
    if((q2>q0*q0))
    {
      
	//double differential cross section
      result = Valencia2020_dsdEdc(E, 3, q0,Eout_lepton,ct, t)*(E-ml)*width_cos;//GeV;

      

      // Why should we bother to do anything, if weight=0??
      if(result>0)
      {
        double phi = 2.0 * Pi * frandom();

        //momentum transfer
        vec qq (cos(phi) * lp * sqrt(1.0 - ct * ct),
                sin(phi) * lp * sqrt(1.0 - ct * ct),
                E - lp * ct);

        //lepton momentum
        vec kprim  (-cos(phi) * lp * sqrt(1.0 - ct * ct),
                    -sin(phi) * lp * sqrt(1.0 - ct * ct),
                    lp * ct);

        //4-momentum transfer
        vect qqq (qq, q0);
        unsigned licz=0;
        vect suma;


        Generate_nucleon_kinematics_3p3h(inc_nucleon_3p3h, out_nucleon_3p3h, t, result, qqq, ile_PB);
        Isospin_model_3p3h_2020Valencia (inc_nucleon_3p3h, out_nucleon_3p3h, t);
        meclep.set_momentum(kprim);
        
      }
      else 
      {
        result = 0;
      }
    }
  }
  return result;
}

////////////////////////////////////////

double Valencia2020_dsdEdc(double E, int outgoing_pair, double q0, double Ep, double ct, nucleus &T)
{
  double result=0;
  double lp=sqrt(Ep*Ep-ml2);
  double vq2=E*E+lp*lp-2*E*lp*ct;
  double q02 = (outgoing_pair == 3) ? q0*q0 : (q0+Bmec)*(q0+Bmec);
  if((vq2>q0*q0)&&(vq2>0)&&(fabs(ct)<=1))
  {
    double q=sqrt(vq2);
    q0 = (outgoing_pair == 3) ? q0-Bmec : q0;
    // interpolating the nuclear tensor elements
    double w00=0;
    double w03=0;
    double w11=0;
    double w12=0;
    double w33=0;

    int n = 0;
    int m = 0;
    int k;
    double y2 = 0;
    double y1 = 0;
    double x2 = 0;
    double x1 = 0;
    //q0-=Bmec;
    if((q0<=q)and(q>=0.01*GeV)and(q0>=0)and(q<=qmax_Nieves))
    {

      /*          q0
       *          | 
       *q0=1200 _ |    __________________________
       *          |   |                         |
       *          |   |                         |
       *          |   |             2D Matrix   |
       *          |   |                         |
       *          |   |                         | 
       *          |   |                         |
       *          |   |             k3 ------k4 |
       *          |   |               |      |  |
       *          |   |               | (m>n)|  |
       *          |   |               |______|  |
       *          |   |               k1     k2 |
       *q0=10 _   |   |_________________________|
       *          |         (n<0)
       *          |_______________________________ q
       *              ,                         ,
       *              q=10MeV                 q=1200MeV
       *
       */
      double spacing = 10*MeV;
      int l = outgoing_pair;
      int n = (q0 > spacing) ? (int)((q0-spacing)/spacing) : -1;    // q0 = 10MeV, q = 10MeV corresponds to 0th row and 0th column of the 2D matrix
      int m = (int)((q - spacing)/spacing);
      double y2 = (n<0) ? spacing : (n+2)*spacing;
      double y1 = (n<0) ? 0 : (n+1)*spacing;
      double x2 = (m+2)*spacing;
      double x1 = (m+1)*spacing;
      int k1;
      int k2;
      int k3;
      int k4;

      // x - direction corresponds to q (in MeV); i.e traversing column in the 2D matrix
      // y - direction corresponds to q0 (in MeV); i.e traversing row in the 2D matrix
      // k is the index number in 1D array corresponding to index [n][m];




      //cout << "\n-----------------\n";
      //cout << "q0 = " << q0 << "\t" << "q = " << q << "\n";
      //cout << "index = [" << n << "," << m << "]\n";
      //cout << "Base point 1: (" << x1 << "," << y1 << ")\n";
      //cout << "Base point 2: (" << x1 << "," << y2 << ")\n";
      //cout << "Base point 3: (" << x2 << "," << y1 << ")\n";
      //cout << "Base point 4: (" << x2 << "," << y2 << ")\n";

      if (n<0) {

        k3 = ConvertCoordinate_2D_to_1D(n+1, m);
        k4 = ConvertCoordinate_2D_to_1D(n+1, m+1);

        w00 = Bilinear_Interpolation(q, q0, x1, x2, y1, y2, 0, 0, T.GetValue_W00(l, k3), T.GetValue_W00(l, k4));
        w03 = Bilinear_Interpolation(q, q0, x1, x2, y1, y2, 0, 0, T.GetValue_W03(l, k3), T.GetValue_W03(l, k4));
        w11 = Bilinear_Interpolation(q, q0, x1, x2, y1, y2, 0, 0, T.GetValue_W11(l, k3), T.GetValue_W11(l, k4));
        w12 = Bilinear_Interpolation(q, q0, x1, x2, y1, y2, 0, 0, T.GetValue_W12(l, k3), T.GetValue_W12(l, k4));
        w33 = Bilinear_Interpolation(q, q0, x1, x2, y1, y2, 0, 0, T.GetValue_W33(l, k3), T.GetValue_W33(l, k4));

        //cout << "\n-----------------\n";
        //cout << "q0 = " << q0 << "\t" << "q = " << q << "\n";
        //cout << "index = [" << n << "," << m << "]\n";
        //cout << "Base point 1: (" << y1 << "," << x1 << ") \t Value = " << 0 << "\n";
        //cout << "Base point 2: (" << y1 << "," << x2 << ")\t Value = " << 0 << "\n";
        //cout << "Base point 3: (" << y2 << "," << x1 << ")\t Value = " << T.GetValue_W00(l, k3) << "\n";
        //cout << "Base point 4: (" << y2 << "," << x2 << ")\t Value = " << T.GetValue_W00(l, k4) << "\n";

      } else {
        if (n==m) {
          k1 = ConvertCoordinate_2D_to_1D(n, m);
          k2 = ConvertCoordinate_2D_to_1D(n,m+1);
          k4 = ConvertCoordinate_2D_to_1D(n+1, m+1);
        
          w00 = Bilinear_Interpolation(q, q0, x1, x2, y1, y2, T.GetValue_W00(l, k1), T.GetValue_W00(l, k2), 0, T.GetValue_W00(l, k4));
          w03 = Bilinear_Interpolation(q, q0, x1, x2, y1, y2, T.GetValue_W03(l, k1), T.GetValue_W03(l, k2), 0, T.GetValue_W03(l, k4));
          w11 = Bilinear_Interpolation(q, q0, x1, x2, y1, y2, T.GetValue_W11(l, k1), T.GetValue_W11(l, k2), 0, T.GetValue_W11(l, k4));
          w12 = Bilinear_Interpolation(q, q0, x1, x2, y1, y2, T.GetValue_W12(l, k1), T.GetValue_W12(l, k2), 0, T.GetValue_W12(l, k4));
          w33 = Bilinear_Interpolation(q, q0, x1, x2, y1, y2, T.GetValue_W33(l, k1), T.GetValue_W33(l, k2), 0, T.GetValue_W33(l, k4));
        } else {
            k1 = ConvertCoordinate_2D_to_1D(n, m);
            k2 = ConvertCoordinate_2D_to_1D(n,m+1);
            k3 = ConvertCoordinate_2D_to_1D(n+1, m);
            k4 = ConvertCoordinate_2D_to_1D(n+1, m+1);
        
            w00 = Bilinear_Interpolation(q, q0, x1, x2, y1, y2, T.GetValue_W00(l, k1), T.GetValue_W00(l, k2), T.GetValue_W00(l, k3), T.GetValue_W00(l, k4));
            w03 = Bilinear_Interpolation(q, q0, x1, x2, y1, y2, T.GetValue_W03(l, k1), T.GetValue_W03(l, k2), T.GetValue_W03(l, k3), T.GetValue_W03(l, k4));
            w11 = Bilinear_Interpolation(q, q0, x1, x2, y1, y2, T.GetValue_W11(l, k1), T.GetValue_W11(l, k2), T.GetValue_W11(l, k3), T.GetValue_W11(l, k4));
            w12 = Bilinear_Interpolation(q, q0, x1, x2, y1, y2, T.GetValue_W12(l, k1), T.GetValue_W12(l, k2), T.GetValue_W12(l, k3), T.GetValue_W12(l, k4));
            w33 = Bilinear_Interpolation(q, q0, x1, x2, y1, y2, T.GetValue_W33(l, k1), T.GetValue_W33(l, k2), T.GetValue_W33(l, k3), T.GetValue_W33(l, k4));
        }
      }

      q0+=Bmec; 

    }
 
    double w1=0.5*w11;
    double w2=0.5*(w00+w11+q02/vq2*(w33-w11)-2*q0/q*w03);
    double w3=(1-2*ap)*w12/q;
    double w4=0.5*(w33-w11)/vq2;
    double w5=(w03-q0/q*(w33-w11))/q;
    double sts=1.0-ct*ct;
    double st2=0.5*(1.0-ct);
    double ct2=0.5*(1.0+ct);
    result=(2*w1*st2+w2*ct2-w3*(E+Ep)*st2);
    result+=ml2/((Ep+lp)*Ep)*(w1*ct-0.5*w2*ct+0.5*w3*(Ep*(1.0-ct)+lp-E*ct)+0.5*w4*(ml2*ct+2*Ep*(Ep+lp)*sts)-0.5*w5*(Ep+lp));
    result*=2*G*G*Ep*lp/Pi/cm2;
    
  }




  return result;
}

void Isospin_model_2p2h_2020Valencia (particle *in_p, particle *out_p, double ratio_pp, double ratio_np, bool &flag_2p2h_pn)
{
  // Variables and flags for anti neutrinos 
  double& ratio_nn_antinu = ratio_pp;
  double& ratio_pn_antinu = ratio_np;
  bool& flag_2p2h_pn_antinu = flag_2p2h_pn;

  if (ap) { // For anti-neutrinos
    if (frandom() < ratio_nn_antinu) { // initial pairs will be np  
      in_p[0].set_neutron(); in_p[1].set_proton();
      out_p[0].set_neutron(); out_p[1].set_neutron();
    } else {
        if (frandom() < ratio_pn_antinu) {// initial pairs will be pp
          in_p[0].set_proton(); in_p[1].set_proton();
          out_p[0].set_neutron(); out_p[1].set_proton();
          flag_2p2h_pn_antinu = true;
        } else {  // Initial pairs will be pp
            in_p[0].set_proton(); in_p[1].set_proton();
            out_p[0].set_neutron(); out_p[1].set_proton();
          }
      }

  } 
  else {
    if (frandom() < ratio_pp) { // Initial pairs will be  np 
      in_p[0].set_neutron(); in_p[1].set_proton();
      out_p[0].set_proton(); out_p[1].set_proton();
    } else {
        if (frandom() < ratio_np) { // Initial pairs will be nn 
          in_p[0].set_neutron(); in_p[1].set_neutron();
          out_p[0].set_neutron(); out_p[1].set_proton();
        } else { // Initial pairs will be nn
            in_p[0].set_neutron(); in_p[1].set_neutron();
            out_p[0].set_neutron(); out_p[1].set_proton();
            flag_2p2h_pn = true;
          }   
      }
  }
}


void Isospin_model_3p3h_2020Valencia (particle *in_p, particle *out_p, nucleus &T_nucleus)
{
  
    float NP = (float)T_nucleus.p;     // <- No of protons in the nucleus
    float NN = (float)T_nucleus.n;    // <- No of neutron in the nucleus

    double ppp = (NP >= 3) ? (NP * (NP - 1) * (NP - 2)) / 6 : 0;      // inital pairs npp
    double nnn = (NN >=3 ) ? (NN * (NN - 1) * (NN - 2)) / 6 : 0;      // inital pairs pnn

    double npp = (NP >= 2 && NN >= 1) ? ((NP * (NP - 1)) / 2) * NN : 0;   // initial pairs nnp 
    double pnn = (NN >= 2 && NP >= 1) ? ((NN * (NN - 1)) / 2) * NP : 0;   // initial pairs ppn

    double nnp = (NN >= 2 && NP >= 1) ? ((NN * (NN - 1)) / 2) * NP : 0;   // initial pairs nnn
    double ppn = (NP >= 2 && NN >= 1) ? ((NP * (NP - 1)) / 2) * NN : 0;   // initial pairs ppp

    double total_pairs_nu = ppp + npp + nnp;
    double total_pairs_antinu = nnn + pnn + ppn;

    if (ap) {
      if (frandom() < ppn/total_pairs_antinu) {
        in_p[0].set_proton(); in_p[1].set_proton(); in_p[2].set_proton();
        out_p[0].set_neutron(); out_p[1].set_proton(); out_p[2].set_proton();
      } else {
          if (frandom() < pnn/(pnn+nnn)) {
            in_p[0].set_proton(); in_p[1].set_proton(); in_p[2].set_neutron();
            out_p[0].set_neutron(); out_p[1].set_neutron(); out_p[2].set_proton();
          } else {
              in_p[0].set_proton(); in_p[1].set_neutron(); in_p[2].set_neutron();
              out_p[0].set_neutron(); out_p[1].set_neutron(); out_p[2].set_neutron();
            }
        }
    } 

    else {
      if (frandom() < nnp/total_pairs_nu) {
        in_p[0].set_neutron(); in_p[1].set_neutron(); in_p[2].set_neutron();
        out_p[0].set_neutron(); out_p[1].set_neutron(); out_p[2].set_proton();
      } else {
          if (frandom() < npp/(npp+ppp)) {
            in_p[0].set_neutron(); in_p[1].set_neutron(); in_p[2].set_proton();
            out_p[0].set_neutron(); out_p[1].set_proton(); out_p[2].set_proton();
          } else {
              in_p[0].set_neutron(); in_p[1].set_proton(); in_p[2].set_proton();
              out_p[0].set_proton(); out_p[1].set_proton(); out_p[2].set_proton();
            }
        }   
    }
}

void Generate_nucleon_kinematics_2p2h(particle *inc_nucleon_mec, particle *out_nucleon_mec, bool &flag_2p2h_pn, nucleus &T, double &xsec_result, vect qqq,double mec_central, double mec_smearing, double binding, int ile_PB, double* sampling, int* strength)
{
  vect suma(vec(0,0,0), 0);   // 4-momentum vector of the whole system ... initial nucleon pairs + 4-momentum energy transfer
  unsigned int licz = 0;      // <-- To keep count of hadronic trials

  // Mass of the nucleons
  double Mass1 = inc_nucleon_mec[0].mass();
  double Mass2 = inc_nucleon_mec[1].mass();
  double Mass3 = out_nucleon_mec[0].mass();
  double Mass4 = out_nucleon_mec[1].mass();

  // Flag for the type of dyanmics decided by mec_do_cc_Nieves ()
  bool pp=false;
  bool np=false;
  bool pn=false;

  if(out_nucleon_mec[0].pdg == 2212)
  {
     pp = true;
     np = false;
     pn = false;
  }
  else {
     pp = false;
     np = true;
     pn = false;
     if (flag_2p2h_pn == true) { pp = false; np = false; pn = true;}
  }
  //Pauli-NoPauli Blocking Version  <- Hemant Prasad 24/09/2023
  //-------------------------------------------------------------------------------------------------
  do{

    particle probe1; // <- dummy nucleon inside the nucleus . Randomly sampled between [0, r_max] weighted by the nuclear density !!
    vec N1 = vec(0,0,0);
    vec N2 = vec(0,0,0);
    double lengthN1 = 0;
    double lengthN2 = 0;
    do {
      probe1 = T.get_nucleon();
      vec pos = vec(probe1.r);
      double pos_magnitude = pos.length();

      N1 = rand_from_ball( T.localkf_(inc_nucleon_mec[0].pdg, pos_magnitude)  );
      N2 = rand_from_ball( T.localkf_(inc_nucleon_mec[1].pdg, pos_magnitude)  );
      
      // Arange such that N1 contains higher momentum and N2 contains lower momentum 
      vec temp_vector = N1;
      if (N1.length() < N2.length()) { N1 = N2; N2 = temp_vector; } 
      double Ek1_in = sqrt ( std::pow(N1.length() ,2) + Mass1 * Mass1 ) - Mass1;        // <-- K.E of nucleon assuming first nucleon has higher momentum 
      double Ek2_in = sqrt ( std::pow(N2.length(), 2) + Mass2 * Mass2 ) - Mass2;        // <-- K.E of nucleon assuming second nucleon has lower momentum 
      // ==> We provide an extra K.E of about ~ 15 MeV to initial nucleon. 0.8 fraction of the extra kinetic energy goes to nucleon with higher momentum <==
      double E_KE = 14.3*MeV;
      double fraction = 0.5;
      double frac_Ek1_in = fraction * E_KE;
      double frac_Ek2_in = (1 - fraction) * E_KE; 
      
      N1 = N1 * sqrt( std::pow( Ek1_in + Mass1 + frac_Ek1_in, 2) - Mass1 * Mass1 ) / N1.length() ;
      N2 = N2 * sqrt( std::pow( Ek2_in + Mass2 + frac_Ek2_in, 2) - Mass2 * Mass2 ) / N2.length() ;

      lengthN1 = N1.length();
      lengthN2 = N2.length();

    } while ( (lengthN1 > 246*MeV) && (lengthN2 > 246*MeV) );

    
    //vec N2 = -N1;                                        // <- legacy : momentum are back to back
    if(mec_smearing > 0.0)  // smearing of back to back
    {
       N2 = -N1;
       N2.x = N2.x * rand_gauss (mec_smearing, 1.0);
       N2.y = N2.y * rand_gauss (mec_smearing, 1.0);
       N2.z = N2.z * rand_gauss (mec_smearing, 1.0);
    }

    if(mec_central > 0.0)   // cm motion
    {
      N2 = -N1;
      vec central = rand_gauss(mec_central);
      inc_nucleon_mec[0].set_momentum(central + N1);
      inc_nucleon_mec[1].set_momentum(central + N2);
    }
    
    inc_nucleon_mec[0].r = probe1.r;     // <- Get a random position of incoming nucleon 1.
    inc_nucleon_mec[1].r = probe1.r;     // <- set the position of initial nucleon 2 at nucleon 1

    inc_nucleon_mec[0].set_momentum(N1);                    // <- set the momentum for initial nucleon 1
    inc_nucleon_mec[1].set_momentum(N2);                    // <- set the momentum for initial nucleon 2
                                                            //
    inc_nucleon_mec[0].set_fermi(T.Ef(inc_nucleon_mec[0])); // local fermi energy is used in energy balance
    //inc_nucleon_mec[1].set_fermi(T.Ef(inc_nucleon_mec[1])); // local fermi energy for initial nucleon 2
    
        //-------------------------------------------------------------------------------------------------------------------------------------------------------------//

    suma = inc_nucleon_mec[0].p4() + inc_nucleon_mec[1].p4() + qqq; // <-- Total 4-momentum sum of the initial system
    double bilans = inc_nucleon_mec[0].Ek() + inc_nucleon_mec[1].Ek()  - 2*(inc_nucleon_mec[0].his_fermi + 0.5*binding);

    if (bilans > 0) { // to avoid perpetuum mobile
      vect roznica (bilans, 0, 0, 0);
      suma -= roznica;
    }

    licz++;     // <- hadronic trials are counted here ...

    if( (suma*suma > pow(Mass3 + Mass4, 2)) && (suma.t > suma.length()) ) { //  <-- Hadronic condition is checked here
      vec dir_cm;                               // <- direction vector which have component || to boost
      vec dir_cm2;                              // <- direction vector which have component anti || to boost
      vec trans = suma.v();                     // <- boost direction
      suma.boost2(trans);                       // <- boost to CM frame . Note! boost(2) boost in the reverese direction
      double Ecm3 = 0.5*(suma.t*suma.t + Mass3*Mass3 - Mass4*Mass4)/suma.t; // Energy of outgoing nucleon 1
      double Ecm4 = 0.5*(suma.t*suma.t + Mass4*Mass4 - Mass3*Mass3)/suma.t; // Energy of outgoing nucleon 2

      if (!PB) {  // <-- A piece of code if the mec_pauli_blocking is turned of in data/params.txt

        // dir_cm is sampled differently for different outgoing pairs /
        // sampling[0] corresponds for sampling direction for pp
        // sampling[1] corresponds for sampling direction for np
        // sampling[2] corresponds for sampling direction for pn;
        dir_cm = rand_direc2(trans, pp*sampling[0] + np*sampling[1] + pn*sampling[2], pp*strength[0] + np*strength[1] + pn*strength[2], 1);
        dir_cm2 = -dir_cm;        // <-- dir_cm2 is opposite of dir_cm



        double cos_theta_hadronic = dir_cm*trans/dir_cm.length()/trans.length();
        if (cos_theta_hadronic < 0) {  // <-- Here the dir_cm is set as the vector which has component || to the direction of boost
                                       // <-- if cos(0)_cm is negetive then pi/2 <= theta <= pi, so this must have anti || component so reverse the // direction to have || component.
          dir_cm *= -1.0;
          dir_cm2 *= -1.0;              // <-- Vice versa for dir_cm2.
        }

        // Assign momentum to the outgoing nucleon pairs.
        if (pp) {
          out_nucleon_mec[0].set_momentum(dir_cm*sqrt(Ecm3*Ecm3 - Mass3*Mass3));      // <- will always be leading proton
          out_nucleon_mec[1].set_momentum(dir_cm2*sqrt(Ecm4*Ecm4 - Mass4*Mass4));     // <- will always be subleading proton
        } else if (np) {
            out_nucleon_mec[0].set_momentum(dir_cm*sqrt(Ecm3*Ecm3 - Mass3*Mass3));    // <- neutron will have higher momentum in np
            out_nucleon_mec[1].set_momentum(dir_cm2*sqrt(Ecm4*Ecm4 - Mass4*Mass4));   // <- proton will have lower momentum in np
          } else if (pn) {
              out_nucleon_mec[0].set_momentum(dir_cm2*sqrt(Ecm3*Ecm3 - Mass3*Mass3)); // <- neutron will have lower momentum in pn
              out_nucleon_mec[1].set_momentum(dir_cm*sqrt(Ecm4*Ecm4 - Mass4*Mass4));  // <- proton will have higher momentum in pn
            }

        out_nucleon_mec[0].boost2(-trans);      // <<- Boost back to LAB frame
        out_nucleon_mec[1].boost2(-trans);      // <<- Boost back to LAB frame

      }
      else {      // <-- A piece of code if the mec_pauli_blocking is turned on in data/params.txt

        int count_pb_trials=0;  // <-- To keep count of pauli blocked trials.... dynamics should be found in a finite number of trials.
        while (PB && count_pb_trials < ile_PB ) {

          //--- Boost parameters--- //
          double beta = trans.length();
          double gamma = 1.0/sqrt( 1 - beta*beta );
          // ----------------------//


          double kappa = (gamma*(np*Ecm4 + pn*Ecm3 + pp*Ecm4) - np*Mass4 - pn*Mass3 -pp*Mass4 - inc_nucleon_mec[0].his_fermi)/beta/gamma/(np*sqrt(Ecm4 * Ecm4 - Mass4*Mass4) + pn*sqrt(Ecm3*Ecm3 - Mass3*Mass3) + pp*sqrt(Ecm4*Ecm4 - Mass4*Mass4)); // <-- Allowed cosine boundary determined by the nucleus which have anti || component to boost

          if (!( gamma*(pp*Ecm4 + np*Ecm4 + pn*Ecm3) <= (pp*Mass4 + np*Mass4 + pn*Mass3 + inc_nucleon_mec[0].his_fermi) )) { // <- Pauli Blocking condiiton is checked
            PB = false;     // <-- Set the PB flag false to get out of the do{ ...  }while(PB && (licz < calls_max)); loop
            // dir_cm is sampled differently for different outgoing pairs /
            // sampling[0] corresponds for sampling direction for pp
            // sampling[1] corresponds for sampling direction for np
            // sampling[2] corresponds for sampling direction for pn;
            dir_cm = rand_direc2(trans, pp*sampling[0] + np*sampling[1] + pn*sampling[2], pp*strength[0] + np*strength[1] + pn*strength[2], TMath::Min(kappa, 1.0));
            dir_cm2 = -dir_cm;    // <-- dir_cm2 is opposite of dir_cm

            //if (np || pn)
            //{

            double cos_theta_hadronic = dir_cm*trans/dir_cm.length()/trans.length();
            if (cos_theta_hadronic < 0) {       // <-- Here the dir_cm is set as the vector which has component || to the direction of boost
                                                // <-- if cos(0)_cm is negetive then pi/2 <= theta <= pi, so this must have anti || component so reverse the // direction to have || component.
              dir_cm *= -1.0;
              dir_cm2 *= -1.0;
              cos_theta_hadronic = dir_cm*trans/dir_cm.length()/trans.length();
            }

            // Assign momentum to the outgoing nucleon pairs.
            if (pp) {
              out_nucleon_mec[0].set_momentum(dir_cm*sqrt(Ecm3*Ecm3 - Mass3*Mass3));      // <- Will always be leading proton.
              out_nucleon_mec[1].set_momentum(dir_cm2*sqrt(Ecm4*Ecm4 - Mass4*Mass4));     // <- Will always be sub-leading proton
            } else if (np) {
                out_nucleon_mec[0].set_momentum(dir_cm*sqrt(Ecm3*Ecm3 - Mass3*Mass3));    // <- Neutron will have higher momentum in np event.
                out_nucleon_mec[1].set_momentum(dir_cm2*sqrt(Ecm4*Ecm4 - Mass4*Mass4));   // <- Proton will have lower momentum in np event.
              } else if (pn) {
                  out_nucleon_mec[0].set_momentum(dir_cm2*sqrt(Ecm3*Ecm3 - Mass3*Mass3));   // <- Neutron will have lower momentum in pn event.
                  out_nucleon_mec[1].set_momentum(dir_cm*sqrt(Ecm4*Ecm4 - Mass4*Mass4));    // <- Proton will have higher momentum in pn event.
                }

            out_nucleon_mec[0].boost2(-trans);    // <<- Boost back to LAB frame
            out_nucleon_mec[1].boost2(-trans);    // <<- Boost back to LAB frame

          }

          count_pb_trials++; // <- Pauli blocking trials are counted here ...
        }

      }

    }
  }while(PB && (licz < calls_max));
  if(licz >= calls_max)
  {
    xsec_result =0;
    vec zero (0, 0, 0);
    out_nucleon_mec[0].set_momentum(zero);
    out_nucleon_mec[1].set_momentum(zero);
    out_nucleon_mec[0].set_energy(Mass3);
    out_nucleon_mec[1].set_energy(Mass4);
  }
  //-------------------------------------------------------------------------------------------------
}

void Generate_nucleon_kinematics_3p3h(particle *inc_nucleon_mec, particle *out_nucleon_mec, nucleus &T, double &xsec_result, vect qqq, int ile_PB)
{
  vect suma(vec(0,0,0), 0);
  unsigned int licz = 0;
  unsigned int pb_trials = 0;
  vec N1, N2, N3;
  vect N11, N22, N33;
  vec maksio_vec;

  double length1, length2, length3;
  double localfermi;
  double Winv;

  vec trans;

  inc_nucleon_mec[0].set_neutron();
  inc_nucleon_mec[1].set_proton();
  inc_nucleon_mec[2].set_proton();

  double mass[3];
  double totalrestmass = 0; // total rest mass of the system -> mass[0] + mass[1] + mass[2];
  for(int mass_index = 0; mass_index < 3; mass_index++)
  {
    mass[mass_index] = inc_nucleon_mec[mass_index].mass();
    totalrestmass += mass[mass_index];
  }

  do
  {
    particle probe1 = T.get_nucleon();
    vec pos (probe1.x, probe1.y, probe1.z);

    localfermi = T.localkf( probe1);

    for(int nucleon_index = 0; nucleon_index < 3; nucleon_index++)
    {
      out_nucleon_mec[nucleon_index].r.x = probe1.r.x;
      out_nucleon_mec[nucleon_index].r.y = probe1.r.y;
      out_nucleon_mec[nucleon_index].r.z = probe1.r.z;
    }

    N1 = rand_from_ball (localfermi);
    N2 = rand_from_ball (localfermi);
    N3 = rand_from_ball (localfermi);

    inc_nucleon_mec[0].set_momentum(N1);
    inc_nucleon_mec[1].set_momentum(N2);
    inc_nucleon_mec[2].set_momentum(N3);

    length1 = sqrt(MN2 + N1.norm2());
    length2 = sqrt(MN2 + N2.norm2());
    length3 = sqrt(MN2 + N3.norm2());

    N11 = vect(N1, length1);
    N22 = vect(N2, length2);
    N33 = vect(N3, length3);

    suma = N11 + N22 + N33 + qqq;

    licz++;
  }while(  ( (suma.t < suma.length()) || (suma * suma < totalrestmass*totalrestmass ) ) && (licz < calls_max) ) ;

  if(  ( (suma.t > suma.length()) && (suma * suma >= totalrestmass*totalrestmass ) )  && (licz < calls_max))
  {
    //std::cerr << suma * suma << "\t" << totalrestmass << "\n";
    out_nucleon_mec[0].set_proton();
    out_nucleon_mec[1].set_proton();
    out_nucleon_mec[2].set_proton();

    trans = suma.v();
    suma.boost2(trans);
    Winv = suma.t;

    vect nuc11, nuc22, nuc33;
    //threebodydecay(Winv, mass[0], mass[1], mass[2], nuc11, nuc22, nuc33);

    if(PB)
    {
      while(PB && pb_trials < ile_PB)
      {
        threebodydecay(Winv, mass[0], mass[1], mass[2], nuc11, nuc22, nuc33);
        if( (nuc11.length() >= localfermi) && (nuc22.length() >= localfermi ) && ( nuc33.length() >= localfermi ) )
        {
          PB = false;
          break;
        }
        pb_trials++;
      }
    }
    else
    {
      threebodydecay(Winv, mass[0], mass[1], mass[2], nuc11, nuc22, nuc33);
    }

    nuc11.boost2(-trans);
    nuc22.boost2(-trans);
    nuc33.boost2(-trans);

    out_nucleon_mec[0].set_momentum(vec(nuc11));
    out_nucleon_mec[1].set_momentum(vec(nuc22));
    out_nucleon_mec[2].set_momentum(vec(nuc33));


  }
  else if ((licz >= calls_max) || (pb_trials >= ile_PB))
  {
    xsec_result = 0;
    vec zero(0,0,0);
    for(int nucleon_index = 0; nucleon_index < 3; nucleon_index++)
    {
      out_nucleon_mec[nucleon_index].set_momentum(zero);
    }
  }

}

void threebodydecay (double W, double m1, double m2, double m3, vect &p1, vect &p2, vect &p3)
{
	if (W < m1+m2+m3)
	{
          std::cerr << "impossible kinematics in 3-body decay" << "\n";

	}

	double W_2 = W*W;
	double m3_2 = m3*m3;
	double m1_2 = m1*m1;
	double m2_2 = m2*m2;

	double m23_2_los, m23_2_min, m23_2_max, E3;

	do
	{
		double m12_2_los = (m1+m2)*(m1+m2) + frandom()* ( (W-m3)*(W-m3) - (m1+m2)*(m1+m2) );
		double m23_test = (m2+m3)*(m2+m3) + frandom()* ( (W-m1)*(W-m1) - (m2+m3)*(m2+m3) );
		m23_2_los = m23_test;

		E3 = (W_2 + m3_2 - m12_2_los)/2.0/W;

		double E2star = (m12_2_los - m1_2 + m2_2)/2.0/sqrt(m12_2_los);
		double E3star = (W_2 - m12_2_los - m3_2)/2.0/sqrt(m12_2_los);

		m23_2_min = (E2star+E3star)*(E2star+E3star)
		- ( sqrt(E2star*E2star - m2_2) + sqrt (E3star*E3star - m3_2) )*
		( sqrt(E2star*E2star - m2_2) + sqrt (E3star*E3star - m3_2) );

		m23_2_max = (E2star+E3star)*(E2star+E3star)
		- ( sqrt(E2star*E2star - m2_2) - sqrt (E3star*E3star - m3_2) )*
		( sqrt(E2star*E2star - m2_2) - sqrt (E3star*E3star - m3_2) );

	}
	while ( (m23_2_los < m23_2_min) || (m23_2_los > m23_2_max) );

	//double m23_2_los = m23_2_min + los()* ( m23_2_max - m23_2_min);

	double E1 = (W_2 + m1_2 - m23_2_los)/2.0/W;
	double E2 = W - E1 - E3;

	double mom1 = sqrt (E1*E1 - m1_2);
	double mom2 = sqrt (E2*E2 - m2_2);
	double mom3 = sqrt (E3*E3 - m3_2);

	vec axis = rand_dir (); // new -- Y axis
	vec proba = rand_dir ();

	vec kier1 = vecprod (axis, proba);
	kier1 = kier1/kier1.length();//unit vector -- new Z axis!

	// unit vetor -- new X axis
	vec kier2 = vecprod (axis, kier1);

	double kosalpha = (mom1*mom1 + mom2*mom2 - mom3*mom3)/2.0/mom1/mom2;
	double sinalpha = sqrt(1-kosalpha*kosalpha);
	double kosbeta = ( mom1 - mom2*kosalpha )/mom3;
	double sinbeta = sqrt(1 - kosbeta*kosbeta);

	vec pp1 = kier1*mom1;
	vec pp2 = mom2*(-kosalpha*kier1 + sinalpha*kier2);
	vec pp3 = mom3*(-kosbeta*kier1 - sinbeta*kier2);

	p1 = vect(pp1,E1);
	p2 = vect(pp2,E2);
	p3 = vect(pp3,E3);

}


