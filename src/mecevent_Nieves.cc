#include "mecevent_Nieves.h"

void mecevent_Nieves (params & p, event & e, nucleus & t, bool cc)
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
  double mc_sampling = p.MEC_cm_direction;

  // sadly, only CC events available so far...
  if(e.flag.nc)
  {
    cerr<<" MEC error: Wrong Settings!\n";
    e.weight = 0;
    return;
  }

  particle meclepton;
  ap=(e.in[0].pdg<0);
  //if(e.in[0].pdg>0) meclepton.pdg = e.in[0].pdg-1;  //all leptons available
  //else {meclepton.pdg = e.in[0].pdg+1;}
  meclepton.pdg = e.in[0].pdg-1+2*ap;
  meclepton.set_mass (PDG::mass (meclepton.pdg)); //set mass coresponding to pdg
  ml=meclepton.mass();
  ml2=ml*ml;

  //nucleus and Pauli blocking
  nucl=(t.p>7);//+(t.p>15);
  PB=p.MEC_pauli_blocking;

  //binding:either carbon or oxygen or calcium neutrino/antineutrino
  Bmec=qvalues_Nieves[2*nucl+(e.in[0].pdg<0)];

  // Qvalue threshold...
  double q0max=e.in[0].E()-ml-eV;
  if(q0max>qmax_Nieves) q0max=qmax_Nieves;
  width_q0=q0max-Bmec;

  particle mecnucleon[4]; //0,1 - in, 2, 3 - out

  double weight=0;

  if(width_q0>0) weight=Nieves_kin_and_weight (e.in[0].E(), meclepton, mecnucleon, t, mec_central,
                                               mec_smearing, binding, ile_pb, mc_sampling);

  e.weight = weight;
  if (weight>0) mec_do_cc (mecnucleon, p.mec_ratio_pp);
  e.in.push_back (mecnucleon[0]);
  e.in.push_back (mecnucleon[1]);
  e.out.push_back (meclepton);
  e.out.push_back (mecnucleon[2]);
  e.out.push_back (mecnucleon[3]);

  // corrections for A>39;
  // combinatoric factors with respect to hypothetical isosymmetric target here approximated by calcium
  // K.N.: but actually oxygen is used...
  int NN=p.nucleus_n;
  int ZZ=p.nucleus_p;

  if ( NN+ZZ>39 )
  {
    double NZ=(NN+ZZ)/2.0;

    // neutrino on neutron-neutron pair
    if (e.in[0].pdg>0  && (mecnucleon[0].pdg +  mecnucleon[1].pdg) ==4224 )
    {
      double renorm = NN*(NN-1)/NZ/(NZ-1);
      weight=weight*renorm;
    }

    //neutrino on neutron-proton pair
    if (e.in[0].pdg>0  && (mecnucleon[0].pdg +  mecnucleon[1].pdg) ==4324 )
    {
      double renorm = NN*ZZ/NZ/NZ;
      weight=weight*renorm;
    }

    //neutrino on proton-proton pair
    if (e.in[0].pdg<0  && (mecnucleon[0].pdg +  mecnucleon[1].pdg) ==4424 )
    {
      double renorm = ZZ*(ZZ-1)/NZ/(NZ-1);
      weight=weight*renorm;
    }

    // neutrino on neutron-proton pair
    if (e.in[0].pdg<0  && (mecnucleon[0].pdg +  mecnucleon[1].pdg) ==4324 )
    {
      double renorm = NN*ZZ/NZ/NZ;
      weight=weight*renorm;
    }
  }

  /*
  if (e.in[0].pdg>0 && p.nucleus_p==18 && p.nucleus_n==22 && (mecnucleon[0].pdg +  mecnucleon[1].pdg) ==4224)
    weight=weight*1.21;
  if (e.in[0].pdg>0 && p.nucleus_p==18 && p.nucleus_n==22 && (mecnucleon[0].pdg +  mecnucleon[1].pdg) ==4324)
    weight=weight*0.99;
  
  if (e.in[0].pdg>0 && p.nucleus_p==26 && p.nucleus_n==30 && (mecnucleon[0].pdg +  mecnucleon[1].pdg) ==4224)
    weight=weight*1.15;
  if (e.in[0].pdg>0 && p.nucleus_p==26 && p.nucleus_n==30 && (mecnucleon[0].pdg +  mecnucleon[1].pdg) ==4324)
    weight=weight*0.995;
  
  if (e.in[0].pdg>0 && p.nucleus_p==82 && p.nucleus_n==126 && (mecnucleon[0].pdg +  mecnucleon[1].pdg) ==4224)
    weight=weight*1.47;
  if (e.in[0].pdg>0 && p.nucleus_p==82 && p.nucleus_n==126 && (mecnucleon[0].pdg +  mecnucleon[1].pdg) ==4324)
    weight=weight*0.955;
  
  if (e.in[0].pdg<0 && p.nucleus_p==18 && p.nucleus_n==22 && (mecnucleon[0].pdg +  mecnucleon[1].pdg) ==4424)
    weight=weight*0.805;
  if (e.in[0].pdg<0 && p.nucleus_p==18 && p.nucleus_n==22 && (mecnucleon[0].pdg +  mecnucleon[1].pdg) ==4324)
    weight=weight*0.99;
  
  if (e.in[0].pdg<0 && p.nucleus_p==26 && p.nucleus_n==30 && (mecnucleon[0].pdg +  mecnucleon[1].pdg) ==4424)
    weight=weight*0.86;
  if (e.in[0].pdg<0 && p.nucleus_p==26 && p.nucleus_n==30 && (mecnucleon[0].pdg +  mecnucleon[1].pdg) ==4324)
    weight=weight*0.995;
  
  if (e.in[0].pdg<0 && p.nucleus_p==82 && p.nucleus_n==126 && (mecnucleon[0].pdg +  mecnucleon[1].pdg) ==4424)
    weight=weight*0.62;
  if (e.in[0].pdg<0 && p.nucleus_p==82 && p.nucleus_n==126 && (mecnucleon[0].pdg +  mecnucleon[1].pdg) ==4324)
    weight=weight*0.955;
  */
}

////////////////////////////////////////

double Nieves_kin_and_weight (double E, particle &meclep, particle *nucleon, nucleus &t,
                              double mec_central, double mec_smearing, double binding,
                              int ile_PB, double sampling)
{
  // it is assumed that neutrino direction is (0,0,1); but transition to other direction in nuwro.cc!

  // here we extrapolate the cross section from Nieves's data file
  double result=0;
  // Bmec=Qvalue we randomize q0
  double q0_nucl=frandom()*width_q0;
  // de-Forest like subtraction of binding energy...
  double q0=q0_nucl+Bmec;
  // squared final lepton momentum
  double lp2=(E-q0)*(E-q0)-ml2;
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
      result=Nieves_dsdEdc(E,q0,ct)*width_q0*width_cos/GeV;

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

        /*
        //no Pauli blocking version
        do
        {
          nucleon[0] = t.get_nucleon ();
          nucleon[1] = t.get_nucleon ();

          suma = nucleon[0].p4() + nucleon[1].p4() + qqq;
          licz++;
        }//to be able to make Lorentz boost and to "decay"
        while ( (( suma.t < suma.length() ) || ( suma*suma < 4.0 * MN2)) &&(licz<calls_max) );
        if(licz<calls_max)
        {
          vec trans = suma.v();
          nucleon[2].set_proton();  //in mec_do* we decide about isospin  
          nucleon[3].set_proton();
          //cerr<<"\n suma"<<suma;
          //cerr<<" q"<<qqq<<"\n";
          suma.boost2 (trans);      //boost to the CM frame
          double Ecm = suma.t / 2.0;  //each nucleon get the same energy in CM

          vec dir_cm = rand_dir () * sqrt(Ecm * Ecm  - MN2);  //radnomly set direction of nucleons in CM
          nucleon[2].set_momentum (dir_cm);
          nucleon[3].set_momentum (-dir_cm);

          nucleon[2].p4().boost2 (-trans);
          nucleon[3].p4().boost2 (-trans);     //nucleons in the LAB frame
        }
        else
        {
          result=0;
          vect zero (0,0,0);
          nucleon[2].set_momentum (zero);
          nucleon[3].set_momentum (zero);
          nucleon[2].set_energy (MN);
          nucleon[3].set_energy (MN);
        }

        meclep.set_momentum(kprim);*/

        do
        {
          //nucleon[0] = t.get_nucleon ();
          //nucleon[1] = t.get_nucleon ();

          //suma = nucleon[0].p4() + nucleon[1].p4() + qqq;

          nucleon[0] = t.get_nucleon ();
          vec N1=spectral_choice (6, 6); // we take SF carbon distribution; should be improved !!!
          nucleon[0].set_momentum(N1);
          vec ped2=-N1;                  // momenta are roughly back to back

          if (mec_smearing>0.0)          // smearing of back to back (does not make sense with cm motion?)
          {
            ped2.x = -N1.x*rand_gauss (mec_smearing, 1.0); 
            ped2.y = -N1.y*rand_gauss (mec_smearing, 1.0); 
            ped2.z = -N1.z*rand_gauss (mec_smearing, 1.0); 
          }

          vec pos ( nucleon[0].r.x, nucleon[0].r.y, nucleon[0].r.z );
          nucleon[1] = t.get_nucleon (pos); // the same position

          nucleon[1].set_momentum(ped2);

          if (mec_central>0.0) // cm motion
          {
            vec central = rand_gauss(mec_central);
            nucleon[0].set_momentum(central+N1);
            nucleon[1].set_momentum(central+ped2);
          }

          nucleon[0].set_fermi( t.Ef(nucleon[0]) ); // local Fermi energy used in energy balance

          vect suma = nucleon[0].p4() + nucleon[1].p4() + qqq;

          double bilans = nucleon[0].Ek() + nucleon[1].Ek() - 2*(nucleon[0].his_fermi + 0.5*binding);

          if (bilans>0) // to avoid perpetuum mobile
          {
            vect roznica (bilans, 0, 0, 0);
            suma-=roznica;
          }

          licz++;
          if(!( ( suma.t < suma.length() ) || ( suma*suma < 4.0 * MN2)))
          {
            vec trans = suma.v();
            nucleon[2].set_proton();   // in mec_do* we decide about isospin
            nucleon[3].set_proton();
          
            suma.boost2 (trans);       // boost to the CM frame
            double Ecm = suma.t / 2.0; // each nucleon get the same energy in CM
            int licz1=0;
            while(PB&&(licz1<ile_PB))
            {
              licz1++;
              vec dir_cm;

              if (sampling<0.01 && sampling >-0.01)
              {dir_cm = rand_dir () * sqrt(Ecm * Ecm  - MN2);} // randomly set direction of nucleons in CM
              else
              {dir_cm = rand_direc (qq, sampling) * sqrt(Ecm * Ecm  - MN2);} // preferred direction wrt momentum transfer

              nucleon[2].set_momentum (dir_cm);
              nucleon[3].set_momentum (-dir_cm);

              nucleon[2].p4().boost2 (-trans);
              nucleon[3].p4().boost2 (-trans);
              //PB=t.pauli_blocking(nucleon[2])+t.pauli_blocking(nucleon[3]);
              if (nucleon[0].his_fermi<nucleon[2].Ek() && nucleon[0].his_fermi<nucleon[3].Ek() )// PB condition is checked
              {//PB=t.pauli_blocking(nucleon[2])+t.pauli_blocking(nucleon[3]);
                PB=false;
              }
            }
            
          }
        } // to be able to make Lorentz boost and to "decay"
        while  ( PB &&(licz<calls_max) );

        if(licz>=calls_max)
        {
          result=0;
          vect zero (0,0,0);
          nucleon[2].set_momentum (zero);
          nucleon[3].set_momentum (zero);
          nucleon[2].set_energy (MN);
          nucleon[3].set_energy (MN);
        }
        //result*=(1-cosmin);
        meclep.set_momentum(kprim);
      }
      //else result=0;
    }
  }
  return result;
}

////////////////////////////////////////

double Nieves_dsdEdc(double E, double q0, double ct)
{
  double result=0;
  double Ep=E-q0;
  double lp=sqrt(Ep*Ep-ml2);
  double vq2=E*E+lp*lp-2*E*lp*ct;
  if((vq2>q0*q0)&&(vq2>0)&&(fabs(ct)<=1))
  {
    double q=sqrt(vq2);

    // interpolating the nuclear tensor elements
    double w00=0;
    double w03=0;
    double w11=0;
    double w12=0;
    double w33=0;
    q0-=Bmec;
    if((q0<=q)and(q>=0.01*GeV)and(q0>=0)and(q<=qmax_Nieves))
    {
      double spacing=0.01*GeV;
      int m=int((q-spacing)/spacing);
      if (m>118) m=int(118);
      int n=int((q0-spacing)/spacing);
      int pos=int(0.5*m*(m+1)+n)*5;
      // q0-B<0.1 GeV extrapolation with 0 at q0-B=0
      if(n<0)
      {
        int a=pos;
        int b=a+5;
        int c=a+5*(m+1);
        int d=c+5;
        double H1=q0-spacing-n*spacing;
        double H2=q-spacing*m-spacing;
        switch(nucl)
        {
          case 0:
          {
            w00=(H2*(H1*Nieves_C12grid[d])/spacing+(spacing-H2)*(H1*Nieves_C12grid[b])/spacing)/spacing;
            w03=(H2*(H1*Nieves_C12grid[d+1])/spacing+(spacing-H2)*(H1*Nieves_C12grid[b+1])/spacing)/spacing;
            w11=(H2*(H1*Nieves_C12grid[d+2])/spacing+(spacing-H2)*(H1*Nieves_C12grid[b+2])/spacing)/spacing;
            w12=(H2*(H1*Nieves_C12grid[d+3])/spacing+(spacing-H2)*(H1*Nieves_C12grid[b+3])/spacing)/spacing;
            w33=(H2*(H1*Nieves_C12grid[d+4])/spacing+(spacing-H2)*(H1*Nieves_C12grid[b+4])/spacing)/spacing;
            break;
          }
          case 1:
          {
            w00=(H2*(H1*Nieves_O16grid[d])/spacing+(spacing-H2)*(H1*Nieves_O16grid[b])/spacing)/spacing;
            w03=(H2*(H1*Nieves_O16grid[d+1])/spacing+(spacing-H2)*(H1*Nieves_O16grid[b+1])/spacing)/spacing;
            w11=(H2*(H1*Nieves_O16grid[d+2])/spacing+(spacing-H2)*(H1*Nieves_O16grid[b+2])/spacing)/spacing;
            w12=(H2*(H1*Nieves_O16grid[d+3])/spacing+(spacing-H2)*(H1*Nieves_O16grid[b+3])/spacing)/spacing;
            w33=(H2*(H1*Nieves_O16grid[d+4])/spacing+(spacing-H2)*(H1*Nieves_O16grid[b+4])/spacing)/spacing;
            break;
          }
          case 2:
          {
            w00=(H2*(H1*Nieves_Ca40grid[d])/spacing+(spacing-H2)*(H1*Nieves_Ca40grid[b])/spacing)/spacing;
            w03=(H2*(H1*Nieves_Ca40grid[d+1])/spacing+(spacing-H2)*(H1*Nieves_Ca40grid[b+1])/spacing)/spacing;
            w11=(H2*(H1*Nieves_Ca40grid[d+2])/spacing+(spacing-H2)*(H1*Nieves_Ca40grid[b+2])/spacing)/spacing;
            w12=(H2*(H1*Nieves_Ca40grid[d+3])/spacing+(spacing-H2)*(H1*Nieves_Ca40grid[b+3])/spacing)/spacing;
            w33=(H2*(H1*Nieves_Ca40grid[d+4])/spacing+(spacing-H2)*(H1*Nieves_Ca40grid[b+4])/spacing)/spacing;
            break;
          }
        }
      } // now we are in the data grid
      else if(n<m)
      {
        int a=pos;
        int b=a+5;
        int c=a+5*(m+1);
        int d=c+5;
        double H1=q0-spacing-n*spacing;
        double H2=q-spacing*m-spacing;
        switch(nucl)
        {
          case 0:
          {
            w00=(H2*(H1*Nieves_C12grid[d]+(spacing-H1)*Nieves_C12grid[c])/spacing+(spacing-H2)*(H1*Nieves_C12grid[b]+(spacing-H1)*Nieves_C12grid[a])/spacing)/spacing;
            w03=(H2*(H1*Nieves_C12grid[d+1]+(spacing-H1)*Nieves_C12grid[c+1])/spacing+(spacing-H2)*(H1*Nieves_C12grid[b+1]+(spacing-H1)*Nieves_C12grid[a+1])/spacing)/spacing;
            w11=(H2*(H1*Nieves_C12grid[d+2]+(spacing-H1)*Nieves_C12grid[c+2])/spacing+(spacing-H2)*(H1*Nieves_C12grid[b+2]+(spacing-H1)*Nieves_C12grid[a+2])/spacing)/spacing;
            w12=(H2*(H1*Nieves_C12grid[d+3]+(spacing-H1)*Nieves_C12grid[c+3])/spacing+(spacing-H2)*(H1*Nieves_C12grid[b+3]+(spacing-H1)*Nieves_C12grid[a+3])/spacing)/spacing;
            w33=(H2*(H1*Nieves_C12grid[d+4]+(spacing-H1)*Nieves_C12grid[c+4])/spacing+(spacing-H2)*(H1*Nieves_C12grid[b+4]+(spacing-H1)*Nieves_C12grid[a+4])/spacing)/spacing;
            break;
          }
          case 1:
          {
            w00=(H2*(H1*Nieves_O16grid[d]+(spacing-H1)*Nieves_O16grid[c])/spacing+(spacing-H2)*(H1*Nieves_O16grid[b]+(spacing-H1)*Nieves_O16grid[a])/spacing)/spacing;
            w03=(H2*(H1*Nieves_O16grid[d+1]+(spacing-H1)*Nieves_O16grid[c+1])/spacing+(spacing-H2)*(H1*Nieves_O16grid[b+1]+(spacing-H1)*Nieves_O16grid[a+1])/spacing)/spacing;
            w11=(H2*(H1*Nieves_O16grid[d+2]+(spacing-H1)*Nieves_O16grid[c+2])/spacing+(spacing-H2)*(H1*Nieves_O16grid[b+2]+(spacing-H1)*Nieves_O16grid[a+2])/spacing)/spacing;
            w12=(H2*(H1*Nieves_O16grid[d+3]+(spacing-H1)*Nieves_O16grid[c+3])/spacing+(spacing-H2)*(H1*Nieves_O16grid[b+3]+(spacing-H1)*Nieves_O16grid[a+3])/spacing)/spacing;
            w33=(H2*(H1*Nieves_O16grid[d+4]+(spacing-H1)*Nieves_O16grid[c+4])/spacing+(spacing-H2)*(H1*Nieves_O16grid[b+4]+(spacing-H1)*Nieves_O16grid[a+4])/spacing)/spacing;
            break;
          }
          case 2:
          {
            w00=(H2*(H1*Nieves_Ca40grid[d]+(spacing-H1)*Nieves_Ca40grid[c])/spacing+(spacing-H2)*(H1*Nieves_Ca40grid[b]+(spacing-H1)*Nieves_Ca40grid[a])/spacing)/spacing;
            w03=(H2*(H1*Nieves_Ca40grid[d+1]+(spacing-H1)*Nieves_Ca40grid[c+1])/spacing+(spacing-H2)*(H1*Nieves_Ca40grid[b+1]+(spacing-H1)*Nieves_Ca40grid[a+1])/spacing)/spacing;
            w11=(H2*(H1*Nieves_Ca40grid[d+2]+(spacing-H1)*Nieves_Ca40grid[c+2])/spacing+(spacing-H2)*(H1*Nieves_Ca40grid[b+2]+(spacing-H1)*Nieves_Ca40grid[a+2])/spacing)/spacing;
            w12=(H2*(H1*Nieves_Ca40grid[d+3]+(spacing-H1)*Nieves_Ca40grid[c+3])/spacing+(spacing-H2)*(H1*Nieves_Ca40grid[b+3]+(spacing-H1)*Nieves_Ca40grid[a+3])/spacing)/spacing;
            w33=(H2*(H1*Nieves_Ca40grid[d+4]+(spacing-H1)*Nieves_Ca40grid[c+4])/spacing+(spacing-H2)*(H1*Nieves_Ca40grid[b+4]+(spacing-H1)*Nieves_Ca40grid[a+4])/spacing)/spacing;
            break;
          }
        }
      } // an now we are in the external triangle near q0-B=q
      else if(n==m)
      {
        int a=pos;
        int c=a+5*(m+1);
        int d=c+5;
        double H=q0+(m-n+1)*spacing-q;
        double H1=((m+2)*spacing-q);
        switch(nucl)
        {
          case 0:
          {
            w00=(H1*(H*Nieves_C12grid[a]+(spacing-H)*Nieves_C12grid[c])/spacing+(H-H1)*(H*Nieves_C12grid[d]+(spacing-H)*Nieves_C12grid[c])/spacing)/H;
            w03=(H1*(H*Nieves_C12grid[a+1]+(spacing-H)*Nieves_C12grid[c+1])/spacing+(H-H1)*(H*Nieves_C12grid[d+1]+(spacing-H)*Nieves_C12grid[c+1])/spacing)/H;
            w11=(H1*(H*Nieves_C12grid[a+2]+(spacing-H)*Nieves_C12grid[c+2])/spacing+(H-H1)*(H*Nieves_C12grid[d+2]+(spacing-H)*Nieves_C12grid[c+2])/spacing)/H;
            w12=(H1*(H*Nieves_C12grid[a+3]+(spacing-H)*Nieves_C12grid[c+3])/spacing+(H-H1)*(H*Nieves_C12grid[d+3]+(spacing-H)*Nieves_C12grid[c+3])/spacing)/H;
            w33=(H1*(H*Nieves_C12grid[a+4]+(spacing-H)*Nieves_C12grid[c+4])/spacing+(H-H1)*(H*Nieves_C12grid[d+4]+(spacing-H)*Nieves_C12grid[c+4])/spacing)/H;
            break;
          }     
          case 1:
          {
            w00=(H1*(H*Nieves_O16grid[a]+(spacing-H)*Nieves_O16grid[c])/spacing+(H-H1)*(H*Nieves_O16grid[d]+(spacing-H)*Nieves_O16grid[c])/spacing)/H;
            w03=(H1*(H*Nieves_O16grid[a+1]+(spacing-H)*Nieves_O16grid[c+1])/spacing+(H-H1)*(H*Nieves_O16grid[d+1]+(spacing-H)*Nieves_O16grid[c+1])/spacing)/H;
            w11=(H1*(H*Nieves_O16grid[a+2]+(spacing-H)*Nieves_O16grid[c+2])/spacing+(H-H1)*(H*Nieves_O16grid[d+2]+(spacing-H)*Nieves_O16grid[c+2])/spacing)/H;
            w12=(H1*(H*Nieves_O16grid[a+3]+(spacing-H)*Nieves_O16grid[c+3])/spacing+(H-H1)*(H*Nieves_O16grid[d+3]+(spacing-H)*Nieves_O16grid[c+3])/spacing)/H;
            w33=(H1*(H*Nieves_O16grid[a+4]+(spacing-H)*Nieves_O16grid[c+4])/spacing+(H-H1)*(H*Nieves_O16grid[d+4]+(spacing-H)*Nieves_O16grid[c+4])/spacing)/H;
            break;
          }
          case 2:
          {
            w00=(H1*(H*Nieves_Ca40grid[a]+(spacing-H)*Nieves_Ca40grid[c])/spacing+(H-H1)*(H*Nieves_Ca40grid[d]+(spacing-H)*Nieves_Ca40grid[c])/spacing)/H;
            w03=(H1*(H*Nieves_Ca40grid[a+1]+(spacing-H)*Nieves_Ca40grid[c+1])/spacing+(H-H1)*(H*Nieves_Ca40grid[d+1]+(spacing-H)*Nieves_Ca40grid[c+1])/spacing)/H;
            w11=(H1*(H*Nieves_Ca40grid[a+2]+(spacing-H)*Nieves_Ca40grid[c+2])/spacing+(H-H1)*(H*Nieves_Ca40grid[d+2]+(spacing-H)*Nieves_Ca40grid[c+2])/spacing)/H;
            w12=(H1*(H*Nieves_Ca40grid[a+3]+(spacing-H)*Nieves_Ca40grid[c+3])/spacing+(H-H1)*(H*Nieves_Ca40grid[d+3]+(spacing-H)*Nieves_Ca40grid[c+3])/spacing)/H;
            w33=(H1*(H*Nieves_Ca40grid[a+4]+(spacing-H)*Nieves_Ca40grid[c+4])/spacing+(H-H1)*(H*Nieves_Ca40grid[d+4]+(spacing-H)*Nieves_Ca40grid[c+4])/spacing)/H;
            break;
          }
        }
      }
    }
    q0+=Bmec;
    double w1=0.5*w11;
    double q02=q0*q0;
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
    if (result<0) { result=0;}
  }
  return result;
}
