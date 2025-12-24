#include "kaskada7.h"
#include "fsi.h"

// Constructor: set up cascade state from params, event and input_data
kaskada::kaskada(params &p, event &e1, input_data *input)
{
  par = p;
  if (par.nucleus_p + par.nucleus_n < 3) par.kaskada_w = 0;  //for free nucleons and deuteron there is no extra binding energy

  e = &e1;                        // Pointer to current event
  max_step = par.step * fermi;    // set maximum step defined in params
  nucl = make_nucleus(par);       // create nucleus defined in params
  radius = nucl->radius();        // calculate radius of the nucleus

  I = new Interaction(input->get_data_container(0), input->get_data_container(1), input->get_data_container(2),par.kaskada_piN_xsec);
  corr_func = input->get_nucl_data_container(0, par.nucleus_p, par.nucleus_n);
}

// Destructor: free heap-allocated nuclear and interaction objects
kaskada::~kaskada()
{
  delete nucl;
  delete I;
}

// Main cascade driver for a single event
int kaskada::kaskadaevent()
{

  // Global counter of events processed by this cascade instance
  static long long totalEvents = 0;
  totalEvents++;

  // Book-keeping for SF QE topology split
  static int count_transp_corr = 0;
  static int count_transp_uncorr = 0;
  static int count_nontransp_corr = 0;
  static int count_nontransp_uncorr = 0;

if(e->flag.qel and par.sf_method != 0) {

if (e->flag.isTransparent && e->flag.isCorrelated) {count_transp_corr++;}
else if (e->flag.isTransparent && !e->flag.isCorrelated) {count_transp_uncorr++;}
else if (!e->flag.isTransparent && e->flag.isCorrelated) {count_nontransp_corr++;}
else if (!e->flag.isTransparent && !e->flag.isCorrelated) {count_nontransp_uncorr++;}

    // Uncomment block below if you want to periodically print SF topology fractions
    if ((totalEvents % 10000) == 0) {
        double total = totalEvents > 0 ? totalEvents : 1;
        /*std::cout << std::fixed << std::setprecision(2);
        std::cout << "\n[DIAG] After " << totalEvents << " events:\n"
                  << "   Transparent + Correlated     : " << 100.0 * count_transp_corr / total << " %\n"
                  << "   Transparent + Uncorrelated   : " << 100.0 * count_transp_uncorr / total << " %\n"
                  << "   Non-transparent + Correlated : " << 100.0 * count_nontransp_corr / total << " %\n"
                  << "   Non-transparent + Uncorr.    : " << 100.0 * count_nontransp_uncorr / total << " %\n"
                  << std::endl;*/
    }
}

  int result = 0;
  bool nucleon_interaction_occurred = false;
  bool firstNucleonInteracted = false;
  bool needs_final_clean = true;

  if (e->weight <= 0) return result;

  prepare_particles(); // Fill 'parts' queue from primary vertex (e->out), apply formation zone

  if (e->in[0].lepton()) {

    for(int i=1;i< e->in.size();i++)
      nucl->remove_nucleon(e->in[i]);

  }

  // FSI OFF: only check if particles leave nucleus or are jailed, no scatters
  if (par.FSI_on == 0) {

    while (parts.size () > 0) {

      particle p1 = parts.front();    // point a particle from a queue
      parts.pop();                    // remove this particle from a temp vector
      p = &p1;                        // working pointer `p` refers to this particle

      leave_nucleus();                // check if the particle is jailed or escapes (and returns on the mass shell)

    }

    return result;

  }

  // Non-QE events Or QE events using LFG / GFG (sf_method == 0):
  //   → run the classical NuWro cascade for all 'parts'
  if ( (e->flag.qel && (par.sf_method == 0)) || (!e->flag.qel) ) {

      while (!parts.empty() && nucl->Ar() > 0) { // while queue has particles and nucleus not exhausted

        particle p1 = parts.front();
        parts.pop();
        p = &p1;

        X = prepare_interaction();
        if (!move_particle()) continue;
        if (X.r >= radius) leave_nucleus();

        else {

          if (max_step < X.freepath || !make_interaction() || !finalize_interaction()) {

            if (nucleon(p->pdg)) e->nod[13]++;
            if (pion(p->pdg))    e->nod[12]++;
            if (hyperon(p->pdg)) e->nod[14]++;
            parts.push(*p);

          }
        }
      }
    }

  // QE events using Spectral Function (sf_method != 0):
  //   → treat transparent vs non-transparent, correlated vs uncorrelated
  // FSI implemented consistently for QE-SF channel by RWIK DHARMAPAL BANERJEE, 2025
  if (e->flag.qel && (par.sf_method != 0)) {

      // Transparent events
      if (e->flag.isTransparent)
      {

        // Transparent & Uncorrelated: just leave nucleus
        if (!e->flag.isCorrelated)
        {

           while (!parts.empty() && nucl->Ar()>0) {

           particle p1 = parts.front();
           parts.pop();
           p = &p1;

           leave_nucleus();

           }
        }

        // Transparent & Correlated: first nucleon leaves, second gets standard cascade
        else if (e->flag.isCorrelated)
        {

          while (!parts.empty() && nucl->Ar()>0)
          {

            particle p1 = parts.front();
            parts.pop();
            p = &p1;

            if (p1.nucleon_id==1) {
              // Primary nucleon: fully transparent, no interactions
              leave_nucleus();
            }

            else {
              // spectator nucleon: standard cascade
              X = prepare_interaction();
              if (!move_particle()) continue;
              if (X.r >= radius) leave_nucleus();

              else if ( max_step < X.freepath || !make_interaction() || !finalize_interaction() ) {

                   if (nucleon(p1.pdg)) e->nod[13]++;
                   if (pion(p1.pdg))    e->nod[12]++;
                   if (hyperon(p1.pdg)) e->nod[14]++;
                   parts.push(*p);

             }
            }
          }
        }
      }

      // Non-Transparent events
      else if (!e->flag.isTransparent)
      {

        // Non-Transparent & Uncorrelated: full cascade with possible re-draw if no interaction
        if (!e->flag.isCorrelated)
        {

          // Standard cascade loop
          while (!parts.empty() && nucl->Ar() > 0) {

            particle p1 = parts.front();
            parts.pop();
            p = &p1;

            X = prepare_interaction();
            if (!move_particle()) continue;
            if (X.r >= radius) leave_nucleus();

            else {

              if ( max_step < X.freepath || !make_interaction() || !finalize_interaction() ) {

                if (nucleon(p1.pdg)) e->nod[13]++;
                if (pion(p1.pdg))    e->nod[12]++;
                if (hyperon(p1.pdg)) e->nod[14]++;
                parts.push(p1);

              }
              else {

                if (nucleon(p1.pdg) and (p1.nucleon_id==1)) {
                  // Flag that the struck nucleon (id==1) did undergo at least one interaction
                  nucleon_interaction_occurred = true;
                }
              }
             }
            }

          // If no struck nucleon interaction → redraw cascade with reduced mean free path
          if (!nucleon_interaction_occurred)
          {

            int    redrawAttempts = 0;
            // Keep track of mean free path scaling for retries
            double original_scale = par.kaskada_NN_mfp_scale;
            double effective_scale = original_scale;

          this->U_evt = (par.U_switch == 1 and par.sf_method != 0) ? e->optical_potential : 0.0;

          nucleon_interaction_occurred = false;
          while(!parts.empty()) parts.pop();

          while (redrawAttempts < nucl->MAX_EVENT_REDRAWS && !nucleon_interaction_occurred) {

          // Clean up post-state: remove non-leptons so we can retry cascade
          e->post.erase(
            std::remove_if(e->post.begin(), e->post.end(),[](particle& pt){ return !lepton(pt.pdg);}),
          e->post.end());

              // Make interactions more probable (shorter mfp) each redraw
              effective_scale *= nucl->effective_mfp_scale;
              par.kaskada_NN_mfp_scale = effective_scale;

              // New copy of N1 from primary vertex
              particle N1 = e->out[1];
              N1.travelled = 0;
              prepare_single_nucleon_for_redraw(N1, 1);

              // Re-run cascade
              while (!parts.empty() && nucl->Ar() > 0) {

              particle pt = parts.front();
              parts.pop();
              p = &pt;

              X = prepare_interaction();
              if (!move_particle()) continue;
              if (X.r >= radius) leave_nucleus();

              else {

                 if (max_step < X.freepath || !make_interaction() || !finalize_interaction()) {

                     if (nucleon(pt.pdg)) e->nod[13]++;
                     if (pion(pt.pdg))    e->nod[12]++;
                     if (hyperon(pt.pdg)) e->nod[14]++;
                     parts.push(pt);

                 }
                 else {

                    if (nucleon(pt.pdg) and (pt.nucleon_id==1)) nucleon_interaction_occurred = true;

                 }
                }
               }

               redrawAttempts++;

              }
              
       // Restore original mfp scale after redraw loop
       par.kaskada_NN_mfp_scale = original_scale;

       }

      needs_final_clean = false;

      }

        // Non-Transparent & Correlated: both nucleons can interact, redraw if N1 fails to interact
        else if (e->flag.isCorrelated)
        {

    // First pass: standard cascade for both nucleons and secondaries
    while (!parts.empty() && nucl->Ar()>0) {

      particle p1 = parts.front();
      parts.pop();
      p = &p1;

      X = prepare_interaction();
      if (!move_particle()) continue;
      if (X.r >= radius) leave_nucleus();

      else {

      if (max_step < X.freepath || !make_interaction() || !finalize_interaction()) {

        if (nucleon(p1.pdg))  e->nod[13]++;
        if (pion (p1.pdg))    e->nod[12]++;
        if (hyperon (p1.pdg)) e->nod[14]++;
        parts.push(p1);

      }
      else {

          if (nucleon(p1.pdg) and (p1.nucleon_id == 1)) {
          
          // Record if N1 (id=1) has interacted at least once
          firstNucleonInteracted = true;
          
          }
        }
      }
    }
            // If no interaction for N1 → clean post and redraw only N1
            if (!firstNucleonInteracted) {

              int redrawAttempts = 0;

              double cosine_N1 = e->out[1].z / e->out[1].momentum();

              double original_scale = par.kaskada_NN_mfp_scale;
              double effective_scale = original_scale;

              this->U_evt = (par.U_switch == 1 and par.sf_method != 0) ? e->optical_potential : 0.0;

              nucleon_interaction_occurred = false;
              while(!parts.empty()) parts.pop();

              // Redraw attempts only for N1 with progressively shorter mfp
              while (redrawAttempts < nucl->MAX_EVENT_REDRAWS && !firstNucleonInteracted) {

                // Optional diagnostic: which nucleons we are about to remove
                std::vector<particle> removedN1; 
                std::copy_if( e->post.begin(), e->post.end(), std::back_inserter(removedN1), [&](particle &pt)
                { return nucleon(pt.pdg) && fabs(pt.z / pt.momentum() - cosine_N1) < nucl->cosine_threshold;});

                // Remove nucleons close in direction to the original N1
                e->post.erase( std::remove_if( e->post.begin(), e->post.end(), [&](particle &pt) {
                if (lepton(pt.pdg) || pion(pt.pdg) || hyperon(pt.pdg)) return false; // keep these
                if (nucleon(pt.pdg)) { return (fabs(pt.z / pt.momentum() - cosine_N1) < nucl->cosine_threshold);} // remove these
                return false;}), e->post.end());

                // Strengthen scattering (shorter mfp)
                effective_scale *= nucl->effective_mfp_scale;
                par.kaskada_NN_mfp_scale = effective_scale;

                // New N1
                particle N1 = e->out[1];
                N1.travelled = 0;
                prepare_single_nucleon_for_redraw(N1, 1);

                while (!parts.empty() && nucl->Ar() > 0) {

                particle pt = parts.front();
                parts.pop();
                p = &pt;

                X = prepare_interaction();
                if (!move_particle()) continue;
                if (X.r >= radius) leave_nucleus();

                else {

                   if (max_step < X.freepath || !make_interaction() || !finalize_interaction()) {

                       if (nucleon(pt.pdg)) e->nod[13]++;
                       if (pion(pt.pdg))    e->nod[12]++;
                       if (hyperon(pt.pdg)) e->nod[14]++;
                       parts.push(pt);

                   }
                   else {

                       if (nucleon(pt.pdg) and (pt.nucleon_id == 1)) firstNucleonInteracted = true;

                     }
                   }
                 }

          redrawAttempts++;

        }

       par.kaskada_NN_mfp_scale = original_scale;

      }

      needs_final_clean = false;

    }
   }
  }

  if (needs_final_clean) clean();
  
  // Update remaining protons/neutrons in nucleus
  e->pr=nucl->Zr();
  e->nr=nucl->Nr();

  return result;

}

// Build initial queue "parts" from event primary vertex e->out
void kaskada::prepare_particles()
{
  static long long total_events_global = 0;
  static long long events_with_1 = 0;
  static long long events_with_2 = 0;
  bool has_id1 = false;
  bool has_id2 = false;

  // Loop over every particle from the primary vertex.
  for (int i = 0; i < e->out.size(); i++)
  {
    particle p1 = e->out[i];

    if (nucleon_or_pion (p1.pdg)) // formation zone for both nucleons and pions
    {
      if (nucleon (p1.pdg))
      {
        p1.primary = true;

        double U = (e->flag.qel and par.U_switch == 1 and par.sf_method != 0)
                    ? ((par.FSI_on == 1)
                        ? e->optical_potential
                        : (p1.pdg == pdg_proton ? -e->averageCE : 0.0))
                    : 0.0;

        double kaskada_w = (e->flag.qel and par.U_switch == 1 and par.sf_method != 0) ? 0.0 : par.kaskada_w;

        // Tag primary and spectator nucleons form SF-event vertex
        if(e->flag.qel and par.sf_method != 0) {
         if (i == 1) {
            p1.nucleon_id = 1;
            has_id1 = true;
         }
         else if (i == 2) {
            p1.nucleon_id = 2;
            has_id2 = true;
         }
      }

        // add energy substracted in the primary vertex for GFG LFG and SF ( U is negative )
        if (e->flag.qel and (par.sf_method != 0 or par.nucleus_target == 2)) {
            p1.set_energy (p1.E() + nucl->Ef(p1) + kaskada_w - U);
        }

        else if (par.nucleus_target == 1 and (e->flag.qel or e->flag.res))
          p1.set_energy (p1.E() + par.nucleus_E_b);

        p1.set_fermi(nucl->Ef(p1));

        // If kinetic energy is below "barrier" = Ef + kaskada_w - U, jail back to nucleus
        if (p1.Ek() <= kaskada_w + p1.his_fermi - U)  
        {
          p1.endproc=jailed;
          nucl->insert_nucleon (p1);
          if(par.kaskada_writeall) e->all.push_back(p1);
          continue;
        }
      }

      double fz = formation_zone(p1, par, *e);        // calculate formation zone
      p1.krok(fz);      // move particle by a distance defined by its formation zone
      
      parts.push (p1);  // put particle to a queue
      
    }

    else if(hyperon (p1.pdg)) // if a hyperon
    {
      // Add BE
      if(par.nucleus_target) p1.set_fermi(nucl->hyp_BE(p1.r.length(),p1.pdg));
      else p1.set_fermi(0);
        // a new piece of code March 3, 2025 - fixing hyperon potential problem
        if (p1.E() + p1.his_fermi < p1.mass())
        {
            p1.endproc=escape;
            e->post.push_back (p1);

            if(par.kaskada_writeall)
                e->all.push_back(p1);
        }
         else
         {
             p1.set_energy(p1.E() + p1.his_fermi);
             parts.push(p1); // add particle to queue
         }
        // the end of the new piece of code

        /* before March 3, 2025 the below two lines were active
      p1.set_energy(p1.E() + p1.his_fermi);

      parts.push(p1); // add particle to queue
        */
    }
    else              // if not a nucleon, pion nor hyperon
    {
      p1.endproc=escape;
      e->post.push_back (p1);
      if(par.kaskada_writeall) e->all.push_back(p1);
    }
  }

  for (int i = 0; i<18; i++)  // number of dynamics defined in proctable.h
     e->nod[i] = 0;

  // Diagnostics: how often we see 1 vs 2 tagged nucleons
  total_events_global++;
  if (has_id1 && has_id2) events_with_2++;
  else if (has_id1 || has_id2) events_with_1++;
  // Uncomment for debugging distribution of 1-tag vs 2-tag events
  if (total_events_global % 10000 == 0) {
      double perc1 = 100.0 * events_with_1 / total_events_global;
      double perc2 = 100.0 * events_with_2 / total_events_global;
  /*  std::cout << std::fixed << std::setprecision(2)
              << "\n[DIAG] After " << total_events_global << " events:\n"
              << "  " << perc1 << "% had 1 tagged nucleon\n"
              << "  " << perc2 << "% had 2 tagged nucleons\n";*/
  }

  e->r_distance = 10;         // new JS ; default (large) value, if unchanged no absorption took place

}

// Prepare local interaction parameters at current particle position:
interaction_parameters kaskada::prepare_interaction()
{
  interaction_parameters res;

  res.pdg = p->pdg;
  res.Ek  = p->Ek();
  res.r   = p->r.length ();

  res.dens = nucl->density (res.r);
  assert(res.dens>=0);

  res.dens_n = res.dens * nucl->frac_neutron ();
  res.dens_p = res.dens * nucl->frac_proton ();
  res.n = 2;

  I->total_cross_sections (*p, *nucl, res); // calculate cross sections xsec_p and xsec_n // rwik

  corr_func->set_input_point( p->travelled );
  double corr_ii = corr_func->get_value( 1 );
  double corr_ij = corr_func->get_value( 2 );
  corr_func->set_input_point( p->r.length());
  double norm_ii = corr_func->get_value( 3 );
  double norm_ij = corr_func->get_value( 4 );

  if( !e->out[0].nucleon() || e->number_of_interactions() || par.beam_placement != 2 )
  // no correlations for incoming nucleons in the scattering mode (kaskada.cc)
  {
    switch (res.pdg) // mean free path modifications: effective density and scaling
    {
      case pdg_neutron:
        res.xsec_n *= corr_ii / norm_ii / par.kaskada_NN_mfp_scale;
        res.xsec_p *= corr_ij / norm_ij / par.kaskada_NN_mfp_scale;
      case pdg_proton:
        res.xsec_n *= corr_ij / norm_ij / par.kaskada_NN_mfp_scale;
        res.xsec_p *= corr_ii / norm_ii / par.kaskada_NN_mfp_scale;
        break;
    }
  }

  res.xsec = res.dens_n*res.xsec_n + res.dens_p*res.xsec_p; // calculate the inverse of the mean free path

  assert(res.xsec>=0);                      // make sure that the cross section is positive

  if (res.xsec != 0)
  {
    res.freepath = -log (frandom ()) / res.xsec; // choose free path according to the mean free path (1/res.xsec)
    res.prob_proton = res.xsec_p * res.dens_p / res.xsec;
  }
  else
    res.freepath = 2.0 * max_step;

  return res;
}

// Move particle along free path or max_step, and possibly jail nucleons
bool kaskada::move_particle()
{
  p->krok (min (max_step, X.freepath));   // propagate by no more than max_step

  if (!(p->nucleon() || hyperon(p->pdg))) // pion can not be jailed
    return true;

  if(p->nucleon())
  {

    double U = (e->flag.qel and par.U_switch == 1 and par.sf_method != 0)
                ? ((par.FSI_on == 1)
                    ? e->optical_potential
                    : (p->pdg == pdg_proton ? -e->averageCE : 0.0))
                : 0.0;

    double kaskada_w = (e->flag.qel and par.U_switch == 1 and par.sf_method != 0) ? 0.0 : par.kaskada_w;

    // jail nucleon if its kinetic energy is lower than "binding" energy
    if (p->Ek() <= kaskada_w + p->his_fermi - U)
    {
      p->endproc=jailed;
      nucl->insert_nucleon (*p);
      if(par.kaskada_writeall) e->all.push_back(*p);
      return false; // nucleon was jailed
    }
    else
      return true;
  }

  if(hyperon(p->pdg))
  {
    // If hyperon KE falls below its potential energy -  hyperon is jailed
    if(p->Ek() < p->his_fermi)
    {
      p->endproc=jailed;
      if(par.kaskada_writeall) e->all.push_back(*p);
      return false; // hyperon was jailed
    }
    return true;
  }

  return true;
}

// Decide if particle escapes or is jailed when reaching nuclear boundary and apply the corresponding energy shift
bool kaskada::leave_nucleus()
{

  // Nucleons: do final "climb out of potential well" and possible jailing
  if (nucleon (p->pdg))
  {

    double U = (e->flag.qel and par.U_switch == 1 and par.sf_method != 0)
                ? ((par.FSI_on == 1)
                    ? e->optical_potential
                    : (p->pdg == pdg_proton ? -e->averageCE : 0.0))
                : 0.0;

    double kaskada_w = (e->flag.qel and par.U_switch == 1 and par.sf_method != 0) ? 0.0 : par.kaskada_w;

    // If KE still below (Ef + kaskada_w - U) at surface, jail nucleon
    if (p->Ek() <= p->his_fermi + kaskada_w - U)
    {
      p->endproc=jailed;
      nucl->insert_nucleon (*p);
      if(par.kaskada_writeall) e->all.push_back(*p);
      return false; // particle did not escape
    }

      // Otherwise, remove Ef and kaskada_w and add U to energy (reverse of entrance)
      p->set_energy(p->E() - p->his_fermi - kaskada_w + U);

  }
  else if(hyperon(p->pdg))
  {
    // Adjust hyperon for remaining binding energy, jail if it cannot escape
    if (p->Ek() < p->his_fermi)
    {
      p->endproc=jailed;
      //TODO: Perhaps add hyperon to nucleus here (or decay it)
      if(par.kaskada_writeall) e->all.push_back(*p);
      return false;
    }
    // Subtract binding energy from hyperon energy and set momentum so it is on shell
    else {
      p->set_energy(p->E() - p->his_fermi);
    }
  }

  p->endproc=escape;
  e->post.push_back (*p);
  if(par.kaskada_writeall) e->all.push_back (*p);
  return true; // particle has escaped
  
}

// Generate scattering kinematics 
bool kaskada::make_interaction()
{
  int loop = 0;
  static int call=0;
  static int rep=0;
  static int procid=0;

  while(++call && 0 == I->particle_scattering (*p, *nucl, X)) // try to generate kinematics
  {
    if(loop==0)
      procid=I->process_id();
    ++rep;++loop;
    double suma=0;
    for(int i=0;i<X.n;i++)
    {
      suma+=X.p[i].mass();
    }
//CJ  cout<<" Interaction: "<<I.process_name()<<" ("<<I.process_id()<<") ";
//    assert(procid==I->process_id());
    if(loop>100)
      return false;   // it was impossible to make kinematics
  }

  for (int i = 0; i < X.n; i++)
    if (!X.p[i].is_valid ())
    {
      cerr << I->process_name()<< "Interaction: error" << X.p[i] << endl;
      delete nucl;
      delete I;
      exit(18);
      return false;
    }

  if(nucl->pauli_blocking (X.p, X.n)) return false;

  // C Thorpe: Check if hyperon can be moved to new value of potential.
  // Ignore interaction if it can't

  if( (PDG::Lambda(p->pdg) && PDG::Sigma(X.p[0].pdg)) || (PDG::Sigma(p->pdg) && PDG::Lambda(X.p[0].pdg)) )
  {
    double V_old = p->his_fermi;
    double V_new = nucl->hyp_BE(p->r.length(),X.p[0].pdg);

    // Check if hyperon may be moved to new value of potential
    if(X.p[0].Ek() < V_old - V_new) return false;

    X.p[0].set_fermi(V_new);
    X.p[0].set_energy(X.p[0].E() - V_old + V_new);
  }
  else if(PDG::hyperon(X.p[0].pdg)) X.p[0].set_fermi(p->his_fermi);

  //JTS - removed assertion below
  //assert(check(*p,X.p2,nucl->spectator,X.n,X.p,I->process_id()));
  //assert(check2(*p,X.p2,nucl->spectator,X.n,X.p,I->process_id()));

  return true;
}

// Finalize the interaction:
bool kaskada::finalize_interaction()
{
  p->endproc=I->process_id();

  double FE = nucl->Ef(X.p2);

  if (!p->nucleon())
  {
    for (int i = 0; i < X.n; i++){
      if (nucleon(X.p[i].pdg))
      {
        X.p[i].set_fermi (FE);
        break;
      }
    }
  }
  else
  {
    double he = (p->his_fermi > FE) ? p->his_fermi : FE;
    double le = p->his_fermi + FE - he;

    int n_he = -1;
    int n_le = -1;

    for (int i = 0; i < X.n; i++)
      if (nucleon(X.p[i].pdg))
      {
        if (n_he < 0)
          n_he = i;
        else if (X.p[i].Ek() < X.p[n_he].Ek())
          n_le = i;
        else
        {
          n_le = n_he;
          n_he = i;
        }
      }
    X.p[n_he].set_fermi (he);
    X.p[n_le].set_fermi (le);
  }

  for (int i = 0; i < X.n; i++)
  {
    X.p[i].r = p->r;
    X.p[i].travelled = 0;

    double U = (e->flag.qel and par.U_switch == 1 and par.sf_method != 0)
                ? ((par.FSI_on == 1)
                    ? e->optical_potential
                    : (X.p[i].pdg == pdg_proton ? -e->averageCE : 0.0))
                : 0.0;

    double kaskada_w = (e->flag.qel and par.U_switch == 1 and par.sf_method != 0) ? 0.0 : par.kaskada_w;

    // jail nucleon if its kinetic energy is lower than work function
    if (nucleon (X.p[i].pdg) and X.p[i].Ek() <= kaskada_w + X.p[i].his_fermi - U)
    {
      X.p[i].endproc=jailed;
      nucl->insert_nucleon (X.p[i]);
      if(par.kaskada_writeall)
        e->all.push_back(X.p[i]);
      continue;
    }
   //hyperon
    else if( hyperon(p->pdg) && hyperon(X.p[i].pdg) )
    {
      if(X.p[i].Ek() < X.p[i].his_fermi)
      {
        X.p[i].endproc=jailed;
        if(par.kaskada_writeall)
          e->all.push_back(X.p[i]);
        continue;
      }
    }

   // Add particle back to queue if not jailed
   parts.push(X.p[i]);

/*
    else
    {
      parts.push (X.p[i]);
      //double fz = formation_zone(X.p[i], par);
      //X.p[i].krok(fz);
    }
*/
    //procinfo(*p,X.p2,X.n,X.p);

    if(par.kaskada_writeall) e->all.push_back (X.p[i]); //bug fixed

  }//loop over p[i]

  int k = kod(I->process_id());
  e->nod[k]++;

  if (k==8)                               // JS absorption
  {
    e->r_distance = p->r.length()/fermi;  // making length of r
  }

  if(!nucl->remove_nucleon (X.p2))
    return false;                         // remove from the nuclear matter
  if(nucl->spectator!=NULL)
    if(!nucl->remove_nucleon (*nucl->spectator))
    {
      nucl->insert_nucleon (X.p2);
      return false;
    }

  return true;
}

// Prepare a single nucleon (typically N1) to be re-run in a redraw scenario:
void kaskada::prepare_single_nucleon_for_redraw(particle N1, int index)
{

  if (!nucleon(N1.pdg)) return;

    particle pN = N1;
    pN.primary = true;

        if (index == 1) pN.nucleon_id = 1;
        else if (index == 2) pN.nucleon_id = 2;
        else pN.nucleon_id = 0;

    // Here U comes from U_evt (set before redraw loop)
    double U = U_evt;
    double kaskada_w = (par.U_switch == 1 and par.sf_method != 0) ? 0.0 : par.kaskada_w;

    // Reconstruct in-medium energy to consider off-shellness
    pN.set_energy (pN.E() + nucl->Ef(pN) + kaskada_w - U);

    pN.set_fermi(nucl->Ef(pN));

    if (pN.Ek() <= kaskada_w + pN.his_fermi - U)
    {
        pN.endproc = jailed;
        nucl->insert_nucleon(pN);
        if (par.kaskada_writeall) e->all.push_back(pN);
        return;
    }

    double fz = formation_zone(pN, par, *e);
    pN.krok(fz);
    parts.push(pN);
}

// Clean remaining 'parts' queue at the end of event 
void kaskada::clean()
{

  while(! parts.empty())
  {
    particle p0 = parts.front();

    if (p0.Ek() >= 0) { p0.endproc = escape; e->post.push_back(p0); }
    else { p0.endproc = jailed; }

    if(par.kaskada_writeall) e->all.push_back(p0);
    parts.pop();
  }

}

// Helper functions

// Consistency check: charge conservation
bool kaskada::check (particle & p1, particle & p2, particle *spect, int n, particle p[],int k)
{
  int ch1=p1.charge()+p2.charge();
  if(spect) ch1+=spect->charge();
  int i=n;
  while(i)
    ch1-=p[--i].charge();
  if(ch1!=0)
  {
    cout<<endl<<"Proc:"<<k<<endl;
    cout<<p1<<endl<<p2<<endl;
    if(spect)
      cout<<(*spect)<<endl;
    cout<<endl;
    while(i<n)
      cout<<p[i++]<<endl;
    cout<<endl;
  }
  return ch1==0;
}

// Consistency check: four-momentum conservation
bool kaskada::check2 (particle & p1, particle & p2, particle *spect, int n, particle p[],int k)
{
  vect p4=p1+p2;
  if(spect) p4+=*spect;
  int i=n;
  while(i)
    p4-=p[--i];
  double prec=0.001*MeV;
  if(abs(p4.t)>prec ||abs(p4.x)>prec ||abs(p4.y)>prec ||abs(p4.z)>prec)
  {
    cout<<endl<<"Proc:"<<k<<" delta p4 = "<<p4<<endl;
    cout<<p1<<endl<<p2<<endl;
    if(spect)
      cout<<(*spect)<<endl;
    cout<<endl;
    while(i<n)
      cout<<p[i++]<<endl;
    cout<<endl;
    return false;
  }
  return true;
}

// Plcae anywhere in kaskadaevent() to get a list of particles at the moment in a vector
void kaskada::printPDGCounts(const std::vector<particle>& vec,const std::string &vecName,const event* e)
{
    std::map<int,int> counts;
    for (const auto &p : vec) {
        counts[p.pdg]++;
    }

    std::cout << "[DEBUG] Vector '" << vecName << "' size = "
              << vec.size() << ". PDG distribution:" << std::endl;

    for (const auto &kv : counts) {
        std::cout << "    PDG " << kv.first << ": " << kv.second << std::endl;
    }

    if (e) {
        std::cout << "  [FLAGS] isTransparent="
                  << e->flag.isTransparent
                  << ", isCorrelated="
                  << e->flag.isCorrelated
                  << std::endl;
    } else {
        std::cout << "  [FLAGS] (no event info provided)" << std::endl;
    }
}
// Plcae anywhere in kaskadaevent() to get a list of particles at the moment in a queue
void kaskada::printPartsPDGCounts(const std::queue<particle>& q,const std::string &name,const event* e)
{
    std::queue<particle> temp = q;  // copy the queue
    std::map<int, int> pdg_counts;
    std::map<int, int> id_counts;
    int total_particles = 0;

    while (!temp.empty()) {
        const particle& p = temp.front();
        temp.pop();
        pdg_counts[p.pdg]++;
        id_counts[p.nucleon_id]++;
        total_particles++;
    }

    std::cout << "[DEBUG] Queue '" << name << "' contains "
              << total_particles << " particles." << std::endl;

    std::cout << "  PDG distribution:" << std::endl;
    for (const auto& kv : pdg_counts)
        std::cout << "    PDG " << kv.first << ": " << kv.second << std::endl;

    std::cout << "  nucleon_id distribution:" << std::endl;
    for (const auto& kv : id_counts)
        std::cout << "    id=" << kv.first << ": " << kv.second << std::endl;

    if (e) {
        std::cout << "  [FLAGS] isTransparent="
                  << e->flag.isTransparent
                  << ", isCorrelated="
                  << e->flag.isCorrelated
                  << std::endl;
    } else {
        std::cout << "  [FLAGS] (no event info provided)" << std::endl;
    }
}
