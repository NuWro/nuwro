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
  shell_dist = nullptr;
}

// Destructor: free heap-allocated nuclear and interaction objects
kaskada::~kaskada()
{
  delete nucl;
  delete I;
}

// Main cascade driver for a single event
int kaskada::kaskadaevent(bool bare_kaskada)
{

  int result = 0;
  bool intercated = false;
  bool needs_final_clean = true;
  const bool isQESF = (e->flag.qel && par.sf_method != 0 && !bare_kaskada);

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

    // Update remaining protons/neutrons in nucleus
    e->pr=nucl->Zr();
    e->nr=nucl->Nr();

    return result;

  }

  // Non-QE events Or QE events using LFG / GFG (sf_method == 0):
  //   → run the classical NuWro cascade for all 'parts'
  if (!isQESF) {

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
  // FSI implemented consistently for QE-SF channel by RWIK DHARMAPAL BANERJEE, 2025-26
  if (isQESF) {

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
        else
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
      else
      {

        // Non-Transparent & Uncorrelated: full cascade with possible re-draw if no interaction
        if (!e->flag.isCorrelated)
        {
          int nod_before_ntuc[18];
          for (int nod_i = 0; nod_i < 18; ++nod_i) {
            nod_before_ntuc[nod_i] = e->nod[nod_i];
          }

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
                  intercated = true;
                }
              }
             }
            }

          // If no struck nucleon interaction → redraw cascade with reduced mean free path
          if (!intercated)
          {

            int redrawAttempts = 0;
            // Keep track of mean free path scaling for retries
            double original_scale = par.kaskada_NN_mfp_scale;
            double effective_scale = original_scale;

          // safeguard
          while(!parts.empty()) parts.pop();

          while (redrawAttempts < nucl->MAX_EVENT_REDRAWS && !intercated) {

          redrawAttempts++;

          for (int nod_i = 0; nod_i < 18; ++nod_i) {
            e->nod[nod_i] = nod_before_ntuc[nod_i];
          }

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

                    if (nucleon(pt.pdg) and (pt.nucleon_id==1)) {
                      intercated = true;
                      par.kaskada_NN_mfp_scale = original_scale;
                    }

                 }
                }
               }

               if(!parts.empty()) // safeguard
               {
                 while (!parts.empty()) parts.pop();
               }

              }
              
              if (!intercated) { // fallback
                for (int nod_i = 0; nod_i < 18; ++nod_i) {
                  e->nod[nod_i] = nod_before_ntuc[nod_i];
                }

                e->post.erase(
                  std::remove_if(e->post.begin(), e->post.end(),[](particle& pt){ return !lepton(pt.pdg);}),
                  e->post.end());

                while (!parts.empty()) parts.pop();

                particle N1 = e->out[1];
                N1.travelled = 0;
                prepare_single_nucleon_for_redraw(N1, 1);

                while (!parts.empty()) {
                  particle p1 = parts.front();
                  parts.pop();
                  p = &p1;
                  leave_nucleus();
                }
              }

       par.kaskada_NN_mfp_scale = original_scale; // fallback

       }

      needs_final_clean = false;

      }

        // Non-Transparent & Correlated: both nucleons can interact,
        // redraw only N1 if N1 fails to interact.
        else
        {
          int nod_n1_failed_firstpass[18];
          for (int nod_i = 0; nod_i < 18; ++nod_i) {
            nod_n1_failed_firstpass[nod_i] = 0;
          }

          // First pass: standard cascade for both nucleons.
          while (!parts.empty() && nucl->Ar() > 0) {

            particle p1 = parts.front();
            parts.pop();
            p = &p1;

            X = prepare_interaction();
            if (!move_particle()) continue;
            if (X.r >= radius) {
              leave_nucleus();
            }
            else {

              int nod_before_failed_step[18];
              for (int nod_i = 0; nod_i < 18; ++nod_i) {
                nod_before_failed_step[nod_i] = e->nod[nod_i];
              }

              if (max_step < X.freepath || !make_interaction() || !finalize_interaction()) {

                if (nucleon(p1.pdg)) e->nod[13]++;
                if (pion(p1.pdg))    e->nod[12]++;
                if (hyperon(p1.pdg)) e->nod[14]++;

                if (p1.nucleon_id == 1) {
                  for (int nod_i = 0; nod_i < 18; ++nod_i) {
                    nod_n1_failed_firstpass[nod_i] +=
                      e->nod[nod_i] - nod_before_failed_step[nod_i];
                  }
                }

                parts.push(p1);

              }
              else {

                if (nucleon(p1.pdg) and (p1.nucleon_id == 1)) {
                  // Record if N1 (id=1) has interacted at least once.
                  intercated = true;

                }
              }
            }
          }

          // If no interaction for N1: then redraw only N1.
          if (!intercated) {

            for (int nod_i = 0; nod_i < 18; ++nod_i) {
              e->nod[nod_i] -= nod_n1_failed_firstpass[nod_i];
            }

            int nod_before_ntc[18];
            for (int nod_i = 0; nod_i < 18; ++nod_i) {
              nod_before_ntc[nod_i] = e->nod[nod_i];
            }

            int redrawAttempts = 0;

            const double N1_px = e->out[1].x;
            const double N1_py = e->out[1].y;
            const double N1_pz = e->out[1].z;
            const double N1_p  = e->out[1].momentum();

            double original_scale = par.kaskada_NN_mfp_scale;
            double effective_scale = original_scale;

            while (!parts.empty()) parts.pop(); // safeguard

            // Redraw attempts only for N1 with progressively shorter mfp.
            while (redrawAttempts < nucl->MAX_EVENT_REDRAWS && !intercated) {

              redrawAttempts++;

              for (int nod_i = 0; nod_i < 18; ++nod_i) {
                e->nod[nod_i] = nod_before_ntc[nod_i];
              }

              // Remove nucleons close to the original N1 : 1 degree acceptance
              e->post.erase(
                std::remove_if(
                  e->post.begin(),
                  e->post.end(),
                  [&](particle &pt)
                  {
                    if (lepton(pt.pdg) || pion(pt.pdg) || hyperon(pt.pdg))
                      return false;

                    if (nucleon(pt.pdg))
                    {
                      const double pt_p = pt.momentum();

                      if (N1_p <= 0.0 || pt_p <= 0.0)
                        return false;

                      double cos_angle =
                        (N1_px * pt.x + N1_py * pt.y + N1_pz * pt.z) / (N1_p * pt_p);

                      if (cos_angle >  1.0) cos_angle =  1.0;
                      if (cos_angle < -1.0) cos_angle = -1.0;

                      return (1.0 - cos_angle < nucl->cosine_threshold);
                    }

                    return false;
                  }
                ),
                e->post.end()
              );

              // Strengthen scattering for the redrawn N1 only.
              effective_scale *= nucl->effective_mfp_scale;
              par.kaskada_NN_mfp_scale = effective_scale;

              // New N1 from the original primary vertex.
              particle N1 = e->out[1];
              N1.travelled = 0;
              prepare_single_nucleon_for_redraw(N1, 1);

              while (!parts.empty() && nucl->Ar() > 0) {

                particle pt = parts.front();
                parts.pop();
                p = &pt;

                X = prepare_interaction();
                if (!move_particle()) continue;
                if (X.r >= radius) {
                  leave_nucleus();
                }
                else {

                  if (max_step < X.freepath || !make_interaction() || !finalize_interaction()) {

                    if (nucleon(pt.pdg)) e->nod[13]++;
                    if (pion(pt.pdg))    e->nod[12]++;
                    if (hyperon(pt.pdg)) e->nod[14]++;

                    parts.push(pt);

                  }
                  else {

                    if (nucleon(pt.pdg) and (pt.nucleon_id == 1)) {
                      intercated = true;
                      par.kaskada_NN_mfp_scale = original_scale;
                    }

                  }
                }
              }

              if (!parts.empty()) // safeguard
              {
                while (!parts.empty()) parts.pop();
              }
            }

            if (!intercated) { // fallback

              for (int nod_i = 0; nod_i < 18; ++nod_i) {
                e->nod[nod_i] = nod_before_ntc[nod_i];
              }

              e->post.erase(
                std::remove_if(
                  e->post.begin(),
                  e->post.end(),
                  [&](particle &pt)
                  {
                    if (lepton(pt.pdg) || pion(pt.pdg) || hyperon(pt.pdg))
                      return false;

                    if (nucleon(pt.pdg))
                    {
                      const double pt_p = pt.momentum();

                      if (N1_p <= 0.0 || pt_p <= 0.0)
                        return false;

                      double cos_angle =
                        (N1_px * pt.x + N1_py * pt.y + N1_pz * pt.z) / (N1_p * pt_p);

                      if (cos_angle >  1.0) cos_angle =  1.0;
                      if (cos_angle < -1.0) cos_angle = -1.0;

                      return (1.0 - cos_angle < nucl->cosine_threshold);
                    }

                    return false;
                  }
                ),
                e->post.end()
              );

              while (!parts.empty()) parts.pop();

              particle N1 = e->out[1];
              N1.travelled = 0;
              prepare_single_nucleon_for_redraw(N1, 1);

              while (!parts.empty()) {
                particle p1 = parts.front();
                parts.pop();
                p = &p1;
                leave_nucleus();
              }
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

void kaskada::set_shell_sampler(shell_sampler *sampler)
{
  shell_dist = sampler;
}

// Build initial queue "parts" from event primary vertex e->out
void kaskada::prepare_particles()
{
  // Loop over every particle from the primary vertex.
  for (int i = 0; i < e->out.size(); i++)
  {
    particle p1 = e->out[i];

    if (nucleon_or_pion (p1.pdg)) // formation zone for both nucleons and pions
    {
      if (nucleon (p1.pdg))
      {
        p1.primary = true;

        double kaskada_w = par.kaskada_w;

        // Tag primary and spectator nucleons form SF-event vertex
        if(e->flag.qel and par.sf_method != 0) {
         if (i == 1) {
            p1.nucleon_id = 1;
         }
         else if (i == 2) {
            p1.nucleon_id = 2;
         }
      }

        // add in-medium energy for GFG LFG and SF
        if (e->flag.qel and (par.sf_method != 0 or par.nucleus_target == 2)) {
            p1.set_energy (p1.E() + nucl->Ef(p1) + kaskada_w);
        }

        else if (par.nucleus_target == 1 and (e->flag.qel or e->flag.res))
          p1.set_energy (p1.E() + par.nucleus_E_b);

        p1.set_fermi(nucl->Ef(p1));


        // If kinetic energy is below "barrier" = Ef + kaskada_w, jail back to nucleus
        if (p1.Ek() <= kaskada_w + p1.his_fermi)
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

  e->r_distance = 10;         // new JS ; default (large) value, if unchanged no absorption took place

}

// Prepare a single nucleon (typically N1) to be re-run in a redraw scenario:
void kaskada::prepare_single_nucleon_for_redraw(particle N1, int index)
{

  if (!nucleon(N1.pdg)) return;

    particle pN = N1;
    pN.primary = true;

        if (index == 1) {
           pN.nucleon_id = 1;
        }   
        else if (index == 2) {
           pN.nucleon_id = 2;
        }   
        else {
           pN.nucleon_id = 0;
        }   

    double kaskada_w = par.kaskada_w;

    pN.set_energy (pN.E() + nucl->Ef(pN) + kaskada_w);

    pN.set_fermi(nucl->Ef(pN));

    if (pN.Ek() <= pN.his_fermi + kaskada_w)
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

// Prepare local interaction parameters at current particle position:
interaction_parameters kaskada::prepare_interaction()
{
  interaction_parameters res;

  res.pdg = p->pdg;
  res.Ek  = p->Ek();
  res.r   = p->r.length ();

  res.use_consistent_target = e->flag.qel && par.sf_method != 0 && p->nucleon();
  res.target_selected = false;
  res.has_neutron_target = false;
  res.has_proton_target = false;

  res.Ekeff_n_forward = 0.0;
  res.Ekeff_n_reverse = 0.0;
  res.Ekeff_p_forward = 0.0;
  res.Ekeff_p_reverse = 0.0;

  res.xsec_n_forward = 0.0;
  res.xsec_n_reverse = 0.0;
  res.xsec_p_forward = 0.0;
  res.xsec_p_reverse = 0.0;
  res.xsec_n_pair = 0.0;
  res.xsec_p_pair = 0.0;

  res.dens = nucl->density (res.r);
  if(shell_dist && nucl->A()>1) // correct for the primary vertex dist.
  {
    res.dens = std::max(0.0, res.dens * nucl->A() - shell_dist->dens(res.r) * nucl->Ar())
             / (nucl->A() - 1);
  }
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

  // no correlations for incoming nucleons in the scattering mode (kaskada.cc)
  bool use_corr = !e->out[0].nucleon() || e->number_of_interactions() || par.beam_placement != 2;
  double s = 1.0 / par.kaskada_NN_mfp_scale;

  switch (res.pdg) // mean free path modifications: effective density and scaling
  {
    case pdg_neutron:
      res.xsec_n *= (use_corr ? corr_ii / norm_ii : 1.0) * s;
      res.xsec_p *= (use_corr ? corr_ij / norm_ij : 1.0) * s;
      break;

    case pdg_proton:
      res.xsec_n *= (use_corr ? corr_ij / norm_ij : 1.0) * s;
      res.xsec_p *= (use_corr ? corr_ii / norm_ii : 1.0) * s;
      break;

    default:
      break;
  }

  res.xsec = res.dens_n*res.xsec_n + res.dens_p*res.xsec_p; // calculate the inverse of the mean free path

  assert(res.xsec>=0);                      // make sure that the cross section is positive

  if (res.xsec != 0)
  {
    res.freepath = -log (frandom ()) / res.xsec; // choose free path according to the mean free path (1/res.xsec)
    res.prob_proton = res.xsec_p * res.dens_p / res.xsec;
  }
  else
  {
    res.freepath = 2.0 * max_step;
    res.prob_proton = 0.0;
  }

  return res;
}

// Move particle along free path or max_step, and possibly jail nucleons
bool kaskada::move_particle()
{
  p->krok (min (max_step, X.freepath));   // propagate by no more than max_step

  if(e->flag.qel && par.sf_method != 0 && p->nucleon())
  {
    const double old_fermi = p->his_fermi;
    double local_fermi = nucl->Ef(*p);
    if(local_fermi < 0.0 && local_fermi > -1.0e-10) local_fermi = 0.0;
    double new_Ek = p->Ek() + local_fermi - old_fermi;
    if(new_Ek < 0.0 && new_Ek > -1.0e-10) new_Ek = 0.0;

    assert(std::isfinite(local_fermi));
    assert(std::isfinite(new_Ek));
    assert(new_Ek >= 0.0);

    p->set_energy(p->mass() + new_Ek);
    p->set_fermi(local_fermi);
  }

  if (!(p->nucleon() || hyperon(p->pdg))) // pion can not be jailed
    return true;

  if(p->nucleon())
  {

    double kaskada_w = par.kaskada_w;

    // jail nucleon if its kinetic energy is lower than "binding" energy
    if (p->Ek() <= kaskada_w + p->his_fermi)
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
  // Nucleons: go through final "climb-out of potential well" and possible jailing
  if (nucleon (p->pdg))
  {
    double kaskada_w = par.kaskada_w;

    // If KE below (Ef + W) at surface, jail nucleon
    if (p->Ek() <= p->his_fermi + kaskada_w)
    {
      p->endproc=jailed;
      nucl->insert_nucleon (*p);
      if(par.kaskada_writeall) e->all.push_back(*p);
      return false; // particle did not escape
    }

      p->set_energy(p->E() - p->his_fermi - kaskada_w);

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
  if(X.use_consistent_target)
  {
    X.pdg = p->pdg;
    X.Ek = p->Ek();
    X.r = p->r.length();

    X.dens = nucl->density(X.r);
    if(shell_dist && nucl->A()>1)
    {
      X.dens = std::max(0.0, X.dens * nucl->A() - shell_dist->dens(X.r) * nucl->Ar())/(nucl->A() - 1);
    }
    assert(X.dens>=0);

    X.dens_n = X.dens * nucl->frac_neutron();
    X.dens_p = X.dens * nucl->frac_proton();
    X.n = 2;

    I->total_cross_sections(*p, *nucl, X);

    corr_func->set_input_point(p->travelled);
    double corr_ii = corr_func->get_value(1);
    double corr_ij = corr_func->get_value(2);
    corr_func->set_input_point(p->r.length());
    double norm_ii = corr_func->get_value(3);
    double norm_ij = corr_func->get_value(4);

    if(!e->out[0].nucleon() || e->number_of_interactions() || par.beam_placement != 2)
    {
      switch(X.pdg)
      {
        case pdg_neutron:
          X.xsec_n *= corr_ii / norm_ii / par.kaskada_NN_mfp_scale;
          X.xsec_p *= corr_ij / norm_ij / par.kaskada_NN_mfp_scale;
          break;

        case pdg_proton:
          X.xsec_n *= corr_ij / norm_ij / par.kaskada_NN_mfp_scale;
          X.xsec_p *= corr_ii / norm_ii / par.kaskada_NN_mfp_scale;
          break;

        default:
          break;
      }
    }

    const double neutron_weight = X.dens_n * X.xsec_n;
    const double proton_weight = X.dens_p * X.xsec_p;
    const double total_weight = neutron_weight + proton_weight;

    X.xsec = total_weight;
    X.prob_proton = total_weight > 0.0 ? proton_weight / total_weight : 0.0;

    if(!(total_weight > 0.0) || !(X.has_neutron_target || X.has_proton_target)) return false;

    bool choose_proton = false;

    if(!X.has_neutron_target) choose_proton = true;
    else if(!X.has_proton_target) choose_proton = false;
    else choose_proton = frandom() < X.prob_proton;

    particle selected_target;
    double xsec_forward = 0.0;
    double xsec_reverse = 0.0;
    double Ekeff_forward = 0.0;
    double Ekeff_reverse = 0.0;

    if(choose_proton)
    {
      selected_target = X.p2_proton;
      xsec_forward = X.xsec_p_forward;
      xsec_reverse = X.xsec_p_reverse;
      Ekeff_forward = X.Ekeff_p_forward;
      Ekeff_reverse = X.Ekeff_p_reverse;
    }
    else
    {
      selected_target = X.p2_neutron;
      xsec_forward = X.xsec_n_forward;
      xsec_reverse = X.xsec_n_reverse;
      Ekeff_forward = X.Ekeff_n_forward;
      Ekeff_reverse = X.Ekeff_n_reverse;
    }

    const double orientation_weight = xsec_forward + xsec_reverse;
    if(!(orientation_weight > 0.0)) return false;

    X.p2 = selected_target;

    if(frandom() * orientation_weight < xsec_reverse)
    {
      X.p2.x *= -1.0;
      X.p2.y *= -1.0;
      X.p2.z *= -1.0;
      X.Ekeff = Ekeff_reverse;
    }
    else X.Ekeff = Ekeff_forward;

    X.p2.r = p->r;
    X.target_selected = true;
  }

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
    const bool preserve_fermi_by_branch = e->flag.qel && par.sf_method != 0 &&
                                          I->process_id() == nucleon_ + elastic_ &&
                                          X.n == 2 && nucleon(X.p[0].pdg) && nucleon(X.p[1].pdg);

    if(preserve_fermi_by_branch)
    {
      X.p[0].set_fermi(p->his_fermi);
      X.p[1].set_fermi(FE);
    }
    else
    {
      double he = (p->his_fermi > FE) ? p->his_fermi : FE;
      double le = p->his_fermi + FE - he;

      int n_he = -1;
      int n_le = -1;

      for(int i=0;i<X.n;i++)
        if(nucleon(X.p[i].pdg))
        {
          if(n_he < 0)
            n_he = i;
          else if(X.p[i].Ek() < X.p[n_he].Ek())
            n_le = i;
          else
          {
            n_le = n_he;
            n_he = i;
          }
        }

      X.p[n_he].set_fermi(he);
      X.p[n_le].set_fermi(le);
    }
  }

  for (int i = 0; i < X.n; i++)
  {
    X.p[i].r = p->r;
    X.p[i].travelled = 0;

    double kaskada_w = par.kaskada_w;

    // jail nucleon if its kinetic energy is lower than work function
    if (nucleon (X.p[i].pdg) and X.p[i].Ek() <= kaskada_w + X.p[i].his_fermi)
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

// Clean remaining 'parts' queue at the end of event
void kaskada::clean()
{

  while(!parts.empty())
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
