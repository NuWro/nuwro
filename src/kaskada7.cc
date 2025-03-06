#include "kaskada7.h"
#include "fsi.h"

////////////////////////////////////////
// Public methods
////////////////////////////////////////

kaskada::kaskada(params &p, event &e1, input_data *input)
{
  par = p;

  if (par.nucleus_p + par.nucleus_n < 3)
    par.kaskada_w = 0;  //for free nucleons and deuteron there is no extra binding energy

  e = &e1;
  max_step = par.step * fermi;    // set maximum step defined in params
  nucl = make_nucleus(par);       // create nucleus defined in params
  radius = nucl->radius();        // calculate radius of the nucleus
  I = new Interaction(input->get_data_container(0), input->get_data_container(1), input->get_data_container(2),
                      par.kaskada_piN_xsec);
  corr_func = input->get_nucl_data_container(0, par.nucleus_p, par.nucleus_n);
}

////////////////////////////////////////

kaskada::~kaskada()
{
  delete nucl;
  delete I;
}

////////////////////////////////////////

int kaskada::kaskadaevent()
{
  int result = 0;

  if (e->weight <= 0)
    return result;

  prepare_particles();    // copy all nucleons and pions from primary vertex (out) to queue
                          // and other particles to output (post)

  if (e->in[0].lepton())  // if lepton scattering -> remove nucleon from primary vertex 
                          // (in pion or nucleon scattering there is no primary vertex just the cascade)
  {
    for(int i=1;i< e->in.size();i++)
      nucl->remove_nucleon(e->in[i]); // ignores nonnucleons
  }

  if (not par.kaskada_on) // skip the cascade if it is turn off in params,
                          // but make sure about the energy balance
  {
    while (parts.size () > 0)
    {
      particle p1 = parts.front();    // point a particle from a queue
      parts.pop();                    // remove this particle from a temp vector
      p = &p1;
      
      leave_nucleus();    // check if the particle is jailed or escapes (and returns on the mass shell)
    }

    return result;
  }

  while (parts.size () > 0 and nucl->Ar() > 0)  // main loop in cascade
  {
    particle p1 = parts.front();                // point a particle from a queue
    parts.pop();                                // remove this particle from a temp vector
    p = &p1;

    X = prepare_interaction();                  // set the density and the total cross section
                                                // calculate free path

    if (!move_particle()) continue;             // propagate particle, returns false if jailed

    if (X.r >= radius)                          // particle leaves nucleus
      leave_nucleus();
    else if( max_step < X.freepath              // no interaction during max_step 
         or !make_interaction()                 // or unable to generate kinematics or Pauli blocking
         or !finalize_interaction())            // or there was problem with finalizing interaction
    {
      if (nucleon (p1.pdg)) e->nod[13]++;
      if (pion (p1.pdg)) e->nod[12]++;
      if (hyperon (p1.pdg)) e->nod[14]++;
      parts.push (*p);                          // interaction did not happend, 
                                                // p should be further propagated 
    }
  }

  clean();  // if nucleus has evaporated the part queue may not be empty

  return result;
}

////////////////////////////////////////
// Private methods
////////////////////////////////////////

void kaskada::prepare_particles()
{ 
  for (int i = 0; i < e->out.size(); i++)
  {
    particle p1 = e->out[i];

    if (nucleon_or_pion (p1.pdg)) // formation zone for both nucleons and pions
    {           
      if (nucleon (p1.pdg))
      { 
        p1.primary = true;

        // add the binding energy substracted in the primary vertex
        // (for Global Fermi Gas Local Fermi Gas and Spectral Function)
        if (e->flag.qel and (par.sf_method != 0 or par.nucleus_target == 2))
          p1.set_energy (p1.E() + nucl->Ef(p1) + par.kaskada_w);

        else if (par.nucleus_target == 1 and (e->flag.qel or e->flag.res))
          p1.set_energy (p1.E() + par.nucleus_E_b);

        p1.set_fermi(nucl->Ef(p1));
 
        if (p1.Ek() <= par.kaskada_w + p1.his_fermi)  // jailed nucleon if its kinetic energy
                                                      // is lower than binding energy
        {
          p1.endproc=jailed;
          nucl->insert_nucleon (p1);
          if(par.kaskada_writeall) 
            e->all.push_back(p1);
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
        ////////// a new piece of code March 3, 2025 - fixing hyperon potential problem
        if (p1.E() +p1.his_fermi< p1.mass())
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
        ////////// the end of the new piece of code
        ///
        /*// before March 3, 2025 the below two lines where active
      p1.set_energy(p1.E() + p1.his_fermi);

      parts.push(p1); // add particle to queue
        */
    }
    else              // if not a nucleon, pion nor hyperon
    {
      p1.endproc=escape;
      e->post.push_back (p1);
      
      if(par.kaskada_writeall) 
        e->all.push_back(p1);
    }
  }
  
  for (int i = 0; i<18; i++)  // number of dynamics defined in proctable.h
     e->nod[i] = 0;

  e->r_distance = 10;         // new JS ; default (large) value, if unchanged no absorption took place
}

////////////////////////////////////////

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

  I->total_cross_sections (*p, *nucl, res); // calculate cross sections xsec_p and xsec_n

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

////////////////////////////////////////

bool kaskada::move_particle()
{
  p->krok (min (max_step, X.freepath));   // propagate by no more than max_step

  if (!(p->nucleon() || hyperon(p->pdg))) // pion can not be jailed
    return true;

  if(p->nucleon())
  {
    // jail nucleon if its kinetic energy is lower than binding energy
    if (p->Ek() <= par.kaskada_w + p->his_fermi)
    {
      p->endproc=jailed;
      nucl->insert_nucleon (*p);
      if(par.kaskada_writeall) 
        e->all.push_back(*p);
      return false;                       // nucleon was jailed  
    }
    else
      return true;                        // nucleon was not jailed
  }

  if(hyperon(p->pdg))
  {
    // If hyperon KE falls below its potential energy -  hyperon is jailed
    if(p->Ek() < p->his_fermi)
    {
      p->endproc=jailed;
      if(par.kaskada_writeall) 
	      e->all.push_back(*p);
      return false;                       // hyperon was jailed  
    }
    return true; // hyperon was not jailed
  }

  return true;
}

////////////////////////////////////////

bool kaskada::leave_nucleus()
{
  if (nucleon (p->pdg))                 // substract fermi energy and work function
  {
    // jail nucleon if its kinetic energy is lower than binding energy
    if (p->Ek() <= p->his_fermi + par.kaskada_w)
    {
      p->endproc=jailed;
      nucl->insert_nucleon (*p);
      if(par.kaskada_writeall) 
        e->all.push_back(*p);
      return false;                     // particle did not escape
    }
    else
      p->set_energy(p->E() - p->his_fermi - par.kaskada_w);
  }
  else if(hyperon(p->pdg))
  {
    // Adjust hyperon for remaining binding energy, jail if it cannot escape
    if (p->Ek() <  p->his_fermi)
    {
      p->endproc=jailed;
      //TODO: Perhaps add hyperon to nucleus here (or decay it)
      if(par.kaskada_writeall) e->all.push_back(*p);
      return false;
    }
    // Subtract binding energy from hyperon energy and set momentum so it is on shell
    else {
      p->set_energy(p->E() -  p->his_fermi);
    }
  }

  p->endproc=escape;
  e->post.push_back (*p);
  
  if(par.kaskada_writeall)
    e->all.push_back (*p);
    
  return true;                          // particle has escaped
}

////////////////////////////////////////

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
    cout<<I->process_id()<<'/'<<loop<<','<<endl;
    cout<<int(rep*1000./call)/10.<<" % repeated"<<endl;
    double suma=0;
    cout<<"["<<nucl->p<<","<<nucl->n<<"/"<<nucl->Zr()<<","<<nucl->Nr()<<"]"<<p->pdg<<' '<<X.p2<<':'<<endl;
    for(int i=0;i<X.n;i++)
    {
      suma+=X.p[i].mass();
      cout<<' '<<X.p[i].pdg;
    }
    cout<<endl<<p->mass()+X.p2.mass()<<' '<<suma<<endl;
//CJ  cout<<" Interaction: "<<I.process_name()<<" ("<<I.process_id()<<") ";
//    assert(procid==I->process_id());
/* KN: Is it important to keep the lines above? */
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

  //TODO check this
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

////////////////////////////////////////

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

    // jail nucleon if its kinetic energy is lower than work function
    if (nucleon (X.p[i].pdg) and X.p[i].Ek() <= par.kaskada_w + X.p[i].his_fermi)
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

    if(par.kaskada_writeall)
      e->all.push_back (*p);

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

////////////////////////////////////////

void kaskada::clean()
{
  while(! parts.empty())
  {
    particle p0 = parts.front();
    p0.endproc = escape;
    e->post.push_back(p0);

    if(par.kaskada_writeall)
      e->all.push_back(p0);
    
    parts.pop();
  }

  e->pr=nucl->Zr();
  e->nr=nucl->Nr();
}

////////////////////////////////////////

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

////////////////////////////////////////

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
