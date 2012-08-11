#include "kaskada7.h"

kaskada::kaskada(params &p, event &e1)
{
	par = p;
	e = &e1;
	nucl = make_nucleus(par); //create nucleus defined in params
	I = new Interaction(par.xsec);
}

int kaskada::kaskadaevent()
{	
	int result = 0;
	
	if(not par.kaskada_on) //skip cascade if it is turn off in params
		return result;
		
	double max_step = par.step * fermi; //set maximum step defined in params
	
	double radius = nucl->radius(); //calculate radius of the nucleus

	if (e->in[0].lepton()) //remove nucleon from primary vertex 
		nucl->remove_nucleon(e->in[1]);
		
	prepare_particles(); //copy all nucleons and pions from primary vertex (out) to queue and other particles to output (post)
	
	reset_stats();
	
	while (parts.size () > 0 and nucl->Ar() > 0) ///main loop in cascade
    {						
		particle p1 = parts.front(); //point a particle from a queue
		parts.pop(); 				 //remove this particle from a temp vector
		p = &p1;		
		
		X = prepare_interaction();
										
		p->krok (min (max_step, X.freepath)); //propagate by no more than max_step
		
		if (X.xsec == 0 or X.r >= radius) 								   //leaving nucleus
			leave_nucleus();
		else if ((max_step < X.freepath) or !make_interaction()) 		   //no interaction during max_step or unable to generate kinematics or Pauli blocking
			parts.push (*p); 			 								   //put to end of queue for further processing
		else 															   //scattering happened
		{
			result = 1;
			run_test();
			finalize_interaction();
		}		
	}
	
	clean(); // if nucleus has evaporated the part quequ may not be empty   

	return result;
}

void kaskada::prepare_particles()
{
	for (int i = 0; i < e->out.size(); i++)
	{
		particle p1 = e->out[i];      
		
		if (nucleon_or_pion (p1.pdg))
		{
			p1.set_new(); //reset history
		  
			double fz = formation_zone(p1, par, *e); //calculate formation zone
			p1.krok(fz); //move particle by a distance defined by its formation zone
			
			if (nucleon (p1.pdg)) //set Fermi energy from primary vertex
				p1.set_fermi(e->in[1].Ek());
		  		  
			parts.push (p1); //put particle to a queue
		}
		else
		{
			p1.endproc=escape;
			e->post.push_back (p1);
			
			if(par.kaskada_writeall) 
				e->all.push_back(p1);
		}
    }
    
   	for (int i = 0; i<12; i++) //number of dynamics defined in proctable.h
		e->nod[i] = 0;
}

interaction_parameters kaskada::prepare_interaction()
{
	interaction_parameters res;
	
	res.pdg = p->pdg;
	res.Ek = p->Ek();  
	res.r = p->r.length ();
    res.dens = nucl->density (res.r); assert(res.dens>=0);
	res.dens_n = res.dens * nucl->frac_neutron ();
	res.dens_p = res.dens * nucl->frac_proton ();
	res.n = 2;

    I->total_cross_sections (res); //calculate cross sections xsec_p and xsec_n

    res.xsec = res.dens_n*res.xsec_n + res.dens_p*res.xsec_p; assert(res.xsec>=0);    
        
    if (res.xsec != 0)
    {
		res.freepath = -log (frandom ()) / res.xsec;
		res.frac_proton = res.xsec_p * res.dens_p / res.xsec;
	}
	               
    return res;
}

bool kaskada::make_interaction()
{
	int loop = 0;
	static int call=0;
	static int rep=0;
	static int procid=0;
			
	while(++call && 0 == I->particle_scattering (*p, *nucl, X))
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
//CJ	cout<<" Interaction: "<<I.process_name()<<" ("<<I.process_id()<<") ";
		assert(procid==I->process_id());
		if(loop>100) return false; //it was impossible to make kinematics
	}
	
	for (int i = 0; i < X.n; i++)
		if (!X.p[i].is_valid ())
		{
			cerr << I->process_name()<< "Interaction: error" << X.p[i] << endl;
			delete nucl;
			delete I;
            exit(0);
		    return false;
		}

	if ((p->pdg == 2112 or p->pdg == 2212) and nucl->pauli_blocking(X.p, X.n)) //check if there was Pauli blocking
		return false;
	
	assert(check(*p,X.p2,nucl->spectator,X.n,X.p,I->process_id()));
	assert(check2(*p,X.p2,nucl->spectator,X.n,X.p,I->process_id()));

	return true;
}

void kaskada::finalize_interaction()
{
	p->endproc=I->process_id();
	
	nucl->remove_nucleon (X.p2); // remove from the nuclear matter
	if(nucl->spectator!=NULL)
		nucl->remove_nucleon (*nucl->spectator);
		
	for (int i = 0; i < X.n; i++)
	{
		X.p[i].r = p->r;
		X.p[i].travelled = 0;
		X.p[i].set_mother(*p);
		
		double fz = formation_zone(X.p[i], par);
		X.p[i].krok(fz);
		
		if (nucleon (X.p[i].pdg))
		{
			if (nucleon (p->pdg))
				X.p[i].set_fermi((X.p2.Ek()+p->his_fermi)/2.0);
			else
				X.p[i].set_fermi(X.p2.Ek());
		}

		parts.push (X.p[i]);	// queue for further processing
		//procinfo(*p,X.p2,X.n,X.p);
		
		if(par.kaskada_writeall)
			e->all.push_back (*p);
	}
	
	int k = kod(I->process_id());
	e->nod[k]++;
}

void kaskada::leave_nucleus()
{	
	//if (nucleon (p->pdg))	// reduce momentum when leaving nucleus or jailed the nucleon
	//{
	//	double V = p->his_fermi + par.kaskada_w;
	//	p->set_momentum (p->p () /p->momentum() * sqrt(p->momentum()*p->momentum() - 2.0*p->E()*V  + V*V));	
	//}
	
	p->endproc=escape;
	e->post.push_back (*p);
	
	if(par.kaskada_writeall)
		e->all.push_back (*p);
		
}

bool kaskada::check (particle & p1, particle & p2, particle *spect, int n, particle p[],int k)
{
	int ch1=p1.charge()+p2.charge();
	if(spect) ch1+=spect->charge();
	int i=n;
	while(i)
	  ch1-=p[--i].charge();
	if(ch1!=0)
	{   cout<<endl<<"Proc:"<<k<<endl;
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

bool kaskada::check2 (particle & p1, particle & p2, particle *spect, int n, particle p[],int k)
{
	vect p4=p1+p2;
	if(spect) p4+=*spect;
	int i=n;
	while(i)
	  p4-=p[--i];
	double prec=0.001*MeV;  
	if(abs(p4.t)>prec ||abs(p4.x)>prec ||abs(p4.y)>prec ||abs(p4.z)>prec)
	{   cout<<endl<<"Proc:"<<k<<" delta p4 = "<<p4<<endl;
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

void kaskada::reset_stats()
{		
	for (int i = 0; i<20; i++) //place of secondary interaction
		for (int j = 0; j < 11; j++)
			e->place[j][i] = 0;
	
	e->nopp = true;
	e->abs = false;
	e->pabsen = 0;
	e->nofi = 0;
	e->absr = 100;
	e->odl = 0;
	
	for (int i = 0; i < 10; i++)
		e->pen[i] = 0;
}

void kaskada::run_test()
{
	int ktory=2*X.r/fermi;
	ktory=min(ktory,19);
	
	if (I->process_id() >= 20 and I->process_id() <= 25)
	{
		vec help = vecprod(p->r, e->in[0].r);
		e->odl = help.length();
	}
	
	switch(I->process_id())
	{
		case 12: case 13: case 22: case 23: case 24: e->nopp = false; break;
		case 25:
			e->pabsen = p->Ek();
			e->absr = p->r.length();
			e->abs = true;
			break;
	}

			
	if (I->process_id() >= 20 and I->process_id() < 25 and !e->abs)
	{
		e->nofi++;
		if (e->nofi < 10) e->pen[e->nofi] = p->Ek();
	}
	
	int k =kod(I->process_id());
		
	if(k<11)
		e->place[k][ktory]++;
		
	static int stat[11]={0,0,0,0,0,0,0,0,0,0,0};
		stat[k]++;
}

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

	delete nucl;
	delete I;
}
