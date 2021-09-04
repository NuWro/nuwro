#ifndef _beam_atmo_h_
#define _beam_atmo_h_
#include "EnergyProfile.h"
#include "params.h"
#include "beam.h"
#include <cmath>

struct bin{
	int pdg;
	double cosZ;
	double E;
	double w[2];
};

class beam_atmo: public beam
{
	
	private:  
		vector<bin>  bins;
		double f;
		double dcos;
	vec dir;

		bin& find_bin(bool dis){
			int a=0,b=bins.size()-1;
			double tot=bins[b].w[dis];
			double x=tot*frandom();
			while(a<b)
			{
				int s=(a+b)/2;
	//      cout<<a<<" "<<s<< " "<<b <<endl;  
				if(x>=bins[s].w[dis])
					a=s+1;
				else
					b=s;
			}
		//  cout<<a<< " ";
			return bins[a];
		}

	public: 
		beam_atmo(params &p):f(1),dir(p.beam_direction.dir()){
			dcos=0.1;
			stringstream in(p.beam_atmo_files);
			int pdg=0;
			string fname;
			while(in >> pdg >> fname)
			{
				if(abs(pdg)==12 || abs(pdg)== 14 || abs(pdg)==16)
				{

				} 
				else
				{
					cerr<< "bad neutrino pdg code "<< pdg <<" in parameter 'atmo_files"<<endl;
					exit( -13);
				}
					
				ifstream atmo;

				if(! open_data_file(atmo,fname))
				{
					cerr<< "cant't read neutrino atmo_file '"<<fname<<"'."<<endl;
					exit( -13);
				}
				bin b;
				b.pdg=pdg;
				string line;
				double s0=0,s1=0,w,E;    
				while(getline(atmo,line))
				if(line[0]!='#')
				{
					stringstream as(line);
					if(as >> E >> b.cosZ >> w)
					{   b.w[0]=s0+=w;
						b.w[1]=s1+=w*E;
						b.E=E*GeV;
						bins.push_back(b);
					}
				}
				atmo.close();
			}

			f = bins.size()>1 ? sqrt(bins[1].E/bins[0].E) : 1;

			cout<< "Atmo beam created with "<<bins.size()<<" total bins"<<endl;
		}    

		particle shoot(bool dis)
		{
			bin  b=find_bin(dis);
			double Cos=b.cosZ-dcos/2+frandom()*dcos;
			double Sin=sqrt(1-Cos*Cos);
			double phi=2*M_PI*frandom();

			double Elo= b.E/f;
			double Ehi= b.E*f;
			double E=0;	
			double z = frandom();
			if(dis)
				E=Elo + z*(Ehi-Elo);
			else 
				E=exp(log(Elo)*(1-z) + z*log(Ehi));

			particle nu(b.pdg,0);
			vec p(E* Sin*cos(phi), E*Sin*sin(phi), E*Cos);
			nu.set_momentum(p.fromZto(dir));
			return nu;
		} 


};


#endif