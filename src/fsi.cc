#include "fsi.h"

#include <string>
#include <iostream>
#include <cstdlib>
#include <TFile.h>
#include <TTree.h>

bool expr_on[nof_expr];
bool fz_on[nof_fz];
bool options_on[nof_opt];
bool xsec_on[nof_xsec];
bool allinone_on = false;
string allinone = "Formation length for nucleons and pions";

char day[10];
char month[10];
char fullm[10];
char hour[10];
char minute[10];
char year[10];
string date;

ofstream logfile;

#ifdef WIN32
#include <conio.h>
#else
#include <termios.h>
int getch (void)
{
        int key;
        struct termios oldSettings, newSettings;
        tcgetattr(STDIN_FILENO, &oldSettings);
        newSettings = oldSettings;
        newSettings.c_lflag &= ~(ICANON | ECHO);
        tcsetattr(STDIN_FILENO, TCSANOW, &newSettings);
        key = getchar();
        tcsetattr(STDIN_FILENO, TCSANOW, &oldSettings);
        return key;
}
#endif

void get_date()
{
	time_t rawtime;
	struct tm * timeinfo;
	
	time ( &rawtime );
	timeinfo = localtime ( &rawtime );
	 
	strftime(day, sizeof(day), "%d", timeinfo);
	strftime(month, sizeof(month), "%m", timeinfo);
	strftime(fullm, sizeof(fullm), "%B", timeinfo);
	strftime(hour, sizeof(hour), "%H", timeinfo);
	strftime(year, sizeof(year), "%Y", timeinfo);
	strftime(minute, sizeof(minute), "%M", timeinfo);
	
	date = string(year) + dot + string(month) + dot + string(day);
}

void run(string com)
{
	FILE *fp;
	fp = popen(com.c_str(), "w");
	pclose(fp);
}

bool noFile(string filename)
{
	fstream plik;
	plik.open(filename.c_str());
	
	if (plik.is_open())
	{
		plik.close();
		return false;
	}
	
	plik.close();
	return true;	
}

string find_last (string name)
{
	string command;
	command = string("ls ") + name + string(" > tmp/help.txt");
	run(command);
	
	ifstream file ("tmp/help.txt");
	string filename = "empty";	
		
	if (file)
	{
		do
		{
			file >> filename;
		}while (file);
	}
	
	file.close();
	return filename;
}

double ran_exp(double p)
{
	return -p*log(1.0 - frandom00());
}

double formation_zone (particle &p, params &par, event &e)
{
	double flength;
			
	if (par.first_step)
	{	
		string fz;
		
		vect p0 = e.in[1].p4();
		vect q = e.q();
	
		for (int i=0; i<par.formation_zone.size();i++)
		{
			if(par.formation_zone[i]!=' ')
			fz += par.formation_zone[i];
		}

		bool qel = e.flag.qel;
		bool mec = e.flag.mec;

		if (strcmp(fz.c_str(), "fz") == 0)
		{
			double W = e.W();
			
			if (qel) flength = p.momentum()/(p.mass2() - p.p4()*p0);
			else if (mec)
			{
				//vec mom = e.in[1].p() + e.in[2].p();
				//vect mom4 = e.in[1].p4() + e.in[2].p4();
				
				//flength = mom.length() / (mom4 * q);
				flength = 0;
			}
			else if (W < 1400)
			{
				double e = q.t + p0.t;
				vec pa;
				
				pa.x = q.x + p0.x;
				pa.y = q.y + p0.y;
				pa.z = q.z + p0.z;
				
				double mass = sqrt(e*e - pa.length()*pa.length());
				
				flength = pa.length()*ran_exp(1.0/120.0)/mass;
			}
			else if (W > 1800)
			{
				vec help;
				help.x = q.x;
				help.y = q.y;
				help.z = q.z;
								
				double pt = help.dir()*p.p();
				pt = sqrt(p.momentum2() - pt*pt);
				
				double tau = ran_exp(par.tau)/200.0;
					
				flength = tau*p.momentum()*p.mass()/(p.mass2() + pt*pt);
			}
			else
			{
				double e = q.t + p0.t;
				vec pa;
				
				pa.x = q.x + p0.x;
				pa.y = q.y + p0.y;
				pa.z = q.z + p0.z;
				
				double mass = sqrt(e*e - pa.length()*pa.length());
				
				double f1 = pa.length()*ran_exp(1.0/120.0)/mass;
				
				vec help;
				help.x = q.x;
				help.y = q.y;
				help.z = q.z;
								
				double pt = help.dir()*p.p();
				pt = sqrt(p.momentum2() - pt*pt);
												
				double tau = ran_exp(par.tau)/200.0;
					
				double f2 = tau*p.momentum()*p.mass()/(p.mass2() + pt*pt);
				
				
				flength = (W - 1400.0)*f2 + (1800.0 - W)*f1;
				flength /= 400.0;
			}
		}
		else if (strcmp(fz.c_str(), "nofz") == 0) flength = 0;
		else if (strcmp(fz.c_str(), "trans") == 0)
		{
			//double theta[5] = {10.58, 13.44, 12.74, 11.53, 9.09};
			
			//for (int i = 0; i < 5; i++) theta[i] *= 3.14159265/180.0;
			
			if (p.pion())
			{
				const double nu[] = {2.831, 3.282, 3.582, 4.344, 4.733};
				const double Q2p[] = {1.1, 2.15, 3.0, 3.91, 4.69};
				const double PTppi[]   = {2.793, 3.187, 3.418, 4.077, 4.412};
				
				int i;
				
				for (i = 0; i < 5; i++)
				{
					double x = p.momentum() - PTppi[i]*1000;
					if (x < 0)
						x = -x;
					if (x < 10)
						break;
				}
				double M = PDG::mass_proton / 1000.0;
				//flength = sqrt(nu[i]*nu[i] + Q2p[i])/(nu[i] + PDG::mass_proton/1000.0) * ran_exp(1.0/120.0);
				//flength = sqrt(nu[i]*nu[i] + Q2p[i])/sqrt(M*M + 2.0*M*nu[i] - Q2p[i]) * ran_exp(1.0/120.0);
				//flength = sqrt(nu[i]*nu[i] + Q2p[i])/(M*M + 2.0*M*nu[i] - Q2p[i]) * 1e-3;
				flength = PTppi[i] / Q2p[i] * 1e-3;				
				
				//~ int k = 0;
				//~ 
				//~ for (int i = 0; i < 5; i++)
				//~ {
					//~ int energy;
			//~ 
					//~ energy = PTppi[i]*1000.0;
					//~ energy = energy*energy + 140.0*140.0;
					//~ energy = sqrt(energy);
					//~ 
					//~ if (p.E() == energy)
					//~ {
						//~ k = i;
						//~ break;
					//~ }
				//~ }
//~ 
				//~ flength = 2e-6*p.momentum()/0.6; 
				//~ flength *= 0.5*(1.0 - 0.5/q2[k]);
			}
			else
			{
				const double Tp[] = {0.35, 0.7, 0.97, 1.8, 0.625, 1.718, 2.742, 3.65};
				const double Q2[] = {0.6, 1.3, 1.8, 3.3, 1.04, 3.06, 5.00, 6.77};
				
				int i;
				
				for (i = 0; i < 8; i++)
				{
					double x = p.Ek() - Tp[i]*1000;
					if (x < 0)
						x = -x;
					if (x < 10)
						break;
				}
					
				flength = 2.0*p.momentum()/Q2[i] * 1e-6;
			}												
		}
		else if (strcmp(fz.c_str(), "skat8") == 0) flength = 1e-6*p.momentum()/0.08;
		else if (strcmp(fz.c_str(), "cohl") == 0) flength = p.momentum()/(p.mass2() - p.p4()*p0); 
		else if (strcmp(fz.c_str(), "stod") == 0)
		{
			vec help;
			help.x = 0;
			help.y = 0; 
			help.z = 1; 
			
			double pt = help.dir()*p.p();
			pt = sqrt(p.momentum2() - pt*pt);
			
			flength = 2.0*p.momentum()/(p.mass2() + pt*pt);
		}
		else if (strcmp(fz.c_str(), "cosyn") == 0)
		{	
			flength = 1e-6*p.momentum()/0.6;

			double param = 0.9;
			if (p.pdg == 211 or p.pdg == -211 or p.pdg == 111) param = 0.5;
			
			double q2 = -q*q/1000000.0;
						
			double factor = (1.0 - param/q2);
			
			if (factor < 0.5) factor = 0.5;
			
			//flength *= (n-1.0)*factor;
			
			if (flength < 0.5/200.0) flength = 0.5/200.0;

			if (qel) flength = p.momentum()/(p.mass2() - p.p4()*p0);
		}
		else if (strcmp(fz.c_str(), "ranft") == 0)
		{
			vec help;
			help.x = q.x;
			help.y = q.y;
			help.z = q.z;
								
			double pt = help.dir()*p.p();
			pt = sqrt(p.momentum2() - pt*pt);
									
			double tau = ran_exp(par.tau)/200.0;
					
			flength = tau*p.momentum()*p.mass()/(p.mass2() + pt*pt);
		}
		else if (strcmp(fz.c_str(), "rl") == 0)
		{
			vec help;
			help.x = q.x;
			help.y = q.y;
			help.z = q.z;
						
			double pt = 0;
			double tau = 0.342/200.0;
			
			flength = tau*p.momentum()*p.mass()/(p.mass2() + pt*pt);
		}
		else if (strcmp(fz.c_str(), "delta") == 0)
		{
		
				double e = q.t + p0.t;
				vec pa;
				
				pa.x = q.x + p0.x;
				pa.y = q.y + p0.y;
				pa.z = q.z + p0.z;
				
				double mass = sqrt(e*e - pa.length()*pa.length());
				
				flength = pa.length()*ran_exp(1.0/120.0)/mass;
		}
		else if (strcmp(fz.c_str(), "const") == 0)
			flength = par.formation_length / 200.0;		
		else flength = 0;
	}
	else flength = 0;
	
	if (flength < 0) flength = -flength;
	
	flength *= 200.0*fermi; //from 1/MeV to fm

	return flength;
}

double formation_zone (particle &p, params &par)
{
	par.first_step = true;
	
	string fz;
	double fl = 0;
	
	for (int i=0; i<par.formation_zone.size();i++)
	{
		if(par.formation_zone[i]!=' ')
		fz += par.formation_zone[i];
	}

	if (strcmp(fz.c_str(), "nofz") != 0)
	{
		double tau = ran_exp(par.tau);
		fl = fermi*tau*p.momentum()/p.mass();
	}
						
	return fl;
}
