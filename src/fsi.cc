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

double formation1 (particle &p, params &par, vect q, bool qel, vect p0, bool res, int n)//, int nofpi)
{
	double flength;
	
	string fz;
	
	for (int i=0; i<par.formation_zone.size();i++)
	{
		if(par.formation_zone[i]!=' ')
		fz += par.formation_zone[i];
	}
	
	/*if (p0.t == 0 and qel)
	{
		particle phelp;
		phelp.set_mass(mass_proton);
		phelp.set_momentum(rand_from_ball(225.0));
		
		p0.t = phelp.p4().t;
		p0.x = phelp.p4().x;
		p0.y = phelp.p4().y;
		p0.z = phelp.p4().z;
		
		q = p.p4() - p0;
	}*/
	
	if (par.first_step)
	{	
		if (strcmp(fz.c_str(), "fz") == 0)
		{
			if (qel) flength = p.momentum()/(p.mass2() - p.p4()*p0);
			else if (res and p0.t != 0)
			{
				double e = q.t + p0.t;
				vec pa;
				
				pa.x = q.x + p0.x;
				pa.y = q.y + p0.y;
				pa.z = q.z + p0.z;
				
				double mass = sqrt(e*e - pa.length()*pa.length());
				
				flength = pa.length()*ran_exp(1.0/120.0)/mass;
			}
			else if (p0.t != 0)
			{
				vec help;
				help.x = q.x;
				help.y = q.y;
				help.z = q.z;
								
				double pt = help.dir()*p.p();
				pt = sqrt(p.momentum2() - pt*pt);
					
				double tau = ran_exp(0.342)/200.0;
					
				flength = tau*p.momentum()*p.mass()/(p.mass2() + pt*pt);
			}
			else
			{
				double param = 1.0;
				if (p.pdg == 211 or p.pdg == -211 or p.pdg == 111) param = 0.7;
				flength = 2e-6*p.momentum()/param;
			}				
		}
		else if (strcmp(fz.c_str(), "nofz") == 0) flength = 0;
		else if (strcmp(fz.c_str(), "trans") == 0)
		{
			//double theta[5] = {10.58, 13.44, 12.74, 11.53, 9.09};
			
			//for (int i = 0; i < 5; i++) theta[i] *= 3.14159265/180.0;
			
			double q2[5] = {1.1, 2.15, 3.0, 3.91, 4.69};
			const double PTppi[5]   = {2.793, 3.187, 3.418, 4.077, 4.412};
			int k = 0;
			
			for (int i = 0; i < 5; i++)
			{
				int energy;
		
				energy = PTppi[i]*1000.0;
				energy = energy*energy + 140.0*140.0;
				energy = sqrt(energy);
				
				if (p.E() == energy)
				{
					k = i;
					break;
				}
			}
			
			//double pt = p.momentum()*sin(theta[k]);
			//double tau = ran_exp(1.0)/200.0;
			
			//flength = tau*p.momentum()*p.mass()/(p.mass2() + pt*pt);
			
			flength = 2e-6*p.momentum()/0.6; 
			flength *= 0.5*(1.0 - 0.5/q2[k]);
						
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
			
			flength *= (n-1.0)*factor;
			
			//if (flength < 0.5/200.0) flength = 0.5/200.0;

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
					
			double tau = ran_exp(8.0)/200.0;
					
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
		else flength = 0;
	}
	else flength = 0;
	
	if (flength < 0) flength = -flength;
	
	flength *= 200.0*fermi; //from 1/MeV to fm
	
	return flength;
}

double formation2 (particle &p, params &par, int n) //, vect q, bool qel, vect p0, bool res)//, int nofpi)
{
	par.first_step = true;
	
	double tau = ran_exp(8.0)/200.0;
	double pt = 0;
			
	double fl = tau*p.momentum()*p.mass()/(p.mass2() + pt*pt);

	//double fl = (n-1.0)*200.0*fermi*1e-6*p.momentum()/0.6;
	//if (fl < 0.5*fermi) fl = 0.5*fermi;

	string fz;
	
	for (int i=0; i<par.formation_zone.size();i++)
	{
		if(par.formation_zone[i]!=' ')
		fz += par.formation_zone[i];
	}

	if (strcmp(fz.c_str(), "nofz") == 0) fl = 0;
						
	return fl;
	
	//return formation1(p, par, q, qel, p0, res);//, nofpi);
}	

double formation3 (double v, int n)
{	
	return v*(n-1)*0.05*fermi/sqrt(1.0-v*v);
}	
