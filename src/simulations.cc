#include <string.h>
#include "data.h"
#include "simulations.h"
#include "dirs.h"

string fzp (int fz)
{
	string result = "-p 'formation_zone = ";
	result += fzwork[fz] + string("' ");
	return result;
}

string xpar (int xs)
{
	//if (xs == 1) xs = 3;
	
	stringstream temp;
	string x;

	temp << xs;
	temp >> x;
	
	string result = "-p 'xsec = ";
	result += x + string("' ");
	return result;
}

void simlog(string command)
{
	logfile<<day<<" "<<fullm<<" "<<year<<" "<<hour<<":"<<minute<<endl<<"Simulation with params: "<<endl<<endl<<command<<endl<<endl;
}

void nuint_sim()
{
	get_date();
	
	string command = "./kaskada -o root_files/Nuint.root -p 'first_step = 0' -p 'formation_zone = cosyn' -p 'beam_energy = 400 1500' -p 'beam_type = 0' " + events1m + piP + surface + iron + xpar(3);

	simlog(command);	
	run(command);	
}

void ks()
{
	string com = "./kaskada -o root_files/kaskada_new.root -p 'beam_type = 0' -p 'beam_energy = 150 3000' -p 'first_step = 0' " + fzp(0) + events1m + carbon + xpar(1) + piP + surface;
	run(com);
}

void krenz_sim ()
{
	string com = get_bin_dir()+"nuwro -o root_files/E2200_carbon_all_wo_fsi.root ";
	com += events100k + carbon + string(" -p 'beam_type = 0' -p 'beam_energy = 2200' ") + fsioff + numu + ALLdyn;
	run(com);
	
	com = get_bin_dir()+"nuwro -o root_files/E2200_proton.root ";
	com += events100k + hydrogen + string(" -p 'beam_type = 0' -p 'beam_energy = 2200' ") + numu + wocoh;
	//run(com);

	com = get_bin_dir()+"nuwro -o root_files/E2200_neutron.root ";
	com += events100k + fneutron + string(" -p 'beam_type = 0' -p 'beam_energy = 2200' ") + numu + wocoh;
	//run(com);
}

void ft_sim ()
{
	string com = get_bin_dir()+"nuwro -o root_files/E500_carbon_ft.root ";
	com += events100k + carbon + string(" -p 'beam_type = 0' -p 'beam_energy = 500' -p 'formation_zone = ranft' ") + numu + resdis;
	run(com);

	com = get_bin_dir()+"nuwro -o root_files/E1000_carbon_ft.root ";
	com += events100k + carbon + string(" -p 'beam_type = 0' -p 'beam_energy = 1000' -p 'formation_zone = ranft' ") + numu + resdis;
	run(com);

//	com = get_bin_dir()+"nuwro -o root_files/E1500_carbon_ft.root ";
//	com += events100k + carbon + string(" -p 'beam_type = 0' -p 'beam_energy = 1500' -p 'formation_zone = ranft' ") + numu + resdis;
//	run(com);

	com = get_bin_dir()+"nuwro -o root_files/E2000_carbon_ft.root ";
	com += events100k + carbon + string(" -p 'beam_type = 0' -p 'beam_energy = 2000' -p 'formation_zone = ranft' ") + numu + resdis;
	run(com);	

	com = get_bin_dir()+"nuwro -o root_files/E3000_carbon_ft.root ";
	com += events100k + carbon + string(" -p 'beam_type = 0' -p 'beam_energy = 3000' -p 'formation_zone = ranft' ") + numu + resdis;
	run(com);	

	com = get_bin_dir()+"nuwro -o root_files/E4000_carbon_ft.root ";
	com += events100k + carbon + string(" -p 'beam_type = 0' -p 'beam_energy = 4000' -p 'formation_zone = ranft' ") + numu + resdis;
	run(com);	

	com = get_bin_dir()+"nuwro -o root_files/E5000_carbon_ft.root ";
	com += events100k + carbon + string(" -p 'beam_type = 0' -p 'beam_energy = 5000' -p 'formation_zone = ranft' ") + numu + resdis;
	run(com);	

	com = get_bin_dir()+"nuwro -o root_files/E500_carbon_ft_res.root ";
	com += events100k + carbon + string(" -p 'beam_type = 0' -p 'beam_energy = 500' -p 'formation_zone = delta' ") + numu + onlyres;
	//run(com);

	com = get_bin_dir()+"nuwro -o root_files/E1000_carbon_ft_res.root ";
	com += events100k + carbon + string(" -p 'beam_type = 0' -p 'beam_energy = 1000' -p 'formation_zone = delta' ") + numu + onlyres;
	//run(com);

	com = get_bin_dir()+"nuwro -o root_files/E1500_carbon_ft_res.root ";
	com += events100k + carbon + string(" -p 'beam_type = 0' -p 'beam_energy = 1500' -p 'formation_zone = delta' ") + numu + onlyres;
	//run(com);

	com = get_bin_dir()+"nuwro -o root_files/E2000_carbon_ft_res.root ";
	com += events100k + carbon + string(" -p 'beam_type = 0' -p 'beam_energy = 2000' -p 'formation_zone = delta' ") + numu + onlyres;
	//run(com);	
}

void PiT2le (int fz, int xs)
{
	get_date();
	
	string com = get_bin_dir()+"nuwro -o root_files/PiTrans_le_carbon_" + fzwork[fz] + sep + xsec[1] + sep + date + ".root ";
	com += events100k + carbon + string(" -p 'beam_type = 0' -p 'beam_energy = 2000' ") + numu + onlydis + fzp(fz) + xpar(xs);
	run(com);
	
	com = get_bin_dir()+"nuwro -o root_files/PiTrans_le_lithium_" + fzwork[fz] + sep + xsec[3] + sep + date + ".root ";
	com += events100k + lithium + string(" -p 'beam_type = 0' -p 'beam_energy = 2000' ") + numu + onlydis + fzp(fz) + xpar(xs);
	run(com);

	com = get_bin_dir()+"nuwro -o root_files/PiTrans_le_aluminium_" + fzwork[fz] + sep + xsec[2] + sep + date + ".root ";
	com += events100k + aluminium + string(" -p 'beam_type = 0' -p 'beam_energy = 2000' ") + numu + onlydis + fzp(fz) + xpar(xs);
	run(com);
		
}

void K2K (int fz, int xs)
{
	get_date();
	
	string k2k = date + rot + events1m + numu + K2Kbeam;
	
	string command (get_bin_dir()+"nuwro -o root_files/K2K_NC_Hydrogen_1m_");
	command += k2k + hydrogen + NCwocoh;
	
	simlog(command);
	//run(command);
		
	command = get_bin_dir()+"nuwro -o root_files/K2K_NC_Oxygen_1mc_";
	command += fzwork[fz] + sep + xsec[xs] + sep + k2k + oxygen + NCdyn + fzp(fz) + xpar(xs);

	simlog(command);
	//run(command);
	
	//if (noFile("root_files/K2K_CC_Hydrogen_0k.root.txt"))
	{
		k2k = events0k + CCwocoh + numu + K2Kbeam;
		command = get_bin_dir()+"nuwro -o root_files/K2K_CC_Hydrogen_0k.root ";
		command += k2k + hydrogen;
		
		simlog(command);
		run(command);
	}
	//if (noFile("root_files/K2K_CC_Oxygen_0k.root.txt"))
	{
		k2k = events0k + CCdyn + numu + K2Kbeam;
		command = get_bin_dir()+"nuwro -o root_files/K2K_CC_Oxygen_0k.root ";
		command += k2k + carbon;
		
		simlog(command);
		run(command);
	}
}

void MB (int fz, int xs, bool anti)
{
	get_date();
	
	string mb = date + rot + events5m + numu + MBbeam;
	
	string command (get_bin_dir()+"nuwro -o root_files/MB_NC_Hydrogen_5m_");
	command += mb + hydrogen + NCwocoh;
	
	simlog(command);
	//run(command);
	
	command = get_bin_dir()+"nuwro -o root_files/MB_NC_Carbon_5m_";
	command += fzwork[fz] + sep + xsec[xs] + sep + mb + carbon + NCdyn + fzp(fz) + xpar(xs);

	simlog(command);
	run(command);

	if(anti)
	{
		string mba = sep + date + rot + events5m + antinumu + MBbeamanti;
		
		command = get_bin_dir()+"nuwro -o root_files/MB_anti_NC_Hydrogen_5m_";
		command +=  mba + hydrogen + NCwocoh;
			
		simlog(command);
		run(command);
		
		command = get_bin_dir()+"nuwro -o root_files/MB_anti_NC_Carbon_5m_";
		command += fzwork[fz] + sep + xsec[xs] + mba + carbon + NCdyn + fzp(fz) + xpar(xs);

		simlog(command);
		run(command);
	}
	
	if (noFile("root_files/MB_CC_Hydrogen_0k.root.txt"))
	{
		mb = events0k + CCwocoh + numu + MBbeam;
		command = get_bin_dir()+"nuwro -o root_files/MB_CC_Hydrogen_0k.root ";
		command += mb + hydrogen;
		
		simlog(command);
		cout << endl << command << endl << endl;
		run(command);
	}
	if (noFile("root_files/MB_CC_Carbon_0k.root.txt"))
	{
		mb = events0k + CCdyn + numu + MBbeam;
		command = get_bin_dir()+"nuwro -o root_files/MB_CC_Carbon_0k.root ";
		command += mb + carbon;
		
		simlog(command);
		run(command);
	}
}

void simMBback (int fz, int xsprob, bool anti)
{
	get_date();
	int xs = xsprob;
	if (xsprob == 1) xs = 1;
	
	string mb = date + rot + events1m + numu + MBbeam;
	
	string command = get_bin_dir()+"nuwro -o root_files/MB_Carbon_1m_";
	command += fzwork[fz] + sep + xsec[xsprob] + sep + mb + carbon + ALLdyn + fzp(fz) + xpar(xs);

	simlog(command);
	run(command);

	if(anti)
	{
		string mba = sep + date + rot + events1m + antinumu + MBbeamanti;
		
		command = get_bin_dir()+"nuwro -o root_files/MB_anti_Carbon_1m_";
		command += fzwork[fz] + sep + xsec[xsprob] + mba + carbon + ALLdyn + fzp(fz) + xpar(xs);

		simlog(command);
		run(command);
	}
}

void MBCC (int fz, int xsprob)
{
	get_date();
	string mb = date + rot + events100k + numu + MBbeam;

	int xs = xsprob;
	if (xsprob == 1) xs = 2;
	
	string command (get_bin_dir()+"nuwro -o root_files/MB_CC_Hydrogen_100k_");
	command += mb + hydrogen + CCwocoh;
	
	simlog(command);
	run(command);
	
	command = get_bin_dir()+"nuwro -o root_files/MB_CC_Carbon_100k_";
	command += fzwork[fz] + sep + xsec[xsprob] + sep + mb + carbon + CCdyn + fzp(fz) + xpar(xs);

	simlog(command);
	run(command);
}

void simMBCCtotal(int fz, int xs)
{
	for (int i = 0; i < 14; i++)
	{
		double en = MBCCnuen[i]*1000.0;
		
		stringstream temp;
		string energy;

		temp << en;
		temp >> energy;
	
		string help = date + rot + events10k + numu + string("-p 'beam_energy = ") + energy + string("' ");
		
		string com = get_bin_dir()+"nuwro -o root_files/E" + energy + string("_MB_CC_Hydrogen_100k_");
		com += help + hydrogen + CCwocoh;
	
		simlog(com);
		run(com);
	
		com = get_bin_dir()+"nuwro -o root_files/E" + energy + string("_MB_CC_Carbon_100k_");
		com += fzwork[fz] + sep + xsec[xs] + sep + help + carbon + CCdyn + fzp(fz) + xpar(xs);

		simlog(com);
		run(com);
	}
}

void simMBCCratio(int fz, int xs)
{
	get_date();
	
	for (int i = 0; i < 13; i++)
	{
		double en = mbraten[i]*1000.0;
		
		stringstream temp;
		string energy;

		temp << en;
		temp >> energy;
	
		string help = date + rot + events10k + numu + string("-p 'beam_energy = ") + energy + string("' ");
		
		string com = get_bin_dir()+"nuwro -o root_files/E" + energy + string("_MB_CC_Hydrogen_100k_");
		com += help + hydrogen + CCwocoh;
	
		simlog(com);
		run(com);
	
		com = get_bin_dir()+"nuwro -o root_files/E" + energy + string("_MB_CC_Carbon_100k_");
		com += fzwork[fz] + sep + xsec[xs] + sep + help + carbon + CCdyn + fzp(fz) + xpar(xs);

		simlog(com);
		run(com);
	}
}

void simSBCCtotal(int fz, int xs)
{
	get_date();
	double en[6] = {380, 620, 870, 1110, 1430, 2470};
	
	for (int i = 0; i < 6; i++)
	{
		stringstream temp;
		string energy;
		
		temp << en[i];
		temp >> energy;
		
		string help = date + rot + events0k + numu + string("-p 'beam_energy = ") + energy + string("' ");

		string com = get_bin_dir()+"nuwro -o root_files/E" + energy + string("_SB_CC_Hydrogen_0k_");
		com += help + hydrogen + CCwocoh;
	
		simlog(com);
		run(com);
	
		com = get_bin_dir()+"nuwro -o root_files/E" + energy + string("_SB_CC_Carbon_0k_");
		com += fzwork[fz] + sep + xsec[xs] + sep + help + carbon + CCdyn + fzp(fz) + xpar(xs) + fsioff;

		simlog(com);
		run(com);		
	}
}

void simNOMADCCtotal(int fz, int xs)
{
	get_date();
	double en[30] = {4.6, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 16.2, 18.7, 21.2, 23.7, 26.2, 28.7, 32.3, 37.3, 42.4, 47.4, 54.6, 64.7, 74.8, 84.8, 94.8, 107, 122, 136.9, 165.9, 228.3};
	
	for (int i = 29; i < 30; i++)
	{
		stringstream temp;
		string energy;
		
		temp << en[i]*1000.0;
		temp >> energy;
		
		string help = date + rot + events0k + numu + string("-p 'beam_energy = ") + energy + string("' ");

		string com = get_bin_dir()+"nuwro -o root_files/E" + energy + string("_NOMAD_CC_Carbon_0k_");
		com += help + carbon + CCdyn + fsioff;
	
		simlog(com);
		run(com);
	}
}

void simsfg()
{
	for (int i = 5; i < 30; i++)
	{
		stringstream temp;
		string energy;
		
		temp << i*50;
		temp >> energy;
		
		string com  = get_bin_dir()+"nuwro -o root_files/E" + energy + string("_carbon_fg.root ") + events0k + numu + string("-p 'beam_energy = ") + energy + string("' ") + carbon + QEL + fsioff;
		string com2 = get_bin_dir()+"nuwro -o root_files/E" + energy + string("_carbon_sf.root ") + events0k + numu + string("-p 'beam_energy = ") + energy + string("' ") + carbon + QEL + fsioff + sf;

		simlog(com);
		run(com);

		simlog(com2);
		run(com2);
	}
}

void simMINOSCCtotal(int fz, int xs)
{
	get_date();
	double en[13] = {3.48, 4.45, 5.89, 7.97, 10.45, 13.43, 16.42, 19.87, 23.88, 27.89, 32.81, 38.87, 45.77};
	double enbar[11] = {6.07, 7.99, 10.43, 13.42, 16.41, 19.82, 23.82, 27.84, 32.72, 38.74, 45.61};
	
	for (int i = 0; i < 13; i++)
	{
		stringstream temp;
		string energy;
		
		temp << en[i]*1000.0;
		temp >> energy;
		
		string help = date + rot + events0k + numu + string("-p 'beam_energy = ") + energy + string("' ");

		string com = get_bin_dir()+"nuwro -o root_files/E" + energy + string("_MINOS_CC_Iron_0k_");
		com += help + iron + CCdyn + fsioff;
	
		simlog(com);
		//run(com);
		
		if (i < 11)
		{
			stringstream tempbar;
			string energybar;
		
			tempbar << enbar[i]*1000.0;
			tempbar >> energybar;
			
			string helpbar = date + rot + events0k + antinumu + string("-p 'beam_energy = ") + energybar + string("' ");

			string combar = get_bin_dir()+"nuwro -o root_files/E" + energybar + string("_MINOS_CC_Iron_0k_anti");
			combar += helpbar + iron + CCdyn + fsioff;

			simlog(combar);
			run(combar);
		}			
	}	
}

void PNS (int fz, int xsprob)
{
	get_date();
	int xs = xsprob;
	
	string pns = fzwork[fz] + sep + xsec[xsprob] + sep + date + rot + events100k + piP + surface + PNSbeam + fzp(fz) + string("-p 'first_step = 0' ");
	
	string command = "./bin/kaskada -o root_files/PionNucleus_Carbon_100k_";
	command +=  pns + carbon + xpar(xs);
	
	simlog(command);
	run(command);

	if (xsprob == 1) xs = 2;
	
	command = "./bin/kaskada -o root_files/PionNucleus_Iron_100k_";
	command += pns + iron + xpar(xs);

	simlog(command);	
	run(command);
}

void PrT (int fz, int xs, char *which)
{
	get_date();
	string prt = fzwork[fz] + sep + xsec[xs] + sep + date + rot + events100k + proton + pos_ran + fzp(fz) + xpar(xs);
			
	if (strcmp(which, "high_energy") == 0)
	{
		string command = "./bin/kaskada -o root_files/ProtonTransparency_he_Carbon_100k_";
		command += prt + PrThe + carbon;
	
		simlog(command);
		run(command);
		
		command = "./bin/kaskada -o root_files/ProtonTransparency_he_Iron_100k_";
		command += prt + PrThe + iron;
	
		simlog(command);
		run(command);
	}
	else if (strcmp(which, "low_energy") == 0)
	{
		string command = "./kaskada -o root_files/ProtonTransparency_le_Lithium_100k_";
		command += prt + PrTle + lithium;
		
		simlog(command);
		run(command);
	
		command = "./kaskada -o root_files/ProtonTransparency_le_Carbon_100k_";
		command += prt + PrTle + carbon;
		
		simlog(command);
		run(command);

		command = "./kaskada -o root_files/ProtonTransparency_le_Aluminium_100k_";
		command += prt + PrTle + aluminium;
	
		simlog(command);
		run(command);		
	}
		
}

void PiT2 (int fz, int xs)
{
	get_date();
	string pit = fzwork[fz] + sep + xsec[xs] + sep + date + rot + events1m + fzp(fz);
	
	string command = get_bin_dir()+"nuwro -o root_files/PiTrans_he_Carbon_1m_";
	command += pit + pitbeam + carbon + onlydis + xpar(1);
	
	run(command);
	
	command = get_bin_dir()+"nuwro -o root_files/PiTrans_he_Aluminium_1m_";
	command += pit + pitbeam + aluminium + onlydis + xpar(3);
	
	run(command);

	command = get_bin_dir()+"nuwro -o root_files/PiTrans_he_Copper_1m_";
	command += pit + pitbeam + copper + onlydis + xpar(2);
	
	run(command);
}	

void PiT3(int fz)
{
	for (int i = 0; i < 5; i++)
	{
		int energy;
		
		energy = PTppi[i]*1000.0;
		energy = energy*energy + 140.0*140.0;
		energy = sqrt(energy);
		
		stringstream temp;
		string en;
		
		temp << energy;
		temp >> en;
		
		string command = "./bin/kaskada -o root_files/PiTrans_e" + en + string("_carbon_100k_") + fzwork[fz] + rot + xpar(1) + carbon + piP + events100k + pos_ran + fzp(fz) + string("-p 'beam_energy = ") + en + string("' ");
		run(command);
		command = "./bin/kaskada -o root_files/PiTrans_e" + en + string("_aluminium_100k_") + fzwork[fz] + rot + xpar(3) + aluminium + piP + events100k + pos_ran + fzp(fz) + string("-p 'beam_energy = ") + en + string("' ");
		run(command);		
		command = "./bin/kaskada -o root_files/PiTrans_e" + en + string("_copper_100k_") + fzwork[fz] + rot + xpar(2) + copper + piP + events100k + pos_ran + fzp(fz) + string("-p 'beam_energy = ") + en + string("' ");
		run(command);			
	}
}

void PiT (int fz, int xs, char *which)
{
	get_date();
	string pit = fzwork[fz] + sep + xsec[xs] + sep + date + rot + events100k + pos_ran + fzp(fz);	
	
	if (strcmp(which, "high_energy") == 0)
	{
		string command = "./bin/kaskada -o root_files/PionTransparency_he_Carbon_100k_";
		command += pit + PiThe + piP + carbon + xpar(1);
		
		simlog(command);
		run(command);

		command = "./bin/kaskada -o root_files/PionTransparency_he_Aluminium_100k_";
		command += pit + PiThe + piP + aluminium + xpar(3);
		
		simlog(command);
		run(command);
		
		command = "./bin/kaskada -o root_files/PionTransparency_he_Copper_100k_";
		command += pit + PiThe + piP + copper + xpar(2);
		
		simlog(command);
		run(command);
	}
	else if (strcmp(which, "low_energy") == 0)
	{
		string command = "./bin/kaskada -o root_files/PionTransparency_le_pip_Lithium_100k_";
		command += pit + PiTle + piP + lithium + xpar(3);
		
		simlog(command);
		run(command);
	
		command = "./bin/kaskada -o root_files/PionTransparency_le_pip_Carbon_100k_";
		command += pit + PiTle + piP + carbon + xpar(1);
		
		simlog(command);
		run(command);

		command = "./bin/kaskada -o root_files/PionTransparency_le_pip_Aluminium_100k_";
		command += pit + PiTle + piP + aluminium + xpar(3);
		
		simlog(command);
		run(command);
		
		command = "./bin/kaskada -o root_files/PionTransparency_le_pi0_Lithium_100k_";
		command += pit + PiTle + pi0 + lithium + xpar(3);
		
		simlog(command);
		run(command);
	
		command = "./bin/kaskada -o root_files/PionTransparency_le_pi0_Carbon_100k_";
		command += pit + PiTle + pi0 + carbon + xpar(2);
		
		simlog(command);
		run(command);

		command = "./bin/kaskada -o root_files/PionTransparency_le_pi0_Aluminium_100k_";
		command += pit + PiTle + pi0 + aluminium + xpar(2);
		
		simlog(command);
		run(command);	
	}
}

void Nomad (int fz, int xsprob)
{
	get_date();
	int xs = xsprob;
//	if (xsprob == 1) xs = 3;
	
	string command = get_bin_dir()+"nuwro -o root_files/Nomad_100k_";
	command += fzwork[fz] + sep + xsec[xsprob] + sep + date + rot + events5m+ CCdyn + "-p '@beam/nomad.txt' " + fzp(fz) + carbon + xpar(xs); // + nomdet;
	
	simlog(command);
	run(command);
}

void AtmNu (int fz, int xs)
{
	get_date();
	string command = get_bin_dir()+"nuwro -o root_files/AtmNu_100k_";
	command += fzwork[fz] + sep + xsec[xs] + sep + date + rot + events100k + CCdyn + numu + Atmbeam + neon + fzp(fz) + xpar(xs);
	
	simlog(command);
	run(command);
}

void Multiplicity (int fz, int xs)
{
	get_date();
	string pns = fzwork[fz] + sep + xsec[xs] + sep + date + rot + events100k + surface + string("-p 'beam_energy = 3000' ") + fzp(fz) + string("-p 'first_step = 0' ") + barium + xpar(xs);
	
	string command = "./kaskada -o root_files/multiplicity_pip_3GeV_barium_";
	command +=  pns + piP;
	
	simlog(command);
	run(command);
	
	//pns = fzwork[fz] + sep + date + rot + events100k + surface + string("-p 'beam_energy = 1020' ") + fzp(fz) + string("-p 'first_step = 0' ") + barium;
	
	//command = "./kaskada -o root_files/multiplicity_proton_1GeV_barium_";
	//command += pns + proton;
	
	//simlog(command);
	//run(command);
}

void roman_sim ()
{
	double en[21] = {200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000};
	
	for (int i = 0; i < 21; i++)
	{
		stringstream temp;
		string energy;
		
		temp << en[i];
		temp >> energy;
			
		string command = get_bin_dir()+"nuwro -o root_files/E" + energy + "_carbon_100k.root -p 'formation_zone = fz' -p 'xsec = 2' -p 'beam_energy = " + energy + string("' ") + events100k + numu + carbon + CCdyn;
		simlog(command);
		if (i == 3 or i == 8 or i == 14) run(command);
	}
}

void oset_sim ()
{
	get_date();
	
	string command = "./kaskada -o root_files/pionE85_1m.root -p 'formation_zone = nofz' -p 'beam_energy = 225' -p 'xsec = 2' -p 'first_step = 0' " + events1m + piP + surface + carbon;
	
	simlog(command);
	run(command);
	
	command = "./kaskada -o root_files/pionE245_1m.root -p 'formation_zone = nofz' -p 'beam_energy = 385' -p 'xsec = 2' -p 'first_step = 0' " + events1m + piP + surface + carbon;

	simlog(command);	
	run(command);

	command = "./kaskada -o root_files/pionE180_1m.root -p 'formation_zone = nofz' -p 'beam_energy = 320' -p 'xsec = 1' -p 'first_step = 0' " + events1m + pi0 + surface + calcium;

	simlog(command);	
	run(command);
}

void t2k_anal_sim()
{
	string com1 = get_bin_dir()+"nuwro -o root_files/T2K_carbon_fg.root -p 'formation_zone = cosyn' -p 'xsec = 1' " + to1500 + events1m + carbon + ALLdyn;
	string com2 = get_bin_dir()+"nuwro -o root_files/T2K_carbon_sf.root -p 'formation_zone = cosyn' -p 'xsec = 1' " + to1500 + events1m + carbon + ALLdyn + sf;
	string com3 = get_bin_dir()+"nuwro -o root_files/T2K_oxygen_fg.root -p 'formation_zone = cosyn' -p 'xsec = 1' " + to1500 + events1m + oxygen + ALLdyn;
	string com4 = get_bin_dir()+"nuwro -o root_files/T2K_oxygen_sf.root -p 'formation_zone = cosyn' -p 'xsec = 1' " + to1500 + events1m + oxygen + ALLdyn + sf;
		
	//run(com1);
	run(com2);
	//run(com3);
	//run(com4);
}

void dens_test_sim()
{
	string com = get_bin_dir()+"nuwro -o root_files/E40k_carbon.root -p 'formation_zone = nofz' -p 'xsec = 1' " + E40k + evtest + carbon + ALLdyn;
	run(com);
}

void hayato_sim1()
{
	run("mkdir -p root_files/hayato");
	
	for (int i = 0; i < 30; i++)
	{
		double en = i*25.0 + 250;
		
		stringstream temp;
		string energy;

		temp << en;
		temp >> energy;
	
		string com1 = get_bin_dir()+"nuwro -o root_files/hayato/E" + energy + "_oxygen_fg.root -p 'beam_type = 0' - p 'beam_energy = " + energy + "' " + events25k + oxygen + fsioff + NCQEL;
		string com2 = get_bin_dir()+"nuwro -o root_files/hayato/E" + energy + "_oxygen_sf.root -p 'beam_type = 0' - p 'beam_energy = " + energy + "' " + events25k + oxygen + fsioff + NCQEL + sf;
		
		run(com1);
		run(com2);
	}
}

void hayato_sim2()
{
	run("mkdir -p root_files/hayato");
	
	for (int i = 0; i < 30; i++)
	{
		double en = i*50.0 + 1000.0;
		
		stringstream temp;
		string energy;

		temp << en;
		temp >> energy;
	
		string com1 = get_bin_dir()+"nuwro -o root_files/hayato/E" + energy + "_oxygen_fg.root -p 'beam_type = 0' - p 'beam_energy = " + energy + "' " + events25k + oxygen + fsioff + NCQEL;
		string com2 = get_bin_dir()+"nuwro -o root_files/hayato/E" + energy + "_oxygen_sf.root -p 'beam_type = 0' - p 'beam_energy = " + energy + "' " + events25k + oxygen + fsioff + NCQEL + sf;
		
		run(com1);
		run(com2);
	}
}

void hayato_sim3()
{
	run("mkdir -p root_files/hayato");
	
	for (int i = 0; i < 30; i++)
	{
		double en = i*250.0 + 2500.0;
		
		stringstream temp;
		string energy;

		temp << en;
		temp >> energy;
	
		string com1 = get_bin_dir()+"nuwro -o root_files/hayato/E" + energy + "_oxygen_fg.root -p 'beam_type = 0' - p 'beam_energy = " + energy + "' " + events25k + oxygen + fsioff + NCQEL;
		string com2 = get_bin_dir()+"nuwro -o root_files/hayato/E" + energy + "_oxygen_sf.root -p 'beam_type = 0' - p 'beam_energy = " + energy + "' " + events25k + oxygen + fsioff + NCQEL + sf;
		
		run(com1);
		run(com2);
	}
}

void hayato_sim4()
{
	run("mkdir -p root_files/hayato");
	
	for (int i = 0; i < 21; i++)
	{
		double en = i*1000.0 + 10000.0;
		
		stringstream temp;
		string energy;

		temp << en;
		temp >> energy;
	
		string com1 = get_bin_dir()+"nuwro -o root_files/hayato/E" + energy + "_oxygen_fg.root -p 'beam_type = 0' - p 'beam_energy = " + energy + "' " + events25k + oxygen + fsioff + NCQEL;
		string com2 = get_bin_dir()+"nuwro -o root_files/hayato/E" + energy + "_oxygen_sf.root -p 'beam_type = 0' - p 'beam_energy = " + energy + "' " + events25k + oxygen + fsioff + NCQEL + sf;
		
		run(com1);
		run(com2);
	}
}

void ptsim()
{
	for (int i = 0; i < 10; i++)
	{
		double en = (i+1.0)*500.0;
		
		stringstream temp;
		string energy;

		temp << en;
		temp >> energy;
			
		string com = get_bin_dir()+"nuwro -o root_files/ptE" + energy + rot + events100k + numu + oxygen + resdis + string("-p 'beam_energy = ") + energy + string("' ");
	
		simlog(com);
		run(com);
	}
}

void angle_test()
{
	for (int i = 0; i < 10; i++)
	{
		double en = (i+1.0)*100.0 + 100;
		
		stringstream temp;
		string energy;

		temp << en;
		temp >> energy;
			
		string com = "./kaskada -o root_files/iiE" + energy + rot + events100k + piP + string("-p 'nucleus_p = 10' -p 'nucleus_n = 0' -p 'beam_type = 0' -p 'beam_energy = ") + energy + string("' ");
	
		simlog(com);
		run(com);

		com = "./kaskada -o root_files/ijE" + energy + rot + events100k + piM + hydrogen + string("-p 'nucleus_p = 10' -p 'nucleus_n = 0' -p 'beam_type = 0' -p 'beam_energy = ") + energy + string("' ");
	
		simlog(com);
		run(com);
	}
}
