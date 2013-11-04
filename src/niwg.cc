#include "fsi.h"
#include "niwg_ccpi.h"
#include "niwg_ccqe.h"
#include "niwg_tech.h"
#include "niwg.h"
#include "dirs.h"

string pam0 = "-p 'qel_cc_axial_mass = 1030' ";
string pam1 = "-p 'qel_cc_axial_mass = 1236' ";
string pam2  = "-p 'qel_cc_axial_mass = 824' ";

string am0 = "_am1030";
string am1 = "_am1236";
string am2 = "_am824";

string am = "1030";

string pvm0 = "";//-p 'qel_cc_vector_mass = 840' ";
string pvm1 = "";//-p 'qel_cc_vector_mass = 940' ";
string pvm2 = "";//-p 'qel_cc_vector_mass = 740' ";
	
string vm0  = "_vm840_";
string vm1  = "_vm940_";
string vm2  = "_vm740_"; 

string ccqedir = "NIWG/ccqe/results/";
string ccpidir = "NIWG/ccpi/results/";


using namespace std;
bool actions_on[3*nofac] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
void set_am(bool change)
{
	if (change)
	{
		pam0 = "-p 'qel_cc_axial_mass = 1210' ";
		pam1 = "-p 'qel_cc_axial_mass = 1410' ";
		pam2 = "-p 'qel_cc_axial_mass = 1010' ";
	
		am0 = "_am1210";
		am1 = "_am1410";
		am2 = "_am1010";
		
		am = "1210";
	}
}

int main(int argc,char** argv)
{
	set_dirs(argv[0]);
	run("mkdir -p NIWG/");
	run("mkdir -p NIWG/ccqe");
	run("mkdir -p NIWG/ccqe/root_files");
	run("mkdir -p NIWG/ccqe/results");
	run("mkdir -p NIWG/ccqe/table");
	run("mkdir -p NIWG/ccqe/plots");
	run("mkdir -p NIWG/ccpi");
	run("mkdir -p NIWG/ccpi/root_files");
	run("mkdir -p NIWG/ccpi/results");
	run("mkdir -p NIWG/ccpi/table");
	run("mkdir -p NIWG/ccpi/plots");
	run("mkdir -p NIWG/tech/root_files");
	run("mkdir -p NIWG/tech/results");
	run("mkdir -p NIWG/tech/results/tex");
	run("mkdir -p NIWG/tech/results/dis/txt");
	run("mkdir -p NIWG/tech/results/dis/plots");
	run("mkdir -p tmp");
		
	char type;
	
	do{	
		if(system("clear"));
		cout<<"NIWG - CCQE cross section on free neutron/carbon/iron"<<endl;
		for (int i = 0; i < 55; i++) cout<<"-";
		cout<<endl<<endl<<endl<<"Select actions you want to do (type '0' to select all of them). Press 'enter' to go to start."<<endl<<endl;
		
		for (int i = 0; i < nofac; i++)
		{		
			cout<<((char)(97+i))<<" - "<<actions[i]; for (int k = 0; k < 50 - actions[i].length(); k++) cout<<" "; if (actions_on[i]) cout<<"*"; cout<<endl;
			for (int k = 0; k < 50 + 5; k++) cout<<"-"; cout<<endl;
		}
		
		cout<<endl<<endl<<endl;
		
		cout<<"NIWG - CCpi cross section on free neutron/carbon/iron"<<endl;
		for (int i = 0; i < 55; i++) cout<<"-";
		cout<<endl<<endl<<endl<<"Select actions you want to do (type '0' to select all of them). Press 'enter' to go to start."<<endl<<endl;
		
		for (int i = 0; i < nofac; i++)
		{		
			cout<<((char)(97+i+3))<<" - "<<actions[i]; for (int k = 0; k < 50 - actions[i].length(); k++) cout<<" "; if (actions_on[i+3]) cout<<"*"; cout<<endl;
			for (int k = 0; k < 50 + 5; k++) cout<<"-"; cout<<endl;
		}
		
		cout<<endl<<endl<<endl;
		
		cout<<"NIWG - technical note"<<endl;
		for (int i = 0; i < 55; i++) cout<<"-";
		cout<<endl<<endl<<endl<<"Select actions you want to do (type '0' to select all of them). Press 'enter' to go to start."<<endl<<endl;
		
		for (int i = 0; i < nofac; i++)
		{		
			cout<<((char)(97+i+6))<<" - "<<actions[i]; for (int k = 0; k < 50 - actions[i].length(); k++) cout<<" "; if (actions_on[i+6]) cout<<"*"; cout<<endl;
			for (int k = 0; k < 50 + 5; k++) cout<<"-"; cout<<endl;
		}
								
		type = getch();
		
		if (type == '0') for (int i = 0; i < 3*nofac; i++) actions_on[i] = 1;
		else for (int i = 0; i < 3*nofac; i++) if (type == ((char)(97+i))) actions_on[i] = !actions_on[i];
	
	}while(type != 10);
	
	set_am(axialmass);
	
	if (actions_on[0]) ccqe_sim();
	if (actions_on[1]) ccqe_calc();		
	if (actions_on[2]) ccqe_plot();
	if (actions_on[3]) ccpi_sim();
	if (actions_on[4]) ccpi_calc();
	if (actions_on[5]) ccpi_plot();

	set_am(!axialmass);

	if (actions_on[0]) ccqe_sim();
	if (actions_on[1]) ccqe_calc();	
	if (actions_on[2]) ccqe_plot();
	if (actions_on[3]) ccpi_sim();
	if (actions_on[4]) ccpi_calc();
	if (actions_on[5]) ccpi_plot();
	
	if (actions_on[6]) tech_sim();
	if (actions_on[7]) tech_calc();
	
	run("rm -r tmp/");
	
	return 1;
}
