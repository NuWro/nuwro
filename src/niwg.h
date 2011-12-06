#ifndef _niwg_h_
#define _niwg_h_
#include <string>
using namespace std;

const string qelcc = "-p 'number_of_test_events = 5000000' -p 'number_of_events = 100' -p 'beam_particle = 14' -p 'qel_vector_ff_set = 2' -p 'dyn_qel_cc = 1' -p 'dyn_qel_nc = 0' -p 'dyn_res_cc = 0' -p 'dyn_res_nc = 0' -p 'dyn_dis_cc = 0' -p 'dyn_dis_nc = 0' -p 'dyn_coh_cc = 0' -p 'dyn_coh_nc = 0' ";
const string ccpi = "-p 'number_of_test_events = 5000000' -p 'number_of_events = 10000' -p 'beam_particle = 14' -p 'qel_vector_ff_set = 2' -p 'dyn_qel_cc = 1' -p 'dyn_qel_nc = 0' -p 'dyn_res_cc = 1' -p 'dyn_res_nc = 0' -p 'dyn_dis_cc = 1' -p 'dyn_dis_nc = 0' -p 'dyn_coh_cc = 1' -p 'dyn_coh_nc = 0' ";
const int nofac = 3;
const string actions[nofac] = {"Make simulations", "Make calculations", "Make plot/table"};
extern bool actions_on[3*nofac];
const string target[4] = {"neutron", "oxygen_FG", "oxygen_SF", "iron"};
const bool axialmass = 0; //0 = 1030, 1 = 1210

extern  string pam0;
extern  string pam1;
extern  string pam2;
	
extern  string am0;
extern  string am1;
extern  string am2;

extern  string am;

extern  string pvm0;
extern  string pvm1;
extern  string pvm2;
	
extern  string vm0;
extern  string vm1;
extern  string vm2; 

extern  string ccqedir;
extern  string ccpidir;

const string kofpi[3] = {"ccpip_nuwro_", "ccpim_nuwro_", "ccpi0_nuwro_"};
void set_am(bool change);
/*{
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
*/

#endif
