#include "nuwro.h"
#include "dis/singlepion.h"
#include "beam.h"
#include "ff.h"
#include <string>
#include "kaskada7.h"

const int nof_dyn = 10; //number of dynamics in NuWro
const int nof_class = 60;

float mb_nce_start = 0;
int mb_nce_nof = 0;
int mb_nce_mode = 0; //0 - chi2, 1 - only qel, 2 - qel2 (not fixed axial mass for hydrogen), 3 - qel on hydrogen (1030MeV), 4 - mec, 5 - res/dis
int mb_nce_he_ma = 1030;

int mb_nce_ma = 0;
double mb_nce_ds = 0;
int mb_nce_ma_step = 0;
int mb_nce_ma_start = 0;
int mb_nce_ma_end = 0;
double mb_nce_ds_step = 0;
double mb_nce_ds_start = 0;
double mb_nce_ds_end = 0;

double norm [nof_dyn];
double xsec [nof_dyn];

const bool active_class[nof_class] =
{
	0,	//mb NCpi0
	0,	//mb CCpi0
	0,	//mb CCpi+
	0,	//NUINT12 A
	0,	//NUINT12 antiA
	0,	//NUINT12 B
	0,	//NUINT12 antiB
	0,	//NUINT12 C1
	0,	//NUINT12 antiC1
	0,	//NUINT12 C2
	0,	//NUINT12 antiC2
	0,	//NUINT12 D
	0,	//NUINT12 E
	0,	//NUINT12 F
	0,	//NUINT12 G1
	0,	//NUINT12 G1_2
	0,	//NUINT12 G1_3
	0,	//NUINT12 G1_4 (src/particle.h:214: void particle::set_energy(double): Assertion `E>=_mass' failed.)
	0,	//NUINT12 G1_5
	0,	//NUINT12 G1_6
	0,	//NUINT12 G2
	0,	//NUINT12 G2_2
	0,	//NUINT12 G2_3
	0,	//NUINT12 G2_4
	0,	//NUINT12 H
	0,	//NUINT12 antiH
	0,  //NUINT12 B2
	0,  //NUINT12 antiB2
	0,	//MEC OLD CC
	0,	//MEC NEW CC
	0,	//MEC 2 OLD CC
	0,	//MEC 2 NEW CC
	0,	//MEC NC
	0,	//MiniBooNE NC Elastic
	0,	//MiniBooNE NC Elastic ratio
	0,	//Phil total	cc pi+ p
	0,	//Phil total	cc pi+ n
	0,	//Phil total	cc pi0 p
	0,	//Phil total	nc pi0 p
	0,	//Phil total	nc pi0 n
	0,	//Phil total	nc pi- p
	0,	//Phil różniczkowy
	0,	//MEC MB NC cherenkov treshold test
	0,	//MB NCEL BOTH
	0,	//MB NCEL TRANSPARENCY TEST
	0,	//MB NCEL MEC MODELS COMPARE
	0,	//NIEVES KINEMATICS IN TEM MODEL TEST
	0,	//energy test
	0,	//NIWG nucleons momentum in Nieves after energy balance update
	0,	//NIWG lepton kinematics for events with protons above treshold (Nieves)
	0,	//Bodek nu - Q2/2M
	0,	//phd1 -> proton out momentum distributions fg vs lfg vs sf
	1,	//phd2 -> qel cross section...
	0,	//backward muon in MEC
	0,	//backward muon in MEC 2
	0,	//phd3 -> total mec
	0,	//phd4 -> nucleon kin mec
	0,	//phd5 -> total COH
	0,	//phd6 -> proton transparency
	0	//phd7 -> pion transparency
};

void run_command (string com)
{
	FILE *fp;
	fp = popen (com.c_str(), "w");
	pclose (fp);
}

void parabola (double &a, double &b, double &c, const double* x, const double& start, const double &step, const int bins);

class histogram
{
	public:
	
		string name;
		string data;
		string title;
		string xlabel;
		string ylabel;
		bool finalized;
		bool saved;
		bool ploted;
		
		int bins_x;
		double width_x;
		double begin_x;
		int *tnorm [nof_dyn]; //use for total cross section
		int norm_type; 
		/* 0 - cross section
		 * 1 - to number of events defined by normalization
		 * 2 - for multicase - each case nomalized separately 
		 * 3 - for multicase - each normalized together
		 * 4 - for multicase - each normalized respect to number of events in 1st case
		 * 5 - events / total
		 * 6 - total in neutrino energy		
		 * 7 - to number of events defined by normalization in neutrino energy
		 * 8 - sum all dynamics
		*/
		
		double normalization;
		
		virtual void save () {};
		virtual void plot () {};
		virtual void finalize () {};
		void put (double x, int dyn);
		
		histogram (string name_of_hist, string data_file, string plot_title, string x_name, string y_name, int nof_bins, double begin, double end, int typnormy, double norma);
		virtual ~histogram ();
}; 

class hist1D : public histogram
{
	public:
	
		double *part_result [nof_dyn];
		double *result;
		
		hist1D (string name_of_hist, string data_file, string plot_title, string x_name, string y_name, int nof_bins, double begin, double end, int typnormy, double norma = 0);
		~hist1D ();
		
		void save ();
		void plot ();
		void finalize ();
		void put (double x);
		void put (double x, int dyn, double weight);
		void put_over (double x, int dyn, double weight);
};

class hist2D : public histogram
{
	public:
	
		string zlabel;
		
		int bins_y;
		double width_y;
		double begin_y;
		double max;
		
		double **part_result [nof_dyn];
		double **result;
		
		hist2D (string name_of_hist, string data_file, string plot_title, string x_name, string y_name, string z_name, int nof_bins_x, double begin_x, double end_x, int nof_bins_y, double begin_y, double end_y, int typnormy, double norma = 0);
		~hist2D ();
		
		void save ();
		void plot ();
		void finalize ();
		void put (double x, double y, int dyn, double weight);
};

class mhist1D : public histogram
{
	public:
	
		double **part_result [nof_dyn];
		double **result;
		int cases;
		string *cnames;
		
		mhist1D (string name_of_hist, string data_file, string plot_title, string x_name, string y_name, int nof_bins, double begin, double end, int nof_cases, int typnormy, double norma = 0);
		~mhist1D ();
		
		void finalize ();
		void save ();
		void plot ();
		void put (double x, int dyn, double weight, int cas);
};	

class pattern
{
	protected:
		
		bool onlyfsi;
		int loop;
		double *E;
		string name;
		NuWro *N;
		kaskada *k;
		params P;
		int events;

		void run ();
		void runfsi ();
		virtual void set_params () {}
		virtual void calculate(event *e) {}
		virtual bool change_params () {return true;}
		
	public:
		
		pattern (string name_of_class, int nof_events, bool fsi = false, int l = 1) : name (name_of_class), events (nof_events), onlyfsi(fsi), loop(l)
		{			
			for (int i  = 0; i < nof_dyn; i++)
			{
				norm [i] = 0;
				xsec [i] = 0;
			}
		}
		
		inline bool fsi() {return onlyfsi;}
		
		virtual ~pattern () {}		
		virtual void start ();
		virtual void startfsi ();
};

class mb_ncpi0 : public pattern
{
	protected:
		
		hist1D *h1;
		hist1D *h2;
		
		void set_params ();
		void calculate (event *e);
		
	public:
	
		mb_ncpi0 () : pattern ("MB ncpi0 production", 5000000)
		{
			run_command ("mkdir -p analysis/mb/");
			
			h1 = new hist1D ("mb/mb_ncpi0_pi_mom", "data/data/mb_ncpi0_pi_mom.txt", "MiniBooNE NC {/Symbol p}^{0} production on CH_{2}", "{/Symbol p}^{0} momentum [GeV]", "d{/Symbol s} / dp_{/Symbol p} [cm^{2} c / GeV / nucleon]", 30, 0, 1.5, 0);
			h2 = new hist1D ("mb/mb_ncpi0_pi_cos", "data/data/mb_ncpi0_pi_cos.txt", "MiniBooNE NC {/Symbol p}^{0} production on CH_{2}", "{/Symbol p}^{0} cos{/Symbol q}", "d{/Symbol s} / dcos{/Symbol q} [cm^{2} / nucleon]", 20, -1, 1, 0);
		}
		
		~mb_ncpi0 ()
		{
			delete h1;
			delete h2;
		}
};

class mb_ccpi0 : public pattern
{
	protected:
		
		hist1D *h1;
		hist1D *h2;
		hist1D *h3;
		hist1D *h4;
		hist1D *h5;
		
		void set_params ();
		void calculate (event *e);

	public:
	
		mb_ccpi0 () : pattern ("MB ccpi0 production", 5000000)
		{
			run_command ("mkdir -p analysis/mb/");
			
			h1  = new hist1D ("mb/mb_ccpi0_total", "data/data/mb_ccpi0_total.txt", "MiniBooNE CC {/Symbol p}^{0} production on CH_{2}", "Neutrino Energy [GeV]", "{/Symbol s} (E_{/Symbol n}) [cm^{2} / CH_{2}]", 15, 0.5, 2, 0);
			h2  = new hist1D ("mb/mb_ccpi0_q2", "data/data/mb_ccpi0_q2.txt", "MiniBooNE CC {/Symbol p}^{0} production on CH_{2}", "Q^{2} [GeV^{2}/c^{4}]", "d{/Symbol s} / dQ^{2} [cm^{2}c^{4} / GeV^{2} / CH_{2}]", 20, 0, 2, 0);
			h3  = new hist1D ("mb/mb_ccpi0_pi_mom", "data/data/mb_ccpi0_pi_mom.txt", "MiniBooNE CC {/Symbol p}^{0} production on CH_{2}", "{/Symbol p}^{0} momentum [GeV]", "d{/Symbol s} / dp_{/Symbol p}) [cm^{2} c / GeV / CH_{2}]", 28, 0, 1.4, 0);
			h4 = new hist1D ("mb/mb_ccpi0_pi_cos", "data/data/mb_ccpi0_pi_cos.txt", "MiniBooNE CC {/Symbol p}^{0} production on CH_{2}", "{/Symbol p} cos{/Symbol q}", "d{/Symbol s} / dcos{/Symbol q}) [cm^{2} / CH_{2}]", 20, -1, 1, 0);
			h5 = new hist1D ("mb/mb_ccpi0_mu_ek", "data/data/mb_ccpi0_mu_ek.txt", "MiniBooNE CC {/Symbol p}^{0} production on CH_{2}", "{/Symbol m} kinetic energy [GeV]", "d{/Symbol s} / d(KE_{/Symbol m}) [cm^{2} / GeV / CH_{2}]", 24, 0, 1.2, 0);
		}
		
		~mb_ccpi0 ()
		{
			delete h1;
			delete h2;
			delete h3;
			delete h4;
			delete h5;
		}
};

class mb_ccpip : public pattern
{
	protected:
		
		hist1D *h1;
		hist1D *h2;
		hist1D *h3;
		hist1D *h4;
		
		void set_params ();
		void calculate (event *e);
	
	public:
	
		mb_ccpip () : pattern ("MB ccpip production", 5000000)
		{
			run_command ("mkdir -p analysis/mb/");
			
			h1 = new hist1D ("mb/mb_ccpip_total", "data/data/mb_ccpip_total.txt", "MiniBooNE CC {/Symbol p}^{+} production on CH_{2}", "Neutrino Energy [MeV]", "{/Symbol s} (E_{/Symbol n}) [cm^{2}]", 30, 500, 2000, 0);
			h2 = new hist1D ("mb/mb_ccpip_q2", "data/data/mb_ccpip_q2.txt", "MiniBooNE CC {/Symbol p}^{+} production on CH_{2}", "Q^{2} [GeV^{2}/c^{4}]", "d{/Symbol s} / dQ^{2} [cm^{2}c^{4} / GeV^{2}]", 28, 0, 1.4, 0);
			h3 = new hist1D ("mb/mb_ccpip_pi_ek", "data/data/mb_ccpip_pi_ek.txt", "MiniBooNE CC {/Symbol p}^{+} production on CH_{2}", "{/Symbol p}^{+} kinetic energy [MeV]", "d{/Symbol s} / d(KE_{/Symbol p}) [cm^{2} / MeV]", 16, 0, 400, 0);
			h4 = new hist1D ("mb/mb_ccpip_mu_ek", "data/data/mb_ccpip_mu_ek.txt", "MiniBooNE CC {/Symbol p}^{+} production on CH_{2}", "{/Symbol m} kinetic energy [MeV]", "d{/Symbol s} / d(KE_{/Symbol m}) [cm^{2} / MeV]", 30, 0, 1500, 0);
			
		}
		
		~mb_ccpip ()
		{
			delete h1;
			delete h2;
			delete h3;
			delete h4;
		}
};

class A : public pattern
{
	protected:
	
		hist1D *h1;
		hist1D *h2;
		hist2D *h3;
		hist2D *h3b;
		hist1D *h4;
		hist1D *h4r;
		hist1D *h4d;
		hist1D *h4c;
		hist1D *h4b;
		
		void set_params ();
		void calculate (event *e);
		
	public:
	
		A () : pattern ("NUINT12 A", 5000000)
		{
			h1 = new hist1D ("A21", "", "{/Symbol n}_{/Symbol m}, E = 2 GeV, CC, on carbon, no mesons in FS", "#proton", "#events", 10, 0, 10, 1, 100);
			h2 = new hist1D ("A23", "", "{/Symbol n}_{/Symbol m}, E = 2 GeV, CC, on carbon, no mesons in FS", "{/Symbol n}_{/Symbol m} reconstructed energy [MeV]", "#events", 40, 0, 2000, 1, 100);
			h3 = new hist2D ("A31b", "", "{/Symbol n}_{/Symbol m}, E = 2 GeV, NC, on carbon, at least 1pi0", "{/Symbol p}^{0} energy [MeV]", "cos{/Symbol q}", "xsec per nucleus", 20, 0, 2000, 20, -1, 1, 0);
			h3b = new hist2D ("A31", "", "{/Symbol n}_{/Symbol m}, E = 2 GeV, NC, on carbon, at least 1pi0", "{/Symbol p}^{0} energy [MeV]", "cos{/Symbol q}", "#events", 20, 0, 2000, 20, -1, 1, 1, 100);
			h4b = new hist1D ("A11b", "", "{/Symbol n}_{/Symbol m}, E = 2 GeV, NC, on carbon", "Total energy in {/Symbol p}^{0}s / energy transfer", "#events", 25, 0, 1, 1, 100);
			h4r = new hist1D ("A11_res", "", "{/Symbol n}_{/Symbol m}, E = 2 GeV, NC, on carbon", "Total energy in {/Symbol p}^{0}s / energy transfer", "#events", 25, 0, 1, 1, 100);
			h4d = new hist1D ("A11_dis", "", "{/Symbol n}_{/Symbol m}, E = 2 GeV, NC, on carbon", "Total energy in {/Symbol p}^{0}s / energy transfer", "#events", 25, 0, 1, 1, 100);
			h4c = new hist1D ("A11_coh", "", "{/Symbol n}_{/Symbol m}, E = 2 GeV, NC, on carbon", "Total energy in {/Symbol p}^{0}s / energy transfer", "#events", 25, 0, 1, 1, 100);
			h4 = new hist1D ("A11", "", "{/Symbol n}_{/Symbol m}, E = 2 GeV, NC, on carbon", "Total energy in {/Symbol p}^{0}s / energy transfer", "#events", 25, 0, 1, 1, 100);			
		}
		
		~A ()
		{
			delete h1;
			delete h2;
			delete h3;
			delete h3b;
			delete h4;
			delete h4r;
			delete h4d;
			delete h4c;
			delete h4b;
		}
};

class antiA : public A
{
	protected:
	
		void set_params ();
	
	public:
		
		antiA ()
		{
			h1 -> name = "A22";
			h1 -> title = "anti {/Symbol n}_{/Symbol m}, E = 2 GeV, CC, on carbon, no mesons in FS";
			h2 -> name = "A24";
			h2 -> xlabel = "anti {/Symbol n}_{/Symbol m} reconstructed energy [MeV]";
			h2 -> title = "anti {/Symbol n}_{/Symbol m}, E = 2 GeV, CC, on carbon, no mesons in FS";
			h3 -> name = "A32";
			h3 -> title = "anti {/Symbol n}_{/Symbol m}, E = 2 GeV, NC, on carbon, at least 1pi0";	
			h4 -> name = "A12";
			h4 -> title = "anti {/Symbol n}_{/Symbol m}, E = 2 GeV, NC, on carbon";	
		}
		
		~antiA () {}	
};

class B : public pattern
{
	protected:
	
		hist1D *h1;
		
		void set_params ();
		void calculate (event *e);
		
	public:
	
		B () : pattern ("NUINT12 B", 5000000)
		{
			h1 = new hist1D ("B11", "", "{/Symbol n}_{/Symbol m}, E = 5 GeV, CC, on iron, at least 1 {/Symbol p}^{0} in FS", "{/Symbol p}^{0} energy / energy transfer", "probabilty", 20, 0, 1, 1, 1);
		}
		
		~B ()
		{
			delete h1;
		}
};	

class antiB : public B
{
	protected:
	
		void set_params ();
		
	public:
	
		antiB ()
		{
			h1 -> name = "B12";
			h1 -> title = "anti {/Symbol n}_{/Symbol m}, E = 5 GeV, CC, on iron, at least 1 {/Symbol p}^{0} in FS";
		}
		
		~antiB () {}
};

class B2 : public pattern
{
	protected:
	
		hist1D *h1;
		
		void set_params ();
		void calculate (event *e);
		
	public:
	
		B2 () : pattern ("NUINT12 B2", 5000000)
		{
			h1 = new hist1D ("B21", "", "{/Symbol n}_{/Symbol m}, E = 3 GeV, CC, on iron", "Reconstructed energy [MeV]", "#events", 20, 1500, 3500, 1, 100);
		}
		
		~B2 ()
		{
			delete h1;
		}
};	

class antiB2 : public B2
{
	protected:
	
		void set_params ();
		
	public:
	
		antiB2 ()
		{
			h1 -> name = "B22";
			h1 -> title = "anti {/Symbol n}_{/Symbol m}, E = 3 GeV, CC, on iron";
		}
		
		~antiB2 () {}
};

class C1 : public pattern
{
	protected:
	
		mhist1D *h1;
		
		void set_params ();
		void calculate (event *e);
		
	public:
	
		C1 () : pattern ("NUINT12 C1", 5000000)
		{
			h1 = new mhist1D ("C11C12", "", "{/Symbol n}_{/Symbol m}, E = 0-10 GeV, CC, on argon", "neutrino energy [GeV]", "total xsec per Argon", 100, 0, 10, 2, 0);
			h1 -> cnames [0] = "all events (C11)";
			h1 -> cnames [1] = "events with no mesons in FS (C12)";
			
		}
		
		~C1 ()
		{
			delete h1;
		}
};

class antiC1 : public C1
{
	protected:
	
		void set_params ();
		
	public:
	
		antiC1 ()
		{
			h1 -> name = "C13C14";
			h1 -> xlabel = "antineutrino energy [GeV]";
			h1 -> title = "anti {/Symbol n}_{/Symbol m}, E = 0-10 GeV, CC, on argon";
			h1 -> cnames [0] = "all events (C13)";
			h1 -> cnames [1] = "events with no mesons in FS (C14)";
		}
		
		~antiC1 () {}
};

class C2 : public pattern
{
	protected:
	
		mhist1D *h1;
		
		void set_params ();
		void calculate (event *e);
		
	public:

		C2 () : pattern ("NUINT12 C2", 5000000)
		{
			h1 = new mhist1D ("C2122", "", "{/Symbol n}_{/Symbol m}, E = 2,5 GeV, CC, on argon", "#protons", "#events", 10, 0, 10, 2, 4, 100);
			h1 -> cnames [0] = "all protons (C21)";
			h1 -> cnames [1] = "only proton with Tk > 50 MeV (C22)"; 
		}
		
		~C2 ()
		{
			delete h1;
		}
};

class antiC2 : public C2
{
	protected:
		
		void set_params ();
		
	public:
	
		antiC2 ()
		{
			h1 -> name = "C2324";
			h1 -> title = "anti {/Symbol n}_{/Symbol m}, E = 2,5 GeV, CC, on argon";
			h1 -> cnames [0] = "all protons (C23)";
			h1 -> cnames [1] = "only proton with Tk > 50 MeV (C24)"; 
		}
	
		~antiC2 () {}
};

class D : public pattern
{
	protected:
	
		mhist1D *h1;
		mhist1D *h2;
		mhist1D *h2b;
		mhist1D *h3;
		hist1D *h4;
		
		void set_params ();
		void calculate (event *e);
		
	public:
	
		D () : pattern ("NUINT12 D", 5000000)
		{
			h1 = new mhist1D ("D1112", "", "{/Symbol n}_{/Symbol m}, E = 3 GeV, CC, on argon, no meson in FS", "#protons", "#events / #total events", 10, 0, 10, 2, 5);
			h1 -> cnames [0] = "all protons (D11)";
			h1 -> cnames [1] = "only protons with Tk > 50 MeV (D12)";
			
			h2 = new mhist1D ("D13", "", "{/Symbol n}_{/Symbol m}, E = 3 GeV, CC, on argon, no meson in FS", "momentum of most energetic proton [MeV]", "#events", 60, 0, 3000, 4, 3, 100);
			h2 -> cnames [0] = "1 proton";
			h2 -> cnames [1] = "2 proton";
			h2 -> cnames [2] = "3 proton";
			h2 -> cnames [3] = "more protons";
			
			h2b = new mhist1D ("D13b", "", "{/Symbol n}_{/Symbol m}, E = 3 GeV, CC, on argon, no meson in FS", "momentum of most energetic proton [MeV]", "#events", 60, 0, 3000, 4, 3, 100);
			h2b -> cnames [0] = "1 proton";
			h2b -> cnames [1] = "2 proton";
			h2b -> cnames [2] = "3 proton";
			h2b -> cnames [3] = "more protons";

			h3 = new mhist1D ("D14", "", "{/Symbol n}_{/Symbol m}, E = 3 GeV, CC, on argon, no meson in FS", "cos{/Symbol q} of most energetic proton", "#events", 20, -1, 1, 4, 3, 100);
			h3 -> cnames [0] = "1 proton";
			h3 -> cnames [1] = "2 proton";
			h3 -> cnames [2] = "3 proton";
			h3 -> cnames [3] = "more protons";
			
			h4 = new hist1D ("D31", "", "{/Symbol n}_{/Symbol m}, E = 3 GeV, CC, on argon", "total visible energy", "#events", 30, 0, 3, 1, 100);
		}
		
		~D()
		{
			delete h1;
			delete h2;
			delete h2b;
			delete h3;
			delete h4;
		}
};

class E : public D
{
	protected:
	
		void set_params ();
		
	public:
	
		E ()
		{
			h1 -> name = "E1112";
			h1 -> title = "{/Symbol n}_{/Symbol m}, E = 1 GeV, CC, on argon, no meson in FS";
			h1 -> cnames [0] = "all protons (E11)";
			h1 -> cnames [1] = "only protons with Tk > 50 MeV (E12)";
			
			h2 -> name = "E13";
			h2 -> title = "{/Symbol n}_{/Symbol m}, E = 1 GeV, CC, on argon, no meson in FS";
			
			h2b -> name = "E13b";
			h2b -> title = "{/Symbol n}_{/Symbol m}, E = 1 GeV, CC, on argon, no meson in FS";
			
			h3 -> name = "E14";
			h3 -> title = "{/Symbol n}_{/Symbol m}, E = 1 GeV, CC, on argon, no meson in FS";
			
			h4 -> name = "E31";
			h4 -> title = "{/Symbol n}_{/Symbol m}, E = 1 GeV, CC, on argon";
		}
		
		~E () {}
};

class F : public pattern
{
	protected:
	
		mhist1D *h1;
		mhist1D *h2;
		mhist1D *h3;
		
		void set_params ();
		void calculate (event *e);
		void ratio ();
		
	public:
	
		F () : pattern ("NUINT12 F", 5000000)
		{
			h1 = new mhist1D ("Fnu", "", "{/Symbol n}_{/Symbol m}, E = 0.2 - 5 GeV, on CH_{2}", "neutrino energy [GeV]", "neutrino xsec [cm^2]", 50, 0, 5, 4, 6);
			h1 -> cnames [0] = "CC with MEC (F11)";
			h1 -> cnames [1] = "NC with MEC (F12)";
			h1 -> cnames [2] = "CC without MEC (F13)";
			h1 -> cnames [3] = "NC without MEC (F14)";

			h2 = new mhist1D ("Fanu", "", "anti {/Symbol n}_{/Symbol m}, E = 0.2 - 5 GeV, on CH_{2}", "neutrino energy [GeV]", "antineutrino xsec [cm^2]", 50, 0, 5, 4, 6);
			h2 -> cnames [0] = "CC with MEC (F11)";
			h2 -> cnames [1] = "NC with MEC (F12)";
			h2 -> cnames [2] = "CC without MEC (F13)";
			h2 -> cnames [3] = "NC without MEC (F14)";

			h3 = new mhist1D ("F", "", "E = 0.2 - 5 GeV, on CH_{2}", "neutrino energy [GeV]", "neutrino xsec / antineutrino xsec", 50, 0, 5, 4, 6);
			h3 -> cnames [0] = "CC with MEC (F11)";
			h3 -> cnames [1] = "NC with MEC (F12)";
			h3 -> cnames [2] = "CC without MEC (F13)";
			h3 -> cnames [3] = "NC without MEC (F14)";
		}
		
		~F ();
};

class G1 : public pattern
{
	protected:
	
		hist2D *h1;
		hist2D *h2;
		
		void set_params ();
		void calculate (event *e);
		
	public:
	
		G1 () : pattern ("NUINT G1", 5000000)
		{
			h1 = new hist2D ("G11", "", "{/Symbol n}_{/Symbol m}, E = 600 MeV, CC QE, on CH_{2}", "{/Symbol m} kinetic energy [MeV]", "cos{/Symbol q}", "d{/Symbol s} / dTk / dcos{/Symbol q}", 16, 0, 800, 20, -1, 1, 1, 100);
			h2 = new hist2D ("G11b", "", "{/Symbol n}_{/Symbol m}, E = 600 MeV, CC QE + MEC, on CH_{2}", "{/Symbol m} kinetic energy [MeV]", "cos{/Symbol q}", "d{/Symbol s} / dTk / dcos{/Symbol q}", 16, 0, 800, 20, -1, 1, 1, 100);
		}
		
		~G1 () 
		{
			delete h1;
			delete h2;
		}
};

class G1_2 : public G1
{
	protected:

		void set_params ();
		
	public:
	
		G1_2 ()
		{
			h1 -> name = "G12";
			h1 -> title = "{/Symbol n}_{/Symbol m}, E = 800 MeV, CC QE, on CH_{2}";
			
			h2 -> name = "G12b";
			h2 -> title = "{/Symbol n}_{/Symbol m}, E = 800 MeV, CC QE + MEC, on CH_{2}";		
		}
		
		~G1_2 () {}
};

class G1_3 : public G1
{
	protected:

		void set_params ();
		
	public:
	
		G1_3 ()
		{
			h1 -> name = "G13";
			h1 -> title = "{/Symbol n}_{/Symbol m}, E = 1200 MeV, CC QE, on CH_{2}";
			
			h2 -> name = "G13b";
			h2 -> title = "{/Symbol n}_{/Symbol m}, E = 1200 MeV, CC QE + MEC, on CH_{2}";			
		}
		
		~G1_3 () {}
};

class G1_4 : public G1
{
	protected:

		void set_params ();
		
	public:
	
		G1_4 ()
		{
			h1 -> name = "G14";
			h1 -> title = "{/Symbol n}_{/Symbol m}, E = 600 MeV, CC QE, on H_{2}O";
			
			h2 -> name = "G14b";
			h2 -> title = "{/Symbol n}_{/Symbol m}, E = 600 MeV, CC QE + MEC, on H_{2}O";			
		}
		
		~G1_4 () {}
};

class G1_5 : public G1
{
	protected:

		void set_params ();
		
	public:
	
		G1_5 ()
		{
			h1 -> name = "G15";
			h1 -> title = "{/Symbol n}_{/Symbol m}, E = 800 MeV, CC QE, on H_{2}O";
			
			h2 -> name = "G15b";
			h2 -> title = "{/Symbol n}_{/Symbol m}, E = 800 MeV, CC QE + MEC, on H_{2}O";			
		}
		
		~G1_5 () {}
};

class G1_6 : public G1
{
	protected:

		void set_params ();
		
	public:
	
		G1_6 ()
		{
			h1 -> name = "G16";
			h1 -> title = "{/Symbol n}_{/Symbol m}, E = 1200 MeV, CC QE, on H_{2}O";
			
			h2 -> name = "G16b";
			h2 -> title = "{/Symbol n}_{/Symbol m}, E = 1200 MeV, CC QE + MEC, on H_{2}O";			
		}
		
		~G1_6 () {}
};

class G2 : public pattern
{
	protected:
	
		mhist1D *h1;
		
		void set_params ();
		void calculate (event *e);
		
	public:
	
		G2 () : pattern ("NUINT G2", 5000000)
		{
			h1 = new mhist1D ("G2121b", "", "{/Symbol n}_{/Symbol m}, E = 800 MeV, CC, on CH_{2}", "total kinetic energy of all protons", "#events", 40, 0, 800, 2, 4, 100);
			h1 -> cnames[1] = "only QE (G21)";
			h1 -> cnames[0] = "QE + MEC (G21b)";
		}
		
		~G2 () 
		{
			delete h1;
		}
};

class G2_2 : public G2
{
	protected:
	
		void set_params ();
		
	public:
	
		G2_2 ()
		{
			h1 -> name = "G2222b";
			h1 -> title = "anti {/Symbol n}_{/Symbol m}, E = 800 MeV, CC, on CH_{2}";
			h1 -> cnames[1] = "only QE (G22)";
			h1 -> cnames[0] = "QE + MEC (G22b)";
		}
		
		~G2_2 () {}
};

class G2_3 : public G2
{
	protected:
	
		void set_params ();
		
	public:
	
		G2_3 ()
		{
			h1 -> name = "G2323b";
			h1 -> title = "{/Symbol n}_{/Symbol m}, E = 800 MeV, CC, on H_{2}O";
			h1 -> cnames[1] = "only QE (G23)";
			h1 -> cnames[0] = "QE + MEC (G23b)";
		}
		
		~G2_3 () {}
};

class G2_4 : public G2
{
	protected:
	
		void set_params ();
		
	public:
	
		G2_4 ()
		{
			h1 -> name = "G2424b";
			h1 -> title = "anti {/Symbol n}_{/Symbol m}, E = 800 MeV, CC, on H_{2}O";
			h1 -> cnames[1] = "only QE (G24)";
			h1 -> cnames[0] = "QE + MEC (G24b)";
		}
		
		~G2_4 () {}
};

class H : public pattern
{
	protected:
	
		hist1D *h1; 
		hist1D *h2; 
		hist2D *h3; 
		hist2D *h4; 
		hist2D *h5; 
		hist2D *h6; 
		hist2D *h7; 
		hist2D *h8; 
		hist1D *h9; 
		hist1D *h10; 
		hist1D *h11; 
		hist2D *h13; 
		hist2D *h19; 
		hist2D *h20;
		hist2D *h14;
		hist2D *h15; 
		hist2D *h16;
		hist2D *h17;
		hist2D *h17a;
		hist2D *h17b;
		hist2D *h17c;
		hist2D *h17d;
		hist2D *h18;
		hist2D *h18a;
		hist2D *h18b;
		hist2D *h18c;
		hist2D *h18d;
		
		void set_params ();
		void calculate (event *e);
		
	public:
	
		H () : pattern ("NUINT H", 5000000)
		{
			h1 = new hist1D ("H11", "", "{/Symbol n}_{/Symbol m}, E = 5 GeV, CC, on carbon, no pions in FS", "leadin proton kinetic energy [GeV]", "#events", 20, 0, 5, 1, 100);
			h2 = new hist1D ("H21", "", "{/Symbol n}_{/Symbol m}, E = 5 GeV, CC, on carbon, no pions in FS", "sum of all proton kinetic energy [GeV]", "#events", 20, 0, 5, 1, 100);
			h3 = new hist2D ("H31", "", "{/Symbol n}_{/Symbol m}, E = 5 GeV, CC, on carbon, {/Symbol p}^{+} in FS", "{/Symbol p}^{+} energy [GeV]", "cos{/Symbol q}", "d{/Symbol s} / dTk / dcos{/Symbol q}", 20, 0, 5, 20, -1, 1, 0);
			h4 = new hist2D ("H41", "", "{/Symbol n}_{/Symbol m}, E = 5 GeV, CC, on carbon, {/Symbol p}^{-} in FS", "{/Symbol p}^{-} energy [GeV]", "cos{/Symbol q}", "d{/Symbol s} / dTk / dcos{/Symbol q}", 20, 0, 5, 20, -1, 1, 0);
			h5 = new hist2D ("H51", "", "{/Symbol n}_{/Symbol m}, E = 5 GeV, CC, on carbon, {/Symbol p}^{0} in FS", "{/Symbol p}^{0} energy [GeV]", "cos{/Symbol q}", "d{/Symbol s} / dTk / dcos{/Symbol q}", 20, 0, 5, 20, -1, 1, 0);
			h6 = new hist2D ("H61", "", "{/Symbol n}_{/Symbol m}, E = 5 GeV, CC, on carbon, at least 2 {/Symbol p} in FS", "leading {/Symbol p}^{+} energy [GeV]", "cos{/Symbol q}", "d{/Symbol s} / dTk / dcos{/Symbol q}", 20, 0, 5, 20, -1, 1, 0);
			h7 = new hist2D ("H71", "", "{/Symbol n}_{/Symbol m}, E = 5 GeV, CC, on carbon, at least 2 {/Symbol p} in FS", "leading {/Symbol p}^{-} energy [GeV]", "cos{/Symbol q}", "d{/Symbol s} / dTk / dcos{/Symbol q}", 20, 0, 5, 20, -1, 1, 0);
			h8 = new hist2D ("H81", "", "{/Symbol n}_{/Symbol m}, E = 5 GeV, CC, on carbon, at least 2 {/Symbol p} in FS", "leading {/Symbol p}^{0} energy [GeV]", "cos{/Symbol q}", "d{/Symbol s} / dTk / dcos{/Symbol q}", 20, 0, 5, 20, -1, 1, 0);
			h9 = new hist1D ("H91", "", "{/Symbol n}_{/Symbol m}, E = 5 GeV, CC, on carbon, at least 2 {/Symbol p} in FS", "sum of all {/Symbol p}^{+} energies [GeV]", "#events", 20, 0, 5, 1, 100);
			h10 = new hist1D ("H101", "", "{/Symbol n}_{/Symbol m}, E = 5 GeV, CC, on carbon, at least 2 {/Symbol p} in FS", "sum of all {/Symbol p}^{-} energies [GeV]", "#events", 20, 0, 5, 1, 100);
			h11 = new hist1D ("H111", "", "{/Symbol n}_{/Symbol m}, E = 5 GeV, CC, on carbon, at least 2 {/Symbol p} in FS", "sum of all {/Symbol p}^{0} energies [GeV]", "#events", 20, 0, 5, 1, 100);
			h13 = new hist2D ("H131", "", "{/Symbol n}_{/Symbol m}, E = 5 GeV, CC, on carbon", "energy transfer [GeV]", "(sum of all {/Symbol p}^{0} energies) / (energy transfer)", "#events", 20, 0, 5, 20, 0, 1, 1, 100);
			h19 = new hist2D ("H191", "", "{/Symbol n}_{/Symbol m}, E = 5 GeV, CC, on carbon", "energy transfer [GeV]", "(sum of all {/Symbol p}^{+} energies) / (energy transfer)", "#events", 20, 0, 5, 20, 0, 1, 1, 100);
			h20 = new hist2D ("H201", "", "{/Symbol n}_{/Symbol m}, E = 5 GeV, CC, on carbon", "energy transfer [GeV]", "(sum of all {/Symbol p}^{-} energies) / (energy transfer)", "#events", 20, 0, 5, 20, 0, 1, 1, 100);
			h14 = new hist2D ("H141", "", "{/Symbol n}_{/Symbol m}, E = 5 GeV, CC, on carbon", "energy transfer [GeV]", "{Symbol p}^{0} multiplicity", "#events", 20, 0, 5, 10, 0, 10, 1, 100);
			h15 = new hist2D ("H151", "", "{/Symbol n}_{/Symbol m}, E = 5 GeV, CC, on carbon", "energy transfer [GeV]", "{Symbol p}^{+} multiplicity", "#events", 20, 0, 5, 10, 0, 10, 1, 100);
			h16 = new hist2D ("H161", "", "{/Symbol n}_{/Symbol m}, E = 5 GeV, CC, on carbon", "energy transfer [GeV]", "{Symbol p}^{-} multiplicity", "#events", 20, 0, 5, 10, 0, 10, 1, 100);
			h17 = new hist2D ("H171", "", "{/Symbol n}_{/Symbol m}, E = 5 GeV, CC, on carbon", "energy transfer [GeV]", "(total kinetic energy in all neutrons)/(energy transfer)", "#events", 20, 0, 5, 20, 0, 1, 1, 100);
			h17a = new hist2D ("H171_nopi", "", "{/Symbol n}_{/Symbol m}, E = 5 GeV, CC, on carbon", "energy transfer [GeV]", "(total kinetic energy in all neutrons)/(energy transfer)", "#events", 20, 0, 5, 20, 0, 1, 1, 100);
			h17b = new hist2D ("H171_1pi", "", "{/Symbol n}_{/Symbol m}, E = 5 GeV, CC, on carbon", "energy transfer [GeV]", "(total kinetic energy in all neutrons)/(energy transfer)", "#events", 20, 0, 5, 20, 0, 1, 1, 100);
			h17c = new hist2D ("H171_2pi", "", "{/Symbol n}_{/Symbol m}, E = 5 GeV, CC, on carbon", "energy transfer [GeV]", "(total kinetic energy in all neutrons)/(energy transfer)", "#events", 20, 0, 5, 20, 0, 1, 1, 100);
			h17d = new hist2D ("H171_morpi", "", "{/Symbol n}_{/Symbol m}, E = 5 GeV, CC, on carbon", "energy transfer [GeV]", "(total kinetic energy in all neutrons)/(energy transfer)", "#events", 20, 0, 5, 20, 0, 1, 1, 100);
			h18 = new hist2D ("H181", "", "{/Symbol n}_{/Symbol m}, E = 5 GeV, CC, on carbon", "energy transfer [GeV]", "(total kinetic energy in all protons)/(energy transfer)", "#events", 20, 0, 5, 20, 0, 1, 1, 100);
			h18a = new hist2D ("H181_nopi", "", "{/Symbol n}_{/Symbol m}, E = 5 GeV, CC, on carbon, no pi in FS", "energy transfer [GeV]", "(total kinetic energy in all protons)/(energy transfer)", "#events", 20, 0, 5, 20, 0, 1, 1, 100);
			h18b = new hist2D ("H181_1pi", "", "{/Symbol n}_{/Symbol m}, E = 5 GeV, CC, on carbon, 1 pi in FS", "energy transfer [GeV]", "(total kinetic energy in all protons)/(energy transfer)", "#events", 20, 0, 5, 20, 0, 1, 1, 100);
			h18c = new hist2D ("H181_2pi", "", "{/Symbol n}_{/Symbol m}, E = 5 GeV, CC, on carbon, 2 pi in FS", "energy transfer [GeV]", "(total kinetic energy in all protons)/(energy transfer)", "#events", 20, 0, 5, 20, 0, 1, 1, 100);
			h18d = new hist2D ("H181_morepi", "", "{/Symbol n}_{/Symbol m}, E = 5 GeV, CC, on carbon, more than 2 pi in FS", "energy transfer [GeV]", "(total kinetic energy in all protons)/(energy transfer)", "#events", 20, 0, 5, 20, 0, 1, 1, 100);			
		}
		
		~ H ()
		{
			delete h1;
			delete h2;
			delete h3;
			delete h4;
			delete h5;
			delete h6;
			delete h7;
			delete h8;
			delete h9;
			delete h10;
			delete h11;
			delete h13;
			delete h14;
			delete h15;
			delete h16;
			delete h17;
			delete h17a;
			delete h17b;
			delete h17c;
			delete h17d;
			delete h18;
			delete h18a;
			delete h18b;
			delete h18c;
			delete h18d;
		}
};

class antiH : public H
{
	protected:
	
		void set_params ();
		
	public:
	
		antiH ()
		{
			h1 -> title = "anti {/Symbol n}_{/Symbol m}, E = 5 GeV, CC, on carbon, no pions in FS";
			h2 -> title = "anti {/Symbol n}_{/Symbol m}, E = 5 GeV, CC, on carbon, no pions in FS";
			h3 -> title = "anti {/Symbol n}_{/Symbol m}, E = 5 GeV, CC, on carbon, {/Symbol p}^{+} in FS";
			h4 -> title = "anti {/Symbol n}_{/Symbol m}, E = 5 GeV, CC, on carbon, {/Symbol p}^{-} in FS";
			h5 -> title = "anti {/Symbol n}_{/Symbol m}, E = 5 GeV, CC, on carbon, {/Symbol p}^{0} in FS";
			h6 -> title = "anti {/Symbol n}_{/Symbol m}, E = 5 GeV, CC, on carbon, at least 2 {/Symbol p} in FS";
			h7 -> title = "anti {/Symbol n}_{/Symbol m}, E = 5 GeV, CC, on carbon, at least 2 {/Symbol p} in FS";
			h8 -> title = "anti {/Symbol n}_{/Symbol m}, E = 5 GeV, CC, on carbon, at least 2 {/Symbol p} in FS";
			h9 -> title = "anti {/Symbol n}_{/Symbol m}, E = 5 GeV, CC, on carbon, at least 2 {/Symbol p} in FS";
			h10 -> title = "anti {/Symbol n}_{/Symbol m}, E = 5 GeV, CC, on carbon, at least 2 {/Symbol p} in FS";
			h11 -> title = "anti {/Symbol n}_{/Symbol m}, E = 5 GeV, CC, on carbon, at least 2 {/Symbol p} in FS";
			h13 -> title = "anti {/Symbol n}_{/Symbol m}, E = 5 GeV, CC, on carbon";
			h14 -> title = "anti {/Symbol n}_{/Symbol m}, E = 5 GeV, CC, on carbon";
			h15 -> title = "anti {/Symbol n}_{/Symbol m}, E = 5 GeV, CC, on carbon";
			h16 -> title = "anti {/Symbol n}_{/Symbol m}, E = 5 GeV, CC, on carbon";
			h17 -> title = "anti {/Symbol n}_{/Symbol m}, E = 5 GeV, CC, on carbon";
			h17a -> title = "anti {/Symbol n}_{/Symbol m}, E = 5 GeV, CC, on carbon";
			h17b -> title = "anti {/Symbol n}_{/Symbol m}, E = 5 GeV, CC, on carbon";
			h17c -> title = "anti {/Symbol n}_{/Symbol m}, E = 5 GeV, CC, on carbon";
			h17d -> title = "anti {/Symbol n}_{/Symbol m}, E = 5 GeV, CC, on carbon";
			h18 -> title = "anti {/Symbol n}_{/Symbol m}, E = 5 GeV, CC, on carbon";
			h18a -> title = "anti {/Symbol n}_{/Symbol m}, E = 5 GeV, CC, on carbon";
			h18b -> title = "anti {/Symbol n}_{/Symbol m}, E = 5 GeV, CC, on carbon";
			h18c -> title = "anti {/Symbol n}_{/Symbol m}, E = 5 GeV, CC, on carbon";
			h18d -> title = "anti {/Symbol n}_{/Symbol m}, E = 5 GeV, CC, on carbon";
			h1 -> name = "H12";
			h2 -> name = "H22";
			h3 -> name = "H32";
			h4 -> name = "H42";
			h5 -> name = "H52";
			h6 -> name = "H62";
			h7 -> name = "H72";
			h8 -> name = "H82";
			h9 -> name = "H92";
			h10 -> name = "H102";
			h11 -> name = "H112";
			h13 -> name = "H132";
			h14 -> name = "H142";
			h15 -> name = "H152";
			h16 -> name = "H162";
			h17 -> name = "H172";
			h17a -> name = "H172_nopi";
			h17b -> name = "H172_1pi";
			h17c -> name = "H172_2pi";
			h17d -> name = "H172_morepi";
			h18 -> name = "H182";
			h18a -> name = "H182_nopi";
			h18b -> name = "H182_1pi";
			h18c -> name = "H182_2pi";
			h18d -> name = "H182_morepi";
			
		}
		
		~antiH () {}
};

class test_mec_old_cc : public pattern
{
	protected:
	
		hist1D *h1;
		hist1D *h2;
		
		void set_params ();
		void calculate (event *e);
		
	public:
	
		test_mec_old_cc () : pattern ("MEC", 5000000)
		{
			h1 = new hist1D ("mec_old_cc", "", "{/Symbol n}_{/Symbol m}, E = 0 - 2 GeV, on C", "neutrino energy [MeV]", "xsec [cm^2]", 40, 0, 2000, 6);
			h2 = new hist1D ("mec_old_cc_anti", "", "anti {/Symbol n}_{/Symbol m}, E = 0 - 2 GeV, on C", "anti neutrino energy [MeV]", "xsec [cm^2]", 40, 0, 2000, 6);
		}
		
		~test_mec_old_cc ()
		{
			delete h1;
			delete h2;
		}
};

class test_mec_new_cc : public test_mec_old_cc
{
	protected:
	
		void set_params ();
		
	public:
	
		test_mec_new_cc ()
		{
			h1 -> name = "mec_new_cc";
			h2 -> name = "mec_new_cc_anti";
		}
		
		~test_mec_new_cc () {}
};

class test2_mec_old_cc : public pattern
{
	protected:
	
		mhist1D *h1;
		mhist1D *h2;
		
		void set_params ();
		void calculate (event *e);
		
	public:
	
		test2_mec_old_cc () : pattern ("MEC2", 1000000)
		{
			h1 = new mhist1D ("mec_old_momentum", "", "{/Symbol n}_{/Symbol m}, E = 1 GeV, on C", "most energetic nucleon momentum [MeV]", "No. of events", 20, 0, 1000, 3, 3, 100);
			h1 -> cnames[0] = "p+p";
			h1 -> cnames[1] = "n+n";
			h1 -> cnames[2] = "p+n";

			h2 = new mhist1D ("mec_old_cos", "", "{/Symbol n}_{/Symbol m}, E = 1 GeV, on C", "most energetic nucleon cosine", "No. of events", 20, -1, 1, 3, 3, 100);
			h2 -> cnames[0] = "p+p";
			h2 -> cnames[1] = "n+n";
			h2 -> cnames[2] = "p+n";
		}
		
		~test2_mec_old_cc ()
		{
			delete h1;
			delete h2;
		}
};

class test2_mec_new_cc : public test2_mec_old_cc
{
	protected:
	
		void set_params ();
		
	public:
	
		test2_mec_new_cc ()
		{
			h1 -> name = "mec_new_momentum";
			h2 -> name = "mec_new_cos";
		}
		
		~test2_mec_new_cc () {}
};

class test_mec_nc : public pattern
{
	protected:
	
		mhist1D *h1;
		mhist1D *h2;
		
		void set_params ();
		void calculate (event *e);
		
	public:
	
		test_mec_nc () : pattern ("MEC", 50000000)
		{			
			h1 = new mhist1D ("mec_nc", "", "{/Symbol n}_{/Symbol m}, E = 0 - 10 GeV, on C", "neutrino energy [GeV]", "xsec [cm^2]", 150, 0, 1500, 3, 6);
			h1 -> cnames[0] = "QE + MEC";
			h1 -> cnames[1] = "only QE";
			h1 -> cnames[2] = "only MEC";
			h2 = new mhist1D ("mec_nc_anti", "", "anti {/Symbol n}_{/Symbol m}, E = 0 - 10 GeV, on C", "anti neutrino energy [GeV]", "xsec [cm^2]", 150, 0, 1500, 3, 6);
			h2 -> cnames[0] = "QE + MEC";
			h2 -> cnames[1] = "only QE";
			h2 -> cnames[2] = "only MEC";			
		}
		
		~test_mec_nc ()
		{
			delete h1;
			delete h2;
		}
};

class mb_nce : public pattern
{
	protected:
		
		string name;
		double res[51];
													
		void set_params ();
		void calculate (event *e);
		void chi2 ();
		double calc_chi (double *x);
		int t2r (double x, int s, int &w);
		int find_bin (int s, int col, double x);
		void zero();
		void save();
		void load (mhist1D& h, mhist1D& part_nomec, mhist1D& part_mec1, mhist1D& part_mec2, string Ma);
		void load (double* x, string fname);
		
	public:
	
		mb_nce ();		
		~mb_nce ();
		
		void start ();	
};

class mb_nce_he : public pattern
{
	protected:
		
		string name;
		double licznik[30];
		double mianownik[30];
						
		void set_params ();
		void calculate (event *e);
		void chi2 ();
		double calc_chi (double *x);
		int t2r (double x, int s, int &w, const double R[5][30][30]);
		int find_bin (int s, int col, double x, const double R[5][30][30]);
		void zero();
		void save();
		void load (mhist1D& h, mhist1D& nomec, mhist1D& mec, string Ma);
		void load (double* x, double *y, string fname);
		
	public:
	
		mb_nce_he ();		
		~mb_nce_he ();
		
		void start ();	
};

class mb_nce_both : public pattern
{
	protected:
		
		string name;
		
		double ma_res[51];
		double ds_licznik[30];
		double ds_mianownik[30];
		
		double true_ma_res[6][51];
		double true_ratio[6][30];

		double rec_ma_res[6][51];
		double rec_licznik[6][30];
		double rec_mianownik[6][30];
													
		void set_params ();
		void calculate (event *e);
		
		void chi2 ();
		double calc_chi (double *x);
		double calc_chi_ratio (double *x);
		
		int t2r (double x, int s, int &w);
		int find_bin (int s, int col, double x);
		int t2r (double x, int s, int &w, const double R[5][30][30]);
		int find_bin (int s, int col, double x, const double R[5][30][30]);
		
		void zero();
		void save();
		void load (mhist1D& h, mhist1D& part_nomec, mhist1D& part_mec1, mhist1D& part_mec2, string Ma, string ds);
		void load (double* x, string fname);
		void load_ratio (mhist1D& h, mhist1D& nomec, mhist1D& mec, mhist1D& mec2, string ma, string ds);
		void load (double* x, double *y, string fname);
		
		//double minimum;
		
	public:
	
		mb_nce_both ();		
		~mb_nce_both ();
		
		void start ();	
};

class phil1a : public pattern
{
	protected:
	
		hist1D *h1;
		
		void set_params ();
		void calculate (event *e);
		
	public:
		
		phil1a () : pattern ("phil1", 5000000)
		{
			h1 = new hist1D ("phil1a", "", "{/Symbol n}_{/Symbol m} + p -> {/Symbol m}^{-} + {/Symbol p}^{+} + p", "Neutrino Energy [MeV]", "Total cross section", 36, 200, 2000, 6);
		}
		
		~phil1a ()
		{
			delete h1;
		}
};

class phil1b : public phil1a
{
	protected:
	
		void set_params ();
		void calculate (event *e);
		
	public:
	
		phil1b ()
		{
			h1 -> name = "phil1b";
			h1 -> title = "{/Symbol n}_{/Symbol m} + n -> {/Symbol m}^{-} + {/Symbol p}^{+} + n";
		}
		
		~phil1b () {};
};

class phil1c : public phil1a
{
	protected:
	
		void set_params ();
		void calculate (event *e);
		
	public:
	
		phil1c ()
		{
			h1 -> name = "phil1c";
			h1 -> title = "{/Symbol n}_{/Symbol m} + n -> {/Symbol m}^{-} + {/Symbol p}^{0} + p";
		}
		
		~phil1c () {};
};

class phil1d : public phil1a
{
	protected:
	
		void set_params ();
		void calculate (event *e);
		
	public:
	
		phil1d ()
		{
			h1 -> name = "phil1d";
			h1 -> title = "{/Symbol n}_{/Symbol m} + p -> {/Symbol n}_{/Symbol m} + {/Symbol p}^{0} + p";
		}
		
		~phil1d () {};
};

class phil1e : public phil1a
{
	protected:
	
		void set_params ();
		void calculate (event *e);
		
	public:
	
		phil1e ()
		{
			h1 -> name = "phil1e";
			h1 -> title = "{/Symbol n}_{/Symbol m} + n -> {/Symbol n}_{/Symbol m} + {/Symbol p}^{0} + n";
		}
		
		~phil1e () {};
};

class phil1f : public phil1a
{
	protected:
	
		void set_params ();
		void calculate (event *e);
		
	public:
	
		phil1f ()
		{
			h1 -> name = "phil1f";
			h1 -> title = "{/Symbol n}_{/Symbol m} + n -> {/Symbol n}_{/Symbol m} + {/Symbol p}^{-} + p";
		}
		
		~phil1f () {};
};

class phil2 : public pattern
{
	protected:
	
		hist1D *h1;
		hist1D *h2;
		
		void set_params ();
		void calculate (event *e);
		
	public:
		
		phil2 () : pattern ("phil2", 5000000)
		{
			h1 = new hist1D ("phil2a", "", "{/Symbol n}_{/Symbol m} + p -> {/Symbol m}^{-} + {/Symbol p}^{+} + p", "{/Symbol p}^{+} momentum [MeV]", "d{/Symbol s} / dp_{/Symbol p} [cm^2 / MeV]", 24, 0, 600, 0);
			h2 = new hist1D ("phil2b", "", "{/Symbol n}_{/Symbol m} + p -> {/Symbol m}^{-} + {/Symbol p}^{+} + p", "Q^{2} [GeV^2]", "d{/Symbol s} / dQ^{2} [cm^2 / GeV^2]", 24, 0, 600000, 0);
		}
		
		~phil2 ()
		{
			delete h1;
			delete h2;
		}
};

class mec_cher : public pattern
{
	protected:
		
		mhist1D *h1;
		
		void set_params ();
		void calculate (event *e);
		
	public:
	
		mec_cher () : pattern ("mec_cher", 5000000)
		{
			h1 = new mhist1D ("mec_cher", "", "Kinetic energy of the most energetic nucleon", "Tk [MeV]", "No. of events", 50, 0, 1000, 2, 2, 100);
			h1 -> cnames [0] = "before FSI";
			h1 -> cnames [1] = "after FSI";
		}
		
		~mec_cher ()
		{
			delete h1;
		}
};

class mb_nce_test : public pattern
{
	protected:
	
		double licznik;
		double mianownik;
				
		void set_params ();
		void calculate (event *e);
		
	public:
	
		mb_nce_test () : pattern ("mb_nce_test", 10000000)
		{
			licznik = 0;
			mianownik = 0;
		}
		
		~mb_nce_test ()
		{
			std::cout << "Transparency = " << 100.0 * licznik / mianownik << endl;
		}
};

class mb_nce_test2 : public pattern
{
	protected:
	
		hist1D *h1;
		hist1D *h2;
					
		void set_params ();
		void calculate (event *e);
		
	public:
	
		mb_nce_test2 () : pattern ("mb_nce_test2", 1000000)
		{
			h1 = new hist1D ("nucleon_total_energy_test", "", "Total nucleons kinetic energy", "T_k [MeV]", "No. of events", 50, 0, 1000, 1, 100);
			h2 = new hist1D ("nucleon_maximum_energy_test", "", "The kinetic energy of the most energetic nucleon", "T_k [MeV]", "No. of events", 50, 0, 1000, 1, 100);
		}
		
		~mb_nce_test2 ()
		{
			delete h1;
			delete h2;
		}
};

class tem_nieves_test : public pattern
{
	protected:
	
		hist2D *h1;
					
		void set_params ();
		void calculate (event *e);
		
	public:
	
		tem_nieves_test () : pattern ("mb_nce_test2", 1000000)
		{
			h1 = new hist2D ("lepton_kin_nie", "", "", "cos{/Symbol q}", "Lepton kinetic energy [MeV]", "d{/Symbol s} / dTk / dcos{/Symbol q}", 20, -1, 1, 50, 0, 1000, 0);
		}
		
		~tem_nieves_test ()
		{
			delete h1;
		}
};

class energy_test : public pattern
{
	protected:
	
		mhist1D *h1;
		
		void set_params ();
		void calculate (event *e);
		
	public:
	
		energy_test () : pattern ("Energy test", 100000)
		{
			h1 = new mhist1D ("total_nucleon_energy_with_fsi_no_interaction", "", "{/Symbol n}_{/Symbol m}, E = 800 MeV, CC TEM, on carbon", "Total nucleon kinetic energy", "No. of events", 40, 0, 800, 2, 1, 100);
			h1 -> cnames[0] = "in out";
			h1 -> cnames[1] = "in post";
		}
		
		~energy_test ()
		{
			delete h1;
		}
};	

class niwg_nieves : public pattern
{
	protected:
	
		hist1D *h1;
		hist1D *h1norm;
		hist1D *h2;
		hist1D *h3;
		hist1D *h4;
		hist1D *h5;
		hist1D *h5norm;
		hist1D *h6norm;
		hist1D *h6;
				
		void set_params ();
		void calculate (event *e);
		
	public:
	
		niwg_nieves () : pattern ("NIWG Nieves", 5000000)
		{			
			h1 = new hist1D ("proton_above_treshold", "", "Nieves model, ^{16}O({/Symbol n}_{/Symbol m}, {/Symbol m}^{-})", "Neutrino energy [MeV]", "% of events with leading protons above the Cherenkov treshold", 10, 500, 1500, 8);
			h1norm = new hist1D ("", "", "", "", "", 10, 500, 1500, 8);
			h1norm -> saved = true;
			h1norm -> ploted = true;
			h5 = new hist1D ("proton_above_treshold2", "", "Nieves model, ^{16}O({/Symbol n}_{/Symbol m}, {/Symbol m}^{-}), Muon kinetic energy > 400 MeV", "Neutrino energy [MeV]", "% of events with leading protons above the Cherenkov treshold", 10, 500, 1500, 8);
			h5norm = new hist1D ("", "", "", "", "", 10, 500, 1500, 8);
			h5norm -> saved = true;
			h5norm -> ploted = true;
			h6 = new hist1D ("proton_above_treshold3", "", "Nieves model, ^{16}O({/Symbol n}_{/Symbol m}, {/Symbol m}^{-}), Muon kinetic energy < 400 MeV", "Neutrino energy [MeV]", "% of events with leading protons above the Cherenkov treshold", 10, 500, 1500, 8);
			h6norm = new hist1D ("", "", "", "", "", 10, 500, 1500, 8);
			h6norm -> saved = true;
			h6norm -> ploted = true;
			
			h2 = new hist1D ("leading_proton_momentum_500_600", "", "Nieves model, ^{16}O({/Symbol n}_{/Symbol m}, {/Symbol m}^{-}), E_{/Symbol n} = 500-600", "Leading proton momentum [MeV]", "No. of events (normalized to 100)", 100, 0, 2000, 1, 100);
			h3 = new hist1D ("leading_proton_momentum_900_1000", "", "Nieves model, ^{16}O({/Symbol n}_{/Symbol m}, {/Symbol m}^{-}), E_{/Symbol n} = 500-600", "Leading proton momentum [MeV]", "No. of events (normalized to 100)", 100, 0, 2000, 1, 100);
			h4 = new hist1D ("leading_proton_momentum_1400_1500", "", "Nieves model, ^{16}O({/Symbol n}_{/Symbol m}, {/Symbol m}^{-}), E_{/Symbol n} = 500-600", "Leading proton momentum [MeV]", "No. of events (normalized to 100)", 100, 0, 2000, 1, 100);
		}
		
		~niwg_nieves ()
		{
			h1 -> finalize();
			h1norm -> finalize();
			h5 -> finalize();
			h5norm -> finalize();
			h6 -> finalize();
			h6norm -> finalize();
			
			for (int i = 0; i < h1->bins_x; i++)
			{
				if (h1norm->result[i] != 0)
					h1->result[i] /= h1norm->result[i];
				if (h5norm->result[i] != 0)
					h5->result[i] /= h5norm->result[i];
				if (h6norm->result[i] != 0)
					h6->result[i] /= h6norm->result[i];
			}
			
			delete h1;
			delete h2;
			delete h3;
			delete h4;
			delete h5;
			delete h6;
		}
};	

class niwg_nieves2 : public pattern
{
	protected:
	
		hist2D *h1;
		hist2D *h2;
		
		void set_params ();
		void calculate (event *e);
		
	public:
	
		niwg_nieves2 () : pattern ("NIWG Nieves 2", 5000000)
		{
			h1 = new hist2D ("lepton_kinematics", "", "Nieves model, ^{16}O({/Symbol n}_{/Symbol m}, {/Symbol m}^{-}), E_{/Symbol n} = 1000 ", "cos{/Symbol q}", "Lepton energy [MeV]", "No. of events", 20, -1, 1, 50, 0, 1000, 0);
			h2 = new hist2D ("lepton_kinematics_above_treshold", "", "Nieves model, ^{16}O({/Symbol n}_{/Symbol m}, {/Symbol m}^{-}), E_{/Symbol n} = 1000 ", "cos{/Symbol q}", "Lepton energy [MeV]", "No. of events", 20, -1, 1, 50, 0, 1000, 0);
		}
		
		~niwg_nieves2 ()
		{
			delete h1;
			delete h2;
		}
};	

class bodek : public pattern
{
	protected:
	
		static string model;
		
		mhist1D *h1;
		hist1D *h2;
		
		void set_params ();
		void calculate (event *e);
		
	public:
	
		bodek () : pattern ("Bodek", 10000000)
		{
			h1 = new mhist1D (model, "", "{/Symbol n}_{/Symbol m}, E = 10 GeV, CC QEL, on carbon", "{/Symbol n} - Q^{2}/2M", "No. of events", 100, -0.5, 0.5, 8, 1, 100);
			h1 -> cnames[0] = "0.05-0.15";
			h1 -> cnames[1] = "0.15-0.45";
			h1 -> cnames[2] = "0.45-0.55";
			h1 -> cnames[3] = "0.55-0.85";
			h1 -> cnames[4] = "0.85-1.15";
			h1 -> cnames[5] = "1.15-1.25";
			h1 -> cnames[6] = "1.25-1.75";
			h1 -> cnames[7] = "1.75-2.25";
			h2 = new hist1D ("nucleon_momentum_" + model, "", "", "Nucleon momentum [MeV]", "Probability", 100, 0, 1000, 1, 1);
		}
		
		~bodek()
		{
			h1 -> finalize();
			
			string fname = "analysis/average_" + model + ".txt";
			
			ofstream favg(fname.c_str());
			
			favg << "#Q2 | average | RMS \n\n";
			
			double Q2[8] = {0.1, 0.3, 0.5, 0.7, 1.0, 1.2, 1.5, 2.0};
						
			for (int i = 0; i < 8; i++)
			{
				double avg = 0;
				double rms = 0;
				double norm = 0;			
	
				for (int j = 0; j < 100; j++)
				{	
					double x = (h1->begin_x + h1 -> width_x * (j + 0.5));
										
					avg += x * h1 -> result[i][j];
					rms += x * x * h1 -> result[i][j];
					
					norm += h1 -> result[i][j];
				}
				
				avg /= norm;
				rms /= norm;
				rms = sqrt(rms);
				
				favg << Q2[i] << " " << avg << " " << rms << endl;
			}
			
			favg.close();
			
			delete h1;
			delete h2;
		}
};

string bodek::model = "spectral_function";

class phd1 : public pattern
{
	protected:
	
		hist1D *h1;
		hist1D *h2;
		
		void set_params ();
		void calculate (event *e);
		
	public:
	
		phd1 () : pattern ("phd1", 10000000)
		{
			h1 = new hist1D ("proton_momentum_gfg", "", "", "Proton momentum [MeV/c]", "No. of events", 30, 0, 1500, 1, 100);
			h2 = new hist1D ("xsec_Q2_gfg", "", "", "Q^{2} [GeV^2]", "d^{2}{/Symbol s} / dQ^{2} [cm^2 / GeV^{2} / nucleon]", 100, 0, 2, 0);
		}
		
		~phd1 ()
		{
			delete h1;
			delete h2;
		}
};	

class phd2: public pattern
{
	protected:
	
		hist1D *h1;
		
		void set_params ();
		void calculate (event *e);
		
	public:
		
		phd2 () : pattern ("phd2", 1000000) //0
		{
			h1 = new hist1D ("qel_xsec_gfg_nopb", "", "", "Neutrino Energy [MeV]", "{/Synbol s}_{QEL}", 18, 200, 2000, 6);
		}
		
		~phd2 ()
		{
			delete h1;
		}
};

class phd3: public pattern
{
	protected:
	
		hist1D *h1;
		
		void set_params ();
		void calculate (event *e);
		
	public:
		
		phd3 () : pattern ("phd3", 5000000)
		{
			//h1 = new hist1D ("mec_xsec_nie", "", "", "Neutrino Energy [MeV]", "{/Synbol s}_{QEL}", 30, 0, 1500, 6);
			h1 = new hist1D ("qel_xsec", "", "", "Neutrino Energy [MeV]", "{/Synbol s}_{QEL}", 30, 0, 1500, 6);
		}
		
		~phd3 ()
		{
			delete h1;
		}
};

class phd4 : public pattern
{
	protected:
		
		mhist1D *h1;
		mhist1D *h2;
		
		void set_params ();
		void calculate (event *e);
		
	public:
	
		phd4 () : pattern ("phd4", 5000000)
		{
			h1 = new mhist1D ("mec_nucl_tem", "", "Kinetic energy of the most energetic nucleon", "Tk [MeV]", "No. of events", 20, 0, 1000, 2, 2, 100);
			h1 -> cnames [0] = "before FSI";
			h1 -> cnames [1] = "after FSI";
			
			h2 = new mhist1D ("mec_nucl2_tem", "", "Angle of the most energetic nucleon", "Tk [MeV]", "No. of events", 20, -1, 1, 2, 2, 100);
			h2 -> cnames [0] = "before FSI";
			h2 -> cnames [1] = "after FSI";
		}
		
		~phd4 ()
		{
			delete h1;
			delete h2;
		}
};

class phd5: public pattern
{
	protected:
	
		hist1D *h1;
		
		void set_params ();
		void calculate (event *e);
		
	public:
		
		phd5 () : pattern ("phd5", 5000000)
		{
			h1 = new hist1D ("coh_xsec", "", "", "Neutrino Energy [GeV]", "{/Synbol s}_{COH}", 50, 0, 200, 6);
		}
		
		~phd5 ()
		{
			delete h1;
		}
};

class phd6 : public pattern
{
	protected:
			
		void set_params ();
		void calculate (event *e);
		bool change_params();
		int count;
		
		enum {N = 8};
		
		static double Tk[N];
		static double Q2[N];
		
		double T[N];
		
	public:
	
		phd6 () : pattern ("phd6", 100000, true,N), count(0)
		{
			for (int i = 0; i < N; i++)
				T[i] = 0;
		}
		
		~phd6 ()
		{
			ofstream file("analysis/proton_trans_Fe.txt");
			for (int i = 0; i < N; i++)
				file << Q2[i] << " " << T[i]/events << endl;
			file.close();
		}
};	

double phd6::Tk[N] = {0.35, 0.7, 0.97, 1.8, 0.625, 1.718, 2.742, 3.65};
double phd6::Q2[N] = {0.6, 1.3, 1.8, 3.3, 1.04, 3.06, 5.00, 6.77};

class phd7 : public pattern
{
	protected:
			
		void set_params ();
		void calculate (event *e);
		bool change_params();
		int count;
		
		enum {N = 5};
		
		static double mom[N];
		static double Q2[N];
		
		double T[N];
		
	public:
	
		phd7 () : pattern ("phd7", 100000, true,N), count(0)
		{
			for (int i = 0; i < N; i++)
				T[i] = 0;
		}
		
		~phd7 ()
		{
			ofstream file("analysis/pion_trans_C_fz.txt");
			for (int i = 0; i < N; i++)
				file << Q2[i] << " " << T[i]/events << endl;
			file.close();
		}
};	

double phd7::mom[N] = {2.793, 3.187, 3.418, 4.077, 4.412};
double phd7::Q2[N] = {1.1, 2.15, 3.0, 3.91, 4.69};

class mec_bu: public pattern
{
	protected:
	
		hist1D *h1;
		hist1D *h2;
		hist1D *h3;
		
		void set_params ();
		void calculate (event *e);
		
	public:
		
		mec_bu () : pattern ("mec_bu", 5000000)
		{
			h1 = new hist1D ("backward_muons_all", "", "", "Neutrino Energy [GeV]", "% of events with backward muon", 40, 0, 2, 6);
			h2 = new hist1D ("backward_muons_backward", "", "", "Neutrino Energy [GeV]", "% of events with backward muon", 40, 0, 2, 6);
			h3 = new hist1D ("backward_muons", "", "", "Neutrino Energy [GeV]", "% of events with backward muon", 40, 0, 2, 6);
			h1 -> ploted = true;
			h1 -> saved = true;
			h2 -> ploted = true;
			h2 -> saved = true;
			h3 -> finalized = true;
		}
		
		~mec_bu ()
		{
			h1 -> finalize();
			h2 -> finalize();
			
			for (int i = 0; i < h3 -> bins_x; i++)
				if (h1 -> result[i] != 0)
					h3 -> result[i] = 100.0 * h2 -> result[i] / h1 -> result[i];
	
			delete h1;
			delete h2;
			delete h3;
		}
};

class mec_bu2: public pattern
{
	protected:
	
		hist1D *h1;
		hist1D *h2;
		hist1D *h3;
		
		void set_params ();
		void calculate (event *e);
		
	public:
		
		mec_bu2 () : pattern ("mec_bu", 5000000)
		{
			h1 = new hist1D ("backward_muons2_all", "", "", "Neutrino Energy [GeV]", "% of events with backward muon", 40, 0, 2, 0);
			h2 = new hist1D ("backward_muons2_backward", "", "", "Neutrino Energy [GeV]", "% of events with backward muon", 40, 0, 2, 0);
			h3 = new hist1D ("backward_muons2", "", "", "Muon kinetic energy [GeV]", "% of backward muon", 40, 0, 2, 0);
			h1 -> ploted = true;
			h1 -> saved = true;
			h2 -> ploted = true;
			h2 -> saved = true;
			h3 -> finalized = true;
		}
		
		~mec_bu2 ()
		{
			h1 -> finalize();
			h2 -> finalize();
			
			for (int i = 0; i < h3 -> bins_x; i++)
				if (h1 -> result[i] != 0)
					h3 -> result[i] = 100.0 * h2 -> result[i] / h1 -> result[i];
	
			delete h1;
			delete h2;
			delete h3;
		}
};

pattern * choose (int x);
void run_command (string com);

static double dystrR[5][51][51] = {
{
{0, 0, 0, 0.0126585, 0.00678409, 0.00213671, 0.000320215, 0.000313668, 0, 0, 0.000370749, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{0, 0, 0, 0.202595, 0.0898293, 0.0449681, 0.0124884, 0.00595971, 0.00195794, 0.000692646, 0.00111225, 0.000409718, 0, 0.00100136, 0, 0, 0, 0.000723131, 0, 0.000850234, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{0, 0.303607, 0.166667, 0.525387, 0.428689, 0.223813, 0.0896695, 0.0342059, 0.0140319, 0.00450222, 0.0048858, 0.00122915, 0.00310303, 0.00150204, 0.00158658, 0.00126159, 0.00066798, 0.00216939, 0.00232748, 0.00340093, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{0, 0.303607, 0.611111, 0.778557, 0.727773, 0.504534, 0.243486, 0.0934975, 0.0323101, 0.0161115, 0.0137838, 0.00696519, 0.00709262, 0.00300408, 0.00317316, 0.00189239, 0.00200394, 0.00506192, 0.0062066, 0.00595163, 0, 0.00105136, 0, 0, 0, 0, 0.00332291, 0, 0, 0, 0, 0.00282425, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{0, 0.721442, 0.666666, 0.892484, 0.900561, 0.706499, 0.447974, 0.196401, 0.0721846, 0.0365446, 0.0204573, 0.0122915, 0.00982421, 0.00350476, 0.00475974, 0.00315398, 0.00333989, 0.00506192, 0.0062066, 0.00595163, 0, 0.00210272, 0, 0.00387327, 0.00144665, 0, 0.00332291, 0, 0, 0, 0, 0.00282425, 0, 0.00381424, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{0, 0.860721, 0.722222, 0.949366, 0.954696, 0.848441, 0.63592, 0.347331, 0.144504, 0.0594019, 0.0323223, 0.0184372, 0.012484, 0.00601799, 0.00634632, 0.00441557, 0.00400787, 0.00578505, 0.00775825, 0.00595163, 0.000911783, 0.00210272, 0, 0.00387327, 0.00144665, 0, 0.00332291, 0, 0, 0, 0, 0.00282425, 0, 0.00381424, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{0, 1, 0.888889, 0.962025, 0.979762, 0.940795, 0.794398, 0.522065, 0.264239, 0.0989166, 0.0497475, 0.0270413, 0.0178034, 0.00902206, 0.00951949, 0.00756954, 0.00534383, 0.00578505, 0.0108616, 0.00680187, 0.00273535, 0.00315408, 0.00113038, 0.00387327, 0.00144665, 0, 0.00332291, 0, 0, 0, 0, 0.00282425, 0, 0.00381424, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{1, 1, 0.944444, 0.974683, 0.989399, 0.97159, 0.905715, 0.698135, 0.419012, 0.193559, 0.0887063, 0.0454785, 0.0275558, 0.0125268, 0.0116349, 0.00883113, 0.00534383, 0.00650818, 0.0116374, 0.00680187, 0.00364713, 0.00315408, 0.00113038, 0.00387327, 0.00144665, 0, 0.00498436, 0, 0, 0, 0, 0.00282425, 0, 0.00381424, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{1, 1, 0.944444, 0.974683, 0.99229, 0.983633, 0.961574, 0.83853, 0.596305, 0.325293, 0.158137, 0.0799757, 0.0435143, 0.0260451, 0.0158658, 0.0132467, 0.00667979, 0.00867757, 0.0124132, 0.00850233, 0.00364713, 0.00420544, 0.00226075, 0.00387327, 0.00144665, 0.00164982, 0.00498436, 0, 0, 0, 0, 0.00282425, 0, 0.00381424, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.021083},
{1, 1, 1, 0.981012, 0.997109, 0.991753, 0.983669, 0.922408, 0.756323, 0.493072, 0.273532, 0.119733, 0.0634623, 0.0390628, 0.0248564, 0.0157698, 0.00868372, 0.0101238, 0.0155165, 0.00935257, 0.00455891, 0.00420544, 0.00226075, 0.00387327, 0.00144665, 0.00164982, 0.00498436, 0.00210174, 0, 0.00421941, 0, 0.00282425, 0, 0.00381424, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0331634},
{1, 1, 1, 0.981012, 0.999036, 0.996581, 0.991354, 0.966437, 0.871654, 0.671143, 0.418834, 0.209938, 0.0984964, 0.0631405, 0.0413064, 0.0208162, 0.0146955, 0.0144626, 0.0170682, 0.011053, 0.00547069, 0.0052568, 0.00339113, 0.00387327, 0.00144665, 0.00164982, 0.00498436, 0.00210174, 0, 0.00421941, 0, 0.00282425, 0, 0.00381424, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.047619, 0, 0, 0, 0, 0.0331634},
{1, 1, 1, 0.987341, 0.999036, 0.998718, 0.995837, 0.985885, 0.938325, 0.817848, 0.589935, 0.345711, 0.16765, 0.0976874, 0.0571722, 0.0334321, 0.0193714, 0.016632, 0.0201715, 0.0127535, 0.00820604, 0.00841087, 0.00565188, 0.00387327, 0.00144665, 0.00329964, 0.00498436, 0.00210174, 0, 0.00421941, 0, 0.00282425, 0.00321798, 0.00381424, 0, 0, 0, 0, 0, 0.0138475, 0, 0, 0, 0, 0, 0.047619, 0, 0, 0, 0, 0.0331634},
{1, 1, 1, 0.993671, 0.999036, 0.998718, 0.997758, 0.993413, 0.973241, 0.915251, 0.747421, 0.507277, 0.282687, 0.150497, 0.0883749, 0.0510943, 0.031395, 0.0245864, 0.0294814, 0.0170047, 0.00911782, 0.011565, 0.00565188, 0.00516437, 0.0028933, 0.00329964, 0.00498436, 0.00420347, 0, 0.00421941, 0, 0.00282425, 0.00321798, 0.00381424, 0, 0, 0, 0, 0, 0.0138475, 0, 0, 0, 0, 0, 0.047619, 0, 0, 0, 0, 0.0331634},
{1, 1, 1, 0.993671, 0.999036, 0.999573, 0.998719, 0.995295, 0.989558, 0.959386, 0.870305, 0.665236, 0.428568, 0.229661, 0.133328, 0.0750645, 0.0414461, 0.0339871, 0.0349122, 0.0238066, 0.0100296, 0.0157704, 0.00565188, 0.00516437, 0.0028933, 0.00329964, 0.00498436, 0.00420347, 0, 0.00421941, 0, 0.00282425, 0.00321798, 0.00381424, 0, 0, 0, 0, 0, 0.0138475, 0, 0, 0, 0, 0, 0.047619, 0, 0, 0, 0, 0.0331634},
{1, 1, 1, 1, 0.999036, 1, 0.99936, 0.998432, 0.994126, 0.979819, 0.943829, 0.804384, 0.589687, 0.376402, 0.214975, 0.123005, 0.0674973, 0.0484497, 0.0481012, 0.0280866, 0.0155003, 0.0178731, 0.00678226, 0.00645546, 0.00433996, 0.00329964, 0.00498436, 0.00420347, 0, 0.00632912, 0, 0.00282425, 0.00321798, 0.00381424, 0, 0, 0, 0, 0, 0.0138475, 0, 0, 0, 0, 0, 0.047619, 0, 0, 0, 0, 0.0331634},
{1, 1, 1, 1, 0.999036, 1, 0.99936, 0.998745, 0.996737, 0.991307, 0.976171, 0.909323, 0.737326, 0.529782, 0.324978, 0.195546, 0.0948844, 0.07159, 0.058187, 0.0433907, 0.0200592, 0.0241813, 0.00932031, 0.00903766, 0.00433996, 0.00329964, 0.00498436, 0.00420347, 0, 0.00632912, 0, 0.00282425, 0.00321798, 0.00381424, 0.00417852, 0, 0, 0, 0, 0.0138475, 0, 0, 0, 0, 0, 0.047619, 0, 0, 0, 0, 0.0331634},
{1, 1, 1, 1, 1, 1, 1, 0.999686, 0.997716, 0.996502, 0.991102, 0.9606, 0.861682, 0.687086, 0.456426, 0.292058, 0.163051, 0.104131, 0.0830134, 0.0569945, 0.0300888, 0.0283867, 0.0115811, 0.0103288, 0.00433996, 0.00494947, 0.00664581, 0.00420347, 0, 0.00632912, 0, 0.00282425, 0.00321798, 0.00381424, 0.00417852, 0, 0, 0, 0, 0.0138475, 0, 0, 0, 0, 0, 0.047619, 0, 0, 0, 0, 0.0331634},
{1, 1, 1, 1, 1, 1, 1, 1, 0.999021, 0.997887, 0.997034, 0.982382, 0.934402, 0.81287, 0.622289, 0.423894, 0.257338, 0.166453, 0.123356, 0.069748, 0.0428538, 0.0378489, 0.0183633, 0.0154931, 0.00433996, 0.00659929, 0.00830727, 0.00420347, 0, 0.00632912, 0, 0.00282425, 0.00321798, 0.00381424, 0.00835704, 0, 0, 0, 0, 0.0138475, 0, 0, 0, 0, 0, 0.047619, 0, 0, 0, 0, 0.0331634},
{1, 1, 1, 1, 1, 1, 1, 1, 0.999021, 0.999273, 0.998146, 0.993445, 0.967205, 0.910195, 0.772151, 0.574898, 0.378324, 0.250423, 0.176112, 0.100356, 0.0601776, 0.0473111, 0.0251456, 0.0154931, 0.00578661, 0.00659929, 0.0132916, 0.00630521, 0, 0.00632912, 0, 0.00282425, 0.00321798, 0.00381424, 0.00835704, 0, 0, 0, 0, 0.0138475, 0, 0, 0, 0, 0, 0.047619, 0, 0, 0, 0, 0.0331634},
{1, 1, 1, 1, 1, 1, 1, 1, 0.999021, 0.999619, 0.998888, 0.997132, 0.987588, 0.959905, 0.881348, 0.729729, 0.541979, 0.363231, 0.249816, 0.135353, 0.0829722, 0.0641328, 0.0398405, 0.0193664, 0.00867991, 0.00989893, 0.0132916, 0.00630521, 0.00196609, 0.00843883, 0, 0.00847274, 0.00321798, 0.00381424, 0.00835704, 0, 0, 0, 0, 0.0138475, 0, 0, 0, 0, 0, 0.047619, 0, 0, 0, 0, 0.0331634},
{1, 1, 1, 1, 1, 1, 1, 1, 0.999347, 1, 1, 0.997951, 0.99601, 0.979472, 0.946056, 0.84795, 0.692414, 0.51589, 0.363863, 0.19748, 0.112266, 0.0799032, 0.0545353, 0.0322774, 0.0130199, 0.0148484, 0.0166145, 0.0105087, 0.00786434, 0.00843883, 0, 0.00847274, 0.00321798, 0.00381424, 0.00835704, 0.00619888, 0, 0, 0, 0.0138475, 0, 0, 0, 0, 0, 0.047619, 0, 0, 0, 0, 0.0331634},
{1, 1, 1, 1, 1, 1, 1, 1, 0.999674, 1, 1, 0.999181, 0.998227, 0.989486, 0.974086, 0.925538, 0.807306, 0.67064, 0.48647, 0.29373, 0.165197, 0.114598, 0.0658391, 0.0490616, 0.0159132, 0.0230975, 0.0199374, 0.0105087, 0.00786434, 0.014768, 0.00233941, 0.00847274, 0.00321798, 0.00381424, 0.0167141, 0.00619888, 0, 0, 0, 0.0138475, 0, 0, 0, 0, 0, 0.047619, 0, 0, 0, 0, 0.0331634},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.99867, 0.994493, 0.990481, 0.96717, 0.892348, 0.813879, 0.626361, 0.406081, 0.250132, 0.149293, 0.0952289, 0.067137, 0.0274864, 0.028047, 0.0265833, 0.0126104, 0.00786434, 0.0189874, 0.00467883, 0.011297, 0.00321798, 0.00762849, 0.0167141, 0.00619888, 0, 0, 0, 0.0138475, 0, 0, 0, 0, 0, 0.047619, 0, 0, 0, 0, 0.0331634},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.999557, 0.996996, 0.996827, 0.985492, 0.962593, 0.890531, 0.762374, 0.561083, 0.347693, 0.216981, 0.147226, 0.10587, 0.0433995, 0.0428954, 0.0348906, 0.0273225, 0.0157287, 0.0232068, 0.00467883, 0.0169455, 0.00321798, 0.0114427, 0.0167141, 0.0123978, 0, 0, 0, 0.0138475, 0, 0, 0.0238095, 0, 0, 0.047619, 0, 0, 0.0666667, 0, 0.0331634},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.998999, 0.998942, 0.996215, 0.985972, 0.949915, 0.875052, 0.692869, 0.507255, 0.309631, 0.207136, 0.144602, 0.0636527, 0.0610434, 0.0465207, 0.0378312, 0.0235931, 0.0316456, 0.0116971, 0.0197697, 0.00965391, 0.0228855, 0.0167141, 0.0123978, 0.00805171, 0, 0, 0.0138475, 0, 0, 0.0238095, 0, 0, 0.047619, 0, 0, 0.0666667, 0, 0.0331634},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.998999, 0.999471, 0.999369, 0.99666, 0.977394, 0.936382, 0.835055, 0.666186, 0.433691, 0.297566, 0.200168, 0.104159, 0.0923901, 0.0581509, 0.0483399, 0.0334235, 0.0421941, 0.0140365, 0.0282425, 0.0160898, 0.0228855, 0.0167141, 0.0123978, 0.00805171, 0, 0, 0.0138475, 0.0175439, 0.0208334, 0.0238095, 0, 0, 0.047619, 0, 0, 0.0666667, 0, 0.0331634},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.998999, 1, 0.999369, 0.998664, 0.990599, 0.973622, 0.920928, 0.792114, 0.573939, 0.395236, 0.281507, 0.149005, 0.123737, 0.0847342, 0.0588486, 0.0412879, 0.0548524, 0.0233942, 0.0367152, 0.0289618, 0.0305139, 0.0208926, 0.0123978, 0.00805171, 0, 0, 0.0138475, 0.0175439, 0.0416668, 0.0238095, 0, 0, 0.047619, 0, 0, 0.0666667, 0, 0.0331634},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.999369, 0.999332, 0.996384, 0.993793, 0.96429, 0.899704, 0.69905, 0.52285, 0.382213, 0.218444, 0.176531, 0.116302, 0.0777642, 0.0570166, 0.0675106, 0.03977, 0.0367152, 0.0386157, 0.0343282, 0.0376066, 0.0123978, 0.0241551, 0, 0, 0.0138475, 0.0175439, 0.0416668, 0.0238095, 0, 0, 0.047619, 0, 0, 0.0666667, 0, 0.0331634},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.999369, 0.999332, 1, 1, 0.985546, 0.95897, 0.819956, 0.654067, 0.480463, 0.303796, 0.244381, 0.146208, 0.111392, 0.088474, 0.086498, 0.0584854, 0.0451879, 0.0482696, 0.0381424, 0.0417852, 0.0185966, 0.0322068, 0, 0, 0.027695, 0.0350877, 0.0416668, 0.0476191, 0, 0.03125, 0.047619, 0.0500001, 0, 0.0666667, 0, 0.0452438},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.993198, 0.981764, 0.909322, 0.784061, 0.594384, 0.418423, 0.338421, 0.189406, 0.143359, 0.119931, 0.111814, 0.0818795, 0.0677819, 0.0514876, 0.0457709, 0.0459637, 0.0185966, 0.0322068, 0.00811287, 0, 0.0692376, 0.0350877, 0.0416668, 0.0714286, 0, 0.03125, 0.047619, 0.0500001, 0, 0.0666667, 0, 0.0577314},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.995749, 0.993618, 0.967408, 0.886925, 0.704127, 0.560194, 0.422562, 0.239249, 0.198004, 0.151389, 0.137131, 0.0935766, 0.0847274, 0.0643595, 0.0610279, 0.0501422, 0.0247955, 0.0402585, 0.0162257, 0, 0.0830851, 0.0350877, 0.0416668, 0.0714286, 0, 0.03125, 0.047619, 0.0500001, 0, 0.0666667, 0, 0.0698117},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.99915, 0.998176, 0.989486, 0.940378, 0.839904, 0.681713, 0.519901, 0.337275, 0.254751, 0.200834, 0.156118, 0.105274, 0.107321, 0.0804494, 0.0800991, 0.0543207, 0.0309944, 0.0402585, 0.0162257, 0.0237493, 0.0830851, 0.0526316, 0.0416668, 0.0952381, 0, 0.03125, 0.047619, 0.1, 0, 0.0666667, 0, 0.0698117},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.997897, 0.983332, 0.919952, 0.796894, 0.637039, 0.430317, 0.320768, 0.255885, 0.196202, 0.135686, 0.129915, 0.106193, 0.0915418, 0.0626777, 0.0371933, 0.056362, 0.0243386, 0.0237493, 0.0830851, 0.0701755, 0.0833334, 0.0952381, 0, 0.03125, 0.047619, 0.1, 0.04, 0.0666667, 0, 0.0818921},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.997897, 0.994636, 0.959976, 0.878254, 0.762425, 0.575077, 0.402736, 0.336494, 0.236287, 0.152062, 0.146861, 0.131937, 0.12587, 0.0919273, 0.0743864, 0.112724, 0.0324515, 0.0593734, 0.138475, 0.0701755, 0.0833334, 0.119048, 0, 0.0625001, 0.047619, 0.15, 0.04, 0.0666667, 0, 0.0939725},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.997897, 0.996897, 0.980634, 0.9566, 0.871314, 0.734797, 0.499917, 0.411314, 0.28059, 0.217565, 0.180752, 0.170553, 0.144941, 0.11282, 0.111579, 0.136879, 0.0405643, 0.0949974, 0.16617, 0.140351, 0.125, 0.166667, 0, 0.0625001, 0.047619, 0.15, 0.08, 0.0666667, 0, 0.106053},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.998949, 0.998027, 0.993545, 0.985534, 0.952155, 0.846421, 0.647038, 0.497822, 0.348101, 0.271372, 0.217467, 0.20595, 0.157212, 0.14207, 0.136375, 0.185189, 0.0811287, 0.106872, 0.22156, 0.157895, 0.145833, 0.190476, 0.0416667, 0.0625001, 0.047619, 0.15, 0.08, 0.0666667, 0, 0.118133},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.998027, 1, 0.998553, 0.980202, 0.923246, 0.807176, 0.629955, 0.43038, 0.341555, 0.265479, 0.255027, 0.180097, 0.175498, 0.167369, 0.217396, 0.121693, 0.118747, 0.304645, 0.192983, 0.145833, 0.190476, 0.0416667, 0.0625001, 0.142857, 0.25, 0.08, 0.0666667, 0, 0.154374},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.998027, 1, 1, 0.995051, 0.973089, 0.916466, 0.760137, 0.554852, 0.390682, 0.316316, 0.293642, 0.222054, 0.200569, 0.185966, 0.257655, 0.178483, 0.17812, 0.360036, 0.210526, 0.208333, 0.238095, 0.125, 0.125, 0.190476, 0.25, 0.12, 0.0666667, 0, 0.190615},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.998027, 1, 1, 0.99835, 0.989704, 0.954297, 0.858442, 0.656118, 0.465544, 0.375625, 0.332258, 0.252568, 0.229818, 0.210761, 0.297914, 0.219047, 0.225619, 0.387731, 0.280702, 0.25, 0.333333, 0.208333, 0.15625, 0.285714, 0.25, 0.12, 0.133333, 0, 0.263098},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.998027, 1, 1, 1, 0.993027, 0.987925, 0.946916, 0.789029, 0.547667, 0.417989, 0.403053, 0.286896, 0.271604, 0.241755, 0.330121, 0.300176, 0.261243, 0.429273, 0.333333, 0.3125, 0.333333, 0.208333, 0.15625, 0.380952, 0.3, 0.28, 0.266667, 0.2, 0.275178},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.998027, 1, 1, 1, 0.998011, 0.99423, 0.986237, 0.888186, 0.664638, 0.482946, 0.448105, 0.328853, 0.30921, 0.291346, 0.354276, 0.365079, 0.332491, 0.512358, 0.403509, 0.4375, 0.428572, 0.25, 0.21875, 0.428571, 0.4, 0.28, 0.4, 0.4, 0.34766},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.998027, 1, 1, 1, 0.998011, 0.996332, 0.994102, 0.947257, 0.781609, 0.576147, 0.503689, 0.382252, 0.325924, 0.323513, 0.402586, 0.397531, 0.37999, 0.540053, 0.438596, 0.479166, 0.452381, 0.25, 0.3125, 0.47619, 0.5, 0.32, 0.466667, 0.4, 0.395982},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.998027, 1, 1, 1, 0.998011, 0.996332, 1, 0.985232, 0.868993, 0.672171, 0.551959, 0.401323, 0.342638, 0.335911, 0.450896, 0.454321, 0.439363, 0.609291, 0.473684, 0.520833, 0.476191, 0.291667, 0.34375, 0.571429, 0.55, 0.32, 0.6, 0.4, 0.456384},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.998027, 1, 1, 1, 0.998011, 1, 1, 0.993671, 0.936836, 0.799478, 0.619536, 0.435652, 0.388602, 0.366905, 0.499207, 0.486772, 0.486862, 0.650834, 0.508772, 0.5625, 0.571429, 0.291667, 0.34375, 0.619048, 0.55, 0.4, 0.666667, 0.4, 0.504705},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.998027, 1, 1, 1, 0.998011, 1, 1, 0.99789, 0.971927, 0.88703, 0.703946, 0.508122, 0.435117, 0.435092, 0.579724, 0.551675, 0.53436, 0.692376, 0.561403, 0.625, 0.595238, 0.458333, 0.5, 0.619048, 0.55, 0.48, 0.666667, 0.6, 0.577187},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.998027, 1, 1, 1, 0.998011, 1, 1, 0.99789, 0.992982, 0.960461, 0.813357, 0.606321, 0.493617, 0.509478, 0.653776, 0.584127, 0.608135, 0.733919, 0.614035, 0.625, 0.666667, 0.541667, 0.53125, 0.714286, 0.65, 0.52, 0.666667, 0.7, 0.66175},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.998027, 1, 1, 1, 0.998011, 1, 1, 0.99789, 1, 0.991527, 0.874499, 0.747448, 0.585544, 0.546672, 0.710138, 0.691711, 0.655634, 0.817004, 0.68421, 0.645833, 0.690476, 0.583333, 0.53125, 0.761905, 0.65, 0.6, 0.666667, 0.7, 0.710071},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.998027, 1, 1, 1, 0.998011, 1, 1, 1, 1, 0.994352, 0.964602, 0.842804, 0.702542, 0.633456, 0.798707, 0.756614, 0.79813, 0.833829, 0.807018, 0.75, 0.761905, 0.833333, 0.65625, 0.761905, 0.75, 0.68, 0.666667, 0.7, 0.758393},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.998011, 1, 1, 1, 1, 1, 0.996782, 0.930531, 0.811184, 0.751234, 0.903379, 0.845856, 0.845629, 0.861524, 0.894737, 0.833333, 0.928571, 0.916667, 0.75, 0.857143, 0.85, 0.8, 0.733333, 0.9, 0.806714},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.998011, 1, 1, 1, 1, 1, 0.996782, 0.9733, 0.937322, 0.850416, 0.943638, 0.926984, 0.940627, 0.903067, 0.947368, 0.979167, 0.97619, 0.958333, 0.84375, 0.857143, 0.95, 0.92, 0.8, 0.9, 0.891277},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}
},
{
{1, 0, 0.00899888, 0.00390071, 0.00350122, 0.00119568, 0.000307736, 0.000147087, 0.000306809, 0, 0.000178738, 0.000193593, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{1, 0, 0.125984, 0.105475, 0.0726457, 0.0321181, 0.0140015, 0.00367719, 0.00138064, 0.00132383, 0.000714954, 0.00116156, 0.00107337, 0.000228124, 0.000499052, 0, 0.0005977, 0.000635424, 0.000355904, 0, 0.000408629, 0, 0, 0.000551879, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{1, 0.444444, 0.350956, 0.36456, 0.320676, 0.184396, 0.0845898, 0.0248869, 0.0108917, 0.00661911, 0.0037535, 0.00309748, 0.00279075, 0.00159687, 0.000748578, 0.00137438, 0.0011954, 0.00158856, 0.00177952, 0.00114441, 0.000408629, 0, 0, 0.00110376, 0, 0.000643646, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{1, 0.666667, 0.675343, 0.628801, 0.594671, 0.412089, 0.218077, 0.079529, 0.027321, 0.0120922, 0.00915268, 0.00696933, 0.00601085, 0.00296561, 0.00249526, 0.00247388, 0.00209195, 0.00158856, 0.00249132, 0.00190735, 0.000408629, 0.00136076, 0.000998408, 0.00110376, 0, 0.00128729, 0.00147378, 0, 0.00190992, 0, 0, 0.00125965, 0, 0.00151299, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{1, 1, 0.820022, 0.798216, 0.770807, 0.627472, 0.389095, 0.165928, 0.0617434, 0.0237364, 0.0164809, 0.0116156, 0.00968793, 0.00593122, 0.00374288, 0.00467289, 0.0023908, 0.00285941, 0.00391493, 0.00228882, 0.00122589, 0.00181435, 0.00199682, 0.0027594, 0.000596455, 0.00128729, 0.00147378, 0, 0.00190992, 0, 0, 0.0025193, 0, 0.00151299, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{1, 1, 0.955006, 0.901182, 0.881723, 0.786939, 0.562071, 0.294359, 0.119629, 0.0455199, 0.0261328, 0.0178105, 0.0126934, 0.0100375, 0.00548956, 0.00690555, 0.00424426, 0.00381255, 0.00533854, 0.00343322, 0.00204315, 0.00181435, 0.00199682, 0.00386316, 0.000596455, 0.00128729, 0.00294756, 0, 0.00282257, 0.00103791, 0, 0.0025193, 0, 0.00151299, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{1, 1, 0.964005, 0.94669, 0.931217, 0.888966, 0.726, 0.459135, 0.212919, 0.0880712, 0.0439542, 0.0257606, 0.0191573, 0.0130305, 0.00802972, 0.00827993, 0.00573851, 0.0050834, 0.00569445, 0.00419616, 0.00204315, 0.00181435, 0.00199682, 0.00441504, 0.000596455, 0.00128729, 0.00294756, 0.00085068, 0.00282257, 0.00103791, 0, 0.0025193, 0.00145904, 0.00151299, 0, 0, 0, 0, 0, 0.00884956, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{1, 1, 0.982002, 0.975296, 0.96846, 0.947645, 0.851397, 0.638007, 0.360555, 0.165059, 0.0765339, 0.0379569, 0.0273149, 0.0180492, 0.011024, 0.0113036, 0.00812931, 0.00635424, 0.00854167, 0.0049591, 0.00204315, 0.00226794, 0.00199682, 0.00441504, 0.00119291, 0.00193094, 0.00294756, 0.00085068, 0.00282257, 0.00103791, 0, 0.0025193, 0.00145904, 0.00151299, 0, 0, 0, 0, 0, 0.00884956, 0, 0, 0, 0, 0, 0, 0, 0.0398954, 0, 0, 0},
{1, 1, 1, 0.986998, 0.984828, 0.973954, 0.925695, 0.782453, 0.539381, 0.2914, 0.13236, 0.061615, 0.0414833, 0.0246648, 0.0160146, 0.0146021, 0.00992241, 0.0101668, 0.0103212, 0.00692457, 0.00245178, 0.00226794, 0.00199682, 0.00441504, 0.00119291, 0.00257458, 0.00294756, 0.00255204, 0.00282257, 0.00103791, 0, 0.0025193, 0.00145904, 0.00151299, 0, 0, 0, 0, 0, 0.00884956, 0, 0, 0, 0, 0, 0, 0, 0.0398954, 0, 0, 0},
{1, 1, 1, 0.996099, 0.992998, 0.988303, 0.96378, 0.877659, 0.706426, 0.441625, 0.222098, 0.106722, 0.0597586, 0.037676, 0.0230013, 0.0179006, 0.0138075, 0.0123908, 0.0113889, 0.00768751, 0.00367767, 0.00272152, 0.00299522, 0.00496691, 0.00119291, 0.00321823, 0.00294756, 0.00255204, 0.00282257, 0.00103791, 0.00122538, 0.0025193, 0.00145904, 0.00151299, 0, 0, 0, 0, 0, 0.00884956, 0, 0, 0.015625, 0, 0, 0, 0, 0.0398954, 0, 0, 0},
{1, 1, 1, 0.9987, 0.994943, 0.995127, 0.985974, 0.941682, 0.832129, 0.605791, 0.362356, 0.179936, 0.0994942, 0.0561884, 0.0319842, 0.0247725, 0.0158994, 0.0155679, 0.0128125, 0.00921339, 0.00449492, 0.00408229, 0.00449284, 0.00496691, 0.00178937, 0.00321823, 0.00442133, 0.00340272, 0.00373523, 0.00103791, 0.00122538, 0.00377896, 0.00145904, 0.00151299, 0, 0, 0, 0, 0, 0.00884956, 0, 0, 0.015625, 0, 0, 0, 0, 0.0398954, 0, 0, 0},
{1, 1, 1, 1, 0.996888, 0.998545, 0.994613, 0.970276, 0.913226, 0.757165, 0.520733, 0.297524, 0.160765, 0.0865655, 0.0490611, 0.0363173, 0.020681, 0.020969, 0.0139455, 0.0115022, 0.00572081, 0.00453588, 0.00449284, 0.00717444, 0.00178937, 0.00321823, 0.00442133, 0.00340272, 0.00373523, 0.00103791, 0.00245076, 0.00377896, 0.00145904, 0.00151299, 0, 0, 0, 0, 0, 0.00884956, 0, 0, 0.015625, 0, 0, 0, 0, 0.0398954, 0, 0, 0},
{1, 1, 1, 1, 0.998444, 0.999201, 0.99723, 0.987205, 0.957928, 0.872202, 0.694886, 0.452362, 0.252359, 0.139542, 0.0765988, 0.0525349, 0.0344398, 0.0260524, 0.0192841, 0.0156984, 0.00653807, 0.00498946, 0.00499204, 0.00772631, 0.00298228, 0.00321823, 0.00515822, 0.00340272, 0.00556054, 0.00103791, 0.00245076, 0.00377896, 0.00145904, 0.00151299, 0, 0, 0, 0, 0, 0.00884956, 0, 0, 0.015625, 0, 0, 0, 0, 0.0398954, 0, 0, 0},
{1, 1, 1, 1, 0.999611, 0.999601, 0.998769, 0.993677, 0.979283, 0.933312, 0.823153, 0.598845, 0.3897, 0.226134, 0.122028, 0.0731877, 0.047888, 0.0355838, 0.0288935, 0.0218019, 0.00858122, 0.00544305, 0.00549125, 0.00993384, 0.00298228, 0.00321823, 0.00589512, 0.00340272, 0.00556054, 0.00103791, 0.00245076, 0.00503861, 0.00145904, 0.00151299, 0, 0, 0, 0, 0, 0.00884956, 0, 0, 0.015625, 0, 0, 0, 0, 0.0398954, 0, 0, 0},
{1, 1, 1, 1, 0.999611, 0.999601, 0.999385, 0.997207, 0.992637, 0.966904, 0.911618, 0.748482, 0.552531, 0.345431, 0.188402, 0.107547, 0.067911, 0.0495465, 0.0364128, 0.026761, 0.0118503, 0.00907176, 0.00698886, 0.0126932, 0.00417519, 0.00450552, 0.00589512, 0.00340272, 0.00556054, 0.00103791, 0.00245076, 0.00503861, 0.00291808, 0.00151299, 0.00187233, 0, 0, 0, 0, 0.00884956, 0, 0, 0.015625, 0, 0, 0, 0, 0.0398954, 0, 0, 0},
{1, 1, 1, 1, 1, 0.999801, 0.999385, 0.998678, 0.997085, 0.984611, 0.959913, 0.861681, 0.707293, 0.484139, 0.294274, 0.17224, 0.100784, 0.0661361, 0.0492253, 0.033246, 0.0183883, 0.0113397, 0.00748807, 0.0143489, 0.00417519, 0.00579281, 0.00589512, 0.00340272, 0.00647319, 0.00103791, 0.00245076, 0.00629826, 0.00291808, 0.00302598, 0.00187233, 0, 0, 0, 0, 0.00884956, 0, 0, 0.015625, 0, 0, 0, 0, 0.0398954, 0, 0, 0},
{1, 1, 1, 1, 1, 0.999801, 0.999692, 0.999265, 0.998619, 0.993877, 0.98247, 0.930638, 0.829177, 0.644645, 0.427816, 0.264207, 0.155175, 0.10303, 0.0727149, 0.045453, 0.0299498, 0.0172363, 0.00998409, 0.0160045, 0.00596455, 0.00643646, 0.00663201, 0.0042534, 0.00738584, 0.00103791, 0.00245076, 0.00629826, 0.00291808, 0.00302598, 0.00187233, 0, 0, 0, 0, 0.00884956, 0, 0, 0.015625, 0, 0, 0, 0, 0.0398954, 0, 0, 0},
{1, 1, 1, 1, 1, 0.999801, 0.999692, 0.999559, 0.999233, 0.997683, 0.9923, 0.967842, 0.908837, 0.782073, 0.588993, 0.392849, 0.234989, 0.155208, 0.105127, 0.0630006, 0.0406364, 0.0222258, 0.0139777, 0.018212, 0.00835037, 0.00643646, 0.0073689, 0.00510408, 0.0082985, 0.00103791, 0.00367614, 0.00629826, 0.00291808, 0.00302598, 0.00187233, 0, 0, 0, 0, 0.00884956, 0, 0, 0.015625, 0, 0, 0, 0, 0.0398954, 0, 0, 0},
{1, 1, 1, 1, 1, 1, 0.999846, 0.999706, 0.999386, 0.998842, 0.996604, 0.985075, 0.956951, 0.888555, 0.737376, 0.539203, 0.341794, 0.228083, 0.152106, 0.0916108, 0.0537126, 0.0303904, 0.0184706, 0.0209714, 0.0125256, 0.00965468, 0.0127199, 0.00595476, 0.0082985, 0.00207581, 0.00367614, 0.00629826, 0.00437712, 0.00302598, 0.00187233, 0, 0, 0, 0, 0.00884956, 0, 0, 0.015625, 0, 0, 0, 0, 0.0398954, 0, 0, 0},
{1, 1, 1, 1, 1, 1, 0.999846, 1, 0.999386, 0.999504, 0.998749, 0.993593, 0.980182, 0.945865, 0.846759, 0.689155, 0.480221, 0.339387, 0.216168, 0.129821, 0.070875, 0.0444516, 0.0249602, 0.0220752, 0.0149114, 0.0115856, 0.0141937, 0.00680545, 0.0101238, 0.00415163, 0.0061269, 0.00881756, 0.00437712, 0.00453897, 0.00561699, 0, 0.00346042, 0, 0, 0.00884956, 0, 0, 0.015625, 0, 0, 0, 0, 0.0398954, 0, 0, 0},
{1, 1, 1, 1, 1, 1, 0.999846, 1, 0.999693, 0.999669, 0.999821, 0.99729, 0.992486, 0.974678, 0.92455, 0.824185, 0.640258, 0.472123, 0.316177, 0.192, 0.109286, 0.0630487, 0.038938, 0.0259384, 0.0208759, 0.0180221, 0.0149306, 0.00935749, 0.0128618, 0.00415163, 0.0061269, 0.0113369, 0.00437712, 0.00453897, 0.00561699, 0, 0.00346042, 0, 0.00585566, 0.00884956, 0, 0, 0.015625, 0, 0, 0, 0, 0.0398954, 0, 0, 0},
{1, 1, 1, 1, 1, 1, 0.999846, 1, 0.999693, 0.999835, 0.999821, 0.999226, 0.996995, 0.989506, 0.965788, 0.911129, 0.7678, 0.61608, 0.445726, 0.273415, 0.151784, 0.095382, 0.058407, 0.0397354, 0.0310157, 0.0238149, 0.0178781, 0.0136109, 0.0146871, 0.00726535, 0.00735228, 0.0125965, 0.00583616, 0.00605196, 0.00561699, 0, 0.00692085, 0, 0.00585566, 0.00884956, 0, 0, 0.015625, 0, 0, 0, 0, 0.0398954, 0, 0, 0.010369},
{1, 1, 1, 1, 1, 1, 0.999846, 1, 0.999693, 0.999835, 1, 0.999806, 0.997853, 0.994981, 0.98575, 0.956775, 0.867476, 0.754456, 0.588956, 0.392471, 0.224242, 0.139431, 0.0873608, 0.0584993, 0.0441377, 0.0328259, 0.0267208, 0.0195657, 0.017425, 0.011417, 0.00857766, 0.0176351, 0.00583616, 0.00756495, 0.00561699, 0, 0.00692085, 0, 0.00585566, 0.00884956, 0, 0, 0.015625, 0, 0, 0, 0, 0.0398954, 0, 0, 0.010369},
{1, 1, 1, 1, 1, 1, 0.999846, 1, 0.999847, 0.999835, 1, 1, 0.999785, 0.997491, 0.993984, 0.983496, 0.937494, 0.865168, 0.731505, 0.542852, 0.324031, 0.201572, 0.118312, 0.0849895, 0.0614348, 0.0450552, 0.0355635, 0.0238191, 0.0210757, 0.0166065, 0.0134792, 0.0176351, 0.0116723, 0.00756495, 0.00748932, 0, 0.00692085, 0.00830545, 0.00585566, 0.00884956, 0, 0, 0.015625, 0, 0, 0, 0, 0.0398954, 0, 0, 0.010369},
{1, 1, 1, 1, 1, 1, 0.999846, 1, 0.999847, 0.999835, 1, 1, 1, 0.999544, 0.997505, 0.995052, 0.971609, 0.933175, 0.837655, 0.686589, 0.466289, 0.295554, 0.176544, 0.123725, 0.0812114, 0.0624336, 0.0488274, 0.0306246, 0.025639, 0.0197202, 0.0171553, 0.0188948, 0.0160494, 0.00907795, 0.00748932, 0.00516112, 0.00692085, 0.0166109, 0.00585566, 0.0176991, 0, 0, 0.015625, 0, 0, 0, 0, 0.0398954, 0, 0, 0.010369},
{1, 1, 1, 1, 1, 1, 0.999846, 1, 0.999847, 0.999835, 1, 1, 1, 0.999544, 0.999002, 0.999175, 0.987747, 0.972359, 0.91921, 0.81133, 0.623744, 0.421291, 0.264903, 0.168979, 0.116999, 0.0863714, 0.0650735, 0.0425341, 0.0365908, 0.0238719, 0.0183807, 0.0201544, 0.0160494, 0.0121039, 0.011234, 0.0129028, 0.00692085, 0.0207636, 0.0117113, 0.0176991, 0, 0, 0.015625, 0, 0, 0, 0, 0.0398954, 0, 0, 0.010369},
{1, 1, 1, 1, 1, 1, 0.999846, 1, 1, 0.999835, 1, 1, 1, 0.999772, 0.999251, 0.99945, 0.995218, 0.988245, 0.963698, 0.902996, 0.758978, 0.560685, 0.362748, 0.241276, 0.156961, 0.115335, 0.0849695, 0.0629504, 0.0521059, 0.0352888, 0.0306345, 0.025193, 0.0175085, 0.0181559, 0.016851, 0.0154833, 0.0138417, 0.0207636, 0.017567, 0.0176991, 0, 0, 0.015625, 0, 0.038321, 0, 0, 0.0398954, 0, 0, 0.010369},
{1, 1, 1, 1, 1, 1, 0.999846, 1, 1, 0.999835, 1, 1, 1, 1, 0.999501, 0.999725, 0.998506, 0.995552, 0.984696, 0.959072, 0.866632, 0.688538, 0.491868, 0.340614, 0.221587, 0.148805, 0.111498, 0.0825161, 0.0794856, 0.0456679, 0.0379868, 0.0289719, 0.0218856, 0.0302598, 0.0243403, 0.0180639, 0.0138417, 0.0332218, 0.0292783, 0.0265487, 0, 0, 0.03125, 0, 0.038321, 0, 0, 0.062755, 0, 0, 0.0180018},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 0.999835, 1, 1, 1, 1, 0.99975, 0.999725, 1, 0.999047, 0.99395, 0.98314, 0.927518, 0.811701, 0.63572, 0.461476, 0.305091, 0.203632, 0.157922, 0.113309, 0.0986514, 0.0633124, 0.0428883, 0.0403088, 0.0277217, 0.0438767, 0.033702, 0.0283861, 0.0173021, 0.0415272, 0.035134, 0.0265487, 0, 0, 0.03125, 0, 0.038321, 0, 0, 0.062755, 0, 0, 0.0180018},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 0.999835, 1, 1, 1, 1, 0.99975, 1, 1, 0.999682, 0.996797, 0.994202, 0.967292, 0.895161, 0.757284, 0.592317, 0.412111, 0.274618, 0.214969, 0.151589, 0.12238, 0.0809568, 0.0551421, 0.0478667, 0.0393941, 0.0484157, 0.0468083, 0.0309667, 0.0207625, 0.04568, 0.0468454, 0.0265487, 0.0095954, 0, 0.03125, 0, 0.038321, 0.0222223, 0, 0.062755, 0, 0, 0.0180018},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.99975, 1, 1, 1, 0.998932, 0.996491, 0.988254, 0.950105, 0.86661, 0.70766, 0.536173, 0.352055, 0.275393, 0.194256, 0.156211, 0.104829, 0.077199, 0.0718001, 0.0452302, 0.0665716, 0.0524253, 0.0387083, 0.0380647, 0.0664436, 0.070268, 0.0353982, 0.0095954, 0.0240948, 0.03125, 0, 0.038321, 0.0444445, 0, 0.062755, 0, 0.0666666, 0.0180018},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.99975, 1, 1, 1, 0.999288, 0.99878, 0.996426, 0.980949, 0.931007, 0.828156, 0.650693, 0.460924, 0.349083, 0.263161, 0.19113, 0.135966, 0.106608, 0.0871387, 0.0554435, 0.0847275, 0.0692763, 0.0541917, 0.0519064, 0.074749, 0.070268, 0.0353982, 0.0095954, 0.0240948, 0.03125, 0.0159014, 0.0574815, 0.0666668, 0, 0.0856146, 0, 0.0666666, 0.0353078},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.99975, 1, 1, 1, 1, 0.999543, 0.997244, 0.991382, 0.965952, 0.917121, 0.765212, 0.581722, 0.440457, 0.330365, 0.233112, 0.169179, 0.12989, 0.106033, 0.0758701, 0.0938055, 0.0879996, 0.069675, 0.0622876, 0.0830545, 0.0819794, 0.0353982, 0.047977, 0.0347821, 0.0625, 0.0159014, 0.0574815, 0.0666668, 0, 0.108474, 0, 0.0666666, 0.0353078},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.99975, 1, 1, 1, 1, 0.999543, 0.998061, 0.998186, 0.989018, 0.958569, 0.874703, 0.71158, 0.556928, 0.431676, 0.281483, 0.201535, 0.154398, 0.121149, 0.0919195, 0.105909, 0.101106, 0.0774166, 0.0726689, 0.120429, 0.105402, 0.0530973, 0.0575724, 0.0347821, 0.109375, 0.0288483, 0.0574815, 0.0666668, 0, 0.131334, 0, 0.0666666, 0.0353078},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.99975, 1, 1, 1, 1, 0.999543, 0.99847, 1, 0.996006, 0.985611, 0.951687, 0.857746, 0.706591, 0.543208, 0.351065, 0.251355, 0.19361, 0.152788, 0.1211, 0.125578, 0.116084, 0.0903194, 0.100352, 0.141193, 0.105402, 0.0619469, 0.0767632, 0.0347821, 0.140625, 0.0417952, 0.0958025, 0.0888891, 0, 0.131334, 0, 0.0666666, 0.0439607},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.99975, 1, 1, 1, 1, 0.999543, 0.99847, 1, 0.999002, 0.996137, 0.984492, 0.939192, 0.839968, 0.689525, 0.454592, 0.337741, 0.251203, 0.191837, 0.157576, 0.142221, 0.13668, 0.116125, 0.124575, 0.18272, 0.140536, 0.079646, 0.12474, 0.0702005, 0.171875, 0.080636, 0.114963, 0.111111, 0, 0.154193, 0.0416667, 0.0666666, 0.0716356},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.99975, 1, 1, 1, 1, 0.999543, 0.998878, 1, 0.999501, 0.999448, 0.996421, 0.983909, 0.932161, 0.811173, 0.586015, 0.405206, 0.308795, 0.243483, 0.201348, 0.170968, 0.16102, 0.147434, 0.1661, 0.220095, 0.187962, 0.115044, 0.134336, 0.123637, 0.21875, 0.158318, 0.191605, 0.177778, 0.0857142, 0.199913, 0.0416667, 0.0666666, 0.0975944},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.999543, 0.998878, 1, 1, 0.999448, 0.999404, 0.994851, 0.981578, 0.912404, 0.723327, 0.518338, 0.388553, 0.28883, 0.249763, 0.210306, 0.19285, 0.201626, 0.204165, 0.253316, 0.246518, 0.185841, 0.182312, 0.134325, 0.265625, 0.197158, 0.210765, 0.2, 0.0857142, 0.268492, 0.125, 0.0666666, 0.143359},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.999543, 0.999287, 1, 1, 1, 1, 0.999356, 0.994842, 0.963445, 0.8584, 0.627318, 0.453498, 0.336697, 0.28361, 0.231487, 0.237786, 0.235173, 0.24569, 0.323913, 0.293364, 0.274336, 0.280345, 0.187761, 0.296875, 0.235999, 0.268247, 0.266667, 0.171428, 0.33707, 0.166667, 0.0666666, 0.186624},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.999543, 0.999287, 1, 1, 1, 1, 1, 0.998526, 0.989816, 0.944328, 0.771587, 0.547852, 0.404718, 0.330299, 0.264773, 0.278977, 0.273881, 0.290676, 0.361287, 0.35192, 0.309735, 0.309132, 0.283948, 0.3125, 0.27484, 0.386865, 0.333333, 0.285714, 0.405649, 0.208333, 0.2, 0.255847},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.999543, 0.999287, 1, 1, 1, 1, 1, 0.999263, 0.995963, 0.976271, 0.871226, 0.674066, 0.473999, 0.372612, 0.310163, 0.293956, 0.312589, 0.318359, 0.423578, 0.387054, 0.362832, 0.405085, 0.31601, 0.328125, 0.381845, 0.482667, 0.4, 0.314286, 0.451369, 0.291667, 0.333333, 0.333724},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.999543, 0.999287, 1, 1, 1, 1, 1, 1, 0.997665, 0.997262, 0.946029, 0.778423, 0.558395, 0.409088, 0.349925, 0.331403, 0.335814, 0.367331, 0.4568, 0.429227, 0.40708, 0.443467, 0.401508, 0.34375, 0.420686, 0.540149, 0.444444, 0.4, 0.565667, 0.375, 0.333333, 0.4116},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.999543, 0.999287, 1, 1, 1, 1, 1, 1, 0.997665, 0.999087, 0.98028, 0.869372, 0.660797, 0.473286, 0.401664, 0.366977, 0.356459, 0.412317, 0.506633, 0.493639, 0.486726, 0.481849, 0.454945, 0.421875, 0.498368, 0.57847, 0.533333, 0.457143, 0.565667, 0.416667, 0.533333, 0.480824},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.999543, 0.999287, 1, 1, 1, 1, 1, 1, 0.998515, 0.999087, 0.992735, 0.943277, 0.777945, 0.549156, 0.450374, 0.404423, 0.392809, 0.467333, 0.544007, 0.534629, 0.513274, 0.549016, 0.497695, 0.46875, 0.611592, 0.635951, 0.555556, 0.514286, 0.634246, 0.5, 0.533333, 0.532741},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.999543, 0.999287, 1, 1, 1, 1, 1, 1, 0.998515, 1, 0.996886, 0.979169, 0.876197, 0.661502, 0.519972, 0.477444, 0.44518, 0.501937, 0.607146, 0.604897, 0.557522, 0.568207, 0.529757, 0.484375, 0.650433, 0.655111, 0.555556, 0.657143, 0.657105, 0.5, 0.6, 0.584659},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.999543, 0.999287, 1, 1, 1, 1, 1, 1, 0.998515, 1, 0.997924, 0.995098, 0.944219, 0.797193, 0.603186, 0.537359, 0.507114, 0.561407, 0.670995, 0.688489, 0.619469, 0.635375, 0.593881, 0.578125, 0.70222, 0.693432, 0.577778, 0.714286, 0.679965, 0.5, 0.6, 0.636576},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.999543, 0.999287, 1, 1, 1, 1, 1, 1, 0.998515, 1, 1, 0.998775, 0.979489, 0.892031, 0.717071, 0.599146, 0.571628, 0.633195, 0.729133, 0.782179, 0.725664, 0.712138, 0.658005, 0.6875, 0.754008, 0.731753, 0.666667, 0.8, 0.748544, 0.583334, 0.6, 0.714453},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.999543, 0.999287, 1, 1, 1, 1, 1, 1, 0.998515, 1, 1, 0.998775, 0.99874, 0.951852, 0.832058, 0.690891, 0.646464, 0.726626, 0.78821, 0.829024, 0.778761, 0.779306, 0.732816, 0.71875, 0.779902, 0.827556, 0.711111, 0.828571, 0.771403, 0.708333, 0.6, 0.800982},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.999543, 0.999287, 1, 1, 1, 1, 1, 1, 0.998515, 1, 1, 0.998775, 1, 0.989787, 0.928889, 0.808847, 0.739364, 0.816597, 0.867113, 0.900454, 0.823009, 0.856069, 0.83969, 0.828125, 0.870531, 0.923358, 0.777778, 0.857143, 0.817123, 0.833333, 0.666667, 0.861553},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.999543, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.998541, 0.981844, 0.913873, 0.878714, 0.892727, 0.925251, 0.95901, 0.867257, 0.913641, 0.925189, 0.9375, 0.909371, 0.98084, 0.866667, 0.942857, 0.908561, 0.916667, 0.933333, 0.930776},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}
},
{
{1, 0, 0, 0, 0.0035834, 0.000904843, 0.00169688, 0, 0, 0, 0, 0.000556102, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{1, 0, 0.086927, 0.0915903, 0.0501286, 0.0435823, 0.0267304, 0.00955828, 0.00590537, 0.00240603, 0.00109436, 0.000556102, 0, 0.000652829, 0.000738226, 0.00148029, 0, 0, 0, 0.00101899, 0.0011115, 0, 0.00134535, 0, 0, 0.00164879, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{1, 0.285714, 0.347708, 0.274107, 0.290464, 0.219175, 0.134479, 0.0614569, 0.0292844, 0.013955, 0.00820769, 0.00333661, 0.00124791, 0.00195849, 0.00221468, 0.00222043, 0, 0, 0.00195836, 0.00203798, 0.0011115, 0, 0.00269069, 0, 0, 0.00329758, 0.00184831, 0, 0, 0, 0, 0, 0, 0, 0.0029933, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{1, 0.285714, 0.586758, 0.474876, 0.544328, 0.439971, 0.315733, 0.16887, 0.0846213, 0.0428328, 0.0186041, 0.0100098, 0.00436769, 0.00326415, 0.00442936, 0.00296058, 0.000840524, 0.000869638, 0.00489591, 0.00407597, 0.002223, 0, 0.00403604, 0.0015241, 0.00282359, 0.00329758, 0.00369662, 0, 0, 0.00208347, 0, 0, 0, 0, 0.0029933, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.00105532},
{1, 0.857143, 0.73888, 0.7009, 0.693038, 0.627325, 0.493523, 0.301735, 0.163841, 0.0810141, 0.0426805, 0.0177953, 0.0124791, 0.00587547, 0.00590581, 0.00444087, 0.00252158, 0.000869638, 0.00881263, 0.00713294, 0.00333451, 0, 0.00672673, 0.00304819, 0.00423539, 0.00494637, 0.00369662, 0, 0, 0.0062504, 0, 0, 0, 0, 0.0029933, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.00105532},
{1, 0.857143, 0.847538, 0.828662, 0.814993, 0.758717, 0.632875, 0.46057, 0.278297, 0.145622, 0.0705867, 0.0345669, 0.0218384, 0.0124038, 0.00812049, 0.0066613, 0.00504315, 0.00260891, 0.00881263, 0.00917093, 0.00555752, 0, 0.00672673, 0.00304819, 0.00705898, 0.00494637, 0.00739323, 0, 0, 0.0062504, 0.00245462, 0, 0, 0, 0.0029933, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.00881162, 0, 0, 0.00105532},
{1, 1, 0.934805, 0.895585, 0.892264, 0.861187, 0.770711, 0.62138, 0.413412, 0.24233, 0.121474, 0.060312, 0.0386852, 0.0195849, 0.0110734, 0.010362, 0.00924579, 0.00434819, 0.010771, 0.00917093, 0.00555752, 0.0012648, 0.00672673, 0.00304819, 0.00847078, 0.00494637, 0.00739323, 0, 0.00401118, 0.0062504, 0.00245462, 0, 0, 0, 0.00598661, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0101693, 0.00881162, 0, 0, 0.00316595},
{1, 1, 0.956537, 0.932088, 0.931682, 0.914573, 0.850071, 0.741682, 0.559343, 0.36985, 0.194382, 0.103127, 0.0642673, 0.0306969, 0.0169792, 0.015543, 0.0100863, 0.00608747, 0.0127294, 0.0112089, 0.00555752, 0.0050592, 0.00672673, 0.00487831, 0.00988258, 0.00494637, 0.00924154, 0, 0.00601678, 0.00833387, 0.00736386, 0, 0, 0, 0.00598661, 0.00332046, 0, 0.00427052, 0.0046267, 0.00473934, 0, 0, 0, 0.0066225, 0, 0, 0.0101693, 0.00881162, 0, 0, 0.00316595},
{1, 1, 1, 0.975664, 0.967749, 0.953853, 0.904899, 0.834953, 0.676201, 0.502817, 0.300635, 0.178757, 0.110467, 0.0613799, 0.0354187, 0.024423, 0.0134484, 0.00782674, 0.0156669, 0.0112089, 0.00778052, 0.006324, 0.00941742, 0.0094506, 0.00988258, 0.00989274, 0.00924154, 0.00580721, 0.00601678, 0.0125008, 0.00736386, 0, 0, 0, 0.00598661, 0.0099614, 0, 0.00427052, 0.0046267, 0.00947867, 0, 0, 0, 0.013245, 0, 0, 0.0101693, 0.00881162, 0, 0, 0.00422127},
{1, 1, 1, 0.993916, 0.985666, 0.975569, 0.946266, 0.90139, 0.799575, 0.646399, 0.420534, 0.266856, 0.165999, 0.0966328, 0.0538744, 0.0310843, 0.0176511, 0.0121749, 0.0195836, 0.0122279, 0.00778052, 0.006324, 0.0107628, 0.0094506, 0.00988258, 0.00989274, 0.00924154, 0.00580721, 0.00601678, 0.0125008, 0.00736386, 0.00278322, 0, 0, 0.0119732, 0.0099614, 0, 0.00427052, 0.0046267, 0.00947867, 0, 0, 0, 0.013245, 0, 0, 0.0101693, 0.00881162, 0, 0, 0.00527659},
{1, 1, 1, 0.993916, 0.996417, 0.985523, 0.967194, 0.939632, 0.874215, 0.769241, 0.564843, 0.38117, 0.244083, 0.146901, 0.0900475, 0.0562492, 0.0345057, 0.0165231, 0.021542, 0.0163039, 0.00889203, 0.006324, 0.0107628, 0.0094506, 0.00988258, 0.0115415, 0.00924154, 0.00580721, 0.010028, 0.0125008, 0.00736386, 0.00278322, 0, 0, 0.0119732, 0.0099614, 0, 0.00427052, 0.00925339, 0.00947867, 0, 0, 0, 0.013245, 0, 0, 0.0101693, 0.00881162, 0, 0, 0.0063319},
{1, 1, 1, 0.993916, 0.996417, 0.993666, 0.981334, 0.963779, 0.923253, 0.859841, 0.696709, 0.527701, 0.339539, 0.228805, 0.143938, 0.0867215, 0.0546783, 0.0271616, 0.0254587, 0.0193609, 0.0122265, 0.0088536, 0.0121081, 0.0094506, 0.00988258, 0.0115415, 0.00924154, 0.00774294, 0.010028, 0.0125008, 0.00736386, 0.00278322, 0, 0, 0.0119732, 0.0099614, 0, 0.00427052, 0.00925339, 0.014218, 0, 0, 0, 0.013245, 0, 0, 0.0203386, 0.00881162, 0, 0, 0.00738722},
{1, 1, 1, 1, 0.996417, 0.997285, 0.991516, 0.978871, 0.952367, 0.915754, 0.80505, 0.65859, 0.467449, 0.320983, 0.210473, 0.142299, 0.0774184, 0.0524314, 0.0352505, 0.0254748, 0.013338, 0.012648, 0.0121081, 0.0094506, 0.00988258, 0.0115415, 0.00924154, 0.00967868, 0.010028, 0.0125008, 0.00981848, 0.00556645, 0.00274504, 0, 0.0119732, 0.0099614, 0, 0.00427052, 0.00925339, 0.0189573, 0, 0, 0, 0.013245, 0, 0, 0.0203386, 0.00881162, 0, 0, 0.00738722},
{1, 1, 1, 1, 0.998208, 0.997285, 0.994909, 0.987423, 0.972313, 0.945107, 0.867429, 0.764384, 0.594981, 0.424794, 0.304277, 0.208229, 0.123772, 0.0819991, 0.0637536, 0.0336268, 0.0191444, 0.0151776, 0.0134535, 0.0094506, 0.00988258, 0.0115415, 0.0129382, 0.00967868, 0.010028, 0.0125008, 0.0122731, 0.00556645, 0.00274504, 0, 0.0119732, 0.0099614, 0.00360861, 0.00427052, 0.00925339, 0.0189573, 0, 0, 0, 0.013245, 0, 0, 0.0203386, 0.00881162, 0, 0, 0.00738722},
{1, 1, 1, 1, 0.998208, 0.999095, 0.997172, 0.993963, 0.982963, 0.97013, 0.921737, 0.84509, 0.725521, 0.541036, 0.425346, 0.2891, 0.196897, 0.122003, 0.0862748, 0.0519686, 0.033594, 0.0215016, 0.0161441, 0.015547, 0.0112944, 0.0131903, 0.0129382, 0.0116144, 0.0135282, 0.0145843, 0.0147277, 0.00556645, 0.00823513, 0, 0.0119732, 0.0099614, 0.00360861, 0.00854104, 0.00925339, 0.0189573, 0, 0, 0.00653597, 0.013245, 0, 0, 0.0203386, 0.00881162, 0, 0, 0.00844254},
{1, 1, 1, 1, 1, 1, 0.998303, 0.995472, 0.990802, 0.984085, 0.958945, 0.901993, 0.826313, 0.671705, 0.5583, 0.3942, 0.281993, 0.187423, 0.130338, 0.0794814, 0.0491583, 0.0279471, 0.0228709, 0.015547, 0.0155298, 0.0131903, 0.0147865, 0.0135502, 0.0135282, 0.0145843, 0.0147277, 0.00556645, 0.00823513, 0, 0.0119732, 0.0099614, 0.00360861, 0.00854104, 0.00925339, 0.0189573, 0, 0, 0.00653597, 0.0198675, 0, 0, 0.0203386, 0.00881162, 0, 0, 0.00949786},
{1, 1, 1, 1, 1, 1, 0.998869, 0.996479, 0.994675, 0.99182, 0.973735, 0.936471, 0.897481, 0.773639, 0.661818, 0.511883, 0.373698, 0.269169, 0.187131, 0.133488, 0.0780574, 0.0418599, 0.0309429, 0.023384, 0.0240006, 0.0148391, 0.0166348, 0.0154859, 0.0135282, 0.0187512, 0.0147277, 0.00556645, 0.0109802, 0, 0.0119732, 0.0099614, 0.00360861, 0.00854104, 0.00925339, 0.0189573, 0.00480769, 0, 0.00653597, 0.0198675, 0, 0.00819674, 0.0203386, 0.00881162, 0, 0, 0.0105532},
{1, 1, 1, 1, 1, 1, 0.998869, 0.997485, 0.995643, 0.997594, 0.98632, 0.958793, 0.937548, 0.866491, 0.770371, 0.635487, 0.498096, 0.357048, 0.266445, 0.178324, 0.130298, 0.0747446, 0.0457417, 0.0371008, 0.0296477, 0.0197855, 0.0184831, 0.0174216, 0.0215506, 0.0208347, 0.0171824, 0.00556645, 0.0109802, 0, 0.0119732, 0.0099614, 0.00360861, 0.0128116, 0.00925339, 0.0189573, 0.00480769, 0, 0.00653597, 0.0198675, 0, 0.00819674, 0.0203386, 0.00881162, 0, 0, 0.0105532},
{1, 1, 1, 1, 1, 1, 0.998869, 0.998491, 0.99758, 0.999038, 0.990698, 0.974364, 0.961315, 0.920024, 0.859059, 0.746566, 0.6158, 0.457057, 0.354571, 0.253918, 0.163643, 0.116483, 0.0645766, 0.0569141, 0.0451775, 0.0395709, 0.0258763, 0.0232288, 0.0235562, 0.0270851, 0.019637, 0.00556645, 0.0109802, 0.00295588, 0.0149665, 0.0099614, 0.00360861, 0.0128116, 0.00925339, 0.0189573, 0.00480769, 0, 0.00653597, 0.0198675, 0, 0.0163935, 0.0305079, 0.00881162, 0, 0, 0.0105532},
{1, 1, 1, 1, 1, 1, 0.998869, 0.999497, 0.99758, 0.999519, 0.995075, 0.982705, 0.980033, 0.948748, 0.919593, 0.835383, 0.72935, 0.568488, 0.45896, 0.338537, 0.248556, 0.175928, 0.0955195, 0.0919682, 0.0635308, 0.05441, 0.0314212, 0.0329075, 0.029573, 0.0291686, 0.0220916, 0.00556645, 0.0109802, 0.00295588, 0.0149665, 0.0099614, 0.00360861, 0.0128116, 0.00925339, 0.0236967, 0.00480769, 0, 0.00653597, 0.0198675, 0, 0.0163935, 0.0305079, 0.00881162, 0, 0, 0.0116085},
{1, 1, 1, 1, 1, 1, 0.998869, 1, 0.998064, 0.999519, 0.99617, 0.986598, 0.986273, 0.961354, 0.951465, 0.896815, 0.826011, 0.687628, 0.557858, 0.439417, 0.337476, 0.244228, 0.172322, 0.142636, 0.0988258, 0.0807907, 0.0425111, 0.0445219, 0.0335842, 0.031252, 0.0220916, 0.00834967, 0.0137252, 0.00295588, 0.0149665, 0.0132819, 0.00360861, 0.0170821, 0.00925339, 0.0236967, 0.00480769, 0, 0.00653597, 0.0198675, 0, 0.0163935, 0.0305079, 0.00881162, 0, 0.0115932, 0.0116085},
{1, 1, 1, 1, 1, 1, 0.999434, 1, 0.998064, 1, 0.996717, 0.988822, 0.989393, 0.969841, 0.965491, 0.928718, 0.894934, 0.785898, 0.680674, 0.54743, 0.439952, 0.32644, 0.239692, 0.191407, 0.135532, 0.118713, 0.0536009, 0.0561363, 0.0516346, 0.0395859, 0.0270008, 0.0111329, 0.0192153, 0.00886765, 0.0209531, 0.0166023, 0.00360861, 0.0170821, 0.00925339, 0.0236967, 0.00480769, 0, 0.00653597, 0.0198675, 0, 0.0163935, 0.0305079, 0.00881162, 0, 0.0115932, 0.0116085},
{1, 1, 1, 1, 1, 1, 0.999434, 1, 0.998064, 1, 0.997811, 0.992215, 0.990017, 0.982245, 0.97435, 0.956943, 0.926874, 0.861556, 0.775655, 0.67693, 0.555867, 0.446596, 0.319095, 0.266482, 0.190593, 0.173123, 0.0905671, 0.0851724, 0.0676793, 0.062504, 0.039274, 0.0111329, 0.0219603, 0.0177353, 0.0239464, 0.0199228, 0.00721721, 0.0213526, 0.0138801, 0.0236967, 0.00961539, 0, 0.0130719, 0.0198675, 0, 0.0163935, 0.0305079, 0.00881162, 0, 0.0115932, 0.0116085},
{1, 1, 1, 1, 1, 1, 0.999434, 1, 0.998064, 1, 0.999453, 0.992771, 0.991265, 0.987468, 0.979518, 0.972486, 0.952931, 0.929388, 0.860843, 0.779848, 0.657014, 0.547902, 0.426734, 0.345735, 0.261182, 0.214343, 0.133206, 0.116144, 0.0921024, 0.0770883, 0.0613656, 0.0194826, 0.0301955, 0.0236471, 0.0269397, 0.0199228, 0.0108258, 0.0213526, 0.0138801, 0.0236967, 0.00961539, 0, 0.0196079, 0.0198675, 0, 0.0163935, 0.0305079, 0.00881162, 0, 0.0115932, 0.0116085},
{1, 1, 1, 1, 1, 1, 0.999434, 1, 0.998064, 1, 0.999453, 0.993327, 0.991889, 0.98812, 0.983947, 0.980756, 0.969741, 0.957216, 0.907844, 0.860398, 0.755039, 0.655453, 0.534466, 0.449374, 0.365899, 0.28689, 0.194138, 0.170366, 0.128203, 0.10834, 0.0833769, 0.0417483, 0.0439206, 0.0295588, 0.029933, 0.0199228, 0.0108258, 0.0213526, 0.0138801, 0.0236967, 0.0144231, 0, 0.0196079, 0.0198675, 0, 0.0163935, 0.0305079, 0.00881162, 0, 0.0115932, 0.0116085},
{1, 1, 1, 1, 1, 1, 0.999434, 1, 0.998064, 1, 1, 0.993327, 0.991889, 0.989426, 0.985424, 0.985937, 0.981508, 0.976348, 0.947991, 0.90931, 0.837437, 0.760801, 0.63133, 0.551488, 0.450819, 0.370978, 0.255132, 0.241989, 0.210433, 0.143759, 0.12756, 0.0695806, 0.068626, 0.0443382, 0.0329263, 0.0365251, 0.0108258, 0.0298937, 0.0138801, 0.028436, 0.0144231, 0, 0.0196079, 0.02649, 0, 0.0163935, 0.0305079, 0.00881162, 0, 0.0115932, 0.0116085},
{1, 1, 1, 1, 1, 1, 0.999434, 1, 0.998064, 1, 1, 0.996107, 0.993137, 0.989426, 0.988376, 0.991858, 0.982349, 0.982436, 0.968554, 0.940898, 0.899837, 0.83933, 0.740303, 0.636838, 0.543998, 0.468256, 0.356789, 0.320133, 0.298679, 0.197929, 0.169288, 0.128028, 0.0933314, 0.0679853, 0.0419062, 0.0464865, 0.0252602, 0.0298937, 0.0138801, 0.0331753, 0.0240385, 0, 0.0261439, 0.02649, 0.00806451, 0.0163935, 0.0406772, 0.00881162, 0.00951076, 0.0115932, 0.0116085},
{1, 1, 1, 1, 1, 1, 0.999434, 1, 0.998064, 1, 1, 0.996107, 0.993137, 0.990079, 0.989853, 0.992599, 0.982349, 0.985045, 0.976387, 0.962297, 0.931087, 0.883598, 0.816988, 0.717614, 0.635765, 0.570481, 0.451053, 0.413081, 0.381293, 0.262517, 0.230654, 0.158644, 0.137252, 0.109424, 0.0688459, 0.0697298, 0.0396946, 0.0384348, 0.0231334, 0.0331753, 0.0240385, 0, 0.0261439, 0.02649, 0.00806451, 0.0163935, 0.0406772, 0.00881162, 0.00951076, 0.0115932, 0.0116085},
{1, 1, 1, 1, 1, 1, 0.999434, 1, 0.998064, 1, 1, 0.996107, 0.994384, 0.991385, 0.991329, 0.993339, 0.98403, 0.988523, 0.977366, 0.971468, 0.959986, 0.925336, 0.876183, 0.801868, 0.727532, 0.654774, 0.549013, 0.505996, 0.454145, 0.352106, 0.301838, 0.223327, 0.197643, 0.156718, 0.107759, 0.0996139, 0.0577377, 0.0555169, 0.0277601, 0.0379147, 0.0288462, 0, 0.0261439, 0.02649, 0.00806451, 0.0163935, 0.0406772, 0.0176232, 0.00951076, 0.0115932, 0.0137191},
{1, 1, 1, 1, 1, 1, 0.999434, 1, 0.998064, 1, 1, 0.996107, 0.995008, 0.99269, 0.992806, 0.994079, 0.985711, 0.990263, 0.982262, 0.974525, 0.969989, 0.946838, 0.916543, 0.878072, 0.789651, 0.725672, 0.63958, 0.595176, 0.536374, 0.425239, 0.380386, 0.31239, 0.269275, 0.227659, 0.158645, 0.126178, 0.0793893, 0.072599, 0.0416402, 0.056872, 0.0432692, 0.0106952, 0.0326798, 0.0397351, 0.00806451, 0.0163935, 0.0406772, 0.0352464, 0.0190215, 0.0115932, 0.0137191},
{1, 1, 1, 1, 1, 1, 0.999434, 1, 0.998064, 1, 1, 0.996107, 0.995632, 0.99269, 0.993544, 0.994819, 0.987392, 0.991132, 0.98422, 0.978601, 0.976658, 0.960791, 0.939414, 0.920747, 0.851761, 0.778433, 0.711664, 0.666798, 0.618604, 0.490081, 0.444206, 0.379762, 0.329666, 0.31338, 0.215518, 0.199228, 0.126301, 0.132386, 0.0601469, 0.0900474, 0.0673077, 0.026738, 0.0522877, 0.0463576, 0.00806451, 0.0163935, 0.0610158, 0.0352464, 0.0190215, 0.0115932, 0.0137191},
{1, 1, 1, 1, 1, 1, 0.999434, 1, 0.998064, 1, 1, 0.996663, 0.996256, 0.993343, 0.993544, 0.995559, 0.987392, 0.991132, 0.98422, 0.980639, 0.979993, 0.96585, 0.959594, 0.949705, 0.90541, 0.841087, 0.778203, 0.748099, 0.702839, 0.563002, 0.530519, 0.460818, 0.403782, 0.393189, 0.28137, 0.242394, 0.165996, 0.187903, 0.0925337, 0.146919, 0.100962, 0.0374332, 0.0849674, 0.0529801, 0.016129, 0.0163935, 0.0711851, 0.0352464, 0.0190215, 0.0115932, 0.0137191},
{1, 1, 1, 1, 1, 1, 0.999434, 1, 0.998064, 1, 1, 0.997219, 0.996256, 0.993996, 0.994282, 0.995559, 0.989073, 0.992871, 0.9852, 0.984715, 0.983327, 0.975969, 0.970357, 0.969518, 0.943528, 0.893848, 0.842894, 0.81585, 0.753312, 0.656758, 0.599249, 0.544315, 0.483388, 0.455262, 0.377156, 0.318764, 0.22415, 0.273314, 0.143427, 0.203791, 0.139423, 0.0588235, 0.137255, 0.0662251, 0.0483871, 0.0245902, 0.101693, 0.044058, 0.0285323, 0.0231864, 0.0147744},
{1, 1, 1, 1, 1, 1, 0.999434, 1, 0.998064, 1, 1, 0.997219, 0.996256, 0.995301, 0.994282, 0.995559, 0.989073, 0.993741, 0.988137, 0.985734, 0.986662, 0.977234, 0.971702, 0.978663, 0.966117, 0.916932, 0.892798, 0.868274, 0.819496, 0.721346, 0.653251, 0.611493, 0.59319, 0.520291, 0.499881, 0.411737, 0.310757, 0.371536, 0.189694, 0.279621, 0.197116, 0.117647, 0.150327, 0.0860927, 0.0725807, 0.0245902, 0.101693, 0.0881161, 0.0285323, 0.0231864, 0.0168851},
{1, 1, 1, 1, 1, 1, 1, 1, 0.998064, 1, 1, 0.997219, 0.99688, 0.996607, 0.99502, 0.997039, 0.989073, 0.993741, 0.988137, 0.988791, 0.987773, 0.979763, 0.974393, 0.987807, 0.974588, 0.943616, 0.920523, 0.916763, 0.879664, 0.785934, 0.741778, 0.68664, 0.692011, 0.611924, 0.586925, 0.491429, 0.419015, 0.465488, 0.291481, 0.36019, 0.254808, 0.197861, 0.196079, 0.152318, 0.177419, 0.0409837, 0.15254, 0.132174, 0.0475538, 0.0231864, 0.0189957},
{1, 1, 1, 1, 1, 1, 1, 1, 0.998548, 1, 1, 0.997776, 0.997504, 0.996607, 0.99502, 0.99778, 0.989914, 0.994611, 0.989229, 0.991848, 0.988885, 0.982293, 0.978429, 0.990855, 0.98447, 0.963402, 0.951944, 0.945799, 0.933815, 0.858855, 0.795779, 0.761787, 0.744167, 0.7006, 0.70067, 0.608185, 0.53449, 0.572251, 0.407148, 0.478673, 0.379808, 0.347593, 0.294118, 0.192053, 0.201613, 0.131148, 0.183048, 0.176232, 0.076086, 0.0579661, 0.0263829},
{1, 1, 1, 1, 1, 1, 1, 1, 0.998548, 1, 1, 0.997776, 0.997504, 0.996607, 0.99502, 0.99778, 0.991595, 0.99635, 0.991187, 0.992867, 0.989996, 0.984822, 0.983811, 0.990855, 0.987294, 0.971971, 0.972275, 0.9729, 0.963899, 0.898441, 0.8596, 0.803853, 0.801813, 0.762673, 0.766523, 0.687876, 0.613879, 0.666203, 0.532068, 0.57346, 0.504808, 0.486631, 0.45098, 0.278145, 0.306452, 0.254098, 0.284741, 0.299594, 0.142661, 0.082872, 0.0411573},
{1, 1, 1, 1, 1, 1, 1, 1, 0.998548, 1, 1, 0.998332, 0.997504, 0.99726, 0.997047, 0.99778, 0.992435, 0.99635, 0.992167, 0.992867, 0.991108, 0.988617, 0.986501, 0.99238, 0.990117, 0.985161, 0.974124, 0.984514, 0.977938, 0.931776, 0.908692, 0.857504, 0.853968, 0.807012, 0.805436, 0.737683, 0.711312, 0.76085, 0.662252, 0.672986, 0.653846, 0.593583, 0.581699, 0.496688, 0.427419, 0.368852, 0.427111, 0.396522, 0.294834, 0.210397, 0.0601531},
{1, 1, 1, 1, 1, 1, 1, 1, 0.998548, 1, 1, 0.998332, 0.997504, 0.99726, 0.997047, 0.99778, 0.993276, 0.99635, 0.994125, 0.993886, 0.992219, 0.991146, 0.989192, 0.993904, 0.991529, 0.988458, 0.979669, 0.988386, 0.991978, 0.956778, 0.928329, 0.885336, 0.889654, 0.854306, 0.841355, 0.78417, 0.772658, 0.829179, 0.736279, 0.772512, 0.730769, 0.700535, 0.699346, 0.609271, 0.548387, 0.45082, 0.569481, 0.519885, 0.351898, 0.361109, 0.105532},
{1, 1, 1, 1, 1, 1, 1, 1, 0.999032, 1, 1, 0.998888, 0.998128, 0.998566, 0.997047, 0.99778, 0.993276, 0.99635, 0.995104, 0.993886, 0.993331, 0.993676, 0.990537, 0.993904, 0.992941, 0.990107, 0.981517, 0.992257, 0.993983, 0.973446, 0.945511, 0.913168, 0.914359, 0.88682, 0.883261, 0.824015, 0.826787, 0.84199, 0.805679, 0.829384, 0.793269, 0.802139, 0.816993, 0.701987, 0.653226, 0.565574, 0.68221, 0.680901, 0.475538, 0.535007, 0.193123},
{1, 1, 1, 1, 1, 1, 1, 1, 0.999032, 1, 1, 0.998888, 0.998128, 0.999218, 0.997785, 0.99778, 0.993276, 0.99722, 0.996083, 0.994905, 0.993331, 0.996206, 0.990537, 0.993904, 0.992941, 0.990107, 0.981517, 0.994193, 0.997994, 0.98178, 0.967603, 0.943784, 0.93083, 0.904556, 0.898228, 0.837297, 0.841221, 0.863343, 0.838066, 0.862559, 0.836538, 0.850267, 0.869281, 0.807947, 0.790322, 0.704918, 0.745767, 0.751393, 0.637221, 0.674126, 0.289157},
{1, 1, 1, 1, 1, 1, 1, 1, 0.999032, 1, 1, 0.999444, 0.998128, 0.999218, 0.997785, 0.99778, 0.994116, 0.99722, 0.997062, 0.995924, 0.994442, 0.996206, 0.990537, 0.995428, 0.994353, 0.995054, 0.985214, 0.998064, 0.997994, 0.985947, 0.982331, 0.9577, 0.955535, 0.922291, 0.916188, 0.850579, 0.855656, 0.867614, 0.842693, 0.876777, 0.865385, 0.887701, 0.888889, 0.84106, 0.822581, 0.786885, 0.827122, 0.857133, 0.800275, 0.778465, 0.372636},
{1, 1, 1, 1, 1, 1, 1, 1, 0.999032, 1, 1, 0.999444, 0.998128, 0.999218, 0.998524, 0.99852, 0.994116, 0.99722, 0.997062, 0.995924, 0.995554, 0.99747, 0.991883, 0.995428, 0.994353, 0.995054, 0.987062, 0.998064, 0.997994, 0.992197, 0.984785, 0.977182, 0.961025, 0.934115, 0.925168, 0.880463, 0.877307, 0.884696, 0.870453, 0.900474, 0.889423, 0.903743, 0.895425, 0.874172, 0.862903, 0.836066, 0.877968, 0.910002, 0.847828, 0.824838, 0.454224},
{1, 1, 1, 1, 1, 1, 1, 1, 0.999032, 1, 1, 0.999444, 0.998752, 0.999218, 0.998524, 0.99852, 0.994957, 0.998089, 0.998042, 0.996943, 0.996665, 0.99747, 0.991883, 0.995428, 0.994353, 0.996702, 0.98891, 0.998064, 0.997994, 0.99428, 0.98724, 0.99165, 0.977495, 0.946794, 0.946121, 0.887104, 0.898959, 0.888966, 0.879706, 0.919431, 0.899038, 0.925134, 0.901961, 0.887417, 0.887097, 0.885246, 0.918646, 0.938319, 0.895382, 0.871211, 0.523875},
{1, 1, 1, 1, 1, 1, 1, 1, 0.999032, 1, 1, 0.999444, 0.998752, 0.999218, 0.998524, 0.99926, 0.997478, 0.998089, 0.998042, 0.997962, 0.996665, 0.99747, 0.993228, 0.995428, 0.995765, 0.996702, 0.98891, 1, 0.997994, 0.99428, 0.98724, 0.994434, 0.98573, 0.970441, 0.961087, 0.907027, 0.913394, 0.91886, 0.90284, 0.947867, 0.908654, 0.930481, 0.915033, 0.907285, 0.903226, 0.893443, 0.918646, 0.964754, 0.914403, 0.917584, 0.607245},
{1, 1, 1, 1, 1, 1, 1, 1, 0.999032, 1, 1, 0.999444, 0.999376, 0.999218, 0.998524, 0.99926, 0.997478, 0.998089, 0.998042, 0.997962, 0.997777, 0.99747, 0.995964, 0.996952, 0.997176, 0.998351, 0.992607, 1, 0.997994, 0.99428, 0.98724, 0.997217, 0.99122, 0.982265, 0.976054, 0.920309, 0.927828, 0.940213, 0.921346, 0.952607, 0.9375, 0.941176, 0.915033, 0.92053, 0.919355, 0.92623, 0.938984, 0.991188, 0.933425, 0.929177, 0.676896},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.999444, 0.999376, 0.999218, 0.998524, 0.99926, 0.998319, 0.998089, 0.998042, 0.997962, 0.998888, 0.99747, 0.997309, 0.998476, 0.997176, 0.998351, 0.996303, 1, 0.997994, 0.99428, 0.989695, 0.997217, 0.996711, 0.988176, 0.979047, 0.933591, 0.953088, 0.944483, 0.939853, 0.976303, 0.951923, 0.946524, 0.947712, 0.933775, 0.951613, 0.942623, 0.949154, 0.991188, 0.952446, 0.952363, 0.745492},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.999444, 0.999376, 1, 0.998524, 0.99926, 0.999159, 0.998089, 0.998042, 0.998981, 0.998888, 0.998735, 0.997309, 0.998476, 0.997176, 0.998351, 0.996303, 1, 0.997994, 0.996364, 0.992149, 0.997217, 0.996711, 0.997044, 0.994013, 0.953513, 0.971131, 0.957295, 0.953733, 0.976303, 0.956731, 0.967914, 0.954248, 0.953642, 0.967742, 0.967213, 0.969492, 1, 0.971468, 0.97555, 0.821705},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.999376, 1, 0.998524, 0.99926, 0.999159, 0.998959, 0.999021, 1, 1, 1, 1, 0.998476, 0.998588, 0.998351, 0.996303, 1, 1, 1, 0.994604, 0.997217, 0.996711, 1, 1, 0.970116, 0.985566, 0.974377, 0.976867, 0.990521, 0.966346, 0.983957, 0.980392, 0.966887, 0.975806, 0.983607, 0.989831, 1, 0.980979, 0.97555, 0.886079},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.999262, 1, 1, 1, 0.999021, 1, 1, 1, 1, 0.998476, 1, 0.998351, 0.998152, 1, 1, 1, 0.994604, 0.997217, 1, 1, 1, 0.993359, 0.996391, 0.987188, 0.981493, 0.995261, 0.990385, 0.994652, 0.993464, 0.986755, 0.983871, 1, 1, 1, 1, 1, 0.941742},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}
},
{
{1, 0, 0, 0.00236109, 0.000909672, 0.00278219, 0.000582108, 0.000763916, 0.000697952, 0, 0.000477779, 0.000127173, 0.000263188, 0.000144636, 0, 0.00016971, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.000450527, 0, 0, 0, 0, 0, 0.000853395, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{1, 0.08, 0.0579483, 0.0436802, 0.0391623, 0.0310015, 0.0238876, 0.0136395, 0.0102534, 0.00741642, 0.00504797, 0.00241629, 0.00184232, 0.001591, 0.00205094, 0.00152738, 0.00145468, 0.000389043, 0.000655063, 0.000941076, 0.000772475, 0.000575569, 0.00090913, 0.00135581, 0.000739645, 0.000819194, 0.00225264, 0, 0.000536339, 0, 0, 0, 0.000853395, 0, 0.000999368, 0, 0, 0, 0, 0, 0.00187412, 0, 0, 0, 0.00286955, 0, 0, 0, 0, 0, 0.000512889},
{1, 0.22, 0.243422, 0.223973, 0.19795, 0.163074, 0.117711, 0.0798254, 0.0548164, 0.0390174, 0.0215313, 0.0160403, 0.0142121, 0.00766574, 0.00522051, 0.00543069, 0.00418221, 0.0035014, 0.00327531, 0.0030585, 0.00411987, 0.00345341, 0.00303044, 0.00338952, 0.00369822, 0.00204799, 0.00315369, 0.00103371, 0.00160902, 0.00349775, 0.00207735, 0.00233778, 0.00170679, 0.000928049, 0.0029981, 0.00253366, 0, 0.00142811, 0, 0.00345154, 0.00187412, 0.00219898, 0, 0, 0.00286955, 0, 0, 0, 0, 0.00458715, 0.00151519},
{1, 0.4, 0.481689, 0.437404, 0.392869, 0.343375, 0.26207, 0.186348, 0.133661, 0.0968256, 0.0584514, 0.0418809, 0.0261872, 0.0188124, 0.0134243, 0.013407, 0.010001, 0.00739183, 0.00855839, 0.0068228, 0.00695228, 0.00661903, 0.00545478, 0.0064401, 0.00739645, 0.00409598, 0.00540633, 0.00206742, 0.00482707, 0.00757847, 0.00643487, 0.00467557, 0.00426697, 0.00278415, 0.00499683, 0.00715018, 0, 0.00571242, 0, 0.00345154, 0.00374824, 0.00463068, 0, 0.00540544, 0.00286955, 0, 0, 0, 0, 0.00458715, 0.00254097},
{1, 0.62, 0.659643, 0.612281, 0.550444, 0.50455, 0.417918, 0.315306, 0.233111, 0.166526, 0.110011, 0.0783841, 0.0493771, 0.0359095, 0.0249411, 0.0210724, 0.015456, 0.0130329, 0.0135805, 0.0110576, 0.0126515, 0.00949687, 0.00848522, 0.0101686, 0.010355, 0.00655356, 0.00720844, 0.0056854, 0.0053634, 0.00874438, 0.00920466, 0.00779261, 0.00597376, 0.00464025, 0.00799493, 0.00945845, 0.00129032, 0.00571242, 0, 0.00345154, 0.00374824, 0.00463068, 0.00721984, 0.00540544, 0.00573911, 0.00286046, 0, 0, 0, 0.00458715, 0.00254097},
{1, 0.78, 0.777685, 0.731221, 0.684914, 0.645622, 0.562935, 0.43723, 0.338607, 0.249204, 0.172515, 0.126017, 0.0812421, 0.0586174, 0.0406149, 0.0326566, 0.0229113, 0.0194522, 0.0203495, 0.0167041, 0.0180588, 0.0158281, 0.0118187, 0.0128802, 0.0140532, 0.00819195, 0.00946108, 0.00878653, 0.00911778, 0.0122421, 0.00920466, 0.00857187, 0.00853394, 0.00464025, 0.00799493, 0.0106126, 0.00129032, 0.00571242, 0.00151133, 0.00345154, 0.00562236, 0.00463068, 0.00721984, 0.00810816, 0.00573911, 0.00286046, 0, 0.0036888, 0, 0.0091743, 0.00356674},
{1, 0.96, 0.868669, 0.821045, 0.78792, 0.748497, 0.681722, 0.560458, 0.451281, 0.338285, 0.251299, 0.183307, 0.124376, 0.0903139, 0.0623933, 0.0493304, 0.0361853, 0.028058, 0.0260267, 0.0211742, 0.0224362, 0.0201449, 0.0178796, 0.0165794, 0.0159024, 0.0123821, 0.0130653, 0.010854, 0.0112631, 0.0169058, 0.011282, 0.0101304, 0.0102407, 0.00652764, 0.0129918, 0.014075, 0.00387096, 0.00714053, 0.00302266, 0.00345154, 0.00749648, 0.0112276, 0.0120331, 0.00810816, 0.00860866, 0.00572093, 0, 0.0036888, 0, 0.0091743, 0.00456904},
{1, 0.98, 0.909233, 0.886668, 0.869013, 0.834308, 0.772518, 0.676779, 0.571643, 0.444965, 0.339955, 0.25703, 0.178013, 0.132407, 0.0902316, 0.0695258, 0.0543688, 0.0380104, 0.0363272, 0.0299375, 0.0301609, 0.0270517, 0.020607, 0.0223416, 0.0181213, 0.0194596, 0.0148674, 0.0144719, 0.0144812, 0.0192376, 0.0141883, 0.0124682, 0.0119475, 0.00838374, 0.0159899, 0.0163832, 0.00516129, 0.00714053, 0.00604533, 0.00517731, 0.0112447, 0.0112276, 0.0120331, 0.00810816, 0.00860866, 0.00572093, 0.00347874, 0.0036888, 0, 0.0091743, 0.00508193},
{1, 0.98, 0.945915, 0.938612, 0.91859, 0.891998, 0.846996, 0.773941, 0.677617, 0.556299, 0.44301, 0.348897, 0.247717, 0.186553, 0.130987, 0.0966793, 0.0751155, 0.0525996, 0.0507386, 0.0405246, 0.0363408, 0.0333829, 0.0281831, 0.0285106, 0.0232988, 0.0260132, 0.0212564, 0.0170562, 0.0176992, 0.0192376, 0.0155732, 0.0171437, 0.0136543, 0.00838374, 0.0199873, 0.0175374, 0.00774193, 0.00714053, 0.00755666, 0.00690308, 0.0131188, 0.0112276, 0.0120331, 0.0135136, 0.00860866, 0.00914527, 0.00695749, 0.0036888, 0, 0.0183486, 0.00610771},
{1, 0.98, 0.974889, 0.962222, 0.952004, 0.930805, 0.901635, 0.848881, 0.77131, 0.666338, 0.552272, 0.440862, 0.326491, 0.255139, 0.183157, 0.133866, 0.105337, 0.0788513, 0.06646, 0.0525233, 0.0430355, 0.0397142, 0.0360622, 0.0335949, 0.0299556, 0.0301091, 0.023509, 0.0227416, 0.0214536, 0.0244842, 0.018343, 0.0187023, 0.0179213, 0.0120959, 0.0219861, 0.0198456, 0.00903225, 0.00856863, 0.00755666, 0.00862885, 0.0131188, 0.0155249, 0.0120331, 0.0135136, 0.00860866, 0.00914527, 0.00695749, 0.0036888, 0, 0.0183486, 0.00610771},
{1, 1, 0.982615, 0.97757, 0.972045, 0.96121, 0.940156, 0.901057, 0.844117, 0.75921, 0.653549, 0.537449, 0.415015, 0.333475, 0.242667, 0.179107, 0.14738, 0.105501, 0.088974, 0.0699333, 0.0559764, 0.05065, 0.0454565, 0.0434245, 0.0358728, 0.0350243, 0.0262122, 0.0268765, 0.0284261, 0.0291479, 0.024575, 0.0257156, 0.0219708, 0.0139521, 0.0239848, 0.023308, 0.0116129, 0.0128529, 0.00755666, 0.0138062, 0.0131188, 0.0198222, 0.0120331, 0.0135136, 0.0143478, 0.0120057, 0.00695749, 0.0036888, 0, 0.0183486, 0.00610771},
{1, 1, 0.98841, 0.989375, 0.98478, 0.978301, 0.966384, 0.937725, 0.898861, 0.836429, 0.74918, 0.632494, 0.518744, 0.418778, 0.320626, 0.239246, 0.201203, 0.142904, 0.113671, 0.0933105, 0.0707045, 0.0667659, 0.05576, 0.0512204, 0.0454882, 0.0436259, 0.0343217, 0.0361799, 0.0332531, 0.0349775, 0.0280373, 0.0303912, 0.024531, 0.0176643, 0.0269829, 0.0267704, 0.0141935, 0.0128529, 0.0107393, 0.0155319, 0.0131188, 0.0198222, 0.0144397, 0.0135136, 0.0143478, 0.0120057, 0.00695749, 0.0036888, 0, 0.0183486, 0.0066206},
{1, 1, 0.994205, 0.994097, 0.990597, 0.987878, 0.981664, 0.961349, 0.936934, 0.889739, 0.826974, 0.727633, 0.62039, 0.507768, 0.400036, 0.310102, 0.264928, 0.18692, 0.151446, 0.120837, 0.0925914, 0.0828818, 0.0678818, 0.0647785, 0.0528846, 0.0526371, 0.0437828, 0.0403147, 0.0386165, 0.0373093, 0.0314995, 0.0335082, 0.0313581, 0.0185923, 0.0319797, 0.0302329, 0.0193548, 0.0142811, 0.0152732, 0.0189835, 0.0149929, 0.0220212, 0.0192529, 0.0135136, 0.0143478, 0.0120057, 0.00695749, 0.0036888, 0, 0.0183486, 0.00877326},
{1, 1, 0.996137, 0.996458, 0.995755, 0.993045, 0.989377, 0.976527, 0.961167, 0.929856, 0.881924, 0.804981, 0.7115, 0.602219, 0.490257, 0.389563, 0.328763, 0.237193, 0.194781, 0.161792, 0.118376, 0.106192, 0.0821248, 0.0769808, 0.0632396, 0.0649748, 0.0523428, 0.0460001, 0.0466616, 0.041973, 0.0363467, 0.0419982, 0.0364785, 0.0241606, 0.0349778, 0.0302329, 0.0206451, 0.0185654, 0.0152732, 0.0241608, 0.0187412, 0.0258589, 0.0192529, 0.0135136, 0.0172173, 0.0120057, 0.0104362, 0.0036888, 0, 0.0183486, 0.00979903},
{1, 1, 0.996137, 0.998819, 0.997877, 0.995231, 0.994033, 0.987244, 0.97492, 0.957096, 0.92294, 0.863066, 0.790651, 0.693354, 0.583849, 0.472411, 0.413171, 0.296749, 0.247689, 0.202796, 0.153653, 0.134683, 0.104011, 0.0960352, 0.0776627, 0.0735764, 0.0627049, 0.0558204, 0.0552431, 0.0507174, 0.0425787, 0.044336, 0.0407454, 0.031585, 0.0389753, 0.0336953, 0.0219355, 0.0274172, 0.0182959, 0.0258866, 0.0187412, 0.0280579, 0.0192529, 0.0135136, 0.0172173, 0.0120057, 0.013915, 0.00737759, 0.00413223, 0.0183486, 0.0102515},
{1, 1, 0.998068, 0.99941, 0.99909, 0.997615, 0.996507, 0.992227, 0.985432, 0.975097, 0.95354, 0.913344, 0.856933, 0.776813, 0.675373, 0.567003, 0.489642, 0.375356, 0.310673, 0.256907, 0.201048, 0.17042, 0.13371, 0.117728, 0.0965812, 0.0907795, 0.0739681, 0.06409, 0.0643609, 0.0565469, 0.0508881, 0.0490115, 0.0424522, 0.0362253, 0.0419734, 0.0360035, 0.0232258, 0.0302734, 0.0198072, 0.0258866, 0.0187412, 0.0302569, 0.0192529, 0.0135136, 0.0200869, 0.0120057, 0.0173937, 0.0110664, 0.00826446, 0.0183486, 0.012303},
{1, 1, 1, 1, 0.999697, 0.999006, 0.998108, 0.996047, 0.991974, 0.985584, 0.974169, 0.946183, 0.908275, 0.841081, 0.757279, 0.649459, 0.577002, 0.458939, 0.388073, 0.322114, 0.254864, 0.21708, 0.17159, 0.141455, 0.122469, 0.110031, 0.0933408, 0.0708092, 0.072406, 0.0664573, 0.0571201, 0.0544663, 0.0501328, 0.0429128, 0.0469702, 0.0394659, 0.0283871, 0.0317015, 0.0213186, 0.0276123, 0.0206153, 0.0346548, 0.0192529, 0.0162163, 0.0229564, 0.0177267, 0.0173937, 0.0110664, 0.00826446, 0.0183486, 0.0133288},
{1, 1, 1, 1, 0.999697, 0.999603, 0.99869, 0.99783, 0.995231, 0.991935, 0.985422, 0.96997, 0.939501, 0.894549, 0.831488, 0.734269, 0.664742, 0.54074, 0.463428, 0.393635, 0.323986, 0.264852, 0.2189, 0.182129, 0.158045, 0.127234, 0.113615, 0.0936196, 0.0809874, 0.0775334, 0.0626597, 0.0622589, 0.0578133, 0.047553, 0.0509677, 0.0463907, 0.0309677, 0.0359859, 0.0258525, 0.0310639, 0.0206153, 0.0368538, 0.0192529, 0.0216218, 0.0292612, 0.0177267, 0.0208725, 0.0110664, 0.0165289, 0.0183486, 0.0133288},
{1, 1, 1, 1, 0.999697, 0.999801, 0.999127, 0.998593, 0.996859, 0.995971, 0.992117, 0.981677, 0.963104, 0.93189, 0.886156, 0.809974, 0.740898, 0.631582, 0.554516, 0.478805, 0.401343, 0.333345, 0.267993, 0.232294, 0.196913, 0.157586, 0.137042, 0.11171, 0.102492, 0.0903585, 0.0778936, 0.0723893, 0.0603735, 0.0596177, 0.0579633, 0.0510073, 0.036129, 0.0402702, 0.0258525, 0.0310639, 0.0299859, 0.0434507, 0.0192529, 0.0243245, 0.0292612, 0.0205871, 0.0319945, 0.0110664, 0.0289256, 0.0229358, 0.0157626},
{1, 1, 1, 1, 1, 0.999801, 0.999272, 0.999357, 0.998139, 0.997352, 0.995461, 0.989323, 0.977204, 0.961482, 0.924349, 0.872427, 0.821146, 0.721174, 0.641069, 0.566366, 0.483741, 0.397521, 0.328299, 0.279738, 0.241662, 0.190353, 0.170831, 0.138586, 0.130382, 0.112051, 0.0958973, 0.0840782, 0.0672006, 0.0698263, 0.0699557, 0.0637027, 0.0399999, 0.0445545, 0.0318979, 0.0362412, 0.0356082, 0.0434507, 0.0240661, 0.0270272, 0.0321308, 0.0205871, 0.0319945, 0.0110664, 0.0330579, 0.0229358, 0.0178141},
{1, 1, 1, 1, 1, 0.999801, 0.999563, 0.999739, 0.998837, 0.998043, 0.997372, 0.993774, 0.986679, 0.976814, 0.950088, 0.916132, 0.878263, 0.796902, 0.715844, 0.649995, 0.563892, 0.464287, 0.400423, 0.344709, 0.298702, 0.24816, 0.209599, 0.172182, 0.157735, 0.134786, 0.119441, 0.100443, 0.082731, 0.083747, 0.0819481, 0.0694734, 0.0425806, 0.0574074, 0.0364319, 0.0414185, 0.0431047, 0.0478487, 0.031286, 0.0270272, 0.0321308, 0.0205871, 0.0354732, 0.0147552, 0.0371901, 0.0229358, 0.0203786},
{1, 1, 1, 1, 1, 0.999801, 0.999709, 1, 0.999069, 0.998734, 0.998686, 0.996191, 0.991825, 0.985811, 0.971281, 0.947071, 0.920813, 0.856542, 0.785588, 0.733582, 0.643595, 0.542398, 0.484139, 0.412886, 0.375674, 0.316153, 0.258707, 0.214047, 0.18777, 0.164517, 0.145754, 0.127717, 0.105087, 0.103236, 0.0939405, 0.0844771, 0.0490322, 0.0674042, 0.0500338, 0.0500474, 0.0487271, 0.0522466, 0.0409124, 0.0378381, 0.0350003, 0.0205871, 0.0424307, 0.0221328, 0.0371901, 0.0275229, 0.0219172},
{1, 1, 1, 1, 1, 0.999801, 0.999709, 1, 0.999651, 0.999079, 0.999403, 0.997717, 0.995263, 0.991172, 0.983948, 0.968182, 0.953749, 0.905581, 0.855343, 0.806432, 0.723179, 0.630248, 0.56967, 0.487172, 0.441164, 0.379672, 0.315473, 0.258601, 0.23175, 0.204158, 0.181888, 0.159667, 0.134956, 0.121797, 0.11093, 0.0983267, 0.0645161, 0.078829, 0.0651471, 0.0638535, 0.0524753, 0.0588436, 0.0457257, 0.0432435, 0.0378699, 0.0205871, 0.0424307, 0.0221328, 0.0371901, 0.0275229, 0.0249946},
{1, 1, 1, 1, 1, 1, 0.999709, 1, 0.999767, 0.999194, 0.999761, 0.99848, 0.997237, 0.994215, 0.991521, 0.983316, 0.971205, 0.94509, 0.905783, 0.865014, 0.79797, 0.706896, 0.654015, 0.576904, 0.519936, 0.453809, 0.375159, 0.321657, 0.281756, 0.248566, 0.213841, 0.191616, 0.160558, 0.146855, 0.132196, 0.111022, 0.0864515, 0.091682, 0.0787491, 0.0690309, 0.0640641, 0.0654405, 0.0553521, 0.0539307, 0.0407394, 0.0291685, 0.0493882, 0.0258216, 0.0413223, 0.0321101, 0.027559},
{1, 1, 1, 1, 1, 1, 0.999854, 1, 0.999767, 0.999424, 0.999881, 0.998607, 0.998289, 0.996673, 0.995307, 0.99282, 0.985235, 0.968865, 0.94334, 0.914993, 0.858794, 0.787537, 0.735534, 0.65302, 0.594732, 0.537776, 0.448596, 0.395607, 0.345682, 0.294036, 0.261682, 0.239252, 0.190427, 0.173768, 0.159178, 0.131797, 0.104516, 0.111675, 0.0847944, 0.0862886, 0.0715606, 0.0764354, 0.062572, 0.0539307, 0.0493481, 0.0348894, 0.0493882, 0.0331992, 0.0454546, 0.0366972, 0.0306536},
{1, 1, 1, 1, 1, 1, 0.999854, 1, 0.999767, 0.999424, 0.999881, 0.998861, 0.998816, 0.99783, 0.997358, 0.995214, 0.991781, 0.98326, 0.966323, 0.946284, 0.904705, 0.852723, 0.80887, 0.736158, 0.67814, 0.616009, 0.524319, 0.475455, 0.418285, 0.367489, 0.324737, 0.290198, 0.232243, 0.214602, 0.185162, 0.164112, 0.128012, 0.134525, 0.10293, 0.100095, 0.0828053, 0.0940272, 0.0842316, 0.0674443, 0.0579568, 0.0406104, 0.0633032, 0.0451124, 0.0495868, 0.0412844, 0.038448},
{1, 1, 1, 1, 1, 1, 0.999854, 1, 0.999767, 0.999424, 0.999881, 0.998988, 0.998947, 0.998264, 0.998305, 0.996402, 0.995964, 0.991636, 0.979861, 0.967735, 0.938694, 0.902293, 0.860691, 0.800559, 0.751853, 0.68982, 0.606766, 0.542129, 0.4864, 0.439257, 0.397066, 0.342435, 0.288567, 0.260077, 0.234435, 0.207969, 0.152727, 0.166225, 0.121066, 0.124466, 0.107169, 0.116017, 0.103485, 0.0863633, 0.069435, 0.0491918, 0.066782, 0.0561788, 0.053719, 0.0458715, 0.0413947},
{1, 1, 1, 1, 1, 1, 0.999854, 1, 0.999884, 0.999424, 1, 0.998988, 0.998947, 0.998988, 0.998463, 0.99759, 0.997418, 0.995137, 0.989276, 0.98138, 0.962432, 0.9383, 0.905068, 0.858519, 0.813983, 0.754614, 0.684128, 0.611508, 0.557842, 0.518203, 0.457491, 0.415686, 0.352793, 0.319708, 0.288401, 0.256443, 0.191436, 0.206212, 0.155827, 0.150353, 0.127784, 0.135808, 0.134771, 0.0944714, 0.0780436, 0.0577731, 0.0946119, 0.0709339, 0.0743802, 0.0504587, 0.0460107},
{1, 1, 1, 1, 1, 1, 0.999854, 1, 0.999884, 0.999424, 1, 0.998988, 0.99921, 0.999132, 0.998936, 0.998269, 0.998146, 0.997471, 0.993425, 0.991262, 0.980241, 0.965927, 0.934463, 0.897528, 0.866867, 0.808681, 0.754951, 0.684878, 0.626537, 0.58466, 0.534424, 0.483667, 0.416798, 0.403232, 0.354359, 0.304918, 0.235307, 0.259052, 0.190588, 0.176239, 0.165266, 0.157798, 0.154024, 0.0998769, 0.0866523, 0.0806568, 0.101569, 0.0820003, 0.0826446, 0.0504587, 0.0501138},
{1, 1, 1, 1, 1, 1, 0.999854, 1, 1, 0.999655, 1, 0.999116, 0.999342, 0.999421, 0.998936, 0.998269, 0.998509, 0.998249, 0.99607, 0.994555, 0.986936, 0.979165, 0.956319, 0.932248, 0.904335, 0.849742, 0.81407, 0.750579, 0.702161, 0.653449, 0.608579, 0.555359, 0.491043, 0.457987, 0.407325, 0.363779, 0.280468, 0.309329, 0.229882, 0.210755, 0.206497, 0.192981, 0.175683, 0.135012, 0.109609, 0.112122, 0.112006, 0.0967555, 0.0950413, 0.0733944, 0.053704},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 0.99977, 1, 0.999116, 0.999342, 0.999421, 0.998936, 0.998439, 0.998691, 0.998638, 0.997161, 0.995732, 0.992086, 0.985784, 0.974199, 0.952585, 0.937249, 0.890835, 0.859179, 0.80248, 0.764377, 0.709996, 0.679916, 0.627178, 0.556754, 0.528817, 0.473284, 0.415715, 0.334662, 0.369309, 0.285802, 0.243544, 0.236483, 0.23916, 0.194936, 0.164742, 0.149782, 0.137866, 0.136357, 0.115199, 0.107438, 0.0871559, 0.0634662},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 0.99977, 1, 0.999116, 0.999342, 0.999421, 0.998936, 0.998439, 0.998691, 0.998833, 0.99738, 0.996438, 0.99363, 0.989813, 0.981472, 0.969194, 0.95648, 0.917458, 0.895299, 0.85003, 0.813814, 0.763186, 0.736004, 0.693415, 0.637827, 0.604917, 0.546237, 0.47928, 0.421113, 0.435406, 0.358346, 0.288872, 0.281462, 0.283139, 0.238255, 0.194472, 0.181348, 0.166471, 0.157229, 0.133643, 0.123967, 0.105504, 0.0680822},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 0.99977, 1, 0.999243, 0.999342, 0.999421, 0.999094, 0.998609, 0.998873, 0.998833, 0.997816, 0.996908, 0.994403, 0.992979, 0.986017, 0.98004, 0.969424, 0.940396, 0.925484, 0.890345, 0.862084, 0.801661, 0.785339, 0.74095, 0.703538, 0.661528, 0.609197, 0.54662, 0.493371, 0.508239, 0.429378, 0.352725, 0.345182, 0.335915, 0.286387, 0.240418, 0.22726, 0.200796, 0.181581, 0.152087, 0.152893, 0.119266, 0.0768013},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 0.99977, 1, 0.999243, 0.999474, 0.999421, 0.999094, 0.998609, 0.998873, 0.998833, 0.998253, 0.997379, 0.994403, 0.994418, 0.990866, 0.989106, 0.980888, 0.959801, 0.948912, 0.92239, 0.899092, 0.84962, 0.834502, 0.797958, 0.759204, 0.725564, 0.68415, 0.615096, 0.566068, 0.559651, 0.494365, 0.42003, 0.427643, 0.38869, 0.34174, 0.283662, 0.281782, 0.249424, 0.223325, 0.188975, 0.177686, 0.146789, 0.0902381},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 0.99977, 1, 0.999243, 0.999474, 0.999421, 0.999409, 0.998778, 0.998873, 0.998833, 0.998472, 0.997614, 0.994403, 0.994994, 0.991169, 0.993173, 0.988654, 0.973318, 0.967383, 0.949267, 0.924836, 0.89412, 0.87545, 0.845493, 0.810571, 0.782175, 0.742995, 0.686652, 0.624132, 0.642481, 0.571443, 0.506319, 0.50823, 0.45246, 0.401905, 0.343121, 0.350652, 0.306633, 0.310178, 0.218486, 0.227273, 0.201835, 0.109269},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 0.99977, 1, 0.999243, 0.999474, 0.999566, 0.999684, 0.998778, 0.999055, 0.998833, 0.999127, 0.997849, 0.995175, 0.995857, 0.992381, 0.994868, 0.991243, 0.983148, 0.977745, 0.967357, 0.948435, 0.926766, 0.90384, 0.885236, 0.853241, 0.825793, 0.800958, 0.766582, 0.706973, 0.713886, 0.654567, 0.616768, 0.594439, 0.573958, 0.512609, 0.432311, 0.422799, 0.358122, 0.372796, 0.299639, 0.301653, 0.256881, 0.137191},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 0.99977, 1, 0.999243, 0.999474, 0.999566, 0.999684, 0.999118, 0.999055, 0.999027, 0.999127, 0.998084, 0.99569, 0.99672, 0.99329, 0.996224, 0.993092, 0.988473, 0.985855, 0.980278, 0.970425, 0.948918, 0.930153, 0.909393, 0.892497, 0.871268, 0.850926, 0.822264, 0.766328, 0.792829, 0.73769, 0.708725, 0.654411, 0.668514, 0.649787, 0.553934, 0.511755, 0.452517, 0.4702, 0.366037, 0.396694, 0.33945, 0.176684},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 0.99977, 1, 0.999243, 0.999605, 0.999711, 0.999684, 0.999118, 0.999055, 0.999027, 0.999345, 0.99832, 0.996205, 0.997296, 0.994503, 0.996224, 0.993092, 0.990931, 0.99036, 0.989581, 0.981552, 0.96774, 0.943401, 0.930563, 0.917317, 0.900146, 0.890127, 0.859196, 0.808908, 0.82996, 0.802677, 0.765675, 0.744748, 0.771866, 0.736965, 0.654208, 0.603581, 0.601261, 0.574562, 0.520457, 0.491735, 0.417432, 0.233101},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 0.99977, 1, 0.999497, 0.999605, 0.999711, 0.999684, 0.999288, 0.999237, 0.999416, 0.999563, 0.998555, 0.996771, 0.997296, 0.995715, 0.996224, 0.994994, 0.99175, 0.992162, 0.992682, 0.989597, 0.978816, 0.956557, 0.947707, 0.940359, 0.923348, 0.915111, 0.894974, 0.850199, 0.868519, 0.846758, 0.819174, 0.802846, 0.833438, 0.823836, 0.751506, 0.698276, 0.667052, 0.647616, 0.605299, 0.590909, 0.536698, 0.305931},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 0.99977, 1, 0.999624, 0.999605, 0.999711, 0.999684, 0.999457, 0.999237, 0.999611, 0.999563, 0.998555, 0.996771, 0.997296, 0.996927, 0.996902, 0.995437, 0.993388, 0.993063, 0.994344, 0.992921, 0.986395, 0.972483, 0.955499, 0.954866, 0.935412, 0.931101, 0.915748, 0.887618, 0.898605, 0.881518, 0.857141, 0.85907, 0.866423, 0.869562, 0.816371, 0.790102, 0.775749, 0.713712, 0.727029, 0.690083, 0.633028, 0.387988},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 0.99977, 1, 0.999751, 0.999605, 0.999855, 0.999684, 0.999457, 0.9996, 0.999611, 0.999563, 0.998555, 0.997029, 0.997584, 0.99723, 0.99758, 0.996989, 0.994208, 0.993964, 0.994861, 0.99453, 0.989892, 0.980616, 0.973549, 0.965107, 0.951189, 0.944093, 0.936523, 0.902008, 0.920026, 0.91779, 0.879866, 0.877811, 0.912601, 0.891222, 0.873128, 0.86184, 0.830098, 0.786766, 0.800805, 0.752066, 0.715596, 0.484411},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 0.99977, 1, 0.999751, 0.999605, 1, 0.999684, 0.999457, 0.9996, 0.999805, 0.999563, 0.99879, 0.997029, 0.997871, 0.99723, 0.99758, 0.997433, 0.994208, 0.994865, 0.994861, 0.995603, 0.991723, 0.987438, 0.979585, 0.975348, 0.958614, 0.952285, 0.94691, 0.917492, 0.934307, 0.932904, 0.909205, 0.907797, 0.925235, 0.917694, 0.900155, 0.919231, 0.870145, 0.845905, 0.84876, 0.818182, 0.784404, 0.584424},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 0.999885, 1, 0.999751, 0.999605, 1, 0.999684, 0.999457, 0.9996, 1, 0.999563, 0.999201, 0.997286, 0.998735, 0.99723, 0.997985, 0.997802, 0.995027, 0.995316, 0.99628, 0.996139, 0.994638, 0.988823, 0.988311, 0.983882, 0.96975, 0.961279, 0.951526, 0.929105, 0.94002, 0.941972, 0.924736, 0.928412, 0.931832, 0.927321, 0.908263, 0.939318, 0.893028, 0.863299, 0.874581, 0.834711, 0.844037, 0.649299},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 0.999885, 1, 0.999751, 0.999605, 1, 0.999842, 0.999457, 0.9996, 1, 0.999563, 0.999483, 0.998251, 0.998735, 0.99723, 0.998916, 0.997802, 0.995436, 0.996217, 0.997313, 0.996676, 0.995804, 0.992285, 0.991428, 0.987295, 0.976247, 0.964277, 0.961914, 0.937809, 0.948588, 0.944994, 0.93164, 0.941531, 0.949423, 0.944167, 0.935823, 0.942188, 0.915912, 0.88765, 0.896714, 0.842975, 0.880734, 0.707363},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 0.999885, 1, 0.999751, 0.999605, 1, 0.999842, 0.999457, 0.9996, 1, 0.999563, 0.999483, 0.998508, 0.999023, 0.997533, 0.999255, 0.998172, 0.996256, 0.997207, 0.99783, 0.997212, 0.995804, 0.99367, 0.992207, 0.989856, 0.982743, 0.971273, 0.967684, 0.952002, 0.952873, 0.955574, 0.944042, 0.950901, 0.962617, 0.961494, 0.946634, 0.947927, 0.927354, 0.912001, 0.915158, 0.876033, 0.90367, 0.762856},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 0.999885, 1, 0.999751, 0.999737, 1, 0.999842, 0.999457, 0.9996, 1, 0.999563, 0.999483, 0.998766, 0.999023, 0.997896, 0.999594, 0.998542, 0.998199, 0.997658, 0.99783, 0.997748, 0.99697, 0.995884, 0.994545, 0.991562, 0.984599, 0.981012, 0.973455, 0.958454, 0.961441, 0.958596, 0.953062, 0.95465, 0.964816, 0.961494, 0.965553, 0.953666, 0.959954, 0.918959, 0.929913, 0.900826, 0.917431, 0.81117},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 0.999885, 1, 0.999751, 0.999737, 1, 0.999842, 0.999627, 0.9996, 1, 0.999563, 1, 0.999074, 0.999367, 0.997896, 0.999594, 0.998542, 0.99869, 0.998109, 0.99783, 0.998821, 0.997553, 0.996577, 0.995324, 0.993269, 0.989239, 0.988008, 0.975763, 0.967486, 0.97001, 0.966153, 0.956513, 0.965894, 0.971413, 0.973527, 0.970959, 0.962275, 0.968535, 0.946789, 0.948357, 0.933884, 0.944954, 0.856817},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 0.999885, 1, 0.999751, 0.999737, 1, 0.999842, 0.999627, 0.999782, 1, 0.999782, 1, 0.999691, 0.999367, 0.997896, 1, 0.998912, 0.99959, 0.998109, 0.998347, 0.998821, 0.998136, 0.997269, 0.996883, 0.994123, 0.992208, 0.992005, 0.98615, 0.977809, 0.97715, 0.972198, 0.963417, 0.973391, 0.984607, 0.97834, 0.98177, 0.970883, 0.971395, 0.960704, 0.963112, 0.950413, 0.958716, 0.901951},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 0.999885, 1, 0.999879, 0.999737, 1, 0.999842, 0.999627, 0.999782, 1, 0.999782, 1, 1, 0.999367, 0.998199, 1, 0.999355, 0.99959, 0.99946, 0.998864, 1, 0.998136, 0.99917, 0.999221, 0.994976, 0.997032, 0.997002, 0.988459, 0.990712, 0.987147, 0.984289, 0.973771, 0.980887, 0.984607, 0.98556, 0.991892, 0.976622, 0.988558, 0.978097, 0.981556, 0.966942, 0.96789, 0.936842},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 0.999885, 1, 0.999879, 1, 1, 0.999842, 1, 0.999782, 1, 0.999782, 1, 1, 0.999712, 0.999031, 1, 0.999355, 1, 0.99946, 0.999483, 1, 0.998834, 0.99917, 0.999221, 0.997488, 0.998888, 0.997002, 0.994229, 0.995873, 0.990003, 0.993357, 0.989645, 0.996252, 0.993403, 0.995187, 0.994595, 0.988522, 0.994279, 1, 0.981556, 0.987603, 0.981651, 0.970794},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}
},
{
{1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.000361987, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0.00307272, 0.00204961, 0, 0, 0.000506001, 0.000877403, 0.000802055, 0, 0.00109366, 0.00152773, 0.000761093, 0.000723974, 0.00151618, 0.000790918, 0.00042351, 0, 0, 0.000472969, 0, 0.000525132, 0.000532796, 0.000595788, 0, 0.000668414, 0, 0, 0.000928225, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0.015313, 0.0085452, 0.00404888, 0, 0.00913452, 0.00913957, 0.00321058, 0.00122339, 0.00357033, 0.00218157, 0.00492455, 0.0030023, 0.00297763, 0.00301443, 0.00334424, 0.00474636, 0.00379044, 0.00279034, 0.00205021, 0.000842537, 0.00277929, 0.00318546, 0.000499606, 0.00268443, 0.00272362, 0.00122492, 0.000664441, 0.00133683, 0.00149004, 0, 0.00180725, 0.000912338, 0, 0, 0, 0.0012852, 0, 0, 0, 0, 0.000662254},
{1, 1, 1, 1, 1, 1, 0, 0, 0.140609, 0.104094, 0.015313, 0.0494595, 0.0287955, 0.0116039, 0.018269, 0.0142357, 0.00715826, 0.00421704, 0.00666298, 0.00528831, 0.00855625, 0.00494697, 0.00597566, 0.00486767, 0.00556701, 0.00734104, 0.00604407, 0.00402098, 0.00410043, 0.00212992, 0.00418155, 0.00500214, 0.00202639, 0.00373469, 0.00381903, 0.00182071, 0.00132888, 0.00271107, 0.00149004, 0.00154558, 0.00180725, 0.00185543, 0, 0.00115691, 0.00121575, 0.0012852, 0.00151902, 0, 0.00573377, 0, 0.00132451},
{1, 1, 1, 1, 1, 1, 1, 0, 0.570305, 0.104094, 0.030626, 0.0903738, 0.0411687, 0.0170984, 0.0242472, 0.0203845, 0.0119741, 0.00956021, 0.00919298, 0.00877403, 0.0113859, 0.00873093, 0.0112222, 0.00892019, 0.00704886, 0.00957375, 0.00866643, 0.00641586, 0.00530362, 0.00341731, 0.00645979, 0.00637091, 0.00355317, 0.00425983, 0.00598003, 0.00182071, 0.00199332, 0.00337948, 0.00149004, 0.00231837, 0.00180725, 0.00185543, 0, 0.00115691, 0.00249953, 0.00264232, 0.00151902, 0, 0.00757628, 0, 0.00177005},
{1, 1, 1, 1, 1, 1, 1, 0, 0.570305, 0.239765, 0.0894433, 0.13938, 0.0701908, 0.0283948, 0.0411889, 0.0275303, 0.0159218, 0.0149358, 0.0159692, 0.0130894, 0.0130124, 0.0136252, 0.0130858, 0.0119141, 0.00965218, 0.01068, 0.0102032, 0.00960165, 0.00655171, 0.00386215, 0.0096897, 0.00773969, 0.00557956, 0.00531009, 0.00710526, 0.0024165, 0.00199332, 0.00542213, 0.00298008, 0.00309116, 0.00180725, 0.00373117, 0.00197012, 0.00115691, 0.00493102, 0.00392752, 0.00303803, 0, 0.00941879, 0.00203626, 0.00266115},
{1, 1, 1, 1, 1, 1, 1, 1, 0.710913, 0.343859, 0.134571, 0.171749, 0.120364, 0.0731192, 0.0611198, 0.0377226, 0.0238609, 0.023143, 0.0211992, 0.0188045, 0.0174687, 0.0162574, 0.0156785, 0.0159666, 0.0107333, 0.0135961, 0.0139111, 0.0100192, 0.00855703, 0.00648407, 0.0128692, 0.0100043, 0.00707838, 0.00583522, 0.00763806, 0.00483299, 0.00388099, 0.00609055, 0.00451181, 0.00545278, 0.00180725, 0.00464351, 0.00295519, 0.00231381, 0.00493102, 0.00656984, 0.00303803, 0.00152724, 0.0112613, 0.00203626, 0.00313729},
{1, 1, 1, 1, 1, 1, 1, 1, 0.859391, 0.516725, 0.253017, 0.278761, 0.169857, 0.107162, 0.0795561, 0.059104, 0.0334052, 0.0309005, 0.0278905, 0.0257759, 0.0203656, 0.0225907, 0.0200531, 0.0211803, 0.0136365, 0.0184032, 0.0161648, 0.0132936, 0.0122788, 0.0112359, 0.0137956, 0.0104772, 0.00810555, 0.00810449, 0.00982888, 0.00546212, 0.00580388, 0.00679637, 0.00678857, 0.00545278, 0.00180725, 0.00464351, 0.00394025, 0.00231381, 0.00493102, 0.00656984, 0.00303803, 0.00458173, 0.0112613, 0.00203626, 0.00532775},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 0.552047, 0.357773, 0.337671, 0.232403, 0.143645, 0.126723, 0.0869703, 0.0525813, 0.0457067, 0.0371684, 0.0367077, 0.0273627, 0.0294391, 0.0275278, 0.0252739, 0.0187624, 0.0231699, 0.0195659, 0.0165015, 0.0147077, 0.0151452, 0.0156737, 0.0118209, 0.010104, 0.00967988, 0.0119899, 0.00668704, 0.00772677, 0.00746479, 0.00981034, 0.00699836, 0.00180725, 0.00651925, 0.00591037, 0.00231381, 0.00614677, 0.00914024, 0.00303803, 0.00763621, 0.0132069, 0.00203626, 0.0067195},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 0.654269, 0.475407, 0.44423, 0.314741, 0.194325, 0.160607, 0.123696, 0.0780038, 0.060448, 0.0522389, 0.0494181, 0.0335801, 0.037315, 0.0345155, 0.0282883, 0.0239489, 0.0283795, 0.0225569, 0.0185452, 0.0184295, 0.0185625, 0.0180276, 0.0150952, 0.0126579, 0.0118392, 0.0136479, 0.00731617, 0.0102789, 0.0081332, 0.00981034, 0.00777115, 0.00361451, 0.00935838, 0.00695057, 0.00347072, 0.00743055, 0.0117106, 0.00768009, 0.00763621, 0.0150494, 0.00203626, 0.00738175},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 0.829007, 0.624479, 0.535509, 0.418002, 0.28362, 0.227788, 0.160087, 0.0996971, 0.0811117, 0.0713291, 0.0616779, 0.0413118, 0.0441012, 0.0426581, 0.0346016, 0.0294554, 0.0350574, 0.0270848, 0.0221928, 0.0233096, 0.0211845, 0.020331, 0.0173849, 0.0162111, 0.0150781, 0.0158089, 0.0103951, 0.0109433, 0.0101384, 0.0120871, 0.0100895, 0.00454273, 0.00935838, 0.00793563, 0.0078036, 0.00743055, 0.0129958, 0.00768009, 0.00763621, 0.0150494, 0.00203626, 0.0094115},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 0.897778, 0.713923, 0.617791, 0.512713, 0.38986, 0.30053, 0.209828, 0.137443, 0.11021, 0.0930343, 0.074245, 0.05486, 0.0531081, 0.0496866, 0.0420146, 0.0360834, 0.0391, 0.030834, 0.0242144, 0.0266079, 0.0229167, 0.0226597, 0.0201475, 0.018793, 0.0187835, 0.0185027, 0.0140699, 0.0129014, 0.0128495, 0.0135771, 0.0114382, 0.00810806, 0.0120031, 0.00793563, 0.0101174, 0.00986205, 0.0129958, 0.00768009, 0.0108616, 0.0150494, 0.00622275, 0.0104398},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 0.931228, 0.775174, 0.734706, 0.595731, 0.488473, 0.394039, 0.29923, 0.190024, 0.154401, 0.121346, 0.0947088, 0.0716388, 0.0641313, 0.0589432, 0.05034, 0.0423308, 0.0469246, 0.0360994, 0.0310478, 0.0286356, 0.0289087, 0.0272667, 0.0228851, 0.0219025, 0.020884, 0.0201309, 0.0158906, 0.0141951, 0.0135553, 0.0135771, 0.0130271, 0.00903628, 0.0129155, 0.00793563, 0.0159667, 0.0123616, 0.0169234, 0.0108031, 0.0108616, 0.0150494, 0.00622275, 0.0122},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 0.931228, 0.821114, 0.80035, 0.694037, 0.584646, 0.507562, 0.368303, 0.24972, 0.214687, 0.157024, 0.120556, 0.0901564, 0.0830303, 0.0716644, 0.0582016, 0.0511815, 0.0528784, 0.0416922, 0.0342115, 0.032335, 0.0332629, 0.0304461, 0.0255975, 0.0255116, 0.0240642, 0.021789, 0.0194136, 0.0154887, 0.0169722, 0.0173022, 0.0162047, 0.00991531, 0.0138278, 0.00897582, 0.0183453, 0.0161178, 0.0169234, 0.0123222, 0.0108616, 0.0150494, 0.00837297, 0.0135366},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.910557, 0.892083, 0.780877, 0.682951, 0.612274, 0.459148, 0.333758, 0.287012, 0.203669, 0.157173, 0.118032, 0.109003, 0.0917787, 0.0735582, 0.0636562, 0.0627939, 0.0514441, 0.0386501, 0.037706, 0.0384124, 0.033676, 0.0288081, 0.0326738, 0.0272443, 0.0228844, 0.0218634, 0.0181113, 0.0183465, 0.0180473, 0.0177935, 0.00991531, 0.0147402, 0.0110011, 0.0183453, 0.0186174, 0.0182086, 0.0123222, 0.0108616, 0.016995, 0.0104092, 0.015304},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.93956, 0.925358, 0.892689, 0.786905, 0.693241, 0.563455, 0.436578, 0.368411, 0.27452, 0.216099, 0.154888, 0.141699, 0.111997, 0.0922138, 0.0784267, 0.0746612, 0.0566269, 0.0466699, 0.0430545, 0.0453178, 0.0400854, 0.0324916, 0.036255, 0.0315924, 0.0261408, 0.0248423, 0.0233043, 0.0211324, 0.0187923, 0.0201552, 0.0126508, 0.0147402, 0.0129712, 0.0219455, 0.0186174, 0.0197471, 0.0154452, 0.0124744, 0.016995, 0.0125595, 0.0164239},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.985499, 0.941995, 0.946458, 0.84552, 0.776946, 0.671415, 0.536662, 0.449893, 0.356759, 0.276234, 0.204824, 0.182272, 0.137682, 0.114942, 0.097169, 0.0891838, 0.0728123, 0.0572174, 0.0533056, 0.0509357, 0.0487982, 0.0402563, 0.0449161, 0.0359404, 0.02993, 0.0285171, 0.0295966, 0.0239183, 0.0225174, 0.0233328, 0.0162653, 0.0156525, 0.0170768, 0.0267026, 0.0186174, 0.0236746, 0.0180962, 0.0156143, 0.0206801, 0.0125595, 0.0175438},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.985499, 0.967178, 0.959058, 0.932529, 0.851852, 0.755386, 0.643255, 0.54993, 0.451566, 0.356644, 0.26746, 0.239976, 0.171165, 0.143741, 0.124341, 0.109478, 0.092092, 0.072444, 0.0615064, 0.0574198, 0.0594395, 0.0438897, 0.050024, 0.0397044, 0.034252, 0.0357666, 0.0308902, 0.0259609, 0.029431, 0.0249217, 0.0171935, 0.0184406, 0.0251227, 0.030238, 0.0199011, 0.0249598, 0.0180962, 0.017227, 0.0226257, 0.0147097, 0.0184921},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.975253, 0.957715, 0.914382, 0.846175, 0.753884, 0.653087, 0.559108, 0.447748, 0.339323, 0.302944, 0.227637, 0.181541, 0.155918, 0.138683, 0.109342, 0.0877591, 0.0709554, 0.0693568, 0.0695546, 0.055338, 0.0576299, 0.0498877, 0.0402619, 0.0406663, 0.0334071, 0.0300836, 0.0332395, 0.0272833, 0.0190008, 0.0222941, 0.0251227, 0.0349952, 0.0199011, 0.0303164, 0.0180962, 0.0188398, 0.0302019, 0.0168599, 0.0207044},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.983804, 0.980462, 0.949342, 0.908381, 0.832799, 0.738707, 0.661955, 0.542501, 0.42315, 0.379115, 0.283443, 0.229547, 0.194734, 0.172654, 0.128088, 0.103954, 0.0841486, 0.0835652, 0.0801707, 0.0631278, 0.0631815, 0.0584955, 0.0468045, 0.0449368, 0.0385818, 0.0355057, 0.0362196, 0.0327793, 0.0216379, 0.0288337, 0.0291532, 0.0373737, 0.0211849, 0.0369581, 0.0196153, 0.023507, 0.0302019, 0.0211603, 0.0227292},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.991902, 0.997253, 0.976913, 0.952197, 0.893276, 0.82912, 0.741744, 0.642212, 0.527012, 0.461127, 0.359829, 0.293457, 0.252783, 0.207752, 0.161648, 0.129861, 0.106188, 0.101963, 0.0999755, 0.0764428, 0.0723422, 0.0665782, 0.0571364, 0.0528821, 0.0431272, 0.0448636, 0.0407314, 0.037589, 0.0261314, 0.0336784, 0.0301934, 0.039817, 0.0224687, 0.0434561, 0.0227383, 0.023507, 0.0323535, 0.0233106, 0.0245114},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.995951, 1, 0.986047, 0.97164, 0.931366, 0.892075, 0.825979, 0.739741, 0.620021, 0.554059, 0.445308, 0.359425, 0.309572, 0.256206, 0.209219, 0.163436, 0.135269, 0.122655, 0.117376, 0.093764, 0.0821702, 0.0771983, 0.0664026, 0.0638063, 0.0489311, 0.0523658, 0.0505417, 0.0447172, 0.0305266, 0.0410793, 0.0331486, 0.0421308, 0.0237525, 0.0513111, 0.0275504, 0.0251197, 0.034196, 0.0296473, 0.0264982},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.989204, 0.98377, 0.965077, 0.93305, 0.894344, 0.828253, 0.729898, 0.654456, 0.539053, 0.447531, 0.386464, 0.324791, 0.261238, 0.201779, 0.176161, 0.151143, 0.142462, 0.117531, 0.0949959, 0.0911456, 0.0735079, 0.0742681, 0.0604332, 0.0544458, 0.057372, 0.0510725, 0.0394152, 0.0475167, 0.0391692, 0.0456016, 0.0311831, 0.0513111, 0.0308897, 0.0299578, 0.0398266, 0.0338338, 0.0278533},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.993938, 0.990972, 0.978569, 0.962533, 0.939523, 0.892868, 0.817202, 0.747059, 0.630995, 0.545025, 0.472786, 0.390904, 0.329729, 0.259635, 0.2216, 0.188453, 0.173258, 0.142642, 0.109844, 0.106231, 0.0865634, 0.0888338, 0.0675307, 0.0647339, 0.0673491, 0.0668743, 0.0440072, 0.0540052, 0.0492404, 0.0515803, 0.0398974, 0.0591662, 0.0371358, 0.034625, 0.0436147, 0.0399425, 0.0303099},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.996927, 0.993963, 0.99128, 0.983326, 0.970001, 0.944607, 0.886547, 0.824304, 0.725611, 0.648, 0.559832, 0.473732, 0.408978, 0.325028, 0.280363, 0.23618, 0.212067, 0.171014, 0.133773, 0.12622, 0.112704, 0.110749, 0.0836135, 0.075653, 0.0749912, 0.0740457, 0.051187, 0.0605958, 0.0594769, 0.0538941, 0.0472599, 0.0670931, 0.0432969, 0.0394632, 0.0474029, 0.0419788, 0.0333923},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.995957, 0.994447, 0.992822, 0.985941, 0.974914, 0.937731, 0.887581, 0.811155, 0.736656, 0.658949, 0.558653, 0.496235, 0.406811, 0.345636, 0.291844, 0.269928, 0.215588, 0.168889, 0.154757, 0.13396, 0.13091, 0.102913, 0.0934807, 0.0857134, 0.0898042, 0.0556806, 0.0709379, 0.0675228, 0.0575591, 0.052327, 0.0710926, 0.0448159, 0.0442158, 0.0531366, 0.0440151, 0.0357227},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.996953, 0.996833, 0.995815, 0.994292, 0.984398, 0.963156, 0.937533, 0.881525, 0.820527, 0.745955, 0.649725, 0.587673, 0.493605, 0.426798, 0.364571, 0.332121, 0.269024, 0.207586, 0.192368, 0.164926, 0.153388, 0.122142, 0.105068, 0.100161, 0.102428, 0.0665241, 0.0819881, 0.076664, 0.06463, 0.0598937, 0.0763053, 0.047854, 0.0505812, 0.0554681, 0.052388, 0.0399965},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.998947, 0.998438, 0.99823, 0.997948, 0.993954, 0.978664, 0.965809, 0.927869, 0.883525, 0.829835, 0.749267, 0.66491, 0.588859, 0.511726, 0.440952, 0.402577, 0.336818, 0.268382, 0.242293, 0.207528, 0.184969, 0.153713, 0.12635, 0.116098, 0.122008, 0.0837176, 0.0922792, 0.090981, 0.0729874, 0.0675653, 0.0842323, 0.0586571, 0.0553339, 0.0612019, 0.052388, 0.0458434},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.998947, 0.998438, 0.999421, 0.998988, 0.997439, 0.990541, 0.984071, 0.96279, 0.929084, 0.887224, 0.824412, 0.756574, 0.672884, 0.600503, 0.523603, 0.476817, 0.40559, 0.336867, 0.29817, 0.261911, 0.234857, 0.184097, 0.154429, 0.152359, 0.14644, 0.107114, 0.110985, 0.108113, 0.0824369, 0.0774954, 0.0948735, 0.0664222, 0.0648393, 0.0706206, 0.0587247, 0.0494532},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.999494, 0.999573, 0.995838, 0.990899, 0.982032, 0.965234, 0.93368, 0.882399, 0.827935, 0.748252, 0.685541, 0.607354, 0.558148, 0.483836, 0.422202, 0.357694, 0.317952, 0.281633, 0.224653, 0.192014, 0.181254, 0.175682, 0.133196, 0.130553, 0.12322, 0.100308, 0.0961398, 0.105515, 0.0822925, 0.0773991, 0.0781969, 0.069134, 0.0543874},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.998267, 0.995074, 0.99135, 0.985294, 0.960953, 0.931374, 0.89075, 0.81626, 0.764004, 0.690843, 0.646131, 0.553969, 0.50598, 0.438227, 0.379291, 0.337184, 0.275594, 0.23829, 0.220042, 0.209663, 0.169095, 0.155799, 0.138106, 0.133479, 0.109989, 0.124937, 0.0914916, 0.0867335, 0.0950888, 0.0775069, 0.0599718},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.999198, 0.998055, 0.995503, 0.992443, 0.977233, 0.960819, 0.933548, 0.879061, 0.823794, 0.768174, 0.718414, 0.638686, 0.589817, 0.51663, 0.45588, 0.413624, 0.342407, 0.293815, 0.268437, 0.258571, 0.198644, 0.199853, 0.168374, 0.15707, 0.132621, 0.142199, 0.114787, 0.103875, 0.108296, 0.0880302, 0.0662216},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.999599, 0.998796, 0.998521, 0.997352, 0.987185, 0.974818, 0.958973, 0.918657, 0.883584, 0.836547, 0.782684, 0.721164, 0.662271, 0.5976, 0.535014, 0.489531, 0.41586, 0.370436, 0.326268, 0.311645, 0.247095, 0.24274, 0.207674, 0.193711, 0.162411, 0.160696, 0.133525, 0.122544, 0.125394, 0.102626, 0.0736125},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.999166, 0.999635, 0.998472, 0.993392, 0.987449, 0.974769, 0.954868, 0.924967, 0.887129, 0.844525, 0.796647, 0.744189, 0.688081, 0.618377, 0.570734, 0.502566, 0.45196, 0.39283, 0.370086, 0.308296, 0.279999, 0.2542, 0.240829, 0.194565, 0.193257, 0.16619, 0.141298, 0.157748, 0.125595, 0.0836385},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.999556, 0.999635, 0.999613, 0.997097, 0.994107, 0.986078, 0.975501, 0.953403, 0.929296, 0.899256, 0.860134, 0.810036, 0.767446, 0.705503, 0.678289, 0.584797, 0.540206, 0.484275, 0.462021, 0.381902, 0.36562, 0.317857, 0.290844, 0.25185, 0.252386, 0.208797, 0.183388, 0.197265, 0.150144, 0.0961259},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.999556, 1, 1, 0.998538, 0.998532, 0.993576, 0.986422, 0.979476, 0.959795, 0.945574, 0.924219, 0.872582, 0.842114, 0.794436, 0.754858, 0.680981, 0.625628, 0.585865, 0.562582, 0.467328, 0.446627, 0.402422, 0.369013, 0.316634, 0.328654, 0.280627, 0.2489, 0.265142, 0.187594, 0.120522},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.999556, 1, 1, 1, 0.999256, 0.997326, 0.993628, 0.992544, 0.984468, 0.973292, 0.953462, 0.925548, 0.902894, 0.865093, 0.829941, 0.775801, 0.732295, 0.685177, 0.660743, 0.567993, 0.560638, 0.510359, 0.456242, 0.415059, 0.419202, 0.358278, 0.328413, 0.317661, 0.248355, 0.159941},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.999556, 1, 1, 1, 1, 0.998821, 0.999209, 0.997469, 0.99348, 0.990317, 0.976283, 0.960304, 0.943744, 0.914118, 0.895653, 0.855551, 0.812033, 0.773309, 0.761953, 0.681309, 0.665118, 0.623057, 0.568803, 0.527674, 0.499613, 0.467244, 0.420828, 0.397315, 0.340344, 0.204482},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.999556, 1, 1, 1, 1, 0.999559, 0.999605, 0.998294, 0.996476, 0.994899, 0.990519, 0.974544, 0.969421, 0.944019, 0.929155, 0.90394, 0.867522, 0.843024, 0.838427, 0.771813, 0.755402, 0.717803, 0.683483, 0.629543, 0.580878, 0.569285, 0.505265, 0.514541, 0.442856, 0.269615},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.999556, 1, 1, 1, 1, 0.999559, 1, 0.998695, 0.998605, 0.997062, 0.994575, 0.985815, 0.979604, 0.962078, 0.946799, 0.934811, 0.91361, 0.889721, 0.89395, 0.846593, 0.82252, 0.813479, 0.762873, 0.729991, 0.693994, 0.66576, 0.601578, 0.617425, 0.553854, 0.35676},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.999556, 1, 1, 1, 1, 0.999559, 1, 0.999096, 0.99905, 0.997513, 0.997313, 0.990839, 0.986578, 0.973949, 0.960836, 0.954252, 0.932293, 0.913317, 0.921646, 0.887035, 0.881478, 0.8727, 0.840079, 0.793899, 0.778404, 0.729655, 0.681348, 0.704759, 0.643237, 0.442496},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.999556, 1, 1, 1, 1, 0.999559, 1, 0.999519, 0.99905, 0.997963, 0.999104, 0.992894, 0.990371, 0.983245, 0.973714, 0.965124, 0.949971, 0.935256, 0.940626, 0.919516, 0.915037, 0.893827, 0.879034, 0.851049, 0.840319, 0.788312, 0.754923, 0.757999, 0.712257, 0.536637},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.999556, 1, 1, 1, 1, 0.999559, 1, 0.999519, 0.999555, 0.998439, 0.999104, 0.99492, 0.993026, 0.98434, 0.979839, 0.973507, 0.958217, 0.946765, 0.95252, 0.934902, 0.934604, 0.91211, 0.905263, 0.87476, 0.875955, 0.835158, 0.802108, 0.803044, 0.768717, 0.601143},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.999556, 1, 1, 1, 1, 0.999559, 1, 0.999519, 1, 0.998439, 1, 0.998501, 0.996236, 0.988662, 0.984043, 0.979904, 0.965163, 0.955085, 0.965787, 0.951948, 0.949559, 0.925246, 0.925319, 0.905562, 0.899808, 0.850518, 0.844284, 0.829767, 0.80809, 0.660932},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.999556, 1, 1, 1, 1, 1, 1, 0.999519, 1, 0.99889, 1, 0.999001, 0.99892, 0.994583, 0.988976, 0.984415, 0.971328, 0.960342, 0.970467, 0.963769, 0.961624, 0.941543, 0.940553, 0.927786, 0.921018, 0.887825, 0.891127, 0.858229, 0.832753, 0.718695},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.999556, 1, 1, 1, 1, 1, 1, 0.999519, 1, 0.99889, 1, 0.9995, 0.999445, 0.995648, 0.99328, 0.987728, 0.977493, 0.969574, 0.975963, 0.97007, 0.972068, 0.952544, 0.953667, 0.941431, 0.939514, 0.917792, 0.918959, 0.884643, 0.863867, 0.774261},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.999556, 1, 1, 1, 1, 1, 1, 1, 1, 0.99889, 1, 0.9995, 0.999445, 0.996773, 0.996359, 0.991574, 0.982285, 0.977811, 0.982405, 0.976371, 0.982308, 0.964806, 0.96183, 0.961359, 0.959152, 0.948767, 0.936101, 0.907372, 0.897374, 0.830506},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.999556, 1, 1, 1, 1, 1, 1, 1, 1, 0.99943, 1, 1, 0.999445, 0.997306, 0.996954, 0.99479, 0.989081, 0.986876, 0.98481, 0.988228, 0.987935, 0.978762, 0.973658, 0.976288, 0.973577, 0.964297, 0.956553, 0.94095, 0.926452, 0.878793},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.999556, 1, 1, 1, 1, 1, 1, 1, 1, 0.99943, 1, 1, 1, 0.998934, 0.999371, 0.998042, 0.993835, 0.989981, 0.987988, 0.990036, 0.988848, 0.982813, 0.982978, 0.986287, 0.992073, 0.978564, 0.976749, 0.969695, 0.951798, 0.921734},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.999556, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.999467, 1, 0.999336, 0.998626, 0.993873, 0.992711, 0.993748, 0.992599, 0.991789, 0.988957, 0.993717, 0.994715, 0.99232, 0.984386, 0.986896, 0.979068, 0.96219},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}
}
};

static double Mrev [51][51] = {
{0.017965736, 0.000410897, -0.000114951, 2.64396e-06, -9.67534e-05, -1.35961e-05, -2.8951e-05, -0.000115159, 4.37729e-05, 9.52452e-05, 4.26502e-06, 4.8405e-05, 8.24648e-05, -9.59416e-05, -4.96514e-05, -0.000136782, -8.05224e-05, 4.96287e-05, 4.14411e-05, -1.24779e-05, 9.83332e-05, 9.50702e-05, -0.000137379, 5.89401e-05, 0.000128391, -0.000207407, -0.000109347, -4.65733e-05, 8.742e-05, 0.000136528, 0.00013447, -0.000143867, -0.000246086, 0.000224706, 5.70979e-05, 0.000138589, -0.000199036, -4.6457e-05, 3.28348e-05, -8.32483e-05, -6.07291e-05, 0.000379313, 0.000308066, 0.000371455, -0.000798599, 4.15331e-05, 4.46331e-05, -0.000273535, -0.000210617, 0.000442616, -0.000256442},
{0.000410897, 0.001069948, -0.0001096, -6.2837e-05, -1.92761e-05, -5.6709e-05, -4.56913e-06, 9.99485e-06, -3.94752e-06, 2.27827e-05, 8.06401e-06, -5.67898e-06, 3.71524e-05, -2.7878e-06, 1.25567e-05, 1.60544e-07, -9.85258e-05, -2.76498e-05, 3.15176e-05, 1.87005e-05, 5.64557e-05, -1.36524e-05, 1.67219e-05, 6.06075e-06, -2.68488e-06, -6.93407e-05, 2.9713e-06, 1.83319e-06, 4.14368e-05, -5.25881e-05, 5.62995e-05, -3.51059e-05, 4.02116e-05, -0.000104181, 2.85296e-05, 1.53898e-05, 7.20497e-05, -5.57023e-05, 1.92546e-05, -2.29e-05, 0.000130883, -9.53455e-05, -2.12528e-05, -4.31274e-05, 2.42076e-05, -4.6987e-05, 3.83722e-05, 5.05434e-05, -2.59119e-05, -5.19021e-05, 3.98843e-05},
{-0.000114951, -0.0001096, 0.000254012, -7.43654e-05, -4.55369e-05, -4.45989e-05, -1.72471e-05, -4.60799e-06, 7.14859e-06, 1.11382e-05, -2.95438e-06, -2.89264e-06, -8.86647e-06, 1.87034e-06, -9.19841e-06, -1.40444e-05, 3.03036e-05, 2.23746e-05, 2.52406e-05, 5.98781e-06, -1.08933e-05, -4.40666e-06, 6.18351e-06, 2.50943e-06, -1.71965e-05, 2.09407e-05, -1.04552e-06, 1.67884e-06, -3.27843e-07, 4.51441e-07, -2.34602e-05, -3.4446e-05, 4.48407e-05, 5.14078e-06, -4.88281e-05, 5.07297e-07, 4.83132e-07, -3.09784e-05, 1.40378e-05, 5.36181e-05, -3.45689e-05, -4.88781e-05, 4.99183e-05, 1.14867e-05, 2.10336e-05, 2.27733e-05, 1.93506e-05, 3.50583e-05, -3.65015e-05, -2.10393e-05, 1.19656e-05},
{2.64396e-06, -6.2837e-05, -7.43654e-05, 0.000164372, -6.10247e-05, -2.3351e-05, -2.43579e-05, -1.18952e-05, -1.60087e-05, -8.70781e-06, 7.41117e-06, 5.40462e-06, 1.68764e-05, 1.12127e-05, 1.92879e-05, 1.88575e-05, -1.62663e-05, 3.31625e-06, -1.69895e-05, 4.52451e-06, -1.12592e-05, 1.49694e-05, 1.36173e-05, 8.45615e-06, 4.91549e-06, 7.11254e-07, -1.07457e-05, -3.13133e-05, 8.60804e-06, 5.8603e-07, -1.13582e-05, -1.71656e-06, -1.60143e-05, -6.68248e-06, 2.16778e-06, 9.64015e-06, -2.42714e-06, 1.82499e-05, 2.64127e-05, -7.90122e-06, -1.07072e-05, 2.28551e-05, -2.77283e-05, 8.72761e-06, 2.71738e-05, -1.27263e-05, -7.13011e-06, -7.39945e-05, 4.12342e-05, -9.94756e-06, 6.76566e-06},
{-9.67534e-05, -1.92761e-05, -4.55369e-05, -6.10247e-05, 0.000194399, -3.69655e-05, -1.93361e-05, -1.23473e-05, -2.58623e-05, -3.43379e-06, -1.4417e-06, -3.28384e-06, 1.29683e-05, -9.40662e-07, 3.03201e-06, 2.11594e-07, 8.31762e-07, 7.35234e-06, -3.02049e-05, -1.6072e-05, -1.23793e-05, 4.95934e-06, 1.49319e-05, -1.85484e-05, 1.89428e-05, 5.5515e-06, 1.56733e-05, 3.72499e-06, 6.68261e-06, 2.4015e-06, -7.38868e-06, 3.428e-05, -3.08518e-05, 1.8012e-06, 5.87077e-06, 2.66374e-05, -2.76056e-06, -1.48868e-05, -2.69686e-05, -3.02662e-06, 2.73859e-05, 2.39657e-05, -6.01345e-06, 4.39415e-06, -5.59135e-05, -2.50076e-05, -2.43832e-05, 3.71298e-05, -3.55831e-05, 8.70169e-06, 6.85417e-06},
{-1.35961e-05, -5.6709e-05, -4.45989e-05, -2.3351e-05, -3.69655e-05, 0.000190876, -1.95865e-05, -1.3834e-05, 3.04587e-06, -2.18543e-05, -1.60509e-05, 3.8566e-06, -1.59568e-05, -1.92188e-05, -2.4563e-06, -7.18411e-06, 5.23961e-06, -3.82753e-07, 1.17201e-05, -8.94913e-06, 1.03397e-05, 6.63466e-06, -3.18107e-05, -6.89713e-06, -6.3317e-07, 7.21664e-06, 1.73091e-06, 4.52834e-05, -2.71892e-05, 8.51469e-06, 4.47548e-05, 1.39364e-05, -4.63575e-05, 1.39783e-05, 5.23457e-05, -4.20599e-05, -3.41902e-05, 2.39656e-05, -3.75899e-05, 3.57418e-05, -3.6647e-05, 7.64827e-05, -1.94896e-05, -1.92106e-05, -1.24018e-05, -3.57964e-05, -8.27999e-06, 9.35097e-07, 1.58773e-06, -3.3201e-06, -6.88777e-06},
{-2.8951e-05, -4.56913e-06, -1.72471e-05, -2.43579e-05, -1.93361e-05, -1.95865e-05, 0.000182745, -3.50557e-05, -2.73976e-05, -1.31446e-05, -3.12772e-06, -1.83045e-05, -1.92915e-05, -1.12805e-05, -3.17057e-06, -6.57161e-06, 3.33912e-06, -4.99148e-06, 1.90378e-05, 2.70299e-06, 1.58143e-05, -1.7681e-05, 1.08238e-05, -9.77305e-06, 6.72064e-06, 1.98433e-06, 2.98847e-07, 2.46277e-05, -3.05242e-06, -1.59714e-05, 1.19545e-05, 3.66448e-06, 1.19733e-05, 6.58864e-06, 9.47405e-06, -2.56439e-05, -3.22205e-05, 9.25331e-06, 2.28009e-05, -1.77938e-05, 1.80542e-05, -1.77771e-05, -1.5413e-05, -4.16087e-05, -4.2324e-06, 3.97745e-05, -3.05147e-05, 1.53638e-05, 1.85335e-05, 7.11922e-05, -5.72175e-05},
{-0.000115159, 9.99485e-06, -4.60799e-06, -1.18952e-05, -1.23473e-05, -1.3834e-05, -3.50557e-05, 0.000169153, -1.76505e-05, -2.81337e-05, -3.36086e-05, -2.11717e-05, -9.39719e-06, 1.18939e-07, -2.1722e-05, 7.04135e-06, 7.25282e-06, -1.74186e-05, 7.09103e-06, 1.55732e-05, 1.57288e-05, 4.73259e-06, 7.41279e-06, 2.21319e-05, -4.17293e-06, 1.05086e-05, -3.09256e-05, -2.99197e-05, -9.22262e-06, 7.16495e-06, -4.88079e-06, -2.80746e-05, 1.76957e-05, -8.67043e-07, 3.23111e-06, -7.33675e-06, 3.57331e-05, 1.26797e-05, -1.41439e-05, -1.74962e-05, 3.8883e-06, -6.99475e-05, 1.18092e-05, 4.70023e-05, 2.02875e-05, 5.35908e-05, 1.33562e-05, 9.20363e-06, -6.80977e-06, 5.0218e-07, -1.47506e-05},
{4.37729e-05, -3.94752e-06, 7.14859e-06, -1.60087e-05, -2.58623e-05, 3.04587e-06, -2.73976e-05, -1.76505e-05, 0.000184364, -4.26599e-05, -2.1239e-05, -7.39278e-06, -1.48761e-05, -5.78756e-06, -2.66538e-05, -1.57331e-05, 1.12331e-05, 3.24431e-06, 2.01295e-06, 7.67912e-06, 1.52482e-05, 1.50273e-05, -2.75276e-05, 1.59785e-05, -1.84626e-05, -5.38954e-05, 3.75933e-05, 6.46011e-07, -1.0056e-05, 7.18154e-06, 7.26395e-06, -4.49004e-06, 2.62839e-05, -4.38441e-05, -1.92121e-05, 3.87405e-05, -8.76526e-06, 1.80292e-05, 1.03667e-05, -5.99335e-05, 2.64572e-05, -1.59471e-05, 5.00802e-05, -1.90105e-05, 7.72115e-07, -1.41089e-05, 2.21541e-05, 2.76807e-05, 2.39646e-06, -1.07954e-05, 2.85081e-05},
{9.52452e-05, 2.27827e-05, 1.11382e-05, -8.70781e-06, -3.43379e-06, -2.18543e-05, -1.31446e-05, -2.81337e-05, -4.26599e-05, 0.000201438, -2.4407e-05, -2.13964e-05, -1.1149e-05, -2.43911e-05, -1.6328e-07, -1.16482e-05, 1.39289e-05, -8.47199e-06, -5.57601e-06, -1.12594e-05, -1.5865e-05, -1.17533e-05, -2.40088e-05, -6.75852e-06, -2.51291e-06, 2.3678e-05, 1.44722e-05, -3.02561e-06, 1.61188e-05, -1.30364e-05, -1.42555e-05, 3.59925e-05, 2.63542e-05, 2.89174e-06, 8.7681e-06, -2.71037e-05, 2.69741e-05, 5.12937e-06, -8.27412e-06, 1.94915e-05, -3.99272e-05, 9.32271e-06, -2.69437e-05, -3.39865e-05, -1.48649e-05, 4.61962e-06, 3.49435e-05, 5.25736e-06, 3.58668e-05, -2.74568e-05, -2.7246e-05},
{4.26502e-06, 8.06401e-06, -2.95438e-06, 7.41117e-06, -1.4417e-06, -1.60509e-05, -3.12772e-06, -3.36086e-05, -2.1239e-05, -2.4407e-05, 0.000191943, -2.67689e-05, -2.42428e-05, 9.45489e-06, -1.02016e-05, -1.0351e-05, -3.85567e-05, -5.16396e-09, -3.58952e-06, 2.69082e-06, -2.77414e-06, -3.77121e-06, 1.04149e-06, -1.70265e-05, 1.61971e-05, -1.20926e-05, -5.06989e-06, 2.20605e-05, 3.20167e-05, -1.59609e-06, -1.40111e-05, -7.37341e-06, 5.88861e-06, 3.13332e-05, -4.42135e-05, 2.29453e-05, -5.61104e-06, -2.50603e-05, 2.74788e-05, -9.35554e-06, -1.71879e-05, 2.72287e-05, 2.48134e-05, 3.42628e-06, 3.31386e-05, -2.25953e-05, -1.80298e-05, -6.12214e-05, 1.47458e-05, 2.02534e-05, -2.39596e-05},
{4.8405e-05, -5.67898e-06, -2.89264e-06, 5.40462e-06, -3.28384e-06, 3.8566e-06, -1.83045e-05, -2.11717e-05, -7.39278e-06, -2.13964e-05, -2.67689e-05, 0.000204118, -2.54669e-05, -4.00688e-05, -2.36567e-05, 1.50114e-05, -2.10849e-05, -6.06248e-06, -7.08194e-06, 1.58062e-06, -5.71151e-06, -1.78441e-05, 8.38275e-06, -1.78244e-05, 1.70877e-05, -6.77771e-07, -1.55195e-05, -1.18373e-05, -1.93539e-05, 3.82601e-05, 4.18783e-06, 9.7825e-06, 5.43101e-06, 4.7564e-05, -1.9353e-05, 1.39919e-05, 1.11658e-05, -1.35123e-05, -2.95182e-05, 3.40333e-05, -2.94257e-05, -4.05163e-05, -3.27176e-05, 2.25839e-05, 1.24799e-06, -3.55016e-05, 3.06077e-05, 6.62654e-06, 1.15545e-05, -3.9069e-05, 7.52978e-05},
{8.24648e-05, 3.71524e-05, -8.86647e-06, 1.68764e-05, 1.29683e-05, -1.59568e-05, -1.92915e-05, -9.39719e-06, -1.48761e-05, -1.1149e-05, -2.42428e-05, -2.54669e-05, 0.000201481, -2.63924e-05, -1.91794e-05, -4.68517e-05, -9.03522e-06, -2.00966e-05, -2.51891e-05, 3.07187e-06, 1.5366e-05, -6.83882e-07, 1.73913e-05, 5.16778e-06, -1.09331e-05, 2.65005e-05, 6.42981e-06, -5.8514e-06, 1.46309e-05, 2.69701e-06, -1.0088e-05, -2.11366e-06, -1.04257e-05, -4.14178e-05, 7.64532e-06, 1.83102e-05, 2.91142e-05, -1.67712e-05, 1.96373e-05, -4.53247e-05, 1.41788e-05, -3.66186e-05, -4.88455e-05, 3.44733e-05, 2.54032e-05, 2.91994e-05, 2.84064e-05, 1.0232e-05, -2.30424e-05, -3.65931e-05, 2.53522e-05},
{-9.59416e-05, -2.7878e-06, 1.87034e-06, 1.12127e-05, -9.40662e-07, -1.92188e-05, -1.12805e-05, 1.18939e-07, -5.78756e-06, -2.43911e-05, 9.45489e-06, -4.00688e-05, -2.63924e-05, 0.000246271, -7.64422e-06, -2.44347e-05, -1.93929e-05, -3.63935e-05, -1.21472e-05, 2.648e-06, -1.1462e-05, 1.23299e-06, -8.37788e-06, -3.51705e-06, -2.2314e-05, -3.15872e-05, 6.55592e-06, 9.24152e-06, 2.14989e-05, -1.57482e-05, 2.30425e-05, -7.39985e-06, -1.471e-05, 4.19753e-05, -5.04887e-06, -2.03391e-05, 5.42106e-05, 4.26132e-06, -4.78245e-06, 8.51441e-07, -4.9273e-06, 5.88168e-05, -3.77891e-05, -2.69112e-06, -2.70719e-07, -4.13033e-05, -8.99194e-06, -5.42392e-05, -2.24913e-05, -3.1003e-05, 4.0911e-05},
{-4.96514e-05, 1.25567e-05, -9.19841e-06, 1.92879e-05, 3.03201e-06, -2.4563e-06, -3.17057e-06, -2.1722e-05, -2.66538e-05, -1.6328e-07, -1.02016e-05, -2.36567e-05, -1.91794e-05, -7.64422e-06, 0.000245658, -2.65155e-05, -2.92907e-05, -4.3287e-06, -1.87371e-05, -3.24697e-05, -1.01446e-05, -6.50853e-06, -1.72434e-05, -2.00568e-06, -3.51489e-05, -6.84559e-06, 2.19644e-05, 1.38936e-05, -1.92699e-05, 9.85547e-06, -1.39883e-05, -3.65142e-06, 4.95438e-06, 6.05063e-06, 4.99058e-05, 2.84766e-05, 4.51223e-06, 3.97353e-06, -2.38957e-07, -3.2449e-05, 2.94463e-05, -4.85033e-05, 2.27698e-06, -1.90494e-05, -3.34833e-05, 3.80504e-05, -4.36868e-05, -1.99472e-05, -1.73078e-05, 2.29425e-05, 3.8822e-05},
{-0.000136782, 1.60544e-07, -1.40444e-05, 1.88575e-05, 2.11594e-07, -7.18411e-06, -6.57161e-06, 7.04135e-06, -1.57331e-05, -1.16482e-05, -1.0351e-05, 1.50114e-05, -4.68517e-05, -2.44347e-05, -2.65155e-05, 0.000268443, -4.01568e-05, -5.64715e-06, -1.96067e-05, -1.84982e-05, -2.2151e-05, -6.14885e-06, 9.93556e-06, 3.24798e-07, -1.28658e-05, -9.04035e-06, -3.40154e-05, -4.00176e-05, -1.72519e-05, 1.34976e-05, -5.58192e-06, 2.03207e-05, 9.43991e-06, 3.65751e-05, -1.12541e-05, -4.55551e-05, 3.04809e-05, -3.32642e-05, 1.64526e-05, 2.89747e-05, 7.51157e-05, -1.15251e-05, 4.75404e-05, -2.38167e-05, -2.29674e-05, -4.04882e-05, 3.52337e-05, -3.28224e-07, 2.67423e-05, -5.22183e-05, -1.03783e-05},
{-8.05224e-05, -9.85258e-05, 3.03036e-05, -1.62663e-05, 8.31762e-07, 5.23961e-06, 3.33912e-06, 7.25282e-06, 1.12331e-05, 1.39289e-05, -3.85567e-05, -2.10849e-05, -9.03522e-06, -1.93929e-05, -2.92907e-05, -4.01568e-05, 0.000274149, -1.44809e-05, -1.14346e-05, -1.32279e-05, -3.50112e-05, -5.16734e-05, -1.37171e-05, 2.95027e-05, -2.26039e-05, 1.91692e-05, -2.16584e-05, -5.69068e-06, -1.91888e-05, -1.59536e-05, -2.52466e-05, 2.50367e-05, -3.78493e-05, 1.75106e-05, 2.74563e-05, 1.14707e-05, 5.38545e-06, -2.58147e-05, 1.68373e-05, 5.09999e-05, -6.72364e-06, 3.76123e-07, 6.43991e-05, 1.89199e-05, -3.87773e-05, -2.72137e-06, -1.47301e-05, 2.36639e-06, -7.12661e-05, 2.85846e-05, -3.52277e-05},
{4.96287e-05, -2.76498e-05, 2.23746e-05, 3.31625e-06, 7.35234e-06, -3.82753e-07, -4.99148e-06, -1.74186e-05, 3.24431e-06, -8.47199e-06, -5.16396e-09, -6.06248e-06, -2.00966e-05, -3.63935e-05, -4.3287e-06, -5.64715e-06, -1.44809e-05, 0.000272374, -2.80918e-05, -2.50053e-05, -5.56334e-05, -4.43181e-05, -1.47801e-05, -1.49808e-05, 1.2417e-05, -4.82939e-05, 1.26374e-05, -1.68175e-05, 1.36662e-06, -2.40448e-05, 4.45802e-05, 6.56978e-06, 5.15393e-05, -2.97021e-05, 1.13258e-06, 7.70405e-06, -3.38137e-05, 1.58366e-05, 1.1837e-06, 9.44941e-06, 3.02957e-05, 8.91337e-06, 1.98677e-05, 5.05199e-05, 2.49054e-05, -2.37229e-05, -3.41884e-05, 1.79437e-05, -1.09977e-06, -5.09367e-05, -6.44441e-05},
{4.14411e-05, 3.15176e-05, 2.52406e-05, -1.69895e-05, -3.02049e-05, 1.17201e-05, 1.90378e-05, 7.09103e-06, 2.01295e-06, -5.57601e-06, -3.58952e-06, -7.08194e-06, -2.51891e-05, -1.21472e-05, -1.87371e-05, -1.96067e-05, -1.14346e-05, -2.80918e-05, 0.000273281, -1.936e-05, -2.71563e-05, -1.55748e-05, -4.00824e-05, -2.91043e-05, 3.32972e-05, 1.45454e-05, -4.21669e-05, -3.55594e-05, 2.31e-05, -3.36462e-05, -4.54685e-06, -4.91961e-05, -6.52424e-06, -1.54489e-06, 1.80926e-05, -4.54862e-05, 1.62556e-05, 4.9988e-05, 3.55678e-05, -3.55601e-05, 6.20217e-06, 7.89225e-06, 2.56848e-07, 5.3288e-05, 4.10883e-05, 3.21559e-05, -1.5127e-05, -3.28324e-05, 1.47764e-05, -1.17977e-05, -2.43583e-05},
{-1.24779e-05, 1.87005e-05, 5.98781e-06, 4.52451e-06, -1.6072e-05, -8.94913e-06, 2.70299e-06, 1.55732e-05, 7.67912e-06, -1.12594e-05, 2.69082e-06, 1.58062e-06, 3.07187e-06, 2.648e-06, -3.24697e-05, -1.84982e-05, -1.32279e-05, -2.50053e-05, -1.936e-05, 0.000296146, -4.39167e-05, -4.49639e-05, -3.62657e-05, -4.12625e-05, -2.57378e-05, -3.66543e-05, -6.85471e-06, 2.01434e-05, -3.18233e-05, -8.52619e-06, -2.98438e-05, 1.14783e-05, 1.28377e-05, -5.19343e-06, 2.40173e-07, 1.15438e-05, -1.81748e-05, 1.628e-05, -1.21822e-05, 5.66836e-05, 2.37294e-05, 6.02725e-05, 1.79475e-05, -4.7207e-05, -5.43678e-05, 5.70485e-05, -1.03375e-05, 3.96174e-05, -3.48245e-05, -1.64202e-05, -1.10137e-05},
{9.83332e-05, 5.64557e-05, -1.08933e-05, -1.12592e-05, -1.23793e-05, 1.03397e-05, 1.58143e-05, 1.57288e-05, 1.52482e-05, -1.5865e-05, -2.77414e-06, -5.71151e-06, 1.5366e-05, -1.1462e-05, -1.01446e-05, -2.2151e-05, -3.50112e-05, -5.56334e-05, -2.71563e-05, -4.39167e-05, 0.000331257, -6.30693e-05, -3.62338e-05, -1.64104e-05, -5.0966e-05, 1.02438e-05, -3.95282e-06, -9.67782e-06, -2.13306e-05, 7.96323e-06, -2.97501e-05, 2.54642e-05, 1.99917e-05, 6.23287e-07, -5.70096e-05, 3.41113e-05, 1.27654e-05, -8.77713e-06, -3.70525e-05, 2.84176e-05, -3.44194e-05, 3.63448e-05, -1.72516e-05, 8.75924e-05, 4.88488e-05, -1.22627e-05, -6.94926e-06, -4.48862e-05, 8.00848e-06, 2.95728e-05, -2.75493e-05},
{9.50702e-05, -1.36524e-05, -4.40666e-06, 1.49694e-05, 4.95934e-06, 6.63466e-06, -1.7681e-05, 4.73259e-06, 1.50273e-05, -1.17533e-05, -3.77121e-06, -1.78441e-05, -6.83882e-07, 1.23299e-06, -6.50853e-06, -6.14885e-06, -5.16734e-05, -4.43181e-05, -1.55748e-05, -4.49639e-05, -6.30693e-05, 0.000362536, 5.20769e-06, -3.11121e-05, -1.41059e-05, 2.3232e-05, -1.50867e-05, -1.60507e-05, -1.83226e-05, 4.07296e-05, -2.2266e-05, -6.45859e-05, -8.21405e-06, -1.91054e-05, -1.88593e-05, 8.38437e-06, -6.38233e-05, 3.88893e-05, -1.24775e-05, -3.64648e-05, 1.61797e-05, 6.5615e-05, -1.62099e-05, -1.91819e-05, 2.4708e-05, -1.73535e-05, 5.41597e-06, 6.03889e-05, 3.96869e-07, 1.81363e-05, 2.5405e-05},
{-0.000137379, 1.67219e-05, 6.18351e-06, 1.36173e-05, 1.49319e-05, -3.18107e-05, 1.08238e-05, 7.41279e-06, -2.75276e-05, -2.40088e-05, 1.04149e-06, 8.38275e-06, 1.73913e-05, -8.37788e-06, -1.72434e-05, 9.93556e-06, -1.37171e-05, -1.47801e-05, -4.00824e-05, -3.62657e-05, -3.62338e-05, 5.20769e-06, 0.000308164, -2.14936e-05, -1.69338e-05, -1.57882e-05, -1.83313e-05, 1.30853e-05, -1.12452e-05, -6.95374e-05, -2.63076e-05, 2.1927e-05, -4.22906e-05, -3.69735e-05, -4.58772e-06, 4.24247e-05, -2.45815e-05, 2.17907e-05, 2.05219e-05, -1.4247e-05, -6.66802e-06, -0.000102364, 6.79784e-05, 5.53774e-05, -4.81313e-05, 3.61512e-05, -2.3427e-05, 1.05551e-05, 3.94782e-05, 6.37169e-05, -3.74779e-05},
{5.89401e-05, 6.06075e-06, 2.50943e-06, 8.45615e-06, -1.85484e-05, -6.89713e-06, -9.77305e-06, 2.21319e-05, 1.59785e-05, -6.75852e-06, -1.70265e-05, -1.78244e-05, 5.16778e-06, -3.51705e-06, -2.00568e-06, 3.24798e-07, 2.95027e-05, -1.49808e-05, -2.91043e-05, -4.12625e-05, -1.64104e-05, -3.11121e-05, -2.14936e-05, 0.000311398, -4.3718e-05, -6.79344e-06, -2.06542e-05, -1.67096e-05, -1.30952e-05, 1.27321e-05, -4.08029e-05, -1.3496e-05, -5.25037e-05, -1.65434e-05, 3.67864e-05, -6.05523e-06, -4.37665e-05, 2.81732e-08, -2.34731e-06, -1.57953e-05, 1.46012e-05, 3.2046e-05, 5.68144e-05, -2.75795e-05, -2.71591e-06, 4.81765e-06, 9.02434e-06, -4.00405e-05, 1.06515e-05, 3.38493e-05, 3.49637e-06},
{0.000128391, -2.68488e-06, -1.71965e-05, 4.91549e-06, 1.89428e-05, -6.3317e-07, 6.72064e-06, -4.17293e-06, -1.84626e-05, -2.51291e-06, 1.61971e-05, 1.70877e-05, -1.09331e-05, -2.2314e-05, -3.51489e-05, -1.28658e-05, -2.26039e-05, 1.2417e-05, 3.32972e-05, -2.57378e-05, -5.0966e-05, -1.41059e-05, -1.69338e-05, -4.3718e-05, 0.000404181, -5.49112e-05, -3.1669e-05, -2.98687e-05, -1.73861e-05, -1.25331e-05, 1.99964e-05, -7.59504e-05, -4.04572e-05, -3.57427e-06, -2.11953e-05, 5.48494e-05, -3.66966e-05, -7.54261e-06, 6.54459e-05, 1.83489e-06, -7.03417e-05, -8.51502e-05, 4.65527e-05, -1.9402e-05, -6.66808e-06, 8.27693e-05, -2.98674e-05, 2.21796e-05, 4.16289e-05, 2.91688e-05, -6.3517e-05},
{-0.000207407, -6.93407e-05, 2.09407e-05, 7.11254e-07, 5.5515e-06, 7.21664e-06, 1.98433e-06, 1.05086e-05, -5.38954e-05, 2.3678e-05, -1.20926e-05, -6.77771e-07, 2.65005e-05, -3.15872e-05, -6.84559e-06, -9.04035e-06, 1.91692e-05, -4.82939e-05, 1.45454e-05, -3.66543e-05, 1.02438e-05, 2.3232e-05, -1.57882e-05, -6.79344e-06, -5.49112e-05, 0.000395498, -8.03679e-06, -4.62846e-05, -1.28778e-05, -3.65649e-05, -5.25931e-05, -4.70438e-05, -1.94425e-05, -9.70338e-05, -6.05623e-06, -1.99968e-05, 4.77177e-06, 5.19695e-05, -2.06415e-05, 2.68356e-05, -4.33838e-05, -3.14663e-05, -5.75999e-05, 9.11676e-07, -2.32221e-05, -1.08767e-05, -7.88163e-06, 7.91406e-05, -2.2897e-06, 0.000125404, 3.92249e-06},
{-0.000109347, 2.9713e-06, -1.04552e-06, -1.07457e-05, 1.56733e-05, 1.73091e-06, 2.98847e-07, -3.09256e-05, 3.75933e-05, 1.44722e-05, -5.06989e-06, -1.55195e-05, 6.42981e-06, 6.55592e-06, 2.19644e-05, -3.40154e-05, -2.16584e-05, 1.26374e-05, -4.21669e-05, -6.85471e-06, -3.95282e-06, -1.50867e-05, -1.83313e-05, -2.06542e-05, -3.1669e-05, -8.03679e-06, 0.000407611, -0.000111735, -3.90941e-06, -9.66967e-05, 2.00722e-05, 5.82821e-05, 2.61171e-05, 1.48553e-05, -6.30111e-05, -5.35932e-05, -2.14075e-05, -2.74576e-05, -8.42956e-06, 1.02135e-05, -5.88626e-05, 7.44351e-05, -6.26844e-05, 1.13517e-05, 6.70452e-05, -0.000111999, 6.42738e-05, -2.23856e-05, -7.7833e-06, 1.92331e-05, 3.12942e-05},
{-4.65733e-05, 1.83319e-06, 1.67884e-06, -3.13133e-05, 3.72499e-06, 4.52834e-05, 2.46277e-05, -2.99197e-05, 6.46011e-07, -3.02561e-06, 2.20605e-05, -1.18373e-05, -5.8514e-06, 9.24152e-06, 1.38936e-05, -4.00176e-05, -5.69068e-06, -1.68175e-05, -3.55594e-05, 2.01434e-05, -9.67782e-06, -1.60507e-05, 1.30853e-05, -1.67096e-05, -2.98687e-05, -4.62846e-05, -0.000111735, 0.000497386, 1.91935e-05, -5.71868e-05, -9.83573e-06, -4.31914e-05, -8.10933e-05, -6.16319e-05, -9.96793e-06, -2.34663e-05, 4.49479e-05, -5.00732e-05, -4.86356e-05, 6.81448e-05, 5.6396e-05, -2.05868e-06, -1.73941e-05, -9.54347e-05, 3.61341e-05, 1.10182e-06, 6.15458e-05, 2.8005e-05, -8.42321e-05, -7.33274e-05, -1.16119e-06},
{8.742e-05, 4.14368e-05, -3.27843e-07, 8.60804e-06, 6.68261e-06, -2.71892e-05, -3.05242e-06, -9.22262e-06, -1.0056e-05, 1.61188e-05, 3.20167e-05, -1.93539e-05, 1.46309e-05, 2.14989e-05, -1.92699e-05, -1.72519e-05, -1.91888e-05, 1.36662e-06, 2.31e-05, -3.18233e-05, -2.13306e-05, -1.83226e-05, -1.12452e-05, -1.30952e-05, -1.73861e-05, -1.28778e-05, -3.90941e-06, 1.91935e-05, 0.000472848, -8.96164e-05, -3.77777e-05, -4.3715e-05, -5.64336e-05, -2.35189e-05, -3.83532e-06, -6.05772e-05, -3.27635e-05, -6.23039e-05, 3.07955e-05, -8.0211e-06, 3.24983e-05, 5.82808e-05, -0.00010193, -6.44285e-05, -1.25073e-05, -2.07023e-05, -3.5898e-05, -4.06471e-05, 9.36283e-05, 6.04558e-05, -2.85723e-06},
{0.000136528, -5.25881e-05, 4.51441e-07, 5.8603e-07, 2.4015e-06, 8.51469e-06, -1.59714e-05, 7.16495e-06, 7.18154e-06, -1.30364e-05, -1.59609e-06, 3.82601e-05, 2.69701e-06, -1.57482e-05, 9.85547e-06, 1.34976e-05, -1.59536e-05, -2.40448e-05, -3.36462e-05, -8.52619e-06, 7.96323e-06, 4.07296e-05, -6.95374e-05, 1.27321e-05, -1.25331e-05, -3.65649e-05, -9.66967e-05, -5.71868e-05, -8.96164e-05, 0.000653022, -4.93492e-05, -7.53611e-05, -5.61713e-05, 8.38705e-05, -6.58077e-05, -5.63845e-06, -2.53981e-05, -8.04695e-06, -2.44879e-05, -1.77416e-05, -1.05194e-05, -9.54376e-06, -8.18876e-05, -1.55916e-05, 1.76824e-05, 3.60357e-05, 2.9641e-05, 9.00892e-06, -3.03229e-05, -7.26594e-05, -3.06024e-07},
{0.00013447, 5.62995e-05, -2.34602e-05, -1.13582e-05, -7.38868e-06, 4.47548e-05, 1.19545e-05, -4.88079e-06, 7.26395e-06, -1.42555e-05, -1.40111e-05, 4.18783e-06, -1.0088e-05, 2.30425e-05, -1.39883e-05, -5.58192e-06, -2.52466e-05, 4.45802e-05, -4.54685e-06, -2.98438e-05, -2.97501e-05, -2.2266e-05, -2.63076e-05, -4.08029e-05, 1.99964e-05, -5.25931e-05, 2.00722e-05, -9.83573e-06, -3.77777e-05, -4.93492e-05, 0.000584803, 3.34398e-05, -4.79584e-05, -3.04908e-05, 3.0731e-06, -6.94391e-05, -7.39889e-05, -5.02918e-05, -1.5087e-05, 1.52721e-05, -8.74744e-05, 7.92185e-05, 3.1313e-07, -2.37774e-05, -0.00015021, -1.01644e-06, 1.13458e-05, -4.73079e-06, -5.21883e-05, 5.9314e-05, -8.60012e-06},
{-0.000143867, -3.51059e-05, -3.4446e-05, -1.71656e-06, 3.428e-05, 1.39364e-05, 3.66448e-06, -2.80746e-05, -4.49004e-06, 3.59925e-05, -7.37341e-06, 9.7825e-06, -2.11366e-06, -7.39985e-06, -3.65142e-06, 2.03207e-05, 2.50367e-05, 6.56978e-06, -4.91961e-05, 1.14783e-05, 2.54642e-05, -6.45859e-05, 2.1927e-05, -1.3496e-05, -7.59504e-05, -4.70438e-05, 5.82821e-05, -4.31914e-05, -4.3715e-05, -7.53611e-05, 3.34398e-05, 0.000650636, -3.39313e-05, -0.000118697, -0.000142244, -4.42081e-05, -8.89032e-06, -8.41943e-06, 2.73693e-05, -9.61505e-05, 0.000133637, -7.94207e-05, -3.30237e-05, -1.31356e-05, -2.06478e-05, -1.43027e-05, -2.49197e-05, -7.07111e-05, 2.38401e-05, -1.39873e-06, 7.50706e-05},
{-0.000246086, 4.02116e-05, 4.48407e-05, -1.60143e-05, -3.08518e-05, -4.63575e-05, 1.19733e-05, 1.76957e-05, 2.62839e-05, 2.63542e-05, 5.88861e-06, 5.43101e-06, -1.04257e-05, -1.471e-05, 4.95438e-06, 9.43991e-06, -3.78493e-05, 5.15393e-05, -6.52424e-06, 1.28377e-05, 1.99917e-05, -8.21405e-06, -4.22906e-05, -5.25037e-05, -4.04572e-05, -1.94425e-05, 2.61171e-05, -8.10933e-05, -5.64336e-05, -5.61713e-05, -4.79584e-05, -3.39313e-05, 0.000739815, -7.6536e-05, -7.52973e-05, -6.61072e-05, -5.09688e-05, -2.32377e-05, -8.84464e-05, -2.55158e-05, -6.08384e-05, 7.20332e-05, 3.42818e-05, 2.20603e-05, 7.38229e-05, 2.9178e-05, -5.19102e-05, 5.4267e-05, 1.18601e-05, 2.51754e-05, 7.66073e-06},
{0.000224706, -0.000104181, 5.14078e-06, -6.68248e-06, 1.8012e-06, 1.39783e-05, 6.58864e-06, -8.67043e-07, -4.38441e-05, 2.89174e-06, 3.13332e-05, 4.7564e-05, -4.14178e-05, 4.19753e-05, 6.05063e-06, 3.65751e-05, 1.75106e-05, -2.97021e-05, -1.54489e-06, -5.19343e-06, 6.23287e-07, -1.91054e-05, -3.69735e-05, -1.65434e-05, -3.57427e-06, -9.70338e-05, 1.48553e-05, -6.16319e-05, -2.35189e-05, 8.38705e-05, -3.04908e-05, -0.000118697, -7.6536e-05, 0.000706136, -6.82204e-05, -4.46964e-05, -3.08758e-05, -0.000149051, -5.38115e-05, -8.60171e-05, 1.31998e-05, 4.52011e-05, -9.01346e-07, 9.16821e-06, 4.18703e-05, 4.1254e-05, -9.77148e-05, 1.77184e-05, -4.02713e-05, 5.15161e-05, 0.000118307},
{5.70979e-05, 2.85296e-05, -4.88281e-05, 2.16778e-06, 5.87077e-06, 5.23457e-05, 9.47405e-06, 3.23111e-06, -1.92121e-05, 8.7681e-06, -4.42135e-05, -1.9353e-05, 7.64532e-06, -5.04887e-06, 4.99058e-05, -1.12541e-05, 2.74563e-05, 1.13258e-06, 1.80926e-05, 2.40173e-07, -5.70096e-05, -1.88593e-05, -4.58772e-06, 3.67864e-05, -2.11953e-05, -6.05623e-06, -6.30111e-05, -9.96793e-06, -3.83532e-06, -6.58077e-05, 3.0731e-06, -0.000142244, -7.52973e-05, -6.82204e-05, 0.000752103, -0.000168384, -2.82992e-05, -3.17981e-05, 1.56367e-05, -0.000129578, 7.34532e-06, -9.36306e-05, 6.37639e-07, 5.63221e-05, -3.47309e-05, 1.80646e-05, 8.71373e-05, 7.6326e-06, -4.10849e-05, -2.85462e-05, 9.96302e-05},
{0.000138589, 1.53898e-05, 5.07297e-07, 9.64015e-06, 2.66374e-05, -4.20599e-05, -2.56439e-05, -7.33675e-06, 3.87405e-05, -2.71037e-05, 2.29453e-05, 1.39919e-05, 1.83102e-05, -2.03391e-05, 2.84766e-05, -4.55551e-05, 1.14707e-05, 7.70405e-06, -4.54862e-05, 1.15438e-05, 3.41113e-05, 8.38437e-06, 4.24247e-05, -6.05523e-06, 5.48494e-05, -1.99968e-05, -5.35932e-05, -2.34663e-05, -6.05772e-05, -5.63845e-06, -6.94391e-05, -4.42081e-05, -6.61072e-05, -4.46964e-05, -0.000168384, 0.000662503, -0.000101439, -8.44116e-05, -1.77296e-05, -4.65132e-06, -8.53603e-05, -7.8432e-05, 1.03366e-05, -2.60476e-05, -1.60437e-05, -3.82555e-05, 0.000208747, -7.96364e-05, 2.89118e-05, 8.76021e-05, -2.748e-07},
{-0.000199036, 7.20497e-05, 4.83132e-07, -2.42714e-06, -2.76056e-06, -3.41902e-05, -3.22205e-05, 3.57331e-05, -8.76526e-06, 2.69741e-05, -5.61104e-06, 1.11658e-05, 2.91142e-05, 5.42106e-05, 4.51223e-06, 3.04809e-05, 5.38545e-06, -3.38137e-05, 1.62556e-05, -1.81748e-05, 1.27654e-05, -6.38233e-05, -2.45815e-05, -4.37665e-05, -3.66966e-05, 4.77177e-06, -2.14075e-05, 4.49479e-05, -3.27635e-05, -2.53981e-05, -7.39889e-05, -8.89032e-06, -5.09688e-05, -3.08758e-05, -2.82992e-05, -0.000101439, 0.000723819, -2.86586e-05, -0.0001706, -4.40645e-05, -7.98947e-05, -6.16325e-05, 6.10066e-05, 6.44661e-06, -0.000108668, -2.62896e-05, 2.77363e-09, 5.19336e-05, 0.000103895, -0.000153928, -7.53714e-06},
{-4.6457e-05, -5.57023e-05, -3.09784e-05, 1.82499e-05, -1.48868e-05, 2.39656e-05, 9.25331e-06, 1.26797e-05, 1.80292e-05, 5.12937e-06, -2.50603e-05, -1.35123e-05, -1.67712e-05, 4.26132e-06, 3.97353e-06, -3.32642e-05, -2.58147e-05, 1.58366e-05, 4.9988e-05, 1.628e-05, -8.77713e-06, 3.88893e-05, 2.17907e-05, 2.81732e-08, -7.54261e-06, 5.19695e-05, -2.74576e-05, -5.00732e-05, -6.23039e-05, -8.04695e-06, -5.02918e-05, -8.41943e-06, -2.32377e-05, -0.000149051, -3.17981e-05, -8.44116e-05, -2.86586e-05, 0.000761575, -0.000136581, -9.03786e-05, -4.94944e-05, -0.000150109, -3.67547e-05, -6.37717e-05, 6.71884e-05, -1.94696e-05, -2.60602e-05, -6.30476e-06, -1.2354e-05, -2.82305e-05, 3.49805e-05},
{3.28348e-05, 1.92546e-05, 1.40378e-05, 2.64127e-05, -2.69686e-05, -3.75899e-05, 2.28009e-05, -1.41439e-05, 1.03667e-05, -8.27412e-06, 2.74788e-05, -2.95182e-05, 1.96373e-05, -4.78245e-06, -2.38957e-07, 1.64526e-05, 1.68373e-05, 1.1837e-06, 3.55678e-05, -1.21822e-05, -3.70525e-05, -1.24775e-05, 2.05219e-05, -2.34731e-06, 6.54459e-05, -2.06415e-05, -8.42956e-06, -4.86356e-05, 3.07955e-05, -2.44879e-05, -1.5087e-05, 2.73693e-05, -8.84464e-05, -5.38115e-05, 1.56367e-05, -1.77296e-05, -0.0001706, -0.000136581, 0.000867691, 3.55748e-06, -0.000173505, -0.000127591, -8.76251e-05, 1.55289e-05, -5.49429e-05, -2.67585e-05, -4.75479e-05, -7.85904e-05, 0.00010545, -0.000192559, -1.19018e-05},
{-8.32483e-05, -2.29e-05, 5.36181e-05, -7.90122e-06, -3.02662e-06, 3.57418e-05, -1.77938e-05, -1.74962e-05, -5.99335e-05, 1.94915e-05, -9.35554e-06, 3.40333e-05, -4.53247e-05, 8.51441e-07, -3.2449e-05, 2.89747e-05, 5.09999e-05, 9.44941e-06, -3.55601e-05, 5.66836e-05, 2.84176e-05, -3.64648e-05, -1.4247e-05, -1.57953e-05, 1.83489e-06, 2.68356e-05, 1.02135e-05, 6.81448e-05, -8.0211e-06, -1.77416e-05, 1.52721e-05, -9.61505e-05, -2.55158e-05, -8.60171e-05, -0.000129578, -4.65132e-06, -4.40645e-05, -9.03786e-05, 3.55748e-06, 0.000795169, -0.000101928, -0.000121483, -0.000116769, -0.000117642, 1.84039e-05, -9.91049e-05, -0.000183181, 8.42376e-05, 2.4464e-05, -2.1807e-05, 2.56706e-05},
{-6.07291e-05, 0.000130883, -3.45689e-05, -1.07072e-05, 2.73859e-05, -3.6647e-05, 1.80542e-05, 3.8883e-06, 2.64572e-05, -3.99272e-05, -1.71879e-05, -2.94257e-05, 1.41788e-05, -4.9273e-06, 2.94463e-05, 7.51157e-05, -6.72364e-06, 3.02957e-05, 6.20217e-06, 2.37294e-05, -3.44194e-05, 1.61797e-05, -6.66802e-06, 1.46012e-05, -7.03417e-05, -4.33838e-05, -5.88626e-05, 5.6396e-05, 3.24983e-05, -1.05194e-05, -8.74744e-05, 0.000133637, -6.08384e-05, 1.31998e-05, 7.34532e-06, -8.53603e-05, -7.98947e-05, -4.94944e-05, -0.000173505, -0.000101928, 0.000838617, -1.02609e-05, -7.02138e-05, 9.81122e-05, -6.25961e-05, 2.29994e-05, -4.27616e-05, -5.08174e-05, -6.78611e-05, -0.000138051, -8.78255e-05},
{0.000379313, -9.53455e-05, -4.88781e-05, 2.28551e-05, 2.39657e-05, 7.64827e-05, -1.77771e-05, -6.99475e-05, -1.59471e-05, 9.32271e-06, 2.72287e-05, -4.05163e-05, -3.66186e-05, 5.88168e-05, -4.85033e-05, -1.15251e-05, 3.76123e-07, 8.91337e-06, 7.89225e-06, 6.02725e-05, 3.63448e-05, 6.5615e-05, -0.000102364, 3.2046e-05, -8.51502e-05, -3.14663e-05, 7.44351e-05, -2.05868e-06, 5.82808e-05, -9.54376e-06, 7.92185e-05, -7.94207e-05, 7.20332e-05, 4.52011e-05, -9.36306e-05, -7.8432e-05, -6.16325e-05, -0.000150109, -0.000127591, -0.000121483, -1.02609e-05, 0.001053649, -1.34393e-05, -0.000102496, 9.35871e-05, -7.35271e-05, -6.09033e-05, -3.77308e-05, -0.000253417, 1.70948e-05, -0.000115023},
{0.000308066, -2.12528e-05, 4.99183e-05, -2.77283e-05, -6.01345e-06, -1.94896e-05, -1.5413e-05, 1.18092e-05, 5.00802e-05, -2.69437e-05, 2.48134e-05, -3.27176e-05, -4.88455e-05, -3.77891e-05, 2.27698e-06, 4.75404e-05, 6.43991e-05, 1.98677e-05, 2.56848e-07, 1.79475e-05, -1.72516e-05, -1.62099e-05, 6.79784e-05, 5.68144e-05, 4.65527e-05, -5.75999e-05, -6.26844e-05, -1.73941e-05, -0.00010193, -8.18876e-05, 3.1313e-07, -3.30237e-05, 3.42818e-05, -9.01346e-07, 6.37639e-07, 1.03366e-05, 6.10066e-05, -3.67547e-05, -8.76251e-05, -0.000116769, -7.02138e-05, -1.34393e-05, 0.001417483, -8.52384e-05, -8.56471e-05, -0.000221356, -5.11963e-07, -7.12241e-05, -0.000223964, -3.08667e-05, -6.50722e-05},
{0.000371455, -4.31274e-05, 1.14867e-05, 8.72761e-06, 4.39415e-06, -1.92106e-05, -4.16087e-05, 4.70023e-05, -1.90105e-05, -3.39865e-05, 3.42628e-06, 2.25839e-05, 3.44733e-05, -2.69112e-06, -1.90494e-05, -2.38167e-05, 1.89199e-05, 5.05199e-05, 5.3288e-05, -4.7207e-05, 8.75924e-05, -1.91819e-05, 5.53774e-05, -2.75795e-05, -1.9402e-05, 9.11676e-07, 1.13517e-05, -9.54347e-05, -6.44285e-05, -1.55916e-05, -2.37774e-05, -1.31356e-05, 2.20603e-05, 9.16821e-06, 5.63221e-05, -2.60476e-05, 6.44661e-06, -6.37717e-05, 1.55289e-05, -0.000117642, 9.81122e-05, -0.000102496, -8.52384e-05, 0.001494416, -0.000149167, -7.30629e-05, -0.000449775, -0.000164758, -6.51556e-05, -5.01632e-05, -0.000125316},
{-0.000798599, 2.42076e-05, 2.10336e-05, 2.71738e-05, -5.59135e-05, -1.24018e-05, -4.2324e-06, 2.02875e-05, 7.72115e-07, -1.48649e-05, 3.31386e-05, 1.24799e-06, 2.54032e-05, -2.70719e-07, -3.34833e-05, -2.29674e-05, -3.87773e-05, 2.49054e-05, 4.10883e-05, -5.43678e-05, 4.88488e-05, 2.4708e-05, -4.81313e-05, -2.71591e-06, -6.66808e-06, -2.32221e-05, 6.70452e-05, 3.61341e-05, -1.25073e-05, 1.76824e-05, -0.00015021, -2.06478e-05, 7.38229e-05, 4.18703e-05, -3.47309e-05, -1.60437e-05, -0.000108668, 6.71884e-05, -5.49429e-05, 1.84039e-05, -6.25961e-05, 9.35871e-05, -8.56471e-05, -0.000149167, 0.001434163, -0.000155613, -0.0002642, -3.10648e-05, -0.000111317, -0.000144208, -0.000236957},
{4.15331e-05, -4.6987e-05, 2.27733e-05, -1.27263e-05, -2.50076e-05, -3.57964e-05, 3.97745e-05, 5.35908e-05, -1.41089e-05, 4.61962e-06, -2.25953e-05, -3.55016e-05, 2.91994e-05, -4.13033e-05, 3.80504e-05, -4.04882e-05, -2.72137e-06, -2.37229e-05, 3.21559e-05, 5.70485e-05, -1.22627e-05, -1.73535e-05, 3.61512e-05, 4.81765e-06, 8.27693e-05, -1.08767e-05, -0.000111999, 1.10182e-06, -2.07023e-05, 3.60357e-05, -1.01644e-06, -1.43027e-05, 2.9178e-05, 4.1254e-05, 1.80646e-05, -3.82555e-05, -2.62896e-05, -1.94696e-05, -2.67585e-05, -9.91049e-05, 2.29994e-05, -7.35271e-05, -0.000221356, -7.30629e-05, -0.000155613, 0.001425162, -1.10112e-05, -8.78351e-05, -0.000214027, -9.89777e-05, -0.000248322},
{4.46331e-05, 3.83722e-05, 1.93506e-05, -7.13011e-06, -2.43832e-05, -8.27999e-06, -3.05147e-05, 1.33562e-05, 2.21541e-05, 3.49435e-05, -1.80298e-05, 3.06077e-05, 2.84064e-05, -8.99194e-06, -4.36868e-05, 3.52337e-05, -1.47301e-05, -3.41884e-05, -1.5127e-05, -1.03375e-05, -6.94926e-06, 5.41597e-06, -2.3427e-05, 9.02434e-06, -2.98674e-05, -7.88163e-06, 6.42738e-05, 6.15458e-05, -3.5898e-05, 2.9641e-05, 1.13458e-05, -2.49197e-05, -5.19102e-05, -9.77148e-05, 8.71373e-05, 0.000208747, 2.77363e-09, -2.60602e-05, -4.75479e-05, -0.000183181, -4.27616e-05, -6.09033e-05, -5.11963e-07, -0.000449775, -0.0002642, -1.10112e-05, 0.001535232, -0.000237284, -9.62566e-05, -3.37707e-05, -0.000118864},
{-0.000273535, 5.05434e-05, 3.50583e-05, -7.39945e-05, 3.71298e-05, 9.35097e-07, 1.53638e-05, 9.20363e-06, 2.76807e-05, 5.25736e-06, -6.12214e-05, 6.62654e-06, 1.0232e-05, -5.42392e-05, -1.99472e-05, -3.28224e-07, 2.36639e-06, 1.79437e-05, -3.28324e-05, 3.96174e-05, -4.48862e-05, 6.03889e-05, 1.05551e-05, -4.00405e-05, 2.21796e-05, 7.91406e-05, -2.23856e-05, 2.8005e-05, -4.06471e-05, 9.00892e-06, -4.73079e-06, -7.07111e-05, 5.4267e-05, 1.77184e-05, 7.6326e-06, -7.96364e-05, 5.19336e-05, -6.30476e-06, -7.85904e-05, 8.42376e-05, -5.08174e-05, -3.77308e-05, -7.12241e-05, -0.000164758, -3.10648e-05, -8.78351e-05, -0.000237284, 0.001222156, -0.000137115, -0.000196991, -0.000105686},
{-0.000210617, -2.59119e-05, -3.65015e-05, 4.12342e-05, -3.55831e-05, 1.58773e-06, 1.85335e-05, -6.80977e-06, 2.39646e-06, 3.58668e-05, 1.47458e-05, 1.15545e-05, -2.30424e-05, -2.24913e-05, -1.73078e-05, 2.67423e-05, -7.12661e-05, -1.09977e-06, 1.47764e-05, -3.48245e-05, 8.00848e-06, 3.96869e-07, 3.94782e-05, 1.06515e-05, 4.16289e-05, -2.2897e-06, -7.7833e-06, -8.42321e-05, 9.36283e-05, -3.03229e-05, -5.21883e-05, 2.38401e-05, 1.18601e-05, -4.02713e-05, -4.10849e-05, 2.89118e-05, 0.000103895, -1.2354e-05, 0.00010545, 2.4464e-05, -6.78611e-05, -0.000253417, -0.000223964, -6.51556e-05, -0.000111317, -0.000214027, -9.62566e-05, -0.000137115, 0.001624962, -0.000223217, -0.000275659},
{0.000442616, -5.19021e-05, -2.10393e-05, -9.94756e-06, 8.70169e-06, -3.3201e-06, 7.11922e-05, 5.0218e-07, -1.07954e-05, -2.74568e-05, 2.02534e-05, -3.9069e-05, -3.65931e-05, -3.1003e-05, 2.29425e-05, -5.22183e-05, 2.85846e-05, -5.09367e-05, -1.17977e-05, -1.64202e-05, 2.95728e-05, 1.81363e-05, 6.37169e-05, 3.38493e-05, 2.91688e-05, 0.000125404, 1.92331e-05, -7.33274e-05, 6.04558e-05, -7.26594e-05, 5.9314e-05, -1.39873e-06, 2.51754e-05, 5.15161e-05, -2.85462e-05, 8.76021e-05, -0.000153928, -2.82305e-05, -0.000192559, -2.1807e-05, -0.000138051, 1.70948e-05, -3.08667e-05, -5.01632e-05, -0.000144208, -9.89777e-05, -3.37707e-05, -0.000196991, -0.000223217, 0.001634523, -0.000231996},
{-0.000256442, 3.98843e-05, 1.19656e-05, 6.76566e-06, 6.85417e-06, -6.88777e-06, -5.72175e-05, -1.47506e-05, 2.85081e-05, -2.7246e-05, -2.39596e-05, 7.52978e-05, 2.53522e-05, 4.0911e-05, 3.8822e-05, -1.03783e-05, -3.52277e-05, -6.44441e-05, -2.43583e-05, -1.10137e-05, -2.75493e-05, 2.5405e-05, -3.74779e-05, 3.49637e-06, -6.3517e-05, 3.92249e-06, 3.12942e-05, -1.16119e-06, -2.85723e-06, -3.06024e-07, -8.60012e-06, 7.50706e-05, 7.66073e-06, 0.000118307, 9.96302e-05, -2.748e-07, -7.53714e-06, 3.49805e-05, -1.19018e-05, 2.56706e-05, -8.78255e-05, -0.000115023, -6.50722e-05, -0.000125316, -0.000236957, -0.000248322, -0.000118864, -0.000105686, -0.000275659, -0.000231996, 0.001515264}
};

static double wagi[5][51] = {
{5.00709e-05, 0.000366214, 0.000985787, 0.00914913, 0.0643969, 0.160703, 0.229911, 0.255508, 0.267845, 0.277876, 0.282614, 0.276661, 0.29006, 0.281068, 0.288378, 0.271844, 0.282618, 0.279982, 0.295438, 0.292009, 0.296905, 0.286422, 0.29334, 0.281695, 0.277097, 0.27389, 0.284714, 0.253069, 0.284654, 0.294046, 0.284418, 0.285836, 0.259134, 0.236962, 0.222328, 0.166245, 0.146801, 0.158464, 0.118372, 0.103481, 0.0978826, 0.0856827, 0.0804085, 0.05031, 0.0827925, 0.057825, 0.0594655, 0.0733733, 0.0498852, 0.037827, 0.0271847},
{0, 0.00174652, 0.00573062, 0.0235582, 0.0816228, 0.17198, 0.235807, 0.265726, 0.275495, 0.27845, 0.278783, 0.278996, 0.279714, 0.284443, 0.282091, 0.28303, 0.278696, 0.288329, 0.283056, 0.291278, 0.294685, 0.285972, 0.289901, 0.291892, 0.285307, 0.298254, 0.282205, 0.264636, 0.277728, 0.264418, 0.257743, 0.280279, 0.251591, 0.268913, 0.238739, 0.194142, 0.158025, 0.13985, 0.113048, 0.0835846, 0.0860063, 0.0808043, 0.0586006, 0.079863, 0.0622959, 0.0623802, 0.0500734, 0.0651953, 0.0413203, 0.0311915, 0.0207489},
{0, 0.00541418, 0.00693397, 0.0177455, 0.0568918, 0.120846, 0.194671, 0.238473, 0.262165, 0.279815, 0.26836, 0.280639, 0.273256, 0.280427, 0.264726, 0.274111, 0.276162, 0.277796, 0.266221, 0.271817, 0.281908, 0.26488, 0.268081, 0.253919, 0.278849, 0.271646, 0.263844, 0.263598, 0.271344, 0.275889, 0.263377, 0.241991, 0.266467, 0.260049, 0.281903, 0.259528, 0.264781, 0.244003, 0.236043, 0.235499, 0.254184, 0.240617, 0.217183, 0.218809, 0.192143, 0.194097, 0.178495, 0.203401, 0.212731, 0.171704, 0.0910499},
{0, 0.00613405, 0.015021, 0.0308908, 0.0623761, 0.102794, 0.151986, 0.189071, 0.225743, 0.246209, 0.259201, 0.266909, 0.277992, 0.279908, 0.276735, 0.28221, 0.285058, 0.289423, 0.285085, 0.286032, 0.284712, 0.280237, 0.284884, 0.281194, 0.280648, 0.280521, 0.278723, 0.261669, 0.274709, 0.277539, 0.258975, 0.254524, 0.247645, 0.250731, 0.249253, 0.242187, 0.230224, 0.232708, 0.229344, 0.223792, 0.22145, 0.214956, 0.218249, 0.199852, 0.198381, 0.221667, 0.199611, 0.202619, 0.19482, 0.188906, 0.11155},
{0, 0, 0, 0, 0, 0, 1, 0.109076, 0.22495, 0.252944, 0.293526, 0.269803, 0.268593, 0.262123, 0.295095, 0.298222, 0.285676, 0.300638, 0.289924, 0.299128, 0.299028, 0.309271, 0.302556, 0.30564, 0.31022, 0.314945, 0.316786, 0.308994, 0.323376, 0.325566, 0.320721, 0.328176, 0.316279, 0.319579, 0.332979, 0.32596, 0.326457, 0.320243, 0.313269, 0.313701, 0.305215, 0.30412, 0.306179, 0.277683, 0.285149, 0.27625, 0.263819, 0.268088, 0.245121, 0.241333, 0.160561}
};

const double mb_nce_bg [51] = {1.468132e+01,  1.751390e+02,  5.432570e+02,  6.452100e+02,  6.288804e+02,  5.888983e+02,  5.721696e+02,  5.599448e+02,  5.320706e+02,  5.089781e+02,  4.717286e+02,  4.501799e+02,  4.296122e+02,  3.620738e+02,  3.576632e+02,  3.276326e+02,  3.137620e+02,  2.857147e+02,  2.655224e+02,  2.647242e+02,  2.455104e+02,  2.285900e+02,  2.313580e+02,  2.262919e+02,  2.283181e+02,  2.124414e+02,  1.952493e+02,  1.962794e+02,  1.876680e+02,  2.057413e+02,  2.050090e+02,  1.850047e+02,  1.591487e+02,  1.803423e+02,  1.927868e+02,  2.144821e+02,  2.260172e+02,  2.518160e+02,  2.556633e+02,  2.555962e+02,  2.480172e+02,  2.373849e+02,  1.840215e+02,  1.800383e+02,  1.857559e+02,  2.038884e+02,  2.060577e+02,  2.033785e+02,  2.086383e+02,  2.217605e+02,  2.130483e+02};
const double mb_nce_data[51] = {4.000000e+01, 4.950000e+02, 1.902000e+03, 2.722000e+03, 2.908000e+03, 2.878000e+03, 2.934000e+03, 3.003000e+03, 2.892000e+03, 2.720000e+03, 2.651000e+03, 2.656000e+03, 2.633000e+03, 2.409000e+03, 2.364000e+03, 2.250000e+03, 2.104000e+03, 2.070000e+03, 1.929000e+03, 1.832000e+03, 1.762000e+03, 1.618000e+03, 1.546000e+03, 1.458000e+03, 1.371000e+03, 1.203000e+03, 1.174000e+03, 1.055000e+03, 9.930000e+02, 9.080000e+02, 8.320000e+02, 7.570000e+02, 6.840000e+02, 6.330000e+02, 6.520000e+02, 7.090000e+02, 7.010000e+02, 6.710000e+02, 5.890000e+02, 5.110000e+02, 4.890000e+02, 3.920000e+02, 3.590000e+02, 3.320000e+02, 3.450000e+02, 3.430000e+02, 3.560000e+02, 3.410000e+02, 3.150000e+02, 3.240000e+02, 3.360000e+02};

const double mb_nce_bg_he_sp[30] = {6.544238e+01, 6.131658e+01, 6.802268e+01, 6.614949e+01, 5.759906e+01, 5.948360e+01, 6.494799e+01, 7.288111e+01, 7.266334e+01, 8.067760e+01, 7.378685e+01, 7.559236e+01, 6.692164e+01, 5.954329e+01, 5.842109e+01, 6.112299e+01, 5.521307e+01, 5.631305e+01, 5.625658e+01, 5.765196e+01, 6.470708e+01, 6.157652e+01, 6.896612e+01, 5.997896e+01, 7.765791e+01, 6.953401e+01, 7.666214e+01, 8.304232e+01, 9.066559e+01, 1.028902e+02};
const double mb_nce_he_sp[30] = {4.600000e+02, 4.570000e+02, 4.390000e+02, 4.310000e+02, 3.800000e+02, 3.830000e+02, 3.530000e+02, 3.630000e+02, 3.710000e+02, 3.960000e+02, 2.750000e+02, 2.760000e+02, 2.220000e+02, 2.030000e+02, 2.170000e+02, 2.240000e+02, 2.000000e+02, 1.840000e+02, 1.970000e+02, 1.590000e+02, 1.550000e+02, 1.620000e+02, 1.370000e+02, 1.520000e+02, 1.260000e+02, 1.180000e+02, 1.330000e+02, 1.170000e+02, 1.380000e+02, 1.230000e+02};

const double mb_nce_bg_he_all[30] = {4.089986e+02, 3.728637e+02, 3.529947e+02, 3.474604e+02, 3.084260e+02, 2.939212e+02, 3.071228e+02, 3.304660e+02, 3.354977e+02, 3.431949e+02, 3.271199e+02, 3.149625e+02, 2.886637e+02, 2.234001e+02, 2.394022e+02, 2.461791e+02, 2.466851e+02, 2.378370e+02, 2.380596e+02, 2.461632e+02, 2.536219e+02, 2.657772e+02, 2.737242e+02, 2.553884e+02, 2.966817e+02, 2.959479e+02, 3.014583e+02, 3.316176e+02, 3.625725e+02, 3.955619e+02};
const double mb_nce_he_all[30] = {2.280000e+03, 1.972000e+03, 1.815000e+03, 1.516000e+03, 1.356000e+03, 1.163000e+03, 1.154000e+03, 1.136000e+03, 1.150000e+03, 9.980000e+02, 7.680000e+02, 7.140000e+02, 5.730000e+02, 5.010000e+02, 4.520000e+02, 4.960000e+02, 4.620000e+02, 4.150000e+02, 4.040000e+02, 4.060000e+02, 3.950000e+02, 3.930000e+02, 3.910000e+02, 4.100000e+02, 3.750000e+02, 3.760000e+02, 4.230000e+02, 4.230000e+02, 4.390000e+02, 4.380000e+02};

const double M_he [30][30] = {
{1.206371e-02, 1.069998e-02, 1.000041e-02, 9.629539e-03, 9.539918e-03, 9.565212e-03, 9.531436e-03, 9.890880e-03, 9.812638e-03, 9.888926e-03, 9.814626e-03, 9.479057e-03, 9.095917e-03, 8.458753e-03, 8.023143e-03, 7.180851e-03, 6.316530e-03, 5.396447e-03, 4.518830e-03, 3.736281e-03, 2.949985e-03, 2.166038e-03, 1.543153e-03, 6.044698e-04, 8.981331e-05, -1.538775e-04, -5.413937e-04, -9.400317e-04, -1.136515e-03, -1.485396e-03}, 
{1.069998e-02, 1.004124e-02, 9.475049e-03, 9.337847e-03, 9.296995e-03, 9.449259e-03, 9.780269e-03, 9.535039e-03, 9.683942e-03, 9.658871e-03, 9.400156e-03, 9.294999e-03, 8.770377e-03, 8.572908e-03, 7.767952e-03, 7.190231e-03, 6.617065e-03, 6.023973e-03, 5.353260e-03, 4.210453e-03, 3.540511e-03, 2.478098e-03, 1.856306e-03, 1.491926e-03, 9.435929e-04, 3.003647e-04, -4.095794e-05, -2.625562e-04, -7.525252e-04, -8.452042e-04}, 
{1.000041e-02, 9.475049e-03, 9.311699e-03, 9.156919e-03, 9.134467e-03, 9.363105e-03, 9.890488e-03, 9.336832e-03, 9.573611e-03, 9.557580e-03, 9.179387e-03, 9.226143e-03, 8.580059e-03, 8.627177e-03, 7.645219e-03, 7.150819e-03, 6.768508e-03, 6.351317e-03, 5.804593e-03, 4.433567e-03, 3.890020e-03, 2.575071e-03, 2.015742e-03, 1.886699e-03, 1.357885e-03, 5.467148e-04, 2.515424e-04, 1.032703e-04, -5.573032e-04, -4.616821e-04}, 
{9.629539e-03, 9.337847e-03, 9.156919e-03, 9.450951e-03, 9.300223e-03, 9.629718e-03, 1.041162e-02, 9.493558e-03, 9.880628e-03, 9.823312e-03, 9.320631e-03, 9.509662e-03, 8.741478e-03, 9.047659e-03, 7.801266e-03, 7.427642e-03, 7.195474e-03, 6.949310e-03, 6.495255e-03, 4.849579e-03, 4.362820e-03, 2.816633e-03, 2.255699e-03, 2.466039e-03, 1.914370e-03, 8.541227e-04, 5.884735e-04, 5.361497e-04, -3.410147e-04, -1.137496e-04}, 
{9.539918e-03, 9.296995e-03, 9.134467e-03, 9.300223e-03, 9.551854e-03, 9.698937e-03, 1.050331e-02, 9.558200e-03, 1.001501e-02, 9.858181e-03, 9.390601e-03, 9.528456e-03, 8.807038e-03, 9.100351e-03, 7.809130e-03, 7.539927e-03, 7.230225e-03, 6.985408e-03, 6.540526e-03, 4.908083e-03, 4.315294e-03, 2.925158e-03, 2.288335e-03, 2.624285e-03, 2.024329e-03, 8.977604e-04, 6.103385e-04, 5.953450e-04, -2.887999e-04, -1.792823e-04}, 
{9.565212e-03, 9.449259e-03, 9.363105e-03, 9.629718e-03, 9.698937e-03, 1.036074e-02, 1.119170e-02, 9.941306e-03, 1.049861e-02, 1.038943e-02, 9.779403e-03, 1.008240e-02, 9.189442e-03, 9.721484e-03, 8.189445e-03, 7.922299e-03, 7.783264e-03, 7.671817e-03, 7.292865e-03, 5.346795e-03, 4.874975e-03, 3.120347e-03, 2.535671e-03, 3.081196e-03, 2.494403e-03, 1.184355e-03, 9.358139e-04, 9.772102e-04, -1.298339e-04, 1.585546e-04}, 
{9.531436e-03, 9.780269e-03, 9.890488e-03, 1.041162e-02, 1.050331e-02, 1.119170e-02, 1.308458e-02, 1.077219e-02, 1.163399e-02, 1.155599e-02, 1.059523e-02, 1.130479e-02, 1.003386e-02, 1.115961e-02, 8.995298e-03, 8.826220e-03, 9.102942e-03, 9.332440e-03, 9.116808e-03, 6.433954e-03, 6.183370e-03, 3.636953e-03, 3.101503e-03, 4.281853e-03, 3.697824e-03, 1.848189e-03, 1.692185e-03, 1.898850e-03, 2.493618e-04, 9.608016e-04}, 
{9.890880e-03, 9.535039e-03, 9.336832e-03, 9.493558e-03, 9.558200e-03, 9.941306e-03, 1.077219e-02, 1.011809e-02, 1.032563e-02, 1.027048e-02, 9.793481e-03, 9.946953e-03, 9.160004e-03, 9.444436e-03, 8.154163e-03, 7.762886e-03, 7.449055e-03, 7.171483e-03, 6.699920e-03, 4.991145e-03, 4.462430e-03, 2.923034e-03, 2.351793e-03, 2.576774e-03, 2.015731e-03, 9.605867e-04, 6.825717e-04, 6.295132e-04, -2.998752e-04, -2.003970e-04}, 
{9.812638e-03, 9.683942e-03, 9.573611e-03, 9.880628e-03, 1.001501e-02, 1.049861e-02, 1.163399e-02, 1.032563e-02, 1.126508e-02, 1.077760e-02, 1.020616e-02, 1.045446e-02, 9.620326e-03, 1.013960e-02, 8.479301e-03, 8.340837e-03, 8.080938e-03, 7.951718e-03, 7.531300e-03, 5.577977e-03, 4.926539e-03, 3.363898e-03, 2.632607e-03, 3.407568e-03, 2.742209e-03, 1.291029e-03, 1.003741e-03, 1.093593e-03, -6.815230e-05, 1.613219e-05}, 
{9.888926e-03, 9.658871e-03, 9.557580e-03, 9.823312e-03, 9.858181e-03, 1.038943e-02, 1.155599e-02, 1.027048e-02, 1.077760e-02, 1.118960e-02, 1.018856e-02, 1.065627e-02, 9.561927e-03, 1.020130e-02, 8.621228e-03, 8.121004e-03, 8.149276e-03, 8.063100e-03, 7.701885e-03, 5.514694e-03, 5.285542e-03, 3.041134e-03, 2.643639e-03, 3.029214e-03, 2.563534e-03, 1.321605e-03, 1.145040e-03, 1.125382e-03, -1.302806e-04, 3.606202e-04}, 
{9.814626e-03, 9.400156e-03, 9.179387e-03, 9.320631e-03, 9.390601e-03, 9.779403e-03, 1.059523e-02, 9.793481e-03, 1.020616e-02, 1.018856e-02, 1.006423e-02, 9.886933e-03, 9.111829e-03, 9.372399e-03, 8.128266e-03, 7.710559e-03, 7.382721e-03, 7.091429e-03, 6.619684e-03, 4.933251e-03, 4.420442e-03, 2.887384e-03, 2.340028e-03, 2.507412e-03, 1.968671e-03, 9.650146e-04, 6.912607e-04, 6.195130e-04, -2.973148e-04, -2.272584e-04}, 
{9.479057e-03, 9.294999e-03, 9.226143e-03, 9.509662e-03, 9.528456e-03, 1.008240e-02, 1.130479e-02, 9.946953e-03, 1.045446e-02, 1.065627e-02, 9.886933e-03, 1.086657e-02, 9.301993e-03, 1.003286e-02, 8.438135e-03, 7.918646e-03, 8.066112e-03, 8.047024e-03, 7.739109e-03, 5.489504e-03, 5.373482e-03, 2.980811e-03, 2.655509e-03, 3.054776e-03, 2.638821e-03, 1.391185e-03, 1.251094e-03, 1.232359e-03, -7.353837e-05, 5.333216e-04}, 
{9.095917e-03, 8.770377e-03, 8.580059e-03, 8.741478e-03, 8.807038e-03, 9.189442e-03, 1.003386e-02, 9.160004e-03, 9.620326e-03, 9.561927e-03, 9.111829e-03, 9.301993e-03, 9.130468e-03, 8.903445e-03, 7.655048e-03, 7.335468e-03, 7.079197e-03, 6.860074e-03, 6.420194e-03, 4.811273e-03, 4.310421e-03, 2.867709e-03, 2.296081e-03, 2.603119e-03, 2.072177e-03, 1.011471e-03, 7.463800e-04, 7.033335e-04, -2.213334e-04, -1.368971e-04}, 
{8.458753e-03, 8.572908e-03, 8.627177e-03, 9.047659e-03, 9.100351e-03, 9.721484e-03, 1.115961e-02, 9.444436e-03, 1.013960e-02, 1.020130e-02, 9.372399e-03, 1.003286e-02, 8.903445e-03, 1.074486e-02, 8.076929e-03, 7.816542e-03, 8.101880e-03, 8.295281e-03, 8.091315e-03, 5.733095e-03, 5.608317e-03, 3.237945e-03, 2.832381e-03, 3.697703e-03, 3.222366e-03, 1.689925e-03, 1.546671e-03, 1.668484e-03, 2.127231e-04, 8.748946e-04}, 
{8.023143e-03, 7.767952e-03, 7.645219e-03, 7.801266e-03, 7.809130e-03, 8.189445e-03, 8.995298e-03, 8.154163e-03, 8.479301e-03, 8.621228e-03, 8.128266e-03, 8.438135e-03, 7.655048e-03, 8.076929e-03, 7.852361e-03, 6.551367e-03, 6.510903e-03, 6.403162e-03, 6.079255e-03, 4.462071e-03, 4.253550e-03, 2.555585e-03, 2.212736e-03, 2.329251e-03, 1.935072e-03, 1.037732e-03, 8.411023e-04, 7.649158e-04, -1.472802e-04, 1.860825e-04}, 
{7.180851e-03, 7.190231e-03, 7.150819e-03, 7.427642e-03, 7.539927e-03, 7.922299e-03, 8.826220e-03, 7.762886e-03, 8.340837e-03, 8.121004e-03, 7.710559e-03, 7.918646e-03, 7.335468e-03, 7.816542e-03, 6.551367e-03, 7.496125e-03, 6.386561e-03, 6.384438e-03, 6.102274e-03, 4.597698e-03, 4.079300e-03, 2.910344e-03, 2.293072e-03, 2.970566e-03, 2.409359e-03, 1.204493e-03, 9.345490e-04, 1.022878e-03, 9.177251e-05, 1.533163e-04}, 
{6.316530e-03, 6.617065e-03, 6.768508e-03, 7.195474e-03, 7.230225e-03, 7.783264e-03, 9.102942e-03, 7.449055e-03, 8.080938e-03, 8.149276e-03, 7.382721e-03, 8.066112e-03, 7.079197e-03, 8.101880e-03, 6.510903e-03, 6.386561e-03, 7.796980e-03, 7.137735e-03, 7.072099e-03, 4.994571e-03, 5.008331e-03, 2.872796e-03, 2.572943e-03, 3.459612e-03, 3.076972e-03, 1.671563e-03, 1.569231e-03, 1.709748e-03, 3.755339e-04, 1.088461e-03}, 
{5.396447e-03, 6.023973e-03, 6.351317e-03, 6.949310e-03, 6.985408e-03, 7.671817e-03, 9.332440e-03, 7.171483e-03, 7.951718e-03, 8.063100e-03, 7.091429e-03, 8.047024e-03, 6.860074e-03, 8.295281e-03, 6.403162e-03, 6.384438e-03, 7.137735e-03, 8.925672e-03, 7.875995e-03, 5.411354e-03, 5.654345e-03, 3.059013e-03, 2.856520e-03, 4.135696e-03, 3.798533e-03, 2.105837e-03, 2.081166e-03, 2.337271e-03, 7.020049e-04, 1.741795e-03}, 
{4.518830e-03, 5.353260e-03, 5.804593e-03, 6.495255e-03, 6.540526e-03, 7.292865e-03, 9.116808e-03, 6.699920e-03, 7.531300e-03, 7.701885e-03, 6.619684e-03, 7.739109e-03, 6.420194e-03, 8.091315e-03, 6.079255e-03, 6.102274e-03, 7.072099e-03, 7.875995e-03, 9.364826e-03, 5.466919e-03, 5.900992e-03, 3.027996e-03, 2.951977e-03, 4.388607e-03, 4.116588e-03, 2.331762e-03, 2.379216e-03, 2.694861e-03, 9.100071e-04, 2.162774e-03}, 
{3.736281e-03, 4.210453e-03, 4.433567e-03, 4.849579e-03, 4.908083e-03, 5.346795e-03, 6.433954e-03, 4.991145e-03, 5.577977e-03, 5.514694e-03, 4.933251e-03, 5.489504e-03, 4.811273e-03, 5.733095e-03, 4.462071e-03, 4.597698e-03, 4.994571e-03, 5.411354e-03, 5.466919e-03, 5.174355e-03, 3.923259e-03, 2.442148e-03, 2.163641e-03, 3.142313e-03, 2.812217e-03, 1.585830e-03, 1.481513e-03, 1.669844e-03, 5.781431e-04, 1.085182e-03}, 
{2.949985e-03, 3.540511e-03, 3.890020e-03, 4.362820e-03, 4.315294e-03, 4.874975e-03, 6.183370e-03, 4.462430e-03, 4.926539e-03, 5.285542e-03, 4.420442e-03, 5.373482e-03, 4.310421e-03, 5.608317e-03, 4.253550e-03, 4.079300e-03, 5.008331e-03, 5.654345e-03, 5.900992e-03, 3.923259e-03, 5.875622e-03, 2.104383e-03, 2.267903e-03, 3.086628e-03, 3.015363e-03, 1.866477e-03, 1.941444e-03, 2.087329e-03, 7.334793e-04, 1.865778e-03}, 
{2.166038e-03, 2.478098e-03, 2.575071e-03, 2.816633e-03, 2.925158e-03, 3.120347e-03, 3.636953e-03, 2.923034e-03, 3.363898e-03, 3.041134e-03, 2.887384e-03, 2.980811e-03, 2.867709e-03, 3.237945e-03, 2.555585e-03, 2.910344e-03, 2.872796e-03, 3.059013e-03, 3.027996e-03, 2.442148e-03, 2.104383e-03, 3.138840e-03, 1.475414e-03, 2.252104e-03, 1.883129e-03, 1.082427e-03, 8.822823e-04, 1.024384e-03, 4.890690e-04, 3.443540e-04}, 
{1.543153e-03, 1.856306e-03, 2.015742e-03, 2.255699e-03, 2.288335e-03, 2.535671e-03, 3.101503e-03, 2.351793e-03, 2.632607e-03, 2.643639e-03, 2.340028e-03, 2.655509e-03, 2.296081e-03, 2.832381e-03, 2.212736e-03, 2.293072e-03, 2.572943e-03, 2.856520e-03, 2.951977e-03, 2.163641e-03, 2.267903e-03, 1.475414e-03, 2.658807e-03, 1.927938e-03, 1.772976e-03, 1.171583e-03, 1.104069e-03, 1.184266e-03, 5.576128e-04, 7.990726e-04}, 
{6.044698e-04, 1.491926e-03, 1.886699e-03, 2.466039e-03, 2.624285e-03, 3.081196e-03, 4.281853e-03, 2.576774e-03, 3.407568e-03, 3.029214e-03, 2.507412e-03, 3.054776e-03, 2.603119e-03, 3.697703e-03, 2.329251e-03, 2.970566e-03, 3.459612e-03, 4.135696e-03, 4.388607e-03, 3.142313e-03, 3.086628e-03, 2.252104e-03, 1.927938e-03, 4.942891e-03, 3.321079e-03, 1.922003e-03, 1.860168e-03, 2.274467e-03, 1.169684e-03, 1.450528e-03}, 
{8.981331e-05, 9.435929e-04, 1.357885e-03, 1.914370e-03, 2.024329e-03, 2.494403e-03, 3.697824e-03, 2.015731e-03, 2.742209e-03, 2.563534e-03, 1.968671e-03, 2.638821e-03, 2.072177e-03, 3.222366e-03, 1.935072e-03, 2.409359e-03, 3.076972e-03, 3.798533e-03, 4.116588e-03, 2.812217e-03, 3.015363e-03, 1.883129e-03, 1.772976e-03, 3.321079e-03, 4.371322e-03, 1.901235e-03, 1.933544e-03, 2.299154e-03, 1.169386e-03, 1.658926e-03}, 
{-1.538775e-04, 3.003647e-04, 5.467148e-04, 8.541227e-04, 8.977604e-04, 1.184355e-03, 1.848189e-03, 9.605867e-04, 1.291029e-03, 1.321605e-03, 9.650146e-04, 1.391185e-03, 1.011471e-03, 1.689925e-03, 1.037732e-03, 1.204493e-03, 1.671563e-03, 2.105837e-03, 2.331762e-03, 1.585830e-03, 1.866477e-03, 1.082427e-03, 1.171583e-03, 1.922003e-03, 1.901235e-03, 2.413497e-03, 1.366661e-03, 1.537661e-03, 8.601952e-04, 1.178926e-03}, 
{-5.413937e-04, -4.095794e-05, 2.515424e-04, 5.884735e-04, 6.103385e-04, 9.358139e-04, 1.692185e-03, 6.825717e-04, 1.003741e-03, 1.145040e-03, 6.912607e-04, 1.251094e-03, 7.463800e-04, 1.546671e-03, 8.411023e-04, 9.345490e-04, 1.569231e-03, 2.081166e-03, 2.379216e-03, 1.481513e-03, 1.941444e-03, 8.822823e-04, 1.104069e-03, 1.860168e-03, 1.933544e-03, 1.366661e-03, 2.600980e-03, 1.675015e-03, 9.226991e-04, 1.454762e-03}, 
{-9.400317e-04, -2.625562e-04, 1.032703e-04, 5.361497e-04, 5.953450e-04, 9.772102e-04, 1.898850e-03, 6.295132e-04, 1.093593e-03, 1.125382e-03, 6.195130e-04, 1.232359e-03, 7.033335e-04, 1.668484e-03, 7.649158e-04, 1.022878e-03, 1.709748e-03, 2.337271e-03, 2.694861e-03, 1.669844e-03, 2.087329e-03, 1.024384e-03, 1.184266e-03, 2.274467e-03, 2.299154e-03, 1.537661e-03, 1.675015e-03, 3.042401e-03, 1.133970e-03, 1.715271e-03}, 
{-1.136515e-03, -7.525252e-04, -5.573032e-04, -3.410147e-04, -2.887999e-04, -1.298339e-04, 2.493618e-04, -2.998752e-04, -6.815230e-05, -1.302806e-04, -2.973148e-04, -7.353837e-05, -2.213334e-04, 2.127231e-04, -1.472802e-04, 9.177251e-05, 3.755339e-04, 7.020049e-04, 9.100071e-04, 5.781431e-04, 7.334793e-04, 4.890690e-04, 5.576128e-04, 1.169684e-03, 1.169386e-03, 8.601952e-04, 9.226991e-04, 1.133970e-03, 1.691497e-03, 1.029720e-03}, 
{-1.485396e-03, -8.452042e-04, -4.616821e-04, -1.137496e-04, -1.792823e-04, 1.585546e-04, 9.608016e-04, -2.003970e-04, 1.613219e-05, 3.606202e-04, -2.272584e-04, 5.333216e-04, -1.368971e-04, 8.748946e-04, 1.860825e-04, 1.533163e-04, 1.088461e-03, 1.741795e-03, 2.162774e-03, 1.085182e-03, 1.865778e-03, 3.443540e-04, 7.990726e-04, 1.450528e-03, 1.658926e-03, 1.178926e-03, 1.454762e-03, 1.715271e-03, 1.029720e-03, 3.021232e-03}
};

const double dystrR_he_sp [5][30][30] = {
{
{0.643117, 0.714286, 0.625, 0.58209, 0.459642, 0.266153, 0.131245, 0.0596301, 0.0468734, 0.016153, 0.0108975, 0.00430735, 0.00479463, 0.00416303, 0.0120413, 0, 0.00305372, 0.00377076, 0.00431742, 0, 0.00653595, 0.00746524, 0.00898433, 0, 0, 0, 0, 0, 0, 0.0109612}, 
{1, 1, 0.875, 0.865672, 0.766507, 0.548243, 0.315116, 0.161879, 0.110152, 0.0497016, 0.0367791, 0.012922, 0.0159821, 0.0124891, 0.0192661, 0, 0.0122149, 0.00377076, 0.00431742, 0, 0.00653595, 0.0149305, 0.0179687, 0, 0.0151515, 0, 0, 0, 0, 0.0109612}, 
{1, 1, 1, 0.940299, 0.921703, 0.750684, 0.506599, 0.292877, 0.187768, 0.0820076, 0.0711197, 0.0344588, 0.0319642, 0.0333042, 0.0240826, 0.00501568, 0.0122149, 0.0113123, 0.00863484, 0, 0.0261438, 0.0149305, 0.026953, 0, 0.0151515, 0, 0.0425532, 0, 0, 0.0222919}, 
{1, 1, 1, 0.955224, 0.967556, 0.89712, 0.695278, 0.500709, 0.292062, 0.15159, 0.10245, 0.0646102, 0.0415534, 0.0416303, 0.0337156, 0.0125392, 0.0152686, 0.0113123, 0.0129523, 0, 0.0261438, 0.0149305, 0.026953, 0, 0.0151515, 0, 0.0425532, 0, 0, 0.0332531}, 
{1, 1, 1, 1, 0.988719, 0.961835, 0.857714, 0.682236, 0.415104, 0.262176, 0.161227, 0.0918901, 0.0559373, 0.0603639, 0.0457569, 0.0225706, 0.0183223, 0.0113123, 0.0129523, 0.00457817, 0.0326797, 0.0149305, 0.026953, 0, 0.0151515, 0, 0.0425532, 0, 0, 0.0332531}, 
{1, 1, 1, 1, 0.995773, 0.991703, 0.940636, 0.833605, 0.58502, 0.386795, 0.243293, 0.153629, 0.0639283, 0.081179, 0.0602064, 0.0401254, 0.0396984, 0.0150831, 0.0215871, 0.00915634, 0.0326797, 0.0149305, 0.0449217, 0, 0.0454545, 0, 0.0425532, 0, 0, 0.0442143}, 
{1, 1, 1, 1, 0.995773, 0.998341, 0.974286, 0.932343, 0.771342, 0.55081, 0.34722, 0.206753, 0.110276, 0.0957496, 0.0866972, 0.0551725, 0.0427521, 0.0263953, 0.0259045, 0.0183127, 0.0392157, 0.0223957, 0.0449217, 0.0162187, 0.0454545, 0, 0.0425532, 0.03125, 0, 0.0442143}, 
{1, 1, 1, 1, 0.995773, 1, 0.994716, 0.978212, 0.911007, 0.739676, 0.492974, 0.268491, 0.163017, 0.131135, 0.115596, 0.0702195, 0.054967, 0.0301661, 0.0474916, 0.0183127, 0.0392157, 0.029861, 0.0449217, 0.0162187, 0.0454545, 0.0242728, 0.0425532, 0.03125, 0, 0.0661367}, 
{1, 1, 1, 1, 0.995773, 1, 0.99712, 0.998853, 0.973114, 0.874032, 0.686405, 0.399443, 0.230142, 0.176929, 0.155223, 0.0833032, 0.0732893, 0.0414784, 0.051809, 0.0228909, 0.0588235, 0.0373262, 0.053906, 0.0162187, 0.0454545, 0.0242728, 0.0638298, 0.03125, 0, 0.0661367}, 
{1, 1, 1, 1, 0.995773, 1, 0.99712, 1, 0.991863, 0.964964, 0.841986, 0.561687, 0.319642, 0.226885, 0.210613, 0.120921, 0.0916116, 0.0452492, 0.0647613, 0.0457817, 0.0588235, 0.0522567, 0.0628903, 0.0298823, 0.0606061, 0.0485457, 0.0638298, 0.03125, 0, 0.0770979}, 
{1, 1, 1, 1, 0.995773, 1, 0.99712, 1, 0.995379, 0.993542, 0.945512, 0.728238, 0.447498, 0.281004, 0.251553, 0.161046, 0.10688, 0.0678737, 0.0863484, 0.0549381, 0.0784314, 0.0597219, 0.0718747, 0.0298823, 0.0757576, 0.0485457, 0.0638298, 0.0625, 0.0476191, 0.0880592}, 
{1, 1, 1, 1, 0.995773, 1, 0.99712, 1, 0.996551, 0.998513, 0.983654, 0.880545, 0.602525, 0.355939, 0.290085, 0.193648, 0.125203, 0.0980398, 0.11657, 0.0778289, 0.0849673, 0.111979, 0.0988277, 0.0572096, 0.0909091, 0.0728185, 0.0638298, 0.09375, 0.0952381, 0.109982}, 
{1, 1, 1, 1, 0.995773, 1, 0.99712, 1, 0.996551, 0.998513, 0.995913, 0.950899, 0.765781, 0.451688, 0.343724, 0.246313, 0.141049, 0.120664, 0.15111, 0.105298, 0.0980392, 0.134033, 0.116796, 0.0845369, 0.121212, 0.0970914, 0.106383, 0.125, 0.0952381, 0.120943}, 
{1, 1, 1, 1, 0.995773, 1, 0.99712, 1, 0.996551, 0.998513, 1, 0.982486, 0.873425, 0.607961, 0.406339, 0.281422, 0.153264, 0.143289, 0.172697, 0.123611, 0.117647, 0.156428, 0.134765, 0.111864, 0.121212, 0.0970914, 0.106383, 0.125, 0.0952381, 0.120943}, 
{1, 1, 1, 1, 0.995773, 1, 0.99712, 1, 0.996551, 0.998513, 1, 0.993972, 0.937353, 0.766458, 0.500817, 0.344118, 0.202123, 0.184767, 0.207236, 0.155658, 0.130719, 0.163893, 0.161718, 0.125528, 0.151515, 0.0970914, 0.106383, 0.15625, 0.0952381, 0.120943}, 
{1, 1, 1, 1, 0.995773, 1, 0.99712, 1, 0.996551, 0.998513, 1, 0.995408, 0.978906, 0.889267, 0.660436, 0.462318, 0.272359, 0.256412, 0.25041, 0.174945, 0.150327, 0.186289, 0.179687, 0.180182, 0.166667, 0.121364, 0.12766, 0.1875, 0.0952381, 0.131904}, 
{1, 1, 1, 1, 0.995773, 1, 0.99712, 1, 0.996551, 0.998513, 1, 0.998279, 0.991692, 0.957957, 0.804931, 0.585202, 0.373588, 0.29789, 0.328124, 0.193257, 0.189542, 0.20122, 0.188671, 0.193846, 0.181818, 0.145637, 0.170213, 0.21875, 0.0952381, 0.153826}, 
{1, 1, 1, 1, 0.995773, 1, 0.99712, 1, 0.996551, 0.998513, 1, 0.998279, 0.99329, 0.987098, 0.927752, 0.723134, 0.480468, 0.369535, 0.405838, 0.253758, 0.222222, 0.283337, 0.261496, 0.234837, 0.257576, 0.16991, 0.212766, 0.28125, 0.0952381, 0.175749}, 
{1, 1, 1, 1, 0.995773, 1, 0.99712, 1, 0.996551, 0.998513, 1, 0.998279, 0.996487, 0.997506, 0.980734, 0.889655, 0.633154, 0.441179, 0.449012, 0.308696, 0.24183, 0.343059, 0.315402, 0.275828, 0.30303, 0.16991, 0.276596, 0.34375, 0.190476, 0.219594}, 
{1, 1, 1, 1, 0.995773, 1, 0.99712, 1, 0.996551, 1, 1, 0.998279, 1, 0.997506, 0.990367, 0.952351, 0.767917, 0.516595, 0.509456, 0.363634, 0.287582, 0.37292, 0.333371, 0.289491, 0.363636, 0.218456, 0.297872, 0.34375, 0.380952, 0.2744}, 
{1, 1, 1, 1, 0.995773, 1, 0.99712, 1, 0.996551, 1, 1, 0.998279, 1, 0.997506, 0.997592, 0.979937, 0.874797, 0.637259, 0.556947, 0.441463, 0.372549, 0.417712, 0.387277, 0.357809, 0.363636, 0.267001, 0.297872, 0.40625, 0.428571, 0.296322}, 
{1, 1, 1, 1, 0.995773, 1, 0.99712, 1, 0.996551, 1, 1, 0.998279, 1, 0.997506, 1, 0.997492, 0.95114, 0.795631, 0.651931, 0.514714, 0.424837, 0.469968, 0.432199, 0.385137, 0.424242, 0.291274, 0.340426, 0.40625, 0.52381, 0.351128}, 
{1, 1, 1, 1, 0.995773, 1, 0.99712, 1, 0.996551, 1, 1, 0.998279, 1, 0.997506, 1, 1, 0.981678, 0.893671, 0.74368, 0.569652, 0.496732, 0.552086, 0.47712, 0.439791, 0.484848, 0.364093, 0.425532, 0.46875, 0.571429, 0.438818}, 
{1, 1, 1, 1, 1, 1, 0.99712, 1, 0.997955, 1, 1, 0.998279, 1, 0.997506, 1, 1, 0.987785, 0.943439, 0.873203, 0.661215, 0.588235, 0.619273, 0.531026, 0.521773, 0.530303, 0.436911, 0.446808, 0.53125, 0.619048, 0.517707}, 
{1, 1, 1, 1, 1, 1, 0.99856, 1, 0.997955, 1, 1, 0.998279, 1, 0.997506, 1, 1, 0.993893, 0.984917, 0.91206, 0.752779, 0.640523, 0.68646, 0.593917, 0.5491, 0.606061, 0.50973, 0.510638, 0.625, 0.619048, 0.539629}, 
{1, 1, 1, 1, 1, 1, 0.99856, 1, 0.997955, 1, 1, 0.998279, 1, 0.997506, 1, 1, 0.996946, 0.996229, 0.959551, 0.844342, 0.718954, 0.738717, 0.674776, 0.658409, 0.666667, 0.606821, 0.574468, 0.625, 0.619048, 0.616358}, 
{1, 1, 1, 1, 1, 1, 0.99856, 1, 0.997955, 1, 1, 0.998279, 1, 0.997506, 1, 1, 1, 1, 0.986192, 0.917593, 0.810458, 0.776043, 0.74665, 0.754055, 0.712121, 0.679639, 0.744681, 0.625, 0.761905, 0.72597}, 
{1, 1, 1, 1, 1, 1, 1, 1, 0.997955, 1, 1, 1, 1, 0.997506, 1, 1, 1, 1, 0.990509, 0.963375, 0.875817, 0.880556, 0.836494, 0.863364, 0.818182, 0.728185, 0.87234, 0.75, 0.904762, 0.835582}, 
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.997506, 1, 1, 1, 1, 0.995683, 0.972531, 0.960784, 0.940278, 0.908368, 0.918018, 0.878788, 0.873822, 0.93617, 0.78125, 0.952381, 0.890388}, 
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}
}, 
{
{0.5, 0.777778, 0.589902, 0.537887, 0.439205, 0.235685, 0.123688, 0.0565784, 0.0310262, 0.0217531, 0.0187123, 0.00535893, 0.00894729, 0.00522144, 0.00395394, 0.00234861, 0.0027109, 0.00342748, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
{0.75, 0.944444, 0.828671, 0.862991, 0.70511, 0.526406, 0.307735, 0.151037, 0.0802111, 0.0535869, 0.0381248, 0.0187563, 0.0162678, 0.0113131, 0.00691939, 0.0105687, 0.0054218, 0.00856871, 0, 0.0100811, 0.00289766, 0.00401036, 0, 0, 0, 0, 0.0161598, 0, 0, 0.0051477}, 
{1, 1, 0.941033, 0.959405, 0.869454, 0.751935, 0.511444, 0.292762, 0.151139, 0.101456, 0.0697019, 0.0341632, 0.0252151, 0.0200155, 0.0128503, 0.0152659, 0.0108436, 0.0102824, 0.00192781, 0.0176419, 0.00579533, 0.0120311, 0, 0, 0, 0, 0.0323195, 0, 0, 0.0051477}, 
{1, 1, 0.955079, 0.984777, 0.952781, 0.885352, 0.713161, 0.469453, 0.25205, 0.156635, 0.103073, 0.0609578, 0.0455498, 0.0356798, 0.0197697, 0.0211375, 0.0216872, 0.0188512, 0.00578344, 0.0229492, 0.00579533, 0.0210725, 0.00454412, 0, 0, 0, 0.0323195, 0, 0.0294118, 0.0051477}, 
{1, 1, 0.955079, 1, 0.983334, 0.963711, 0.867447, 0.637068, 0.391615, 0.246374, 0.158626, 0.0890922, 0.0675114, 0.0549791, 0.0296545, 0.0317062, 0.0311753, 0.0257061, 0.0115669, 0.0229492, 0.00579533, 0.0210725, 0.00454412, 0.00735292, 0, 0, 0.0484793, 0, 0.0294118, 0.0051477}, 
{1, 1, 0.969124, 1, 0.994445, 0.991626, 0.949568, 0.796718, 0.567282, 0.368934, 0.218271, 0.125382, 0.0870327, 0.0723839, 0.0444818, 0.0422749, 0.0365971, 0.0291336, 0.0154225, 0.0279897, 0.00579533, 0.0290932, 0.00908824, 0.00735292, 0, 0, 0.0484793, 0, 0.0294118, 0.0121337}, 
{1, 1, 0.983169, 1, 0.998611, 0.999302, 0.986221, 0.921603, 0.767103, 0.538849, 0.312481, 0.177802, 0.120382, 0.0950101, 0.059309, 0.0493207, 0.0447298, 0.0359886, 0.0250616, 0.0279897, 0.00869299, 0.0290932, 0.0227206, 0.0133397, 0, 0, 0.0484793, 0, 0.0294118, 0.0238051}, 
{1, 1, 0.983169, 1, 0.998611, 1, 0.998886, 0.978641, 0.908509, 0.725077, 0.461697, 0.275424, 0.169999, 0.119377, 0.0790787, 0.0598895, 0.0528625, 0.046271, 0.0366284, 0.0405911, 0.0144883, 0.0371139, 0.0318088, 0.0372867, 0.0102041, 0, 0.0484793, 0, 0.0294118, 0.0238051}, 
{1, 1, 0.983169, 1, 0.998611, 1, 0.999393, 0.991521, 0.979858, 0.881154, 0.64131, 0.391311, 0.226122, 0.161148, 0.109722, 0.0786783, 0.0679521, 0.0582672, 0.0520509, 0.0506722, 0.017386, 0.0451346, 0.0363529, 0.0672206, 0.0204082, 0.0224719, 0.0484793, 0, 0.0294118, 0.0238051}, 
{1, 1, 0.983169, 1, 1, 1, 0.999393, 0.99658, 0.995609, 0.962861, 0.814984, 0.552249, 0.318107, 0.220484, 0.144319, 0.10921, 0.0896393, 0.0736909, 0.0636178, 0.0632736, 0.0376696, 0.0491449, 0.0408971, 0.0791941, 0.0204082, 0.0337079, 0.062734, 0.0198445, 0.0588235, 0.041312}, 
{1, 1, 0.983169, 1, 1, 1, 0.999393, 0.99658, 0.998951, 0.989919, 0.935444, 0.735224, 0.433609, 0.292714, 0.183858, 0.136219, 0.116748, 0.0959695, 0.0751847, 0.0783953, 0.0585811, 0.0772174, 0.0545294, 0.0911676, 0.0204082, 0.0674157, 0.0950535, 0.0198445, 0.0588235, 0.041312}, 
{1, 1, 0.983169, 1, 1, 1, 0.999393, 0.99704, 0.998951, 0.999469, 0.979886, 0.867187, 0.595473, 0.373646, 0.235259, 0.170274, 0.139791, 0.114821, 0.100246, 0.103598, 0.0788648, 0.0852381, 0.0590735, 0.128674, 0.0408163, 0.0786517, 0.0950535, 0.0198445, 0.0588235, 0.058819}, 
{1, 1, 0.983169, 1, 1, 1, 0.999393, 0.997843, 0.998951, 1, 0.992166, 0.954317, 0.741324, 0.482426, 0.304453, 0.208181, 0.169611, 0.149356, 0.117597, 0.133841, 0.0962507, 0.117321, 0.0733616, 0.152621, 0.0510204, 0.0898876, 0.111213, 0.0595334, 0.0588235, 0.0821617}, 
{1, 1, 0.983169, 1, 1, 1, 0.999393, 0.997843, 0.999428, 1, 0.996259, 0.983121, 0.869461, 0.630622, 0.371864, 0.263373, 0.200786, 0.173348, 0.150369, 0.154004, 0.104944, 0.125342, 0.0960822, 0.158608, 0.0714286, 0.101124, 0.111213, 0.0595334, 0.0882353, 0.0996687}, 
{1, 1, 0.983169, 1, 1, 1, 0.999393, 0.997843, 0.999428, 1, 0.998014, 0.995179, 0.943715, 0.779434, 0.486528, 0.341106, 0.261103, 0.205883, 0.177359, 0.186767, 0.113637, 0.137373, 0.119599, 0.182555, 0.102041, 0.101124, 0.143533, 0.0793778, 0.0882353, 0.105504}, 
{1, 1, 0.983169, 1, 1, 1, 0.999393, 0.997843, 0.999428, 1, 0.998599, 0.997858, 0.981944, 0.898737, 0.662478, 0.433876, 0.335653, 0.249045, 0.214388, 0.21701, 0.142613, 0.165445, 0.160496, 0.206502, 0.122449, 0.123596, 0.159693, 0.0793778, 0.147059, 0.134683}, 
{1, 1, 0.983169, 1, 1, 1, 0.999393, 0.997843, 0.999428, 1, 0.998599, 0.998528, 0.991705, 0.965814, 0.814972, 0.543086, 0.429179, 0.31074, 0.258728, 0.257335, 0.183181, 0.201538, 0.201393, 0.224462, 0.142857, 0.146067, 0.208172, 0.119067, 0.176471, 0.158025}, 
{1, 1, 0.983169, 1, 1, 1, 0.999393, 0.997843, 0.999428, 1, 0.998599, 0.999198, 0.996585, 0.992168, 0.927659, 0.703966, 0.533548, 0.396427, 0.307359, 0.300623, 0.201052, 0.241642, 0.215025, 0.242423, 0.183673, 0.168539, 0.224332, 0.158756, 0.205882, 0.187204}, 
{1, 1, 0.983169, 1, 1, 1, 0.999393, 0.997843, 0.999428, 1, 0.998599, 0.999198, 0.997399, 0.99913, 0.964402, 0.854386, 0.678759, 0.490683, 0.370977, 0.353549, 0.238721, 0.281745, 0.274099, 0.290317, 0.244898, 0.191011, 0.272811, 0.198445, 0.264706, 0.216382}, 
{1, 1, 0.983169, 1, 1, 1, 0.999393, 0.997843, 0.999428, 1, 0.998599, 0.999198, 0.998212, 0.99913, 0.993332, 0.938936, 0.804815, 0.598848, 0.434594, 0.406475, 0.282186, 0.32987, 0.31954, 0.350185, 0.265306, 0.235955, 0.33745, 0.297667, 0.264706, 0.24556}, 
{1, 1, 0.983169, 1, 1, 1, 0.999393, 0.997843, 0.999428, 1, 0.998599, 0.999198, 0.999187, 0.99913, 0.998275, 0.984734, 0.903763, 0.739948, 0.515562, 0.474522, 0.337242, 0.369973, 0.355893, 0.398079, 0.326531, 0.314607, 0.418249, 0.341322, 0.323529, 0.28641}, 
{1, 1, 0.983169, 1, 1, 1, 0.999393, 0.998646, 0.999428, 1, 0.998599, 0.999198, 0.999187, 0.99913, 0.998275, 0.99178, 0.966114, 0.851341, 0.611953, 0.557691, 0.415479, 0.406066, 0.437687, 0.428013, 0.367347, 0.382022, 0.515207, 0.4207, 0.323529, 0.338931}, 
{1, 1, 0.983169, 1, 1, 1, 0.999393, 0.998646, 0.999428, 1, 0.998599, 0.999198, 0.999187, 0.99913, 0.998275, 0.996477, 0.994578, 0.931887, 0.741361, 0.615658, 0.497082, 0.482263, 0.504691, 0.48788, 0.428571, 0.483146, 0.612166, 0.440544, 0.352941, 0.350602}, 
{1, 1, 0.983169, 1, 1, 1, 0.999393, 0.998646, 0.999428, 1, 0.998599, 0.999198, 0.999187, 1, 0.998275, 0.997651, 0.998645, 0.973017, 0.864741, 0.72479, 0.564303, 0.550439, 0.568309, 0.535774, 0.530612, 0.516854, 0.676805, 0.619144, 0.352941, 0.432301}, 
{1, 1, 1, 1, 1, 1, 0.999393, 0.998646, 0.999428, 1, 0.998599, 0.999198, 0.999187, 1, 0.998275, 1, 1, 0.993582, 0.93607, 0.81552, 0.659926, 0.610996, 0.641015, 0.61479, 0.591837, 0.617977, 0.725284, 0.698522, 0.529412, 0.496494}, 
{1, 1, 1, 1, 1, 1, 1, 0.998646, 0.999428, 1, 0.998599, 0.999198, 0.999187, 1, 1, 1, 1, 0.995295, 0.974938, 0.896169, 0.758447, 0.687192, 0.704632, 0.693495, 0.683673, 0.662921, 0.773763, 0.738211, 0.588235, 0.597339}, 
{1, 1, 1, 1, 1, 1, 1, 0.998646, 0.999428, 1, 0.998599, 0.999198, 0.999187, 1, 1, 1, 1, 0.997009, 0.988433, 0.939014, 0.836684, 0.755368, 0.781882, 0.741389, 0.755102, 0.775281, 0.806083, 0.758056, 0.676471, 0.679039}, 
{1, 1, 1, 1, 1, 1, 1, 0.998646, 0.999428, 1, 0.998599, 0.999198, 1, 1, 1, 1, 1, 0.997009, 0.994217, 0.976818, 0.909125, 0.827555, 0.840956, 0.807243, 0.826531, 0.876404, 0.854562, 0.857278, 0.735294, 0.772409}, 
{1, 1, 1, 1, 1, 1, 1, 0.999449, 0.999428, 1, 0.999299, 0.999198, 1, 1, 1, 1, 1, 0.997009, 1, 0.984379, 0.965228, 0.911772, 0.92275, 0.885071, 0.897959, 0.932584, 0.886882, 0.980156, 0.823529, 0.889122}, 
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}
}, 
{
{0.176842, 0.0617348, 0.05584, 0.390659, 0.413974, 0.353296, 0.230769, 0.146864, 0.136251, 0.079267, 0.0769039, 0.0349822, 0.037566, 0.0559495, 0.0135807, 0.00831827, 0.0236326, 0, 0, 0, 0.0266872, 0, 0, 0, 0, 0, 0, 0, 0, 0.00485399}, 
{0.198947, 0.185204, 0.11168, 0.494835, 0.581946, 0.604021, 0.435897, 0.308128, 0.249793, 0.189255, 0.153808, 0.0849569, 0.0845234, 0.0808159, 0.0543226, 0.0417499, 0.0551426, 0.018657, 0, 0, 0.0266872, 0, 0.0185185, 0, 0, 0, 0.0277778, 0, 0, 0.00485399}, 
{0.243158, 0.370409, 0.16752, 0.572967, 0.747536, 0.729384, 0.621795, 0.454735, 0.417835, 0.303642, 0.247191, 0.154921, 0.126785, 0.118116, 0.101855, 0.0833412, 0.0787752, 0.0466425, 0, 0.0131579, 0.0266872, 0, 0.0185185, 0, 0.0245428, 0, 0.0555556, 0, 0, 0.00485399}, 
{0.265263, 0.370409, 0.16752, 0.599011, 0.851029, 0.820557, 0.788462, 0.661354, 0.522294, 0.426829, 0.329588, 0.220396, 0.173743, 0.161632, 0.142597, 0.124933, 0.0866527, 0.074628, 0.0117647, 0.0131579, 0.0800615, 0.016093, 0.0185185, 0, 0.0245428, 0, 0.0555556, 0, 0, 0.00970799}, 
{0.353684, 0.493878, 0.22336, 0.625055, 0.851029, 0.866144, 0.884615, 0.80246, 0.649461, 0.536817, 0.417479, 0.31035, 0.258072, 0.205148, 0.176548, 0.18316, 0.110285, 0.102614, 0.0235294, 0.0131579, 0.0934051, 0.032186, 0.037037, 0, 0.0490855, 0.1, 0.0555556, 0, 0, 0.019416}, 
{0.375789, 0.493878, 0.2792, 0.625055, 0.851029, 0.877541, 0.929487, 0.893171, 0.800166, 0.6908, 0.549314, 0.385312, 0.314421, 0.298397, 0.21729, 0.224752, 0.133918, 0.158585, 0.0352941, 0.0394737, 0.106749, 0.032186, 0.037037, 0.0239873, 0.0490855, 0.1, 0.0555556, 0.05, 0, 0.0388319}, 
{0.42, 0.493878, 0.2792, 0.677143, 0.851029, 0.877541, 0.942308, 0.918369, 0.863749, 0.779152, 0.681149, 0.505251, 0.380162, 0.36678, 0.271613, 0.266343, 0.165428, 0.205227, 0.0705882, 0.0921053, 0.160123, 0.0643719, 0.0555556, 0.0239873, 0.0736283, 0.1, 0.0555556, 0.05, 0, 0.0436859}, 
{0.46421, 0.555613, 0.341463, 0.703187, 0.851029, 0.877541, 0.948718, 0.938527, 0.931875, 0.858343, 0.829464, 0.650177, 0.483468, 0.422729, 0.319145, 0.332889, 0.204815, 0.233213, 0.105882, 0.118421, 0.213497, 0.0965579, 0.0555556, 0.0239873, 0.098171, 0.1, 0.111111, 0.05, 0, 0.0533939}, 
{0.486316, 0.617348, 0.397303, 0.72923, 0.851029, 0.877541, 0.961538, 0.953645, 0.963667, 0.933135, 0.917913, 0.750127, 0.633732, 0.510472, 0.421475, 0.366162, 0.244203, 0.289184, 0.2, 0.144737, 0.213497, 0.112651, 0.0555556, 0.0239873, 0.098171, 0.166667, 0.138889, 0.15, 0.0588235, 0.0631019}, 
{0.530526, 0.679083, 0.397303, 0.72923, 0.851029, 0.888937, 0.967949, 0.963724, 0.981833, 0.959532, 0.972844, 0.865068, 0.774604, 0.622371, 0.536911, 0.465981, 0.291468, 0.345155, 0.270588, 0.210526, 0.253528, 0.144837, 0.111111, 0.0719618, 0.098171, 0.233333, 0.138889, 0.15, 0.0588235, 0.0776639}, 
{0.552631, 0.679083, 0.397303, 0.72923, 0.871728, 0.91173, 0.967949, 0.968764, 0.981833, 0.97273, 0.978338, 0.935033, 0.854432, 0.721836, 0.672718, 0.540846, 0.338733, 0.447768, 0.294118, 0.276316, 0.293559, 0.193116, 0.12963, 0.0959491, 0.122714, 0.233333, 0.194444, 0.15, 0.0588235, 0.0922259}, 
{0.596842, 0.679083, 0.397303, 0.72923, 0.871728, 0.923127, 0.967949, 0.968764, 0.981833, 0.97273, 0.983831, 0.945028, 0.920172, 0.827519, 0.747411, 0.590755, 0.417508, 0.494411, 0.329412, 0.302632, 0.320246, 0.241395, 0.166667, 0.143924, 0.171799, 0.233333, 0.333333, 0.15, 0.117647, 0.111642}, 
{0.663158, 0.740818, 0.397303, 0.72923, 0.892426, 0.934524, 0.967949, 0.968764, 0.990917, 0.97713, 0.983831, 0.970015, 0.953043, 0.889685, 0.801734, 0.632347, 0.464774, 0.541053, 0.341176, 0.328947, 0.33359, 0.273581, 0.185185, 0.167911, 0.245428, 0.333333, 0.388889, 0.25, 0.176471, 0.151538}, 
{0.663158, 0.740818, 0.397303, 0.755274, 0.892426, 0.94592, 0.974359, 0.968764, 0.990917, 0.97713, 0.983831, 0.985008, 0.962434, 0.939418, 0.862847, 0.690575, 0.543549, 0.587696, 0.4, 0.407895, 0.386964, 0.354046, 0.222222, 0.167911, 0.294513, 0.366667, 0.416667, 0.25, 0.176471, 0.170954}, 
{0.685263, 0.864287, 0.453143, 0.807362, 0.892426, 0.94592, 0.974359, 0.968764, 0.990917, 0.97713, 0.983831, 0.985008, 0.971826, 0.970501, 0.917169, 0.790394, 0.622324, 0.643667, 0.447059, 0.460526, 0.453682, 0.370139, 0.240741, 0.191898, 0.319056, 0.366667, 0.444444, 0.25, 0.176471, 0.195224}, 
{0.685263, 0.864287, 0.453143, 0.807362, 0.892426, 0.957317, 0.980769, 0.973803, 0.995458, 0.97713, 0.983831, 0.985008, 0.971826, 0.976717, 0.944331, 0.85694, 0.708977, 0.699638, 0.541176, 0.552632, 0.496343, 0.386232, 0.259259, 0.215885, 0.343599, 0.4, 0.527778, 0.35, 0.176471, 0.248618}, 
{0.685263, 0.864287, 0.453143, 0.807362, 0.913125, 0.957317, 0.987179, 0.973803, 0.995458, 0.981529, 0.983831, 0.985008, 0.981217, 0.982934, 0.957911, 0.906849, 0.779874, 0.774266, 0.588235, 0.592105, 0.563061, 0.466697, 0.388889, 0.287847, 0.392684, 0.466667, 0.555556, 0.35, 0.235294, 0.321428}, 
{0.685263, 0.864287, 0.453143, 0.833406, 0.933824, 0.968714, 0.987179, 0.973803, 0.995458, 0.985929, 0.983831, 0.985008, 0.990609, 0.982934, 0.964702, 0.965077, 0.842894, 0.820908, 0.694118, 0.671053, 0.576405, 0.48279, 0.481481, 0.335822, 0.44177, 0.533333, 0.583333, 0.35, 0.235294, 0.355406}, 
{0.733857, 0.926022, 0.508983, 0.85945, 0.933824, 0.968714, 1, 0.973803, 0.995458, 0.985929, 0.983831, 0.985008, 0.995304, 0.982934, 0.971492, 0.981714, 0.913792, 0.914193, 0.752941, 0.723684, 0.643123, 0.547161, 0.5, 0.431771, 0.44177, 0.533333, 0.611111, 0.45, 0.294118, 0.403946}, 
{0.733857, 0.926022, 0.508983, 0.85945, 0.933824, 0.968714, 1, 0.973803, 0.995458, 0.985929, 0.983831, 0.985008, 0.995304, 0.989151, 0.978282, 0.981714, 0.945302, 0.93285, 0.847059, 0.789474, 0.683153, 0.614677, 0.574074, 0.455758, 0.490855, 0.566667, 0.611111, 0.45, 0.294118, 0.463186}, 
{0.733857, 0.926022, 0.676503, 0.885494, 0.954522, 0.968714, 1, 0.983882, 1, 0.990328, 0.983831, 0.995003, 0.995304, 0.989151, 0.978282, 0.981714, 0.945302, 0.960836, 0.917647, 0.868421, 0.763215, 0.646863, 0.592593, 0.455758, 0.515398, 0.633333, 0.638889, 0.55, 0.294118, 0.49231}, 
{0.755963, 0.926022, 0.676503, 0.885494, 0.975221, 0.968714, 1, 0.988922, 1, 0.990328, 0.983831, 0.995003, 0.995304, 0.989151, 0.991863, 0.981714, 0.982683, 0.970164, 0.941176, 0.894737, 0.776559, 0.679049, 0.703704, 0.503733, 0.539941, 0.633333, 0.666667, 0.55, 0.352941, 0.56512}, 
{0.755963, 0.926022, 0.676503, 0.911538, 0.975221, 0.968714, 1, 0.988922, 1, 0.990328, 0.983831, 0.995003, 0.995304, 0.989151, 0.991863, 0.981714, 0.99056, 0.970164, 0.976471, 0.973684, 0.789902, 0.759514, 0.759259, 0.617558, 0.564483, 0.633333, 0.666667, 0.7, 0.470588, 0.608806}, 
{0.804557, 0.926022, 0.676503, 0.937582, 0.975221, 0.968714, 1, 0.988922, 1, 0.990328, 0.983831, 0.995003, 0.995304, 0.989151, 0.991863, 0.981714, 1, 0.979493, 0.976471, 0.973684, 0.829933, 0.823886, 0.777778, 0.68952, 0.638112, 0.666667, 0.722222, 0.7, 0.470588, 0.663435}, 
{0.843136, 0.926022, 0.676503, 0.937582, 0.975221, 0.968714, 1, 0.988922, 1, 0.990328, 0.983831, 0.995003, 0.995304, 0.989151, 0.991863, 0.990032, 1, 0.979493, 0.988235, 0.986842, 0.923338, 0.932437, 0.814815, 0.761481, 0.785368, 0.733333, 0.75, 0.75, 0.588235, 0.745952}, 
{0.869625, 0.926022, 0.676503, 0.968791, 0.975221, 0.968714, 1, 0.993961, 1, 0.994728, 0.983831, 0.995003, 0.995304, 0.989151, 0.991863, 0.990032, 1, 0.979493, 1, 0.986842, 0.976712, 0.980715, 0.87037, 0.809456, 0.828201, 0.733333, 0.833333, 0.85, 0.705882, 0.814072}, 
{0.896114, 0.926022, 0.810332, 0.968791, 0.975221, 0.968714, 1, 0.993961, 1, 0.994728, 0.983831, 1, 0.995304, 0.989151, 0.991863, 0.990032, 1, 0.988821, 1, 0.986842, 1, 0.980715, 0.888889, 0.851319, 0.877286, 0.833333, 0.888889, 0.85, 0.705882, 0.843196}, 
{0.949066, 0.926022, 0.866172, 0.968791, 0.975221, 0.968714, 1, 1, 1, 0.994728, 0.983831, 1, 0.995304, 0.989151, 0.991863, 1, 1, 0.988821, 1, 1, 1, 0.980715, 0.962963, 0.899294, 0.901829, 0.866667, 0.916667, 0.85, 0.705882, 0.893212}, 
{1, 1, 0.866172, 0.968791, 1, 0.968714, 1, 1, 1, 0.994728, 1, 1, 0.995304, 1, 0.991863, 1, 1, 1, 1, 1, 1, 0.980715, 0.962963, 0.947268, 0.975457, 0.866667, 0.944444, 0.9, 0.823529, 0.941752}, 
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}
}, 
{
{0.128526, 0.251841, 0.314558, 0.383763, 0.335815, 0.310907, 0.252346, 0.185428, 0.134248, 0.107669, 0.0753831, 0.0757714, 0.0663059, 0.0414456, 0.0435524, 0.0397066, 0.0373619, 0.0443565, 0.0217596, 0.0322711, 0.0323194, 0.026631, 0.0267284, 0, 0.0192412, 0, 0.0132982, 0.0160776, 0.0224078, 0.00581763}, 
{0.208855, 0.4757, 0.609313, 0.703566, 0.621258, 0.583568, 0.458208, 0.389399, 0.286921, 0.23777, 0.176679, 0.147308, 0.125634, 0.103461, 0.125002, 0.0893398, 0.114105, 0.0725834, 0.0598389, 0.0645421, 0.0840305, 0.0532621, 0.0890947, 0.0101631, 0.048103, 0.0106197, 0.0398945, 0.0321552, 0.0224078, 0.00872645}, 
{0.273117, 0.559647, 0.699187, 0.831488, 0.738793, 0.782033, 0.645002, 0.56688, 0.426434, 0.334523, 0.29211, 0.231994, 0.19704, 0.159571, 0.190331, 0.138973, 0.162574, 0.104843, 0.0707187, 0.112949, 0.109886, 0.0798931, 0.124733, 0.0304892, 0.0769649, 0.0106197, 0.0531927, 0.0321552, 0.0224078, 0.031997}, 
{0.305249, 0.587629, 0.766592, 0.863468, 0.789165, 0.835279, 0.777816, 0.704843, 0.552784, 0.457894, 0.381627, 0.327822, 0.286298, 0.257024, 0.264992, 0.185297, 0.194887, 0.157264, 0.135997, 0.123706, 0.142205, 0.0865508, 0.142552, 0.0406523, 0.115447, 0.0212393, 0.0664909, 0.0643104, 0.0224078, 0.0487088}, 
{0.321315, 0.587629, 0.766592, 0.879458, 0.797561, 0.888526, 0.870785, 0.779014, 0.681767, 0.554347, 0.485279, 0.405822, 0.355154, 0.31018, 0.311655, 0.228313, 0.239317, 0.177426, 0.184957, 0.150598, 0.174525, 0.126497, 0.196008, 0.0508154, 0.132237, 0.031859, 0.0930872, 0.080388, 0.0672234, 0.0632529}, 
{0.369512, 0.615611, 0.766592, 0.879458, 0.805956, 0.903048, 0.910629, 0.863781, 0.77653, 0.673516, 0.581864, 0.481594, 0.434211, 0.363336, 0.352097, 0.268019, 0.287787, 0.205653, 0.239356, 0.209762, 0.219772, 0.153128, 0.213827, 0.0914677, 0.151479, 0.031859, 0.0930872, 0.0964656, 0.0672234, 0.0807058}, 
{0.385578, 0.615611, 0.789061, 0.879458, 0.814351, 0.912729, 0.937192, 0.922059, 0.859371, 0.76997, 0.694938, 0.593022, 0.503068, 0.428305, 0.409134, 0.334197, 0.340295, 0.25001, 0.266555, 0.225897, 0.226236, 0.173102, 0.267284, 0.101631, 0.199582, 0.0633643, 0.119684, 0.112543, 0.0896312, 0.092341}, 
{0.401643, 0.615611, 0.789061, 0.879458, 0.814351, 0.912729, 0.940512, 0.94325, 0.90412, 0.859694, 0.784456, 0.706679, 0.587225, 0.543477, 0.496239, 0.380521, 0.417038, 0.314528, 0.299194, 0.279683, 0.258555, 0.193075, 0.320741, 0.121957, 0.218823, 0.116463, 0.132982, 0.176854, 0.112039, 0.107459}, 
{0.433775, 0.615611, 0.811529, 0.879458, 0.814351, 0.91757, 0.947153, 0.948548, 0.940972, 0.909043, 0.866906, 0.807604, 0.673933, 0.605493, 0.60512, 0.446699, 0.489742, 0.362917, 0.369913, 0.311954, 0.323194, 0.233022, 0.374198, 0.193098, 0.257305, 0.127082, 0.172876, 0.209009, 0.134447, 0.133638}, 
{0.433775, 0.615611, 0.811529, 0.879458, 0.814351, 0.91757, 0.947153, 0.956495, 0.946237, 0.934143, 0.930511, 0.885604, 0.745673, 0.691726, 0.673559, 0.519883, 0.574563, 0.415339, 0.407992, 0.39801, 0.368441, 0.292941, 0.419608, 0.213424, 0.295788, 0.18018, 0.212771, 0.241164, 0.179262, 0.145273}, 
{0.433775, 0.643594, 0.811529, 0.879458, 0.814351, 0.932091, 0.947153, 0.959144, 0.954134, 0.936386, 0.959284, 0.927946, 0.837481, 0.750789, 0.729555, 0.586061, 0.639189, 0.508084, 0.451512, 0.451795, 0.446008, 0.352861, 0.446336, 0.294729, 0.33427, 0.254518, 0.265963, 0.257242, 0.246486, 0.180179}, 
{0.449841, 0.671576, 0.811529, 0.879458, 0.822747, 0.932091, 0.947153, 0.959144, 0.954134, 0.938629, 0.963996, 0.950232, 0.898686, 0.83691, 0.797994, 0.672547, 0.695736, 0.576635, 0.52767, 0.521716, 0.472951, 0.412781, 0.481974, 0.355708, 0.401614, 0.296997, 0.319156, 0.321552, 0.246486, 0.206358}, 
{0.465906, 0.727541, 0.811529, 0.879458, 0.822747, 0.932091, 0.947153, 0.959144, 0.954134, 0.940873, 0.968707, 0.963603, 0.929289, 0.872348, 0.847769, 0.708945, 0.748245, 0.616959, 0.587509, 0.570122, 0.511734, 0.471005, 0.506433, 0.39636, 0.468959, 0.381954, 0.385647, 0.385863, 0.291301, 0.238355}, 
{0.465906, 0.727541, 0.811529, 0.879458, 0.847933, 0.932091, 0.947153, 0.961793, 0.956766, 0.940873, 0.968707, 0.963603, 0.949691, 0.93141, 0.894432, 0.748651, 0.788636, 0.661316, 0.614708, 0.623907, 0.563445, 0.517609, 0.515342, 0.457338, 0.4882, 0.413813, 0.43884, 0.40194, 0.336117, 0.279079}, 
{0.481972, 0.727541, 0.811529, 0.879458, 0.856328, 0.932091, 0.947153, 0.964442, 0.959398, 0.945359, 0.968707, 0.963603, 0.957342, 0.946176, 0.919319, 0.78174, 0.81691, 0.721802, 0.652788, 0.661557, 0.642218, 0.557555, 0.568799, 0.52848, 0.517062, 0.477531, 0.492032, 0.450173, 0.358525, 0.305258}, 
{0.498038, 0.727541, 0.811529, 0.879458, 0.856328, 0.932091, 0.947153, 0.967091, 0.964663, 0.947602, 0.968707, 0.965832, 0.959892, 0.952082, 0.947317, 0.8413, 0.845183, 0.758094, 0.701747, 0.705804, 0.66161, 0.590844, 0.588397, 0.599621, 0.555544, 0.48815, 0.518629, 0.530561, 0.420039, 0.32562}, 
{0.530169, 0.727541, 0.838454, 0.879458, 0.856328, 0.932091, 0.947153, 0.967091, 0.964663, 0.947602, 0.973419, 0.965832, 0.962442, 0.957989, 0.956649, 0.877698, 0.865379, 0.80245, 0.728946, 0.732697, 0.700393, 0.644106, 0.659673, 0.6606, 0.594027, 0.562488, 0.585119, 0.643104, 0.487262, 0.389614}, 
{0.530169, 0.755523, 0.838454, 0.898619, 0.866389, 0.932091, 0.952948, 0.967091, 0.967817, 0.947602, 0.973419, 0.965832, 0.962442, 0.960942, 0.965982, 0.906635, 0.889613, 0.854872, 0.788785, 0.786482, 0.752104, 0.690711, 0.695311, 0.691089, 0.632509, 0.594347, 0.611716, 0.675259, 0.50967, 0.456516}, 
{0.558208, 0.755523, 0.838454, 0.933771, 0.866389, 0.932091, 0.952948, 0.970266, 0.972411, 0.949845, 0.973419, 0.968061, 0.968048, 0.960942, 0.969093, 0.936415, 0.90577, 0.907293, 0.832304, 0.834888, 0.810279, 0.737315, 0.730949, 0.711415, 0.699853, 0.647445, 0.638312, 0.723492, 0.576894, 0.500149}, 
{0.558208, 0.755523, 0.838454, 0.933771, 0.866389, 0.942733, 0.952948, 0.970266, 0.972411, 0.949845, 0.97753, 0.97296, 0.968048, 0.964477, 0.976549, 0.939724, 0.942122, 0.92019, 0.870384, 0.873598, 0.836135, 0.777261, 0.766587, 0.741904, 0.709474, 0.700544, 0.678207, 0.771725, 0.599301, 0.535054}, 
{0.593525, 0.755523, 0.838454, 0.933771, 0.88104, 0.947573, 0.952948, 0.974889, 0.972411, 0.952533, 0.979886, 0.975188, 0.968048, 0.964477, 0.976549, 0.946453, 0.954239, 0.953249, 0.903023, 0.900491, 0.895592, 0.843534, 0.784406, 0.813046, 0.738336, 0.743022, 0.691505, 0.787803, 0.621709, 0.590189}, 
{0.612777, 0.755523, 0.838454, 0.933771, 0.891101, 0.956021, 0.952948, 0.974889, 0.975565, 0.955221, 0.979886, 0.977859, 0.968048, 0.970969, 0.979659, 0.946453, 0.958278, 0.973411, 0.930222, 0.932762, 0.90852, 0.876823, 0.802225, 0.845538, 0.767197, 0.764262, 0.771294, 0.80388, 0.644117, 0.657665}, 
{0.640815, 0.817037, 0.838454, 0.933771, 0.901161, 0.956021, 0.956927, 0.974889, 0.975565, 0.961824, 0.979886, 0.980529, 0.975555, 0.973922, 0.979659, 0.946453, 0.967157, 0.981476, 0.947621, 0.932762, 0.914983, 0.896796, 0.837863, 0.876028, 0.80568, 0.81736, 0.811188, 0.819958, 0.666525, 0.692571}, 
{0.660067, 0.817037, 0.838454, 0.933771, 0.901161, 0.956021, 0.960906, 0.974889, 0.975565, 0.961824, 0.985064, 0.980529, 0.978611, 0.973922, 0.979659, 0.949762, 0.971624, 0.981476, 0.958501, 0.943519, 0.934375, 0.910112, 0.855681, 0.906517, 0.834542, 0.883184, 0.837785, 0.852113, 0.666525, 0.736203}, 
{0.695385, 0.817037, 0.946151, 0.933771, 0.901161, 0.956021, 0.968863, 0.974889, 0.978198, 0.964512, 0.985064, 0.982758, 0.978611, 0.982615, 0.979659, 0.953727, 0.976464, 0.981476, 0.981339, 0.972538, 0.960231, 0.937447, 0.891319, 0.947169, 0.855691, 0.893803, 0.904276, 0.919612, 0.71134, 0.799456}, 
{0.733889, 0.966468, 0.973076, 0.933771, 0.901161, 0.956021, 0.974658, 0.982686, 0.978198, 0.969885, 0.985064, 0.982758, 0.984723, 0.982615, 0.979659, 0.963467, 0.976464, 0.985508, 0.987858, 0.983864, 0.960231, 0.950763, 0.909138, 0.947169, 0.874932, 0.915043, 0.904276, 0.93569, 0.823379, 0.845997}, 
{0.733889, 1, 0.973076, 0.952932, 0.911222, 0.956021, 0.978637, 0.98586, 0.985946, 0.976488, 0.985064, 0.985428, 0.984723, 0.991307, 0.979659, 0.963467, 0.981304, 0.985508, 0.987858, 0.983864, 0.979622, 0.964078, 0.944776, 0.957332, 0.923035, 0.936282, 0.917574, 0.967845, 0.845787, 0.890206}, 
{0.847702, 1, 1, 0.972094, 0.940525, 0.972463, 0.984432, 0.993658, 0.993694, 0.988466, 0.985064, 0.988099, 0.990832, 0.994846, 0.987115, 0.980481, 0.986144, 0.985508, 0.993481, 1, 0.98579, 0.970736, 0.953686, 0.969511, 0.961518, 0.946902, 0.960105, 0.967845, 0.890603, 0.922203}, 
{0.886187, 1, 1, 0.972094, 0.98994, 0.977303, 0.990226, 0.996829, 1, 0.995069, 0.990243, 0.99077, 0.996944, 1, 0.987115, 0.988105, 0.990985, 0.985508, 0.993481, 1, 1, 0.992022, 0.973272, 0.989837, 0.971138, 0.968141, 0.973404, 1, 0.955184, 0.951291}, 
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}
}, 
{
{0.486389, 0.345543, 0.870717, 0.376381, 0.539017, 0.409833, 0.411368, 0.318519, 0.234332, 0.174507, 0.14426, 0.0995766, 0.0822829, 0.0644396, 0.0421148, 0.0397232, 0.03253, 0.0286028, 0.0150777, 0.0176842, 0.0277036, 0.0166319, 0.0220832, 0.0132468, 0.0278243, 0.0142897, 0.0112106, 0.00881029, 0.00889144, 0.00756906}, 
{0.486389, 1, 1, 0.645641, 0.674513, 0.676133, 0.673055, 0.602613, 0.489101, 0.351854, 0.279812, 0.221429, 0.169289, 0.142186, 0.103504, 0.0904124, 0.0846827, 0.0873415, 0.0332308, 0.0430459, 0.0525458, 0.030442, 0.034579, 0.0165876, 0.0312308, 0.0283823, 0.0224211, 0.0171536, 0.0223496, 0.0138678}, 
{1, 1, 1, 0.752762, 0.892493, 0.901768, 0.887443, 0.788041, 0.627579, 0.532358, 0.390142, 0.309997, 0.251705, 0.214955, 0.192328, 0.143879, 0.116942, 0.136589, 0.074563, 0.0859538, 0.071817, 0.0492975, 0.0378323, 0.0294196, 0.0414504, 0.0393461, 0.0298276, 0.0384789, 0.0223496, 0.0214368}, 
{1, 1, 1, 0.909069, 0.946988, 0.971602, 0.973135, 0.872618, 0.786688, 0.668521, 0.522386, 0.425791, 0.333062, 0.27678, 0.270052, 0.213359, 0.158742, 0.169808, 0.118268, 0.105967, 0.0885302, 0.0659294, 0.0441664, 0.0554377, 0.062271, 0.0540302, 0.0446406, 0.0514609, 0.0353237, 0.0270028}, 
{1, 1, 1, 0.909069, 1, 0.986188, 0.992332, 0.922224, 0.879266, 0.775651, 0.622061, 0.526458, 0.440796, 0.354401, 0.334253, 0.267206, 0.198393, 0.216748, 0.159178, 0.128724, 0.105395, 0.0876067, 0.0850796, 0.0554377, 0.0726813, 0.0722374, 0.0484446, 0.0556326, 0.0531066, 0.0326338}, 
{1, 1, 1, 0.909069, 1, 1, 0.992332, 0.93686, 0.960653, 0.856113, 0.727895, 0.626847, 0.537117, 0.452172, 0.416185, 0.315753, 0.264654, 0.254842, 0.19269, 0.174237, 0.144089, 0.114778, 0.110071, 0.0814559, 0.0794943, 0.0900503, 0.0522487, 0.0725527, 0.0576733, 0.0387927}, 
{1, 1, 1, 0.909069, 1, 1, 0.992332, 0.956921, 0.992864, 0.94289, 0.818861, 0.742918, 0.673853, 0.557285, 0.490713, 0.410323, 0.316265, 0.320376, 0.254266, 0.204118, 0.180225, 0.147593, 0.141569, 0.094465, 0.104293, 0.121956, 0.0818747, 0.0811296, 0.061998, 0.0483649}, 
{1, 1, 1, 0.909069, 1, 1, 0.992332, 0.982134, 0.996529, 0.980123, 0.913435, 0.836298, 0.779886, 0.662899, 0.570735, 0.498212, 0.390459, 0.377581, 0.300624, 0.277046, 0.218919, 0.182931, 0.185908, 0.110815, 0.121707, 0.155044, 0.0890796, 0.0897064, 0.0708895, 0.059947}, 
{1, 1, 1, 0.909069, 1, 1, 0.992332, 0.982134, 1, 0.991481, 0.96552, 0.935182, 0.895498, 0.761672, 0.665052, 0.595814, 0.474059, 0.471718, 0.351727, 0.352993, 0.288481, 0.253906, 0.226649, 0.146147, 0.153129, 0.20262, 0.110896, 0.102221, 0.101888, 0.0751047}, 
{1, 1, 1, 0.909069, 1, 1, 0.992332, 0.982134, 1, 1, 0.982032, 0.958012, 0.948382, 0.870887, 0.786288, 0.687987, 0.575133, 0.573038, 0.431316, 0.426335, 0.381976, 0.349231, 0.292899, 0.188692, 0.194961, 0.220827, 0.144326, 0.136062, 0.12448, 0.0868295}, 
{1, 1, 1, 0.909069, 1, 1, 0.992332, 0.982134, 1, 1, 0.990288, 0.977952, 0.975023, 0.931835, 0.859033, 0.79719, 0.646638, 0.627031, 0.543854, 0.512151, 0.45169, 0.438911, 0.36531, 0.256724, 0.236411, 0.246475, 0.174153, 0.187289, 0.164613, 0.0979143}, 
{1, 1, 1, 0.909069, 1, 1, 0.992332, 0.982134, 1, 1, 0.990288, 0.98799, 0.986964, 0.970896, 0.908167, 0.843595, 0.732521, 0.70719, 0.632107, 0.580697, 0.537359, 0.515679, 0.45946, 0.366593, 0.32746, 0.319706, 0.248218, 0.225768, 0.18672, 0.118357}, 
{1, 1, 1, 0.909069, 1, 1, 0.992332, 0.982134, 1, 1, 0.990288, 0.990604, 0.989326, 0.988806, 0.929482, 0.887984, 0.812095, 0.773112, 0.698288, 0.648828, 0.623331, 0.584132, 0.532216, 0.434802, 0.397497, 0.403696, 0.325684, 0.298788, 0.239101, 0.150441}, 
{1, 1, 1, 0.909069, 1, 1, 0.992332, 0.982134, 1, 1, 0.990288, 0.990604, 0.994181, 1, 0.959598, 0.913455, 0.864112, 0.829413, 0.734875, 0.696532, 0.673167, 0.646641, 0.586142, 0.502835, 0.473013, 0.454006, 0.39655, 0.363464, 0.283316, 0.176449}, 
{1, 1, 1, 0.909069, 1, 1, 0.992332, 0.982134, 1, 1, 0.990288, 0.990604, 0.994181, 1, 0.985635, 0.948131, 0.908602, 0.869428, 0.781233, 0.741631, 0.731284, 0.718064, 0.642632, 0.541685, 0.528089, 0.504119, 0.4556, 0.41052, 0.335938, 0.207169}, 
{1, 1, 1, 0.909069, 1, 1, 0.992332, 0.982134, 1, 1, 0.990288, 0.990604, 0.994181, 1, 0.994951, 0.968937, 0.948253, 0.902389, 0.830667, 0.789611, 0.772384, 0.742863, 0.696386, 0.612703, 0.576353, 0.536813, 0.49563, 0.457342, 0.393128, 0.238061}, 
{1, 1, 1, 0.909069, 1, 1, 0.992332, 0.982134, 1, 1, 0.990288, 0.990604, 0.994181, 1, 0.997248, 0.980411, 0.977823, 0.932525, 0.879397, 0.845131, 0.819511, 0.814286, 0.743978, 0.6737, 0.666829, 0.619028, 0.576699, 0.499759, 0.450559, 0.283315}, 
{1, 1, 1, 0.909069, 1, 1, 0.992332, 0.982134, 1, 1, 0.990288, 0.990604, 0.994181, 1, 0.997248, 0.989743, 0.987499, 0.963048, 0.915282, 0.880499, 0.852482, 0.846952, 0.775304, 0.722218, 0.726838, 0.680499, 0.646758, 0.576017, 0.499341, 0.330102}, 
{1, 1, 1, 0.909069, 1, 1, 0.992332, 0.982134, 1, 1, 0.990288, 0.990604, 0.994181, 1, 0.997248, 0.992139, 0.995027, 0.974847, 0.953538, 0.910795, 0.871905, 0.874273, 0.815872, 0.769531, 0.76867, 0.7314, 0.694195, 0.635588, 0.553173, 0.379326}, 
{1, 1, 1, 0.909069, 1, 1, 0.992332, 0.982134, 1, 1, 0.990288, 0.990604, 0.994181, 1, 0.997248, 0.994409, 1, 0.984339, 0.974345, 0.951374, 0.907738, 0.906939, 0.860383, 0.808736, 0.810501, 0.760571, 0.745435, 0.686815, 0.615656, 0.440878}, 
{1, 1, 1, 0.909069, 1, 1, 0.992332, 0.982134, 1, 1, 0.990288, 0.990604, 0.994181, 1, 0.997248, 0.994409, 1, 0.989084, 0.979371, 0.953841, 0.938606, 0.926243, 0.876132, 0.854268, 0.831513, 0.778581, 0.77486, 0.738043, 0.686545, 0.496366}, 
{1, 1, 1, 0.909069, 1, 1, 0.992332, 0.982134, 1, 1, 0.990288, 0.990604, 0.994181, 1, 0.997248, 0.994409, 1, 0.994313, 0.987191, 0.971387, 0.960587, 0.945397, 0.911229, 0.895927, 0.862934, 0.814996, 0.819097, 0.780226, 0.748785, 0.558948}, 
{1, 1, 1, 0.909069, 1, 1, 0.992332, 0.982134, 1, 1, 0.990288, 0.990604, 0.994181, 1, 0.997248, 0.994409, 1, 0.994313, 0.987191, 0.981532, 0.971426, 0.958908, 0.923897, 0.918428, 0.894355, 0.840841, 0.844919, 0.809473, 0.788918, 0.616419}, 
{1, 1, 1, 0.909069, 1, 1, 0.992332, 0.982134, 1, 1, 0.996763, 0.990604, 0.994181, 1, 0.997248, 0.994409, 1, 0.994313, 0.989845, 0.986742, 0.985429, 0.961581, 0.936565, 0.938118, 0.912151, 0.869421, 0.866937, 0.851657, 0.82698, 0.669365}, 
{1, 1, 1, 0.909069, 1, 1, 0.992332, 0.982134, 1, 1, 0.996763, 0.990604, 0.994181, 1, 0.997248, 0.994409, 1, 0.997234, 0.994476, 0.98921, 0.993861, 0.972569, 0.961902, 0.964491, 0.950385, 0.916799, 0.90092, 0.885731, 0.858463, 0.718316}, 
{1, 1, 1, 0.909069, 1, 1, 1, 0.988307, 1, 1, 0.996763, 0.996868, 0.994181, 1, 0.997248, 0.994409, 1, 0.997234, 0.996989, 0.992331, 0.993861, 0.986379, 0.977823, 0.984004, 0.967609, 0.945773, 0.922937, 0.894307, 0.884653, 0.772717}, 
{1, 1, 1, 0.909069, 1, 1, 1, 0.988307, 1, 1, 0.996763, 1, 0.994181, 1, 0.997248, 0.994409, 1, 0.997234, 1, 0.992331, 0.993861, 0.992403, 0.987066, 0.990332, 0.978504, 0.96026, 0.955763, 0.911228, 0.933435, 0.838038}, 
{1, 1, 1, 0.909069, 1, 1, 1, 0.994154, 1, 1, 0.996763, 1, 0.994181, 1, 0.997248, 0.994409, 1, 0.997234, 1, 0.995453, 0.99729, 1, 0.990146, 0.993495, 0.982101, 0.974746, 0.974178, 0.944834, 0.951218, 0.891313}, 
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.996763, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.996836, 0.985698, 0.982187, 0.988789, 0.974503, 0.969001, 0.933852}, 
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}
}
};
const double wagi_he_sp [5][30] = {
{1.53277e-05, 0.00116772, 0.00306204, 0.0141346, 0.0668621, 0.161602, 0.244934, 0.29114, 0.322217, 0.336695, 0.333732, 0.345346, 0.351968, 0.320019, 0.297136, 0.301623, 0.277044, 0.256726, 0.252079, 0.253534, 0.200828, 0.195963, 0.174863, 0.12946, 0.145562, 0.0977367, 0.118892, 0.0865945, 0.0634849, 0.0299603}, 
{1.23742e-05, 0.00135496, 0.00596889, 0.0186842, 0.0737575, 0.168028, 0.255962, 0.311085, 0.338727, 0.337627, 0.338846, 0.338461, 0.321049, 0.337977, 0.316884, 0.30687, 0.296033, 0.263323, 0.254957, 0.223758, 0.222257, 0.178464, 0.169613, 0.144059, 0.101851, 0.103778, 0.0741747, 0.0713196, 0.0575695, 0.0307658}, 
{0.000424941, 0.0032811, 0.00390886, 0.00903655, 0.0130492, 0.0256413, 0.0497544, 0.0660861, 0.0835523, 0.0943314, 0.0809751, 0.0939374, 0.113351, 0.0909093, 0.0893751, 0.0824097, 0.0916395, 0.0855717, 0.0740682, 0.0705289, 0.0729222, 0.0666715, 0.0618411, 0.0503419, 0.0533647, 0.0415776, 0.0537977, 0.0329344, 0.0288509, 0.0197954}, 
{0.000117112, 0.00167338, 0.00230415, 0.00358008, 0.00757969, 0.0143996, 0.0236757, 0.032522, 0.0369691, 0.0477057, 0.050218, 0.059336, 0.0584924, 0.0564083, 0.0582992, 0.0613737, 0.0558138, 0.0601222, 0.0512278, 0.0553726, 0.0512078, 0.0573987, 0.0474546, 0.0450097, 0.0506736, 0.0514097, 0.0463652, 0.0405205, 0.0321288, 0.0196686}, 
{0.000176631, 0.000474435, 0.00100881, 0.00203761, 0.00372719, 0.00683947, 0.0145476, 0.0192749, 0.0277843, 0.035056, 0.0384888, 0.0462484, 0.0510015, 0.0550801, 0.0584018, 0.063287, 0.0657064, 0.0740125, 0.0722455, 0.079515, 0.0781702, 0.0854761, 0.0812141, 0.0852834, 0.0854817, 0.0892814, 0.0957286, 0.0858351, 0.0941735, 0.0521107}
};

const double dystrR_he_all [5][30][30] = {
{
{0.710065, 0.697841, 0.650334, 0.595794, 0.398951, 0.217805, 0.129088, 0.0616899, 0.0482679, 0.0243102, 0.0138978, 0.0132294, 0.0135411, 0.00691787, 0.0141362, 0.0020831, 0.0024912, 0.00291014, 0.00662513, 0, 0.00458716, 0.00554738, 0.00663563, 0, 0, 0, 0, 0, 0, 0.0041986}, 
{0.825306, 0.946678, 0.925072, 0.852804, 0.697503, 0.481716, 0.29312, 0.160184, 0.125024, 0.0676456, 0.0440096, 0.0300667, 0.0311445, 0.0155652, 0.0201946, 0.0083324, 0.0174384, 0.00582028, 0.00662513, 0, 0.00458716, 0.0110948, 0.0132713, 0, 0.010989, 0, 0, 0, 0, 0.0041986}, 
{0.934547, 1, 0.987512, 0.953271, 0.876401, 0.702799, 0.469208, 0.293609, 0.212317, 0.116266, 0.0801552, 0.0577281, 0.0487479, 0.0363188, 0.026253, 0.0166648, 0.0174384, 0.0145507, 0.0132503, 0.0035994, 0.0229358, 0.0110948, 0.0199069, 0, 0.010989, 0, 0.0253165, 0, 0, 0.0127373}, 
{0.934547, 1, 1, 0.978972, 0.963527, 0.870494, 0.660998, 0.481794, 0.309799, 0.180741, 0.127639, 0.0926055, 0.0690595, 0.0570724, 0.0383698, 0.0312465, 0.0199296, 0.0174608, 0.0198754, 0.00719879, 0.0275229, 0.0110948, 0.0199069, 0, 0.010989, 0, 0.0253165, 0, 0, 0.0253332}, 
{0.934547, 1, 1, 0.997664, 0.991407, 0.944585, 0.829373, 0.657298, 0.434728, 0.284388, 0.193826, 0.125078, 0.0893711, 0.0830144, 0.0585644, 0.0458282, 0.0323856, 0.0174608, 0.0198754, 0.0143976, 0.0366972, 0.0221895, 0.0199069, 0.00841214, 0.010989, 0.0155767, 0.0253165, 0, 0, 0.0253332}, 
{0.934547, 1, 1, 1, 0.997216, 0.985797, 0.927158, 0.813149, 0.588996, 0.410478, 0.28213, 0.191224, 0.115099, 0.110686, 0.078759, 0.0749916, 0.0523152, 0.0320115, 0.0364382, 0.0215964, 0.0458716, 0.0332843, 0.0331782, 0.00841214, 0.043956, 0.0311535, 0.0379747, 0.0148813, 0, 0.0337304}, 
{0.934547, 1, 1, 1, 0.997216, 0.993688, 0.969217, 0.914392, 0.768817, 0.567048, 0.390178, 0.259776, 0.169263, 0.141816, 0.11309, 0.102072, 0.0697536, 0.0582028, 0.0463759, 0.0431927, 0.0688073, 0.0443791, 0.0398138, 0.0183973, 0.043956, 0.0311535, 0.0506329, 0.0297627, 0, 0.037929}, 
{0.934547, 1, 1, 1, 0.997216, 0.998072, 0.991929, 0.971752, 0.902438, 0.743561, 0.530314, 0.330734, 0.230198, 0.190241, 0.153479, 0.124986, 0.099648, 0.0814839, 0.0728765, 0.0539909, 0.100917, 0.072116, 0.0530851, 0.0268095, 0.043956, 0.0467302, 0.0632911, 0.0297627, 0, 0.0505248}, 
{0.934547, 1, 1, 1, 0.997216, 0.998949, 0.996976, 0.994104, 0.965849, 0.872764, 0.70426, 0.458465, 0.300612, 0.250773, 0.210942, 0.150435, 0.129542, 0.110585, 0.099377, 0.0683885, 0.133028, 0.0887581, 0.0663563, 0.0436338, 0.0659341, 0.0623069, 0.101266, 0.0297627, 0, 0.0673192}, 
{0.934547, 1, 1, 1, 0.997216, 0.998949, 0.996976, 0.995892, 0.988563, 0.959627, 0.84696, 0.602785, 0.387275, 0.311304, 0.279603, 0.202513, 0.169402, 0.125136, 0.135815, 0.111581, 0.151376, 0.110948, 0.0862632, 0.0688702, 0.0989011, 0.0778837, 0.126582, 0.0297627, 0, 0.0841136}, 
{0.934547, 1, 1, 1, 0.997216, 0.998949, 0.996976, 0.995892, 0.992349, 0.989222, 0.944245, 0.74951, 0.505224, 0.368376, 0.332109, 0.25459, 0.20677, 0.157148, 0.168941, 0.122379, 0.197248, 0.12759, 0.119441, 0.0856945, 0.10989, 0.0934604, 0.202532, 0.089288, 0.0255547, 0.113504}, 
{0.934547, 1, 1, 1, 0.997216, 0.998949, 0.996976, 0.995892, 0.993295, 0.996621, 0.980147, 0.884305, 0.641988, 0.435826, 0.384101, 0.298336, 0.239155, 0.186249, 0.215317, 0.165572, 0.215596, 0.199706, 0.145984, 0.119343, 0.131868, 0.124614, 0.227848, 0.119051, 0.102219, 0.142894}, 
{0.934547, 1, 1, 1, 0.997216, 0.998949, 0.996976, 0.995892, 0.993295, 0.998735, 0.990571, 0.952857, 0.784372, 0.51884, 0.435139, 0.350413, 0.259556, 0.233483, 0.265005, 0.190768, 0.243119, 0.244456, 0.165891, 0.152992, 0.164835, 0.186921, 0.253165, 0.163695, 0.127774, 0.168086}, 
{0.934547, 1, 1, 1, 0.997216, 0.998949, 0.996976, 0.995892, 0.993295, 0.998735, 0.997224, 0.981721, 0.880991, 0.653871, 0.491684, 0.38166, 0.279486, 0.268404, 0.298131, 0.212364, 0.266055, 0.27774, 0.199069, 0.203464, 0.197802, 0.218074, 0.278481, 0.193457, 0.127774, 0.193277}, 
{0.934547, 1, 1, 1, 0.997216, 0.998949, 0.996976, 0.996786, 0.993295, 0.998735, 0.997224, 0.992545, 0.939218, 0.792479, 0.572928, 0.439986, 0.324327, 0.323697, 0.331257, 0.244759, 0.288991, 0.29993, 0.238883, 0.262349, 0.252747, 0.233651, 0.303797, 0.252983, 0.204438, 0.235263}, 
{0.934547, 1, 1, 1, 0.997216, 0.998949, 0.996976, 0.996786, 0.993295, 0.998735, 0.997224, 0.996153, 0.975779, 0.897976, 0.708797, 0.544417, 0.389098, 0.405754, 0.37432, 0.292317, 0.330275, 0.322119, 0.272061, 0.30441, 0.285714, 0.264804, 0.329114, 0.297626, 0.229992, 0.281448}, 
{0.934547, 1, 1, 1, 0.997216, 0.998949, 0.996976, 0.996786, 0.994428, 0.998735, 0.997224, 0.998559, 0.986612, 0.957121, 0.831984, 0.650655, 0.474171, 0.440676, 0.454686, 0.324711, 0.380734, 0.344309, 0.291968, 0.329647, 0.296703, 0.295958, 0.379747, 0.327389, 0.229992, 0.331831}, 
{0.934547, 1, 1, 1, 0.997216, 0.998949, 0.996976, 0.996786, 0.994428, 0.998735, 0.997224, 0.998559, 0.987966, 0.981333, 0.936996, 0.765225, 0.566346, 0.507609, 0.530875, 0.383076, 0.408257, 0.423072, 0.378933, 0.396944, 0.373626, 0.327111, 0.417722, 0.401796, 0.281102, 0.348626}, 
{0.934547, 1, 1, 1, 0.997216, 0.998949, 0.996976, 0.996786, 0.994428, 0.998735, 0.997224, 0.998559, 0.990674, 0.989981, 0.981424, 0.903544, 0.690906, 0.565812, 0.570626, 0.426268, 0.431193, 0.495188, 0.418747, 0.42218, 0.428571, 0.342688, 0.481013, 0.476202, 0.357766, 0.407406}, 
{0.934547, 1, 1, 1, 0.997216, 0.998949, 0.996976, 0.996786, 0.994428, 1, 0.997224, 0.998559, 0.993651, 0.989981, 0.989502, 0.957925, 0.803336, 0.624015, 0.620314, 0.47666, 0.472477, 0.522925, 0.465197, 0.439004, 0.505494, 0.373842, 0.506329, 0.520846, 0.459985, 0.457789}, 
{0.934547, 1, 1, 1, 0.997216, 0.998949, 0.996976, 0.996786, 0.994428, 1, 0.997224, 0.998559, 0.993651, 0.989981, 0.995561, 0.980839, 0.890528, 0.72005, 0.660065, 0.548648, 0.541284, 0.561757, 0.511646, 0.497889, 0.527473, 0.420572, 0.518987, 0.56549, 0.48554, 0.48718}, 
{1, 1, 1, 1, 0.997216, 0.998949, 0.996976, 0.996786, 0.994428, 1, 0.997224, 0.998559, 0.993651, 0.989981, 0.99758, 0.995421, 0.952808, 0.842275, 0.732942, 0.606238, 0.582569, 0.606136, 0.55146, 0.53995, 0.582418, 0.467302, 0.556962, 0.610134, 0.562204, 0.54596}, 
{1, 1, 1, 1, 0.997216, 0.998949, 0.996976, 0.996786, 0.994428, 1, 0.997224, 0.998559, 0.993651, 0.989981, 0.99758, 0.997504, 0.97772, 0.917939, 0.803337, 0.657344, 0.637615, 0.667157, 0.59259, 0.590423, 0.626374, 0.576339, 0.632911, 0.639897, 0.689977, 0.613138}, 
{1, 1, 1, 1, 0.998608, 0.998949, 0.996976, 0.996786, 0.995562, 1, 0.998612, 0.998559, 0.993651, 0.99171, 0.99758, 0.997504, 0.985687, 0.956348, 0.902714, 0.729332, 0.701835, 0.717083, 0.639039, 0.640896, 0.659341, 0.623069, 0.670886, 0.699422, 0.741087, 0.681143}, 
{1, 1, 1, 1, 0.998608, 0.998949, 0.997984, 0.996786, 0.995562, 1, 1, 0.998559, 0.993651, 0.993783, 0.99758, 0.997504, 0.99067, 0.988359, 0.932527, 0.801319, 0.743119, 0.76701, 0.685488, 0.674544, 0.714286, 0.669799, 0.708861, 0.803591, 0.770008, 0.706335}, 
{1, 1, 1, 1, 0.998608, 0.998949, 0.997984, 0.996786, 0.995562, 1, 1, 0.998559, 0.997637, 0.993783, 0.99758, 0.997504, 0.993161, 0.99709, 0.968966, 0.873307, 0.798165, 0.805842, 0.759796, 0.764935, 0.758242, 0.732106, 0.746835, 0.803591, 0.770008, 0.756718}, 
{1, 1, 1, 1, 0.998608, 1, 0.997984, 0.998929, 0.997214, 1, 1, 0.998559, 0.997637, 0.993783, 0.99758, 0.997504, 1, 1, 0.989405, 0.930898, 0.862385, 0.833578, 0.812881, 0.8339, 0.791209, 0.794413, 0.848101, 0.803591, 0.872226, 0.84069}, 
{1, 1, 1, 1, 1, 1, 1, 0.998929, 0.997214, 1, 1, 1, 1, 0.993783, 0.99758, 1, 1, 1, 0.992718, 0.971205, 0.912844, 0.911242, 0.879238, 0.901197, 0.868132, 0.825567, 0.924051, 0.863117, 0.948891, 0.899471}, 
{1, 1, 1, 1, 1, 1, 1, 1, 0.998866, 1, 1, 1, 1, 0.995855, 0.99758, 1, 1, 1, 0.996687, 0.978404, 0.972477, 0.955621, 0.932323, 0.949527, 0.912088, 0.919027, 0.962025, 0.877998, 0.974445, 0.937258}, 
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}
}, 
{
{0.561833, 0.697024, 0.634817, 0.529525, 0.378461, 0.210499, 0.117868, 0.0601169, 0.0407954, 0.0238673, 0.0239483, 0.0123844, 0.0124929, 0.00662781, 0.00656495, 0.0067254, 0.00657844, 0.0041296, 0, 0.00378612, 0, 0, 0, 0, 0.0112115, 0, 0, 0.0103743, 0, 0}, 
{0.824022, 0.907281, 0.878081, 0.844657, 0.663092, 0.477332, 0.285619, 0.155127, 0.0933754, 0.0616946, 0.0469971, 0.030961, 0.0249858, 0.0154649, 0.0139505, 0.0163331, 0.0109641, 0.0123888, 0.00445441, 0.0132514, 0.00218867, 0.00563591, 0.00334048, 0, 0.0112115, 0, 0.00834153, 0.0103743, 0, 0.00227971}, 
{0.917661, 0.966584, 0.961196, 0.954415, 0.84324, 0.693795, 0.475772, 0.287016, 0.163039, 0.111782, 0.0869109, 0.053478, 0.0367847, 0.0265112, 0.0246185, 0.0269016, 0.0230245, 0.0178949, 0.00890882, 0.0246098, 0.00875467, 0.0112718, 0.00334048, 0, 0.0112115, 0.00694445, 0.0166831, 0.0103743, 0, 0.00227971}, 
{0.917661, 0.99354, 0.983495, 0.9838, 0.945927, 0.850551, 0.665921, 0.456909, 0.262947, 0.180302, 0.127857, 0.0845699, 0.0659348, 0.043449, 0.0352866, 0.0414555, 0.035085, 0.0302837, 0.0193024, 0.0361685, 0.00875467, 0.017625, 0.00668097, 0, 0.0112115, 0.0138889, 0.0166831, 0.0103743, 0.0298507, 0.00227971}, 
{0.936389, 0.99354, 0.987549, 0.995035, 0.982376, 0.941367, 0.827076, 0.622553, 0.396862, 0.272231, 0.193714, 0.118492, 0.0936968, 0.0693541, 0.0508783, 0.0616317, 0.0493383, 0.0399195, 0.0267265, 0.0361685, 0.00875467, 0.0204429, 0.00668097, 0.00498202, 0.022423, 0.0138889, 0.0333661, 0.0103743, 0.0298507, 0.00744848}, 
{0.936389, 0.99354, 0.989576, 0.998492, 0.995295, 0.979205, 0.927201, 0.775407, 0.562547, 0.39472, 0.255082, 0.165876, 0.123541, 0.0951289, 0.0722144, 0.0837295, 0.0613988, 0.0509318, 0.0326657, 0.0456338, 0.0153207, 0.0373507, 0.0200429, 0.00498202, 0.022423, 0.0138889, 0.0500492, 0.0207486, 0.0298507, 0.0105423}, 
{0.936389, 0.99354, 0.991604, 0.998492, 0.998986, 0.996105, 0.971818, 0.908233, 0.752549, 0.553799, 0.350929, 0.228504, 0.16796, 0.126145, 0.100936, 0.0971803, 0.0778449, 0.0660736, 0.0519681, 0.0494199, 0.0218867, 0.0429866, 0.0367453, 0.0130948, 0.0392403, 0.0138889, 0.0500492, 0.0207486, 0.0298507, 0.0182955}, 
{0.936389, 0.99354, 0.991604, 0.998492, 0.998986, 0.998353, 0.99175, 0.971636, 0.8924, 0.733479, 0.493209, 0.3354, 0.235977, 0.171803, 0.141146, 0.123121, 0.094291, 0.0867216, 0.0742402, 0.0721366, 0.0350187, 0.0664152, 0.0601287, 0.037433, 0.0560576, 0.0277778, 0.0667322, 0.031123, 0.0298507, 0.0260486}, 
{0.936389, 0.99354, 0.991604, 0.998492, 0.998986, 0.999102, 0.996377, 0.989381, 0.96947, 0.879462, 0.661423, 0.449674, 0.301912, 0.233663, 0.18628, 0.159631, 0.124039, 0.123888, 0.103936, 0.10072, 0.0481507, 0.0861409, 0.0734906, 0.0820531, 0.0896922, 0.0694444, 0.0834153, 0.0518716, 0.0447761, 0.0415549}, 
{0.955116, 0.99354, 0.991604, 0.998492, 0.999447, 0.999102, 0.996733, 0.995176, 0.991044, 0.955175, 0.822802, 0.598992, 0.394976, 0.297866, 0.230756, 0.197101, 0.176667, 0.154172, 0.135117, 0.142367, 0.0875467, 0.105867, 0.0935335, 0.110448, 0.100904, 0.0833333, 0.107456, 0.0726202, 0.0597015, 0.0552235}, 
{0.955116, 0.99354, 0.991604, 0.998492, 0.999447, 0.999102, 0.997089, 0.99561, 0.995752, 0.983996, 0.933563, 0.763452, 0.50186, 0.365617, 0.277532, 0.23767, 0.213945, 0.188585, 0.170752, 0.166977, 0.129606, 0.134046, 0.130279, 0.134786, 0.152425, 0.125, 0.165847, 0.124492, 0.0746269, 0.0784829}, 
{0.955116, 0.99354, 0.991604, 0.998492, 0.999447, 0.999102, 0.997089, 0.995972, 0.996144, 0.995254, 0.974587, 0.879415, 0.64414, 0.442205, 0.330872, 0.27514, 0.247933, 0.218869, 0.206388, 0.206731, 0.166813, 0.153772, 0.13696, 0.204819, 0.202877, 0.159722, 0.199213, 0.124492, 0.134328, 0.11208}, 
{0.955116, 0.99354, 0.991604, 0.998492, 0.999447, 0.999102, 0.997089, 0.996604, 0.996144, 0.997056, 0.988058, 0.957701, 0.771368, 0.543095, 0.399804, 0.318644, 0.285211, 0.255894, 0.237569, 0.246868, 0.193077, 0.201677, 0.160825, 0.233214, 0.236511, 0.194444, 0.232579, 0.197112, 0.134328, 0.158599}, 
{0.955116, 0.99354, 0.991604, 0.998492, 0.999447, 0.999102, 0.997089, 0.996604, 0.996536, 0.997056, 0.99205, 0.983596, 0.882092, 0.66924, 0.463973, 0.379173, 0.318103, 0.290307, 0.276173, 0.280461, 0.219341, 0.221403, 0.190889, 0.2768, 0.281358, 0.208333, 0.249262, 0.207486, 0.179104, 0.19478}, 
{0.955116, 0.99354, 0.991604, 0.998492, 0.999447, 0.999102, 0.997089, 0.996604, 0.996536, 0.997056, 0.994045, 0.994292, 0.946146, 0.802533, 0.563268, 0.446614, 0.375664, 0.331582, 0.315082, 0.320215, 0.243416, 0.255218, 0.218198, 0.317364, 0.314992, 0.243056, 0.29097, 0.269732, 0.19403, 0.225793}, 
{0.955116, 0.99354, 0.991604, 0.998492, 0.999447, 0.999102, 0.997089, 0.996604, 0.996536, 0.997056, 0.995043, 0.996543, 0.981543, 0.907173, 0.713441, 0.523476, 0.440352, 0.377264, 0.357213, 0.354811, 0.282812, 0.28058, 0.254944, 0.349815, 0.337415, 0.263889, 0.315995, 0.290481, 0.238806, 0.254221}, 
{0.955116, 0.99354, 0.991604, 0.998492, 0.999447, 0.999102, 0.997089, 0.996604, 0.996536, 0.997056, 0.995542, 0.997106, 0.991259, 0.966145, 0.841679, 0.61571, 0.523795, 0.435078, 0.403242, 0.411603, 0.328774, 0.320031, 0.295383, 0.37821, 0.365444, 0.305556, 0.366044, 0.321604, 0.253731, 0.30074}, 
{0.955116, 0.99354, 0.991604, 0.998492, 0.999447, 0.999102, 0.997089, 0.996604, 0.996536, 0.997056, 0.995542, 0.997669, 0.995424, 0.988446, 0.936871, 0.748296, 0.610411, 0.508035, 0.445152, 0.451691, 0.351027, 0.367937, 0.322107, 0.414717, 0.415896, 0.340278, 0.391068, 0.394224, 0.328358, 0.339506}, 
{0.955116, 0.99354, 0.991604, 0.998492, 0.999447, 0.999102, 0.997089, 0.996604, 0.996536, 0.997056, 0.995542, 0.997669, 0.996118, 0.995811, 0.968195, 0.873287, 0.730347, 0.58512, 0.498605, 0.493338, 0.386046, 0.421478, 0.392257, 0.459337, 0.471953, 0.375, 0.432776, 0.435721, 0.402985, 0.378272}, 
{0.955116, 0.99354, 0.991604, 0.998492, 0.999447, 0.999102, 0.997089, 0.996604, 0.997221, 0.997056, 0.995542, 0.997669, 0.996812, 0.995811, 0.993033, 0.944384, 0.836699, 0.673379, 0.552058, 0.534985, 0.425442, 0.469383, 0.432343, 0.516126, 0.499982, 0.430556, 0.50785, 0.52909, 0.447761, 0.409284}, 
{0.955116, 0.99354, 0.991604, 0.998492, 0.999447, 0.999102, 0.997089, 0.996604, 0.997221, 0.997841, 0.995542, 0.997669, 0.998475, 0.998381, 0.997136, 0.983531, 0.916736, 0.788091, 0.619167, 0.589884, 0.47797, 0.517288, 0.47623, 0.55669, 0.544828, 0.486111, 0.61629, 0.572661, 0.522388, 0.468725}, 
{0.955116, 0.99354, 0.991604, 0.998492, 0.999447, 0.999102, 0.99771, 0.997236, 0.997905, 0.997841, 0.99614, 0.997669, 0.998475, 0.998381, 0.997136, 0.989295, 0.967171, 0.878943, 0.693407, 0.652355, 0.537064, 0.562376, 0.558891, 0.589141, 0.578463, 0.5625, 0.683022, 0.634907, 0.552239, 0.516255}, 
{0.977558, 0.99354, 0.991604, 0.998492, 0.999447, 0.999102, 0.99771, 0.997236, 0.998375, 0.998381, 0.99614, 0.997669, 0.998475, 0.998381, 0.997136, 0.993138, 0.991508, 0.94364, 0.793077, 0.695895, 0.598701, 0.618735, 0.614829, 0.629705, 0.63452, 0.631944, 0.749754, 0.645281, 0.567164, 0.547355}, 
{0.977558, 0.99354, 0.991604, 0.998492, 0.999447, 0.999102, 0.99771, 0.997868, 0.998845, 0.998921, 0.99614, 0.997669, 0.998475, 0.999118, 0.997136, 0.994099, 0.994797, 0.976677, 0.88959, 0.777868, 0.658663, 0.672276, 0.67197, 0.674325, 0.713001, 0.652778, 0.791462, 0.749024, 0.567164, 0.61502}, 
{1, 0.99354, 0.994033, 0.998492, 0.999447, 0.999102, 0.99771, 0.9985, 0.998845, 0.998921, 0.997011, 0.998344, 0.998475, 1, 0.997136, 0.996021, 0.998087, 0.993195, 0.944527, 0.846018, 0.730889, 0.717645, 0.732099, 0.727862, 0.757847, 0.729167, 0.824828, 0.800896, 0.686567, 0.680369}, 
{1, 1, 0.996462, 0.998492, 0.999447, 0.999102, 0.998758, 0.9985, 0.99953, 0.99946, 0.997011, 0.998344, 0.998475, 1, 0.998568, 0.996021, 0.998087, 0.994572, 0.975949, 0.910758, 0.807927, 0.776822, 0.778866, 0.785246, 0.813905, 0.770833, 0.866536, 0.832019, 0.746269, 0.758627}, 
{1, 1, 0.996462, 0.998492, 1, 0.999102, 0.999379, 0.9985, 0.99953, 0.99946, 0.997011, 0.998344, 0.998475, 1, 1, 0.996021, 1, 0.995948, 0.988122, 0.948994, 0.871832, 0.828105, 0.835654, 0.817697, 0.859856, 0.840278, 0.89156, 0.842393, 0.820896, 0.820652}, 
{1, 1, 0.996462, 0.998492, 1, 0.999102, 0.999379, 0.998934, 0.99953, 0.99946, 0.997609, 0.998344, 0.999169, 1, 1, 0.997172, 1, 0.997598, 0.995546, 0.980694, 0.928738, 0.878828, 0.87908, 0.862317, 0.904702, 0.902778, 0.916585, 0.915013, 0.850746, 0.883703}, 
{1, 1, 0.996462, 0.998492, 1, 1, 1, 0.999566, 0.99953, 0.99946, 0.998804, 0.998344, 1, 1, 1, 0.998849, 1, 0.997598, 1, 0.988266, 0.971113, 0.938005, 0.939209, 0.91505, 0.943942, 0.944444, 0.941609, 0.989626, 0.895522, 0.945728}, 
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}
}, 
{
{0.201441, 0.283988, 0.333903, 0.445461, 0.436683, 0.367036, 0.250211, 0.19522, 0.156124, 0.155904, 0.118329, 0.0917868, 0.0638283, 0.0652522, 0.0369757, 0.0377011, 0.0201163, 0.0295885, 0.00477228, 0.00498508, 0.0103858, 0.0116927, 0.00647902, 0.00668425, 0.00788273, 0.00404123, 0, 0.00978661, 0.0059326, 0.00270377}, 
{0.263957, 0.383383, 0.492506, 0.635338, 0.680456, 0.621884, 0.487182, 0.402764, 0.325917, 0.294421, 0.243675, 0.207773, 0.164009, 0.122929, 0.0980661, 0.0915942, 0.0585201, 0.0528365, 0.0238614, 0.0204317, 0.0233681, 0.0146159, 0.0194371, 0.0100264, 0.00788273, 0.00808245, 0.00519269, 0.00978661, 0.0059326, 0.00592354}, 
{0.347312, 0.496978, 0.559287, 0.735273, 0.80058, 0.759891, 0.667789, 0.57682, 0.509573, 0.444903, 0.384027, 0.328825, 0.256851, 0.208871, 0.196132, 0.154429, 0.120698, 0.0803115, 0.0548812, 0.0353869, 0.0415433, 0.0263086, 0.0226766, 0.0233949, 0.0197068, 0.0242474, 0.0103854, 0.0195732, 0.0177978, 0.00862731}, 
{0.416774, 0.553776, 0.626067, 0.780244, 0.860643, 0.841901, 0.793745, 0.723646, 0.630071, 0.572937, 0.507625, 0.444854, 0.373081, 0.307874, 0.279883, 0.244194, 0.201163, 0.118354, 0.10499, 0.0777601, 0.0752972, 0.0555404, 0.0421136, 0.0334213, 0.0275896, 0.0484947, 0.0103854, 0.0342531, 0.0177978, 0.0124126}, 
{0.493182, 0.596374, 0.65111, 0.820219, 0.889093, 0.887185, 0.887362, 0.835183, 0.7492, 0.694304, 0.633999, 0.560585, 0.48327, 0.421394, 0.392417, 0.346724, 0.267391, 0.184115, 0.174188, 0.117641, 0.122033, 0.0906186, 0.0647902, 0.0467898, 0.043355, 0.0774549, 0.0259634, 0.0342531, 0.0237304, 0.0167386}, 
{0.53486, 0.624773, 0.6845, 0.840206, 0.898577, 0.915218, 0.92827, 0.901541, 0.842563, 0.789287, 0.741162, 0.652619, 0.57321, 0.54575, 0.47923, 0.445465, 0.373459, 0.27288, 0.288723, 0.187432, 0.179155, 0.143236, 0.103664, 0.0735268, 0.0709446, 0.0895786, 0.0467342, 0.0587196, 0.0415282, 0.0253907}, 
{0.583483, 0.624773, 0.701195, 0.855196, 0.90806, 0.930312, 0.946993, 0.934014, 0.883894, 0.854037, 0.813571, 0.753718, 0.677656, 0.642832, 0.598195, 0.562301, 0.468554, 0.393873, 0.388941, 0.279656, 0.259645, 0.245547, 0.152257, 0.127001, 0.134007, 0.105744, 0.0778903, 0.0880795, 0.0533934, 0.0318797}, 
{0.625161, 0.653172, 0.718851, 0.870186, 0.923866, 0.936781, 0.958908, 0.952368, 0.924973, 0.917424, 0.889952, 0.834863, 0.756086, 0.717633, 0.6834, 0.686177, 0.58772, 0.537589, 0.522663, 0.399298, 0.394661, 0.359551, 0.252682, 0.20387, 0.197068, 0.178486, 0.160973, 0.127226, 0.0889889, 0.0410725}, 
{0.645999, 0.68157, 0.752241, 0.890173, 0.933349, 0.94325, 0.967418, 0.962251, 0.946882, 0.953043, 0.934533, 0.894724, 0.825717, 0.776701, 0.755857, 0.752602, 0.684834, 0.641148, 0.647131, 0.523925, 0.542659, 0.496941, 0.388741, 0.324186, 0.283778, 0.295681, 0.254442, 0.220199, 0.155127, 0.0610804}, 
{0.680731, 0.69577, 0.768936, 0.89517, 0.939672, 0.951876, 0.974227, 0.966487, 0.966052, 0.962277, 0.960918, 0.938622, 0.892447, 0.832845, 0.813732, 0.815648, 0.736039, 0.7278, 0.754836, 0.641074, 0.677674, 0.587559, 0.566914, 0.47124, 0.41374, 0.400753, 0.430993, 0.347425, 0.250049, 0.0902811}, 
{0.701569, 0.69577, 0.777284, 0.89517, 0.945994, 0.958345, 0.977631, 0.969311, 0.972899, 0.968873, 0.966473, 0.961236, 0.927262, 0.882182, 0.878092, 0.85694, 0.778101, 0.803885, 0.793014, 0.733298, 0.742586, 0.684024, 0.654381, 0.592442, 0.539864, 0.498087, 0.550425, 0.471094, 0.416162, 0.146626}, 
{0.722408, 0.709969, 0.793979, 0.900167, 0.949155, 0.964814, 0.977631, 0.971775, 0.972899, 0.970193, 0.973417, 0.967887, 0.956275, 0.917195, 0.905421, 0.878483, 0.814676, 0.833473, 0.833578, 0.795612, 0.79503, 0.757104, 0.72565, 0.686022, 0.665988, 0.627406, 0.645191, 0.56896, 0.564477, 0.220709}, 
{0.743247, 0.738368, 0.810674, 0.905164, 0.952316, 0.969126, 0.979333, 0.973186, 0.978376, 0.97415, 0.978972, 0.974538, 0.972232, 0.936293, 0.929536, 0.896436, 0.834792, 0.852494, 0.843123, 0.810567, 0.815801, 0.803875, 0.751566, 0.749522, 0.717226, 0.74056, 0.738659, 0.726251, 0.665331, 0.293513}, 
{0.757139, 0.752567, 0.810674, 0.91016, 0.955477, 0.971283, 0.982737, 0.974598, 0.983853, 0.976789, 0.980361, 0.98385, 0.975133, 0.955391, 0.944005, 0.910799, 0.854908, 0.873629, 0.862212, 0.840477, 0.849555, 0.830183, 0.777482, 0.779601, 0.764522, 0.785014, 0.770924, 0.755611, 0.730589, 0.352455}, 
{0.764085, 0.795165, 0.827369, 0.920154, 0.955477, 0.971283, 0.98444, 0.97601, 0.985222, 0.976789, 0.981749, 0.98518, 0.979485, 0.968124, 0.958474, 0.945429, 0.884169, 0.890537, 0.874143, 0.855433, 0.867731, 0.838953, 0.787201, 0.803381, 0.799994, 0.821385, 0.812465, 0.789864, 0.795848, 0.42005}, 
{0.764085, 0.795165, 0.827369, 0.925151, 0.958639, 0.975668, 0.986142, 0.977422, 0.986592, 0.979427, 0.983138, 0.98651, 0.979485, 0.976081, 0.968119, 0.966973, 0.907942, 0.909558, 0.898004, 0.887836, 0.883821, 0.862338, 0.800159, 0.820091, 0.831525, 0.83755, 0.843622, 0.824117, 0.825511, 0.493233}, 
{0.786293, 0.795165, 0.827369, 0.925151, 0.964961, 0.977825, 0.989546, 0.978834, 0.986592, 0.983385, 0.984527, 0.98651, 0.982386, 0.980856, 0.97455, 0.977744, 0.928059, 0.932806, 0.907549, 0.897806, 0.8994, 0.879877, 0.835793, 0.842635, 0.866998, 0.861797, 0.869585, 0.853477, 0.861106, 0.561368}, 
{0.786293, 0.795165, 0.827369, 0.935144, 0.968122, 0.979981, 0.989546, 0.978834, 0.986592, 0.984704, 0.984527, 0.98651, 0.990721, 0.980856, 0.976158, 0.990311, 0.948175, 0.948022, 0.929024, 0.917746, 0.901996, 0.888647, 0.861709, 0.859345, 0.882763, 0.886044, 0.87997, 0.882837, 0.878904, 0.620428}, 
{0.801563, 0.823564, 0.844064, 0.940141, 0.971283, 0.979981, 0.99295, 0.980246, 0.986592, 0.986023, 0.985916, 0.98651, 0.992171, 0.980856, 0.977765, 0.993902, 0.966462, 0.97127, 0.940954, 0.935194, 0.914979, 0.903263, 0.871428, 0.88274, 0.886704, 0.910292, 0.890954, 0.892624, 0.891416, 0.678151}, 
{0.801563, 0.823564, 0.852412, 0.940141, 0.971283, 0.979981, 0.99295, 0.980246, 0.989602, 0.988661, 0.985916, 0.98651, 0.995073, 0.982447, 0.980981, 0.993902, 0.975606, 0.975497, 0.96243, 0.950149, 0.925364, 0.921239, 0.887625, 0.892766, 0.894587, 0.914333, 0.890954, 0.90241, 0.891416, 0.726414}, 
{0.818211, 0.823564, 0.877455, 0.945138, 0.974444, 0.979981, 0.99295, 0.984481, 0.990971, 0.992619, 0.990082, 0.989171, 0.995073, 0.982447, 0.980981, 0.993902, 0.977435, 0.983951, 0.976747, 0.967597, 0.940943, 0.932932, 0.897344, 0.896109, 0.906411, 0.930498, 0.901339, 0.91709, 0.903281, 0.756753}, 
{0.832104, 0.823564, 0.877455, 0.945138, 0.977606, 0.983744, 0.99295, 0.987305, 0.990971, 0.993938, 0.991471, 0.990501, 0.996811, 0.984039, 0.984196, 0.993902, 0.987942, 0.986064, 0.981519, 0.975075, 0.946136, 0.938778, 0.925674, 0.912819, 0.914294, 0.930498, 0.916917, 0.921983, 0.909214, 0.795687}, 
{0.840427, 0.823564, 0.89415, 0.960657, 0.977606, 0.986328, 0.99295, 0.987305, 0.99234, 0.993938, 0.991471, 0.991831, 0.996811, 0.984039, 0.984196, 0.993902, 0.98977, 0.986064, 0.988677, 0.99003, 0.951329, 0.956317, 0.938632, 0.932021, 0.922177, 0.930498, 0.92211, 0.94645, 0.928181, 0.825486}, 
{0.879283, 0.848345, 0.89415, 0.971641, 0.977606, 0.986328, 0.99499, 0.987305, 0.99234, 0.993938, 0.994523, 0.991831, 0.996811, 0.985946, 0.987002, 0.993902, 0.993791, 0.990708, 0.988677, 0.992522, 0.959118, 0.96801, 0.948351, 0.946052, 0.934001, 0.934539, 0.937688, 0.94645, 0.928181, 0.854498}, 
{0.891406, 0.897908, 0.89415, 0.971641, 0.977606, 0.988912, 0.99499, 0.987305, 0.995351, 0.993938, 0.995912, 0.991831, 0.996811, 0.985946, 0.987002, 0.995697, 0.993791, 0.992821, 0.99345, 0.997507, 0.977293, 0.987728, 0.95483, 0.95942, 0.957649, 0.942621, 0.94288, 0.951343, 0.940046, 0.887591}, 
{0.899729, 0.931938, 0.908718, 0.977629, 0.983123, 0.988912, 0.99499, 0.990409, 0.995351, 0.996838, 0.995912, 0.991831, 0.996811, 0.988723, 0.988928, 0.995697, 0.993791, 0.992821, 0.995836, 0.997507, 0.990276, 0.996497, 0.964548, 0.966105, 0.964528, 0.942621, 0.958459, 0.96113, 0.964404, 0.91681}, 
{0.936815, 0.948954, 0.938727, 0.977629, 0.983123, 0.991496, 0.99499, 0.993513, 0.99672, 0.998419, 0.995912, 0.993162, 0.996811, 0.988723, 0.990854, 0.995697, 0.993791, 0.997467, 0.995836, 0.997507, 0.994807, 0.996497, 0.967788, 0.971937, 0.976352, 0.954745, 0.968844, 0.966994, 0.970337, 0.93844}, 
{0.968725, 0.965969, 0.955422, 0.982037, 0.98864, 0.991496, 0.99499, 0.996896, 0.998361, 0.998419, 0.995912, 0.993162, 0.996811, 0.990631, 0.990854, 0.997849, 0.993791, 0.997467, 0.995836, 1, 0.997404, 0.996497, 0.9864, 0.978622, 0.980293, 0.971711, 0.974037, 0.971887, 0.970337, 0.959474}, 
{1, 0.982985, 0.979994, 0.988025, 1, 0.991496, 0.997029, 0.998308, 1, 0.998419, 1, 1, 0.998549, 0.998093, 0.995268, 1, 0.997809, 1, 1, 1, 0.997404, 0.996497, 0.993521, 0.988648, 0.992117, 0.975753, 0.989615, 0.981673, 0.982202, 0.979885}, 
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}
}, 
{
{0.223107, 0.488257, 0.419136, 0.429306, 0.413677, 0.35753, 0.303865, 0.267405, 0.210886, 0.192604, 0.15914, 0.143452, 0.126001, 0.097834, 0.081552, 0.0781217, 0.0717993, 0.0501364, 0.0474308, 0.0364846, 0.0495664, 0.0332603, 0.0231157, 0.0248106, 0.0174661, 0.0207476, 0.00833279, 0.0114473, 0.0177526, 0.00987805}, 
{0.364256, 0.660814, 0.668296, 0.666524, 0.687208, 0.603719, 0.554834, 0.476014, 0.406154, 0.369825, 0.309716, 0.290976, 0.256525, 0.237891, 0.183993, 0.185489, 0.15712, 0.109828, 0.106959, 0.0835306, 0.0967489, 0.0677525, 0.059234, 0.0387666, 0.0365199, 0.0456652, 0.0437471, 0.0462434, 0.0329692, 0.017655}, 
{0.459873, 0.732217, 0.76068, 0.792481, 0.797586, 0.760276, 0.714354, 0.657284, 0.555273, 0.522176, 0.464589, 0.418691, 0.387411, 0.356752, 0.290393, 0.281208, 0.237871, 0.175166, 0.172351, 0.139218, 0.141685, 0.114563, 0.105618, 0.0744319, 0.0698642, 0.0686293, 0.0583295, 0.0668486, 0.0729005, 0.0281141}, 
{0.478086, 0.756018, 0.804979, 0.840764, 0.865969, 0.847764, 0.815688, 0.771943, 0.689372, 0.64653, 0.589844, 0.539018, 0.512205, 0.468358, 0.39743, 0.38535, 0.333095, 0.257444, 0.249832, 0.210267, 0.193361, 0.157679, 0.144626, 0.1163, 0.120675, 0.107492, 0.0979102, 0.0805855, 0.0855809, 0.0427492}, 
{0.496299, 0.767918, 0.811794, 0.857558, 0.891296, 0.890126, 0.875847, 0.847273, 0.77584, 0.737277, 0.696937, 0.64364, 0.610657, 0.571742, 0.509564, 0.484487, 0.434472, 0.353632, 0.353092, 0.292837, 0.243088, 0.221735, 0.205305, 0.176776, 0.164729, 0.14812, 0.131241, 0.10577, 0.103334, 0.0543714}, 
{0.528171, 0.773868, 0.822017, 0.868055, 0.898894, 0.912228, 0.910328, 0.898138, 0.843092, 0.812472, 0.780673, 0.723089, 0.696558, 0.656555, 0.605133, 0.571318, 0.522079, 0.452539, 0.445691, 0.379442, 0.337674, 0.293184, 0.267428, 0.227948, 0.228242, 0.195815, 0.164573, 0.137822, 0.128963, 0.0724163}, 
{0.546384, 0.773868, 0.828832, 0.874353, 0.903959, 0.923279, 0.927936, 0.928399, 0.892142, 0.871083, 0.839332, 0.806775, 0.764001, 0.746499, 0.695818, 0.653362, 0.610699, 0.548686, 0.526952, 0.484096, 0.428669, 0.378182, 0.357286, 0.28067, 0.294931, 0.259349, 0.212486, 0.183612, 0.16954, 0.0911336}, 
{0.564597, 0.773868, 0.842463, 0.876452, 0.906492, 0.929726, 0.938207, 0.941276, 0.920364, 0.912111, 0.890387, 0.8643, 0.825733, 0.81212, 0.771758, 0.729466, 0.71002, 0.639999, 0.62522, 0.59451, 0.547748, 0.486897, 0.469975, 0.392318, 0.368197, 0.345906, 0.285398, 0.27519, 0.250695, 0.126395}, 
{0.578257, 0.773868, 0.849278, 0.884849, 0.909024, 0.93433, 0.944076, 0.947715, 0.938378, 0.938752, 0.924605, 0.907199, 0.868713, 0.855455, 0.839496, 0.799358, 0.787818, 0.715115, 0.723751, 0.681194, 0.64773, 0.596532, 0.589888, 0.525832, 0.474581, 0.466027, 0.383308, 0.36448, 0.329314, 0.169825}, 
{0.587363, 0.779819, 0.852685, 0.889048, 0.91409, 0.935251, 0.947744, 0.949646, 0.943275, 0.948995, 0.946874, 0.935697, 0.893364, 0.88839, 0.874538, 0.842511, 0.843429, 0.77158, 0.783279, 0.758879, 0.725244, 0.697795, 0.698707, 0.640581, 0.604782, 0.563182, 0.520511, 0.501848, 0.415541, 0.231631}, 
{0.596469, 0.791719, 0.856093, 0.893246, 0.915356, 0.939856, 0.949945, 0.951669, 0.947479, 0.952725, 0.957171, 0.949419, 0.925911, 0.903867, 0.903335, 0.868492, 0.876948, 0.827238, 0.824855, 0.801125, 0.799388, 0.763084, 0.765165, 0.7307, 0.684173, 0.660338, 0.618422, 0.614032, 0.549954, 0.304117}, 
{0.619235, 0.797669, 0.856093, 0.893246, 0.920422, 0.942619, 0.952146, 0.954244, 0.950481, 0.954323, 0.959887, 0.957335, 0.9442, 0.926254, 0.92576, 0.8932, 0.899802, 0.856277, 0.861913, 0.84337, 0.827662, 0.807431, 0.814285, 0.791176, 0.769916, 0.723931, 0.701749, 0.703321, 0.626036, 0.393659}, 
{0.628342, 0.80957, 0.856093, 0.895345, 0.921688, 0.94354, 0.952146, 0.954244, 0.950481, 0.955922, 0.963689, 0.96427, 0.954197, 0.93662, 0.942962, 0.906325, 0.914276, 0.871604, 0.878921, 0.869293, 0.857707, 0.849, 0.850036, 0.820639, 0.822314, 0.783991, 0.755913, 0.769716, 0.712263, 0.482049}, 
{0.628342, 0.80957, 0.856093, 0.897861, 0.927697, 0.945382, 0.952146, 0.954888, 0.951081, 0.955922, 0.966405, 0.965853, 0.960625, 0.950859, 0.952519, 0.924101, 0.923417, 0.885317, 0.889315, 0.884655, 0.879052, 0.865259, 0.860149, 0.847, 0.841368, 0.812255, 0.79966, 0.788031, 0.750305, 0.543308}, 
{0.632895, 0.80957, 0.856093, 0.897861, 0.928963, 0.947223, 0.952146, 0.956655, 0.952883, 0.958053, 0.966948, 0.966381, 0.962911, 0.953954, 0.958253, 0.932306, 0.931797, 0.902051, 0.898764, 0.897137, 0.898491, 0.884969, 0.877485, 0.869015, 0.855971, 0.842285, 0.830908, 0.813216, 0.79849, 0.60246}, 
{0.642001, 0.81552, 0.86204, 0.897861, 0.93023, 0.947223, 0.953758, 0.957943, 0.955732, 0.958586, 0.966948, 0.967964, 0.965311, 0.956553, 0.963987, 0.94598, 0.940177, 0.91415, 0.910102, 0.910986, 0.901861, 0.893592, 0.887887, 0.896927, 0.87185, 0.858533, 0.849656, 0.83611, 0.821169, 0.657843}, 
{0.66202, 0.81552, 0.866123, 0.902476, 0.931747, 0.947223, 0.953758, 0.958587, 0.957171, 0.96029, 0.968577, 0.967964, 0.966454, 0.95903, 0.965898, 0.954184, 0.946271, 0.924637, 0.921441, 0.921548, 0.913095, 0.907142, 0.908113, 0.909333, 0.879789, 0.877965, 0.872572, 0.872742, 0.856674, 0.712224}, 
{0.66202, 0.82147, 0.869531, 0.908655, 0.9363, 0.948327, 0.957676, 0.958587, 0.957891, 0.96029, 0.96912, 0.967964, 0.967025, 0.960268, 0.968418, 0.960163, 0.952365, 0.93851, 0.933725, 0.932109, 0.927699, 0.916997, 0.919671, 0.919792, 0.900026, 0.897396, 0.889237, 0.897926, 0.866819, 0.755609}, 
{0.67452, 0.82147, 0.869531, 0.91327, 0.9363, 0.949431, 0.958409, 0.959359, 0.959539, 0.962633, 0.96912, 0.969124, 0.969279, 0.960887, 0.970966, 0.966317, 0.956174, 0.950609, 0.944118, 0.943631, 0.940056, 0.926852, 0.926895, 0.922893, 0.911141, 0.918594, 0.899653, 0.916242, 0.890146, 0.794275}, 
{0.68998, 0.82147, 0.873614, 0.91327, 0.941545, 0.953062, 0.959143, 0.961254, 0.961698, 0.964201, 0.970719, 0.970916, 0.969964, 0.962989, 0.973605, 0.967, 0.965316, 0.954597, 0.952722, 0.951691, 0.945673, 0.935475, 0.935563, 0.929096, 0.918675, 0.929192, 0.910069, 0.929978, 0.902181, 0.819416}, 
{0.710902, 0.828601, 0.894544, 0.915786, 0.946539, 0.955087, 0.959143, 0.962377, 0.961698, 0.965478, 0.971913, 0.971444, 0.969964, 0.962989, 0.974368, 0.970403, 0.971713, 0.964391, 0.960468, 0.957452, 0.95713, 0.951119, 0.939897, 0.94336, 0.926614, 0.936258, 0.918813, 0.939136, 0.914861, 0.8424}, 
{0.738183, 0.835731, 0.898627, 0.920817, 0.95084, 0.960508, 0.961302, 0.96392, 0.964184, 0.967685, 0.973866, 0.972077, 0.969964, 0.966913, 0.978407, 0.971087, 0.975327, 0.970639, 0.966137, 0.964173, 0.961624, 0.957278, 0.942787, 0.951726, 0.932966, 0.943324, 0.935891, 0.943715, 0.930078, 0.869069}, 
{0.767948, 0.855942, 0.90271, 0.923332, 0.953875, 0.961612, 0.965698, 0.966939, 0.965623, 0.972099, 0.975818, 0.976159, 0.972903, 0.967532, 0.980631, 0.972964, 0.978332, 0.97656, 0.969159, 0.966809, 0.96387, 0.962206, 0.948566, 0.956378, 0.941218, 0.957804, 0.944224, 0.946005, 0.93515, 0.889042}, 
{0.784317, 0.894847, 0.90271, 0.923332, 0.95691, 0.964323, 0.967456, 0.966939, 0.967062, 0.972737, 0.979913, 0.976792, 0.974273, 0.970838, 0.982158, 0.97597, 0.980504, 0.978493, 0.971049, 0.970707, 0.969201, 0.966146, 0.95468, 0.96754, 0.945981, 0.97052, 0.95297, 0.952873, 0.93515, 0.904335}, 
{0.823189, 0.894847, 0.929074, 0.937059, 0.962727, 0.96653, 0.974011, 0.968834, 0.973356, 0.974944, 0.982516, 0.980138, 0.97664, 0.977111, 0.983685, 0.977608, 0.982746, 0.979901, 0.97728, 0.979483, 0.976387, 0.971204, 0.964425, 0.973743, 0.951375, 0.972286, 0.963386, 0.965229, 0.942758, 0.922238}, 
{0.855926, 0.933753, 0.937241, 0.937059, 0.971007, 0.972551, 0.976572, 0.973044, 0.977282, 0.980237, 0.987963, 0.982956, 0.97938, 0.982239, 0.985212, 0.983582, 0.983659, 0.980707, 0.978413, 0.983615, 0.979694, 0.97785, 0.970777, 0.980866, 0.957321, 0.977936, 0.969518, 0.967518, 0.960511, 0.940289}, 
{0.876848, 0.940883, 0.955439, 0.955597, 0.972525, 0.976966, 0.983527, 0.976833, 0.981208, 0.986219, 0.989454, 0.98683, 0.982745, 0.984803, 0.985212, 0.984775, 0.987727, 0.986874, 0.982941, 0.983615, 0.985532, 0.985974, 0.978287, 0.982417, 0.966848, 0.981469, 0.977733, 0.979871, 0.965583, 0.956588}, 
{0.914139, 0.962274, 0.971772, 0.969323, 0.978461, 0.984508, 0.987846, 0.982519, 0.986902, 0.989065, 0.990105, 0.992466, 0.985484, 0.989528, 0.988266, 0.989484, 0.991227, 0.988807, 0.988963, 0.990998, 0.991649, 0.990157, 0.983985, 0.99039, 0.977873, 0.983235, 0.986479, 0.985358, 0.970655, 0.969922}, 
{0.943905, 0.99287, 0.981803, 0.981681, 0.994062, 0.988738, 0.993044, 0.990938, 0.994698, 0.996275, 0.994496, 0.99402, 0.99327, 0.994995, 0.991668, 0.995084, 0.995144, 0.99376, 0.99269, 0.992673, 1, 0.994095, 0.987161, 0.993492, 0.983267, 0.991618, 0.992198, 0.99268, 0.987966, 0.982576}, 
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}
}, 
{
{0.401107, 0.654132, 0.684835, 0.517853, 0.541612, 0.461949, 0.400936, 0.303939, 0.253896, 0.200538, 0.173462, 0.141282, 0.110771, 0.0814875, 0.0679642, 0.0581383, 0.0416593, 0.0375804, 0.0310621, 0.0326306, 0.0271854, 0.0206447, 0.0237633, 0.0239188, 0.0222048, 0.0133942, 0.00970684, 0.00957539, 0.0129231, 0.00690365}, 
{0.796679, 0.828251, 0.840762, 0.806971, 0.812756, 0.724235, 0.681406, 0.58216, 0.504244, 0.396112, 0.347372, 0.299307, 0.227046, 0.193445, 0.161818, 0.134505, 0.104678, 0.0889609, 0.0662833, 0.0691079, 0.0611584, 0.0457732, 0.0473102, 0.0378521, 0.0365896, 0.0292519, 0.0290724, 0.0225928, 0.0248034, 0.0142841}, 
{1, 0.957655, 0.929517, 0.934638, 0.929102, 0.880759, 0.845407, 0.773811, 0.687914, 0.591843, 0.514977, 0.453678, 0.360348, 0.306351, 0.265499, 0.220063, 0.171873, 0.153392, 0.118111, 0.11286, 0.100148, 0.0826679, 0.0709232, 0.0596696, 0.0585989, 0.0524131, 0.0486787, 0.0420172, 0.0375514, 0.0233508}, 
{1, 1, 0.972687, 0.976312, 0.972037, 0.948497, 0.941177, 0.892692, 0.826708, 0.746197, 0.663887, 0.591556, 0.504467, 0.430342, 0.382828, 0.326576, 0.250063, 0.232032, 0.177321, 0.165327, 0.14373, 0.12552, 0.0968674, 0.0975818, 0.0821618, 0.0713379, 0.0690974, 0.0596255, 0.054471, 0.03382}, 
{1, 1, 0.972687, 0.976312, 0.984599, 0.981764, 0.972817, 0.941758, 0.908068, 0.855927, 0.783803, 0.711321, 0.637355, 0.555401, 0.500273, 0.436579, 0.345484, 0.321857, 0.256223, 0.23593, 0.196464, 0.175922, 0.148429, 0.13157, 0.11512, 0.102394, 0.09554, 0.0828853, 0.0726668, 0.0442437}, 
{1, 1, 0.972687, 0.982219, 0.994538, 0.991495, 0.988819, 0.971698, 0.956693, 0.91451, 0.880087, 0.817742, 0.745461, 0.677277, 0.62028, 0.549757, 0.461838, 0.420913, 0.364109, 0.325345, 0.265312, 0.241186, 0.205145, 0.180103, 0.155064, 0.138085, 0.118348, 0.112783, 0.0962522, 0.0578757}, 
{1, 1, 0.972687, 0.982219, 0.997023, 0.992852, 0.992874, 0.986018, 0.981723, 0.957642, 0.932204, 0.894617, 0.839688, 0.786087, 0.723864, 0.673891, 0.582346, 0.538748, 0.465207, 0.424285, 0.362493, 0.320932, 0.281176, 0.243757, 0.216397, 0.191537, 0.164968, 0.155823, 0.136582, 0.0757145}, 
{1, 1, 0.972687, 0.982219, 0.997023, 0.992852, 0.996164, 0.99261, 0.992802, 0.982754, 0.967856, 0.947944, 0.913949, 0.8865, 0.818989, 0.784659, 0.704804, 0.653878, 0.581505, 0.53979, 0.479738, 0.415248, 0.388172, 0.328041, 0.297257, 0.273898, 0.239277, 0.218135, 0.185664, 0.0999022}, 
{1, 1, 0.972687, 0.982219, 0.997023, 0.994286, 0.996975, 0.99498, 0.996025, 0.994037, 0.985064, 0.981659, 0.963969, 0.936155, 0.891174, 0.868789, 0.797145, 0.770511, 0.704969, 0.66095, 0.602357, 0.543158, 0.507878, 0.444678, 0.396869, 0.375657, 0.326253, 0.286701, 0.261692, 0.138264}, 
{1, 1, 0.972687, 0.982219, 0.997023, 0.994286, 0.996975, 0.99498, 0.996932, 0.997241, 0.992602, 0.990223, 0.984653, 0.972292, 0.939382, 0.91885, 0.872135, 0.856291, 0.804585, 0.762221, 0.710687, 0.679353, 0.629258, 0.574255, 0.518064, 0.476769, 0.431018, 0.39287, 0.364319, 0.185853}, 
{1, 1, 0.972687, 0.982219, 0.997023, 0.994286, 0.996975, 0.99498, 0.996932, 0.998053, 0.995127, 0.994495, 0.992177, 0.987509, 0.966173, 0.953799, 0.912715, 0.900693, 0.874175, 0.838581, 0.795788, 0.777354, 0.738387, 0.677008, 0.630251, 0.589937, 0.533822, 0.505486, 0.475288, 0.257878}, 
{1, 1, 0.972687, 0.982219, 0.997023, 0.994286, 0.996975, 0.99498, 0.996932, 0.998053, 0.99548, 0.996282, 0.995715, 0.99388, 0.978914, 0.969265, 0.939815, 0.931226, 0.909922, 0.882021, 0.855634, 0.83339, 0.81507, 0.766082, 0.714566, 0.682299, 0.629837, 0.603923, 0.569279, 0.335753}, 
{1, 1, 0.972687, 0.982219, 0.997023, 0.994286, 0.996975, 0.99498, 0.996932, 0.998053, 0.99548, 0.99702, 0.996414, 0.996485, 0.983926, 0.978654, 0.960201, 0.948731, 0.930805, 0.909889, 0.893752, 0.875181, 0.854381, 0.819151, 0.780945, 0.751913, 0.708063, 0.675925, 0.63644, 0.414927}, 
{1, 1, 0.972687, 0.982219, 0.997023, 0.994286, 0.996975, 0.99498, 0.996932, 0.998053, 0.99548, 0.99702, 0.997483, 0.998113, 0.990336, 0.983339, 0.971941, 0.961707, 0.944975, 0.92763, 0.917311, 0.902314, 0.880398, 0.848961, 0.82151, 0.786223, 0.758643, 0.726331, 0.695316, 0.474633}, 
{1, 1, 0.972687, 0.982219, 0.997023, 0.994286, 0.996975, 0.99498, 0.996932, 0.998053, 0.99548, 0.99702, 0.997833, 0.998113, 0.994269, 0.989434, 0.979336, 0.969618, 0.956373, 0.93862, 0.935333, 0.920746, 0.907023, 0.872219, 0.857726, 0.817204, 0.806028, 0.771034, 0.728356, 0.533048}, 
{1, 1, 0.972687, 0.982219, 0.997023, 0.994286, 0.996975, 0.99498, 0.996932, 0.998053, 0.99548, 0.99702, 0.997833, 0.998113, 0.995676, 0.992669, 0.986371, 0.976302, 0.965716, 0.950658, 0.947849, 0.931187, 0.925809, 0.894721, 0.881328, 0.849354, 0.841124, 0.809178, 0.762429, 0.5893}, 
{1, 1, 0.972687, 0.982219, 0.997023, 0.994286, 0.996975, 0.99498, 0.996932, 0.998053, 0.99548, 0.99702, 0.997833, 0.998113, 0.996022, 0.994805, 0.991033, 0.982054, 0.974485, 0.961727, 0.958919, 0.948057, 0.939175, 0.917835, 0.911296, 0.880983, 0.874499, 0.837988, 0.802992, 0.643224}, 
{1, 1, 0.972687, 0.982219, 0.997023, 0.994286, 0.996975, 0.99498, 0.996932, 0.998053, 0.99548, 0.99702, 0.997833, 0.998113, 0.996022, 0.996255, 0.992558, 0.987896, 0.981056, 0.969446, 0.967371, 0.959051, 0.948764, 0.933712, 0.928572, 0.900082, 0.899124, 0.867706, 0.833936, 0.687996}, 
{1, 1, 0.972687, 0.982219, 0.997023, 0.994286, 0.996975, 0.99498, 0.996932, 0.998053, 0.99548, 0.99702, 0.998276, 0.998113, 0.996022, 0.996628, 0.993745, 0.989988, 0.987603, 0.976721, 0.970889, 0.965008, 0.95851, 0.945382, 0.942158, 0.917712, 0.913946, 0.890966, 0.861751, 0.731718}, 
{1, 1, 0.972687, 0.982219, 1, 0.994286, 0.996975, 0.99498, 0.996932, 0.998053, 0.995926, 0.99702, 0.998276, 0.998113, 0.996022, 0.996981, 0.994529, 0.99167, 0.992183, 0.984938, 0.978387, 0.973089, 0.967674, 0.955051, 0.953571, 0.93059, 0.931398, 0.910441, 0.889625, 0.770985}, 
{1, 1, 0.972687, 0.982219, 1, 0.994286, 0.996975, 0.995681, 0.996932, 0.998527, 0.996348, 0.99702, 0.998276, 0.998504, 0.996022, 0.997427, 0.994529, 0.992512, 0.993559, 0.985409, 0.984468, 0.976923, 0.972501, 0.965024, 0.960753, 0.939449, 0.939336, 0.924607, 0.907821, 0.802589}, 
{1, 1, 0.972687, 0.989692, 1, 0.994286, 0.996975, 0.995681, 0.996932, 0.998527, 0.996348, 0.99702, 0.998276, 0.998504, 0.996022, 0.997427, 0.995393, 0.993439, 0.995412, 0.988759, 0.988448, 0.981819, 0.979234, 0.974204, 0.967895, 0.950475, 0.952485, 0.941832, 0.928277, 0.833324}, 
{1, 1, 0.972687, 0.989692, 1, 0.996655, 0.996975, 0.995681, 0.996932, 0.998527, 0.996348, 0.99702, 0.998695, 0.998504, 0.996877, 0.997427, 0.99585, 0.994153, 0.995412, 0.991988, 0.991888, 0.985173, 0.981665, 0.978812, 0.975776, 0.957001, 0.959562, 0.950922, 0.940099, 0.861449}, 
{1, 1, 0.972687, 0.989692, 1, 0.996655, 0.996975, 0.995681, 0.996932, 0.998527, 0.997192, 0.99702, 0.998695, 0.998894, 0.996877, 0.997427, 0.996307, 0.995357, 0.995867, 0.992983, 0.994424, 0.98634, 0.984616, 0.982845, 0.980262, 0.965563, 0.965682, 0.963787, 0.951422, 0.88429}, 
{1, 1, 0.972687, 0.989692, 1, 0.998282, 0.996975, 0.99716, 0.997767, 0.999027, 0.997614, 0.99702, 0.998695, 0.999284, 0.997483, 0.997873, 0.997729, 0.996629, 0.996659, 0.993453, 0.995951, 0.988523, 0.989477, 0.989023, 0.988104, 0.976021, 0.9738, 0.972112, 0.962407, 0.905167}, 
{1, 1, 1, 0.989692, 1, 0.998282, 0.998974, 0.997899, 0.997767, 0.999027, 0.997614, 0.998673, 0.998695, 0.99961, 0.997898, 0.998285, 0.998877, 0.997343, 0.998426, 0.995739, 0.996856, 0.991797, 0.992531, 0.993019, 0.99252, 0.9834, 0.980089, 0.974886, 0.969824, 0.925894}, 
{1, 1, 1, 0.989692, 1, 0.998282, 1, 0.998599, 0.999457, 0.999527, 0.99806, 0.999558, 0.998695, 1, 0.99873, 0.998285, 0.999334, 0.997833, 0.998941, 0.997668, 0.996856, 0.993629, 0.994305, 0.994999, 0.995591, 0.986598, 0.987931, 0.979528, 0.98263, 0.946904}, 
{1, 1, 1, 0.989692, 1, 0.998282, 1, 0.9993, 1, 1, 0.999156, 0.999558, 0.999138, 1, 0.999169, 0.998708, 0.999334, 0.998265, 0.999456, 0.998264, 0.998065, 0.995138, 0.994896, 0.997243, 0.996329, 0.989796, 0.99233, 0.986843, 0.986918, 0.963682}, 
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.999578, 0.999558, 1, 1, 1, 1, 0.999334, 0.999246, 1, 1, 0.999412, 0.999022, 0.998583, 0.998576, 0.997067, 0.994711, 0.997322, 0.99445, 0.991206, 0.981428}, 
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}
}
};
const double wagi_he_all [5][30] = {
{9.02367e-05, 0.00938543, 0.0306496, 0.0902926, 0.203012, 0.305818, 0.349929, 0.373422, 0.398957, 0.39581, 0.392528, 0.412284, 0.415416, 0.385162, 0.354342, 0.363123, 0.339601, 0.332648, 0.328547, 0.322477, 0.286147, 0.263711, 0.236757, 0.210279, 0.200699, 0.152301, 0.19984, 0.181844, 0.118299, 0.0782168}, 
{0.000165185, 0.0139627, 0.041355, 0.109704, 0.222017, 0.312993, 0.364312, 0.395127, 0.41218, 0.397784, 0.397144, 0.402759, 0.376253, 0.399391, 0.381706, 0.375072, 0.365975, 0.327829, 0.331026, 0.297894, 0.294255, 0.253979, 0.230727, 0.212615, 0.185398, 0.167911, 0.143696, 0.136423, 0.113446, 0.0694706}, 
{0.00135231, 0.0142653, 0.0261478, 0.0470999, 0.0854436, 0.135518, 0.187378, 0.235887, 0.277126, 0.314592, 0.320297, 0.352904, 0.366919, 0.3551, 0.377501, 0.381835, 0.394745, 0.377701, 0.365189, 0.372316, 0.374759, 0.367046, 0.353512, 0.361316, 0.3323, 0.342946, 0.287786, 0.336525, 0.286066, 0.17769}, 
{0.000413224, 0.00786946, 0.0151926, 0.0272694, 0.0502502, 0.0756879, 0.107152, 0.133804, 0.162063, 0.200833, 0.217807, 0.250563, 0.260988, 0.269083, 0.284657, 0.297026, 0.29593, 0.300551, 0.294927, 0.310192, 0.294644, 0.310218, 0.292647, 0.294994, 0.307031, 0.309064, 0.295974, 0.284552, 0.283876, 0.187062}, 
{0.000868729, 0.0036663, 0.00572203, 0.0179749, 0.0397606, 0.069585, 0.114765, 0.16091, 0.212802, 0.251982, 0.295101, 0.327512, 0.344281, 0.378636, 0.386682, 0.407135, 0.416779, 0.417526, 0.422171, 0.416552, 0.431648, 0.430345, 0.423333, 0.416413, 0.416814, 0.404488, 0.400741, 0.394343, 0.390526, 0.291449}
};

const double Mrev_he [30][30] =
{
{3296.549105, -2886.565814, -963.945346, 2.006991703, 345.0547419, 404.8243074, 644.8525768, -160.1357024, 41.99376179, -327.6292682, -372.7449309, -280.4141071, -271.2800492, -27.7776367, -133.4901853, 23.80658281, 54.04264886, 130.5473967, 210.8220553, 100.5236602, 75.6099038, 62.54538006, 32.51437964, 132.7080437, 80.42952877, -39.25391127, -51.64013669, 12.90654019, -8.071476809, 96.19175383},
{-2886.565814, 6283.217399, -1537.473099, -984.6538157, -618.5972926, -331.831257, -86.70209592, -58.86244424, -17.52512879, 177.9747748, 115.634592, 132.3465153, -16.56437174, 8.081366663, -2.364100489, -22.13428826, -28.85544286, -15.40465296, 11.08806086, -45.16708344, -14.33110434, -84.25090222, -22.35073022, -28.42232196, 10.13803707, 48.19616368, 83.00221337, 74.89304945, 21.01594607, -15.15500768},
{-963.945346, -1537.473099, 5411.006344, -1284.181431, -889.1258786, -558.4316427, -207.8651974, -185.0064629, 192.7701978, 36.88257043, 130.9137206, 50.34320704, 88.88376086, 53.69021983, -22.5346889, 30.71792485, -40.01989933, -36.51892953, -55.35565665, -56.70582835, -99.40138547, -29.1312497, -57.63993364, 52.41592979, 55.21138356, 38.43874231, 50.6294679, 86.50019766, 57.0530813, -99.30217291},
{2.006991703, -984.6538157, -1284.181431, 4839.452697, -947.4099115, -700.3164256, -522.0307721, -276.0545539, -46.72632469, -33.39573894, 98.89286026, 42.84418326, 83.46444875, 15.71728697, 6.34757166, -0.263031595, -45.99969384, -59.83285674, -57.31941883, -53.15002436, -66.97106955, -14.54062834, -24.86986457, -0.293598188, 19.00631212, 38.34723402, 41.02321636, 50.19357724, 40.77179039, -83.92877562},
{345.0547419, -618.5972926, -889.1258786, -947.4099115, 4197.103808, -723.1945376, -384.4433598, -450.8448994, -477.0167111, 24.8630212, -129.6467811, 118.7619755, -16.15497118, 66.16910972, 30.93638084, -97.31433934, 20.46548729, 36.750398, 9.601593194, -29.42336294, 73.1647504, -84.59473872, -22.32737374, -45.50942842, 6.632338021, 32.85949182, 54.79611847, 35.41516866, -27.95848342, 70.34365035},
{404.8243074, -331.831257, -558.4316427, -700.3164256, -723.1945376, 3978.772929, -618.6431373, -458.839191, -469.6761701, -217.9787909, -173.3578633, -58.84445679, -22.79217988, -27.44648355, 14.83900139, -35.77015187, -1.592267813, -10.84927208, -30.96509166, -1.437988386, 22.13597317, -2.521540369, 4.195513932, -37.761331, -15.30383875, 9.082834924, 13.48631762, -17.54001584, -20.49571625, -8.770406643},
{644.8525768, -86.70209592, -207.8651974, -522.0307721, -384.4433598, -618.6431373, 2796.562973, -135.4077947, -682.0826396, -440.3376917, 46.62462616, -285.9616379, -7.817394977, -206.8722264, 101.8411055, 51.02362085, -94.51621646, -150.1833331, -159.2295476, 3.731262196, -25.62172912, 128.5617005, 134.4912034, -104.8765405, -132.5691769, 34.25528277, -28.24192699, -109.1500501, 40.23578247, -115.9228519},
{-160.1357024, -58.86244424, -185.0064629, -276.0545539, -450.8448994, -458.839191, -135.4077947, 3757.608168, -504.6316079, -471.5899595, -586.5108573, -216.4414003, -193.0783124, -40.04339554, -122.9163822, -78.03856004, 32.06311294, 60.55609769, 62.67933249, 35.19260254, 43.56322183, -14.87425909, -38.99585278, 46.67342426, 56.34913839, -40.21297781, -14.40394942, -1.799860686, -49.47370967, 91.51318268},
{41.99376179, -17.52512879, 192.7701978, -46.72632469, -477.0167111, -469.6761701, -682.0826396, -504.6316079, 2988.744961, -214.2806968, -407.2544245, -68.33677242, -276.1782414, -110.0611431, 36.13002056, -185.8269163, 3.902257919, 34.12520895, 66.74424172, 2.273847318, 195.1113673, -92.21180563, 66.26518344, -177.8188522, -110.2583706, 14.28096183, 41.03939339, -43.13303703, -75.12052328, 214.3239783},
{-327.6292682, 177.9747748, 36.88257043, -33.39573894, 24.8630212, -217.9787909, -440.3376917, -471.5899595, -214.2806968, 2842.573856, -414.373675, -538.7922541, -158.1651437, -154.8760844, -137.8344373, 44.21543897, -59.80830659, -55.7035842, -72.4427624, 36.864462, -93.9916599, 137.8037087, 11.99853964, 106.2373207, 34.1691453, -35.2861938, -87.17836948, -47.90429464, 44.33164374, -86.01391689},
{-372.7449309, 115.634592, 130.9137206, 98.89286026, -129.6467811, -173.3578633, 46.62462616, -586.5108573, -407.2544245, -414.373675, 2505.263351, -254.6151805, -235.509119, -74.96407745, -154.5467609, -122.5947367, -13.72025982, 22.41160649, 18.70814539, 5.970086647, 26.43167784, -24.3802005, -39.02744481, 43.28964328, 50.40349263, -24.62768047, -7.058567594, 6.175687865, -40.96109116, 80.56290676},
{-280.4141071, 132.3465153, 50.34320704, 42.84418326, 118.7619755, -58.84445679, -285.9616379, -216.4414003, -68.33677242, -538.7922541, -254.6151805, 1940.642632, -130.0894069, -141.6305641, -122.8252433, 20.60849077, -86.91047987, -78.92610044, -96.70179585, -0.618455746, -104.8347779, 98.04531753, 3.723283504, 82.45875477, 24.9624799, -15.9973904, -62.0529797, -21.95980073, 53.51561925, -84.32915614},
{-271.2800492, -16.56437174, 88.88376086, 83.46444875, -16.15497118, -22.79217988, -7.817394977, -193.0783124, -276.1782414, -158.1651437, -235.509119, -130.0894069, 1653.965486, -77.9301673, -88.93618844, -106.3680007, -57.33091795, -29.88725321, -3.600737493, -39.56823244, 4.878224169, -59.97804463, -17.64969383, -25.81676895, -11.51030405, -3.789581251, 9.904062685, 16.63340379, 2.40380335, 70.95042097},
{-27.7776367, 8.081366663, 53.69021983, 15.71728697, 66.16910972, -27.44648355, -206.8722264, -40.04339554, -110.0611431, -154.8760844, -74.96407745, -141.6305641, -77.9301673, 1072.813671, -53.68363515, -37.05635493, -73.57148889, -74.99916222, -67.03632852, -34.57736733, -47.95052253, 1.76185288, 1.989204016, -18.95592006, -19.62276138, 4.045829835, -0.13488037, -8.886075954, 17.73445427, -28.27580897},
{-133.4901853, -2.364100489, -22.5346889, 6.34757166, 30.93638084, 14.83900139, 101.8411055, -122.9163822, 36.13002056, -137.8344373, -154.5467609, -122.8252433, -88.93618844, -53.68363515, 1039.118749, -60.18403868, -60.02786442, -52.23429094, -49.36653553, -42.97704423, -61.51106089, -26.23315401, -43.53056214, 26.14813477, 22.13947385, -10.4586414, -2.669773384, 18.85768009, 12.45594352, -11.38422195},
{23.80658281, -22.13428826, 30.71792485, -0.263031595, -97.31433934, -35.77015187, 51.02362085, -78.03856004, -185.8269163, 44.21543897, -122.5947367, 20.60849077, -106.3680007, -37.05635493, -60.18403868, 934.6355888, -54.5190689, -45.36062418, -32.08485867, -66.11326569, 5.746300925, -93.66935984, -33.27639307, -56.47095828, -22.25694866, 2.681422459, 32.57769619, 17.35924222, -19.6211986, 57.8093309},
{54.04264886, -28.85544286, -40.01989933, -45.99969384, 20.46548729, -1.592267813, -94.51621646, 32.06311294, 3.902257919, -59.80830659, -13.72025982, -86.91047987, -57.33091795, -73.57148889, -60.02786442, -54.5190689, 924.0583674, -88.90738795, -83.19028749, -62.22819658, -70.12885281, -33.34028343, -27.69861354, -26.93074039, -26.0806914, -6.918327318, -4.541278569, 0.113796202, 20.28337453, -29.42544221},
{130.5473967, -15.40465296, -36.51892953, -59.83285674, 36.750398, -10.84927208, -150.1833331, 60.55609769, 34.12520895, -55.7035842, 22.41160649, -78.92610044, -29.88725321, -74.99916222, -52.23429094, -45.36062418, -88.90738795, 754.996976, -105.0213349, -64.40815641, -83.3451773, -22.75185045, -22.42518144, -30.17538903, -35.81483204, -6.201763462, -9.12245676, -8.93570029, 22.56722374, -52.04508404},
{210.8220553, 11.08806086, -55.35565665, -57.31941883, 9.601593194, -30.96509166, -159.2295476, 62.67933249, 66.74424172, -72.4427624, 18.70814539, -96.70179585, -3.600737493, -67.03632852, -49.36653553, -32.08485867, -83.19028749, -105.0213349, 693.4611271, -61.14562602, -95.30697626, -10.17829159, -29.02599845, -16.00403609, -29.00837829, -9.698948386, -21.01659883, -18.3156169, 18.95974886, -70.98033492},
{100.5236602, -45.16708344, -56.70582835, -53.15002436, -29.42336294, -1.437988386, 3.731262196, 35.19260254, 2.273847318, 36.864462, 5.970086647, -0.618455746, -39.56823244, -34.57736733, -42.97704423, -66.11326569, -62.22819658, -64.40815641, -61.14562602, 728.2543032, -54.63256511, -70.79407308, -48.79199439, -54.27339474, -44.32510324, -27.3252393, -9.296434096, -7.669338321, -3.506215714, 3.50687471},
{75.6099038, -14.33110434, -99.40138547, -66.97106955, 73.1647504, 22.13597317, -25.62172912, 43.56322183, 195.1113673, -93.9916599, 26.43167784, -104.8347779, 4.878224169, -47.95052253, -61.51106089, 5.746300925, -70.12885281, -83.3451773, -95.30697626, -54.63256511, 627.7420746, -10.44733007, -56.73609939, -1.276570936, -26.08709554, -44.65655586, -48.77079476, -21.48639199, 21.43630291, -71.61907498},
{62.54538006, -84.25090222, -29.1312497, -14.54062834, -84.59473872, -2.521540369, 128.5617005, -14.87425909, -92.21180563, 137.8037087, -24.3802005, 98.04531753, -59.97804463, 1.76185288, -26.23315401, -93.66935984, -33.34028343, -22.75185045, -10.17829159, -70.79407308, -10.44733007, 678.6967497, -74.2998535, -99.81563685, -66.59305479, -51.87672109, -14.56992317, -15.74406422, -43.21381453, 65.46983536},
{32.51437964, -22.35073022, -57.63993364, -24.86986457, -22.32737374, 4.195513932, 134.4912034, -38.99585278, 66.26518344, 11.99853964, -39.02744481, 3.723283504, -17.64969383, 1.989204016, -43.53056214, -33.27639307, -27.69861354, -22.42518144, -29.02599845, -48.79199439, -56.73609939, -74.2998535, 717.3498155, -48.96568859, -44.64524944, -74.46712857, -54.22594629, -33.43372554, -38.67378718, 1.044959767},
{132.7080437, -28.42232196, 52.41592979, -0.293598188, -45.50942842, -37.761331, -104.8765405, 46.67342426, -177.8188522, 106.2373207, 43.28964328, 82.45875477, -25.81676895, -18.95592006, 26.14813477, -56.47095828, -26.93074039, -30.17538903, -16.00403609, -54.27339474, -1.276570936, -99.81563685, -48.96568859, 635.358373, -113.3800116, -67.24641012, -41.63034839, -63.85983875, -63.8828023, 28.26029597},
{80.42952877, 10.13803707, 55.21138356, 19.00631212, 6.632338021, -15.30383875, -132.5691769, 56.34913839, -110.2583706, 34.1691453, 50.40349263, 24.9624799, -11.51030405, -19.62276138, 22.13947385, -22.25694866, -26.0806914, -35.81483204, -29.00837829, -44.32510324, -26.08709554, -66.59305479, -44.64524944, -113.3800116, 705.148319, -80.92288587, -71.13705908, -79.5173994, -56.27949392, -0.910286466},
{-39.25391127, 48.19616368, 38.43874231, 38.34723402, 32.85949182, 9.082834924, 34.25528277, -40.21297781, 14.28096183, -35.2861938, -24.62768047, -15.9973904, -3.789581251, 4.045829835, -10.4586414, 2.681422459, -6.918327318, -6.201763462, -9.698948386, -27.3252393, -44.65655586, -51.87672109, -74.46712857, -67.24641012, -80.92288587, 806.2133205, -104.7139744, -96.28723896, -87.72617306, -33.19165768},
{-51.64013669, 83.00221337, 50.6294679, 41.02321636, 54.79611847, 13.48631762, -28.24192699, -14.40394942, 41.03939339, -87.17836948, -7.058567594, -62.0529797, 9.904062685, -0.13488037, -2.669773384, 32.57769619, -4.541278569, -9.12245676, -21.01659883, -9.296434096, -48.77079476, -14.56992317, -54.22594629, -41.63034839, -71.13705908, -104.7139744, 785.9793239, -105.2230183, -85.72599743, -80.901673},
{12.90654019, 74.89304945, 86.50019766, 50.19357724, 35.41516866, -17.54001584, -109.1500501, -1.799860686, -43.13303703, -47.90429464, 6.175687865, -21.95980073, 16.63340379, -8.886075954, 18.85768009, 17.35924222, 0.113796202, -8.93570029, -18.3156169, -7.669338321, -21.48639199, -15.74406422, -33.43372554, -63.85983875, -79.5173994, -96.28723896, -105.2230183, 816.2217626, -127.1982082, -121.0110877},
{-8.071476809, 21.01594607, 57.0530813, 40.77179039, -27.95848342, -20.49571625, 40.23578247, -49.47370967, -75.12052328, 44.33164374, -40.96109116, 53.51561925, 2.40380335, 17.73445427, 12.45594352, -19.6211986, 20.28337453, 22.56722374, 18.95974886, -3.506215714, 21.43630291, -43.21381453, -38.67378718, -63.8828023, -56.27949392, -87.72617306, -85.72599743, -127.1982082, 955.4251598, -173.7644282},
{96.19175383, -15.15500768, -99.30217291, -83.92877562, 70.34365035, -8.770406643, -115.9228519, 91.51318268, 214.3239783, -86.01391689, 80.56290676, -84.32915614, 70.95042097, -28.27580897, -11.38422195, 57.8093309, -29.42544221, -52.04508404, -70.98033492, 3.50687471, -71.61907498, 65.46983536, 1.044959767, 28.26029597, -0.910286466, -33.19165768, -80.901673, -121.0110877, -173.7644282, 735.9323265}
};
