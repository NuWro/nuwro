#include "nuwro.h"
#include "dis/singlepion.h"
#include "beam.h"
#include <string>

const int nof_dyn = 10; //number of dynamics in NuWro
const int nof_class = 33;

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
	1	//MEC NC
};

void run_command (string com)
{
	FILE *fp;
	fp = popen (com.c_str(), "w");
	pclose (fp);
}

class histogram
{
	public:
	
		string name;
		string data;
		string title;
		string xlabel;
		string ylabel;
		
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
		void put (double x, int dyn, double weight);
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
	
		string name;
		NuWro *N;
		params P;
		int events;

		void run ();
		virtual void set_params () {}
		virtual void calculate(event *e) {}
		
	public:
		
		pattern (string name_of_class, int nof_events) : name (name_of_class), events (nof_events)
		{
			for (int i  = 0; i < nof_dyn; i++)
			{
				norm [i] = 0;
				xsec [i] = 0;
			}
		}
		
		virtual ~pattern () {}		
		void start ();
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
	
		test_mec_nc () : pattern ("MEC", 5000000)
		{
			h1 = new mhist1D ("mec_nc", "", "{/Symbol n}_{/Symbol m}, E = 0 - 10 GeV, on C", "neutrino energy [GeV]", "xsec [cm^2]", 100, 0, 10, 3, 6);
			h1 -> cnames[0] = "QE + MEC";
			h1 -> cnames[1] = "only QE";
			h1 -> cnames[2] = "only MEC";
			h2 = new mhist1D ("mec_nc_anti", "", "anti {/Symbol n}_{/Symbol m}, E = 0 - 10 GeV, on C", "anti neutrino energy [GeV]", "xsec [cm^2]", 100, 0, 10, 3, 6);
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

pattern * choose (int x);
void run_command (string com);
