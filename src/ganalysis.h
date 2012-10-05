#include "nuwro.h"
#include "dis/singlepion.h"
#include "beam.h"
#include <string>

const int nof_dyn = 10; //number of dynamics in NuWro
const int nof_class = 24;

double norm [nof_dyn];

const bool active_class[nof_class] =
{
	0,	//mb_ncpi0
	1,	//NUINT12 A
	1,	//NUINT12 antiA
	1,	//NUINT12 B
	1,	//NUINT12 antiB
	1,	//NUINT12 C1
	1,	//NUINT12 antiC1
	1,	//NUINT12 C2
	1,	//NUINT12 antiC2
	1,	//NUINT12 D
	1,	//NUINT12 E
	1,	//NUINT12 F
	1,	//NUINT12 G1
	1,	//NUINT12 G1_2
	1,	//NUINT12 G1_3
	1,	//NUINT12 G1_4 (src/particle.h:214: void particle::set_energy(double): Assertion `E>=_mass' failed.)
	1,	//NUINT12 G1_5
	1,	//NUINT12 G1_6
	1,	//NUINT12 G2
	1,	//NUINT12 G2_2
	1,	//NUINT12 G2_3
	1,	//NUINT12 G2_4
	1,	//NUINT12 H
	1	//NUINT12 antiH
};

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
		double normalization; // 0 - cross section		
		
		virtual void save () {};
		virtual void plot () {};
		virtual void finalize () {};
		
		histogram (string name_of_hist, string data_file, string plot_title, string x_name, string y_name, int nof_bins, double begin, double end, double norma)
			: name (name_of_hist), data (data_file), title (plot_title), xlabel (x_name), ylabel (y_name), bins_x (nof_bins), normalization (norma)
		{
			width_x = (end - begin) / bins_x;
		}
		
		virtual ~histogram () {}
}; 

class hist1D : public histogram
{
	public:
	
		double *part_result [nof_dyn];
		double *result;
		
		hist1D (string name_of_hist, string data_file, string plot_title, string x_name, string y_name, int nof_bins, double begin, double end, double norma);
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
		
		double **part_result [nof_dyn];
		double **result;
		
		hist2D (string name_of_hist, string data_file, string plot_title, string x_name, string y_name, string z_name, int nof_bins_x, double begin_x, double end_x, int nof_bins_y, double begin_y, double end_y, double norma);
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
		
		mhist1D (string name_of_hist, string data_file, string plot_title, string x_name, string y_name, int nof_bins, double begin, double end, double norma, int nof_cases);
		~mhist1D ();
		
		void finalize ();
		void save ();
		void plot ();
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
				norm [i] = 0;
		}
		
		virtual ~pattern () {}		
		void start ();
};

class mb_ncpi0 : public pattern
{
	protected:
		
		hist1D *h1;
		hist2D *h2;
		
		void set_params ();
		void calculate (event *e);
		
	public:
	
		mb_ncpi0 () : pattern ("MB NC pi0 production", 1000000)
		{
			h1 = new hist1D ("mb_ncpi0", "data/data/mb_ncpi0.txt", "MiniBooNE NC {/Symbol p}^{0} production on CH_{2}", "{/Symbol p}^{0} momentum [MeV]", "d{/Symbol s} / dp_{{/Symbol p}^{0}} [cm^{2} / MeV]", 30, 0, 1500, 0);
			h2 = new hist2D ("mb_ncpi0_3d", "", "", "{/Symbol p}^{0} momentum [MeV]", "cos{/Symbol q}", "d{/Symbol s} / dp_{{/Symbol p}^{0}} / dcos{/Symbol q} [cm^{2} / MeV]", 30, 0, 1500, 20, -1, 1, 0);
		}
		
		~mb_ncpi0 ()
		{
			delete h1;
			delete h2;
		}
};

class A : public pattern
{
	protected:
	
		hist1D *h1;
		hist1D *h2;
		hist2D *h3;
		
		void set_params ();
		void calculate (event *e);
		
	public:
	
		A () : pattern ("NUINT12 A", 5000000)
		{
			h1 = new hist1D ("A21", "", "{/Symbol n}_{/Symbol m}, E = 2 GeV, CC, on oxygen, no mesons in FS", "#proton", "#events", 10, 0, 10, 100);
			h2 = new hist1D ("A23", "", "{/Symbol n}_{/Symbol m}, E = 2 GeV, CC, on oxygen, no mesons in FS", "{/Symbol n}_{/Symbol m} reconstructed energy [MeV]", "#events", 40, 0, 2000, 100);
			h3 = new hist2D ("A31", "", "{/Symbol n}_{/Symbol m}, E = 2 GeV, NC, on oxygen, no mesons in FS", "{/Symbol p}^{0} energy [MeV]", "cos{/Symbol q}", "#events", 40, 0, 2000, 20, -1, 1, 100);
		}
		
		~A ()
		{
			delete h1;
			delete h2;
			delete h3;
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
			h1 -> title = "anti {/Symbol n}_{/Symbol m}, E = 2 GeV, CC, on oxygen, no mesons in FS";
			h2 -> name = "A24";
			h2 -> xlabel = "anti {/Symbol n}_{/Symbol m} reconstructed energy [MeV]";
			h2 -> title = "anti {/Symbol n}_{/Symbol m}, E = 2 GeV, CC, on oxygen, no mesons in FS";
			h3 -> name = "A32";
			h3 -> title = "anti {/Symbol n}_{/Symbol m}, E = 2 GeV, NC, on oxygen, no mesons in FS";	
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
			h1 = new hist1D ("B11", "", "{/Symbol n}_{/Symbol m}, E = 5 GeV, CC, on iron, at least 1 {/Symbol p}^{0} in FS", "{/Symbol p}^{0} energy / energy transfer", "probabilty", 20, 0, 1, 1);
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

class C1 : public pattern
{
	protected:
	
		mhist1D *h1;
		
		void set_params ();
		void calculate (event *e);
		
	public:
	
		C1 () : pattern ("NUINT12 C1", 5000000)
		{
			h1 = new mhist1D ("C11C12", "", "{/Symbol n}_{/Symbol m}, E = 0-10 GeV, CC, on argon", "neutrino energy [GeV]", "total xsec per Argon", 10, 0, 10, 100, 2);
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
			h1 = new mhist1D ("C2122", "", "{/Symbol n}_{/Symbol m}, E = 2,5 GeV, CC, on argon", "#protons", "#events", 10, 0, 10, 100, 2);
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
		mhist1D *h3;
		
		void set_params ();
		void calculate (event *e);
		
	public:
	
		D () : pattern ("NUINT12 D", 5000000)
		{
			h1 = new mhist1D ("D1112", "", "{/Symbol n}_{/Symbol m}, E = 3 GeV, CC, on argon, no meson in FS", "#protons", "#events", 10, 0, 10, 100, 2);
			h1 -> cnames [0] = "all protons (D11)";
			h1 -> cnames [1] = "only protons with Tk > 50 MeV (D12)";
			
			h2 = new mhist1D ("D13", "", "{/Symbol n}_{/Symbol m}, E = 3 GeV, CC, on argon, no meson in FS", "momentum of most energetic proton [MeV]", "d{/Symbol s}/dp", 30, 0, 3000, 0, 5);
			h2 -> cnames [0] = "1 proton";
			h2 -> cnames [1] = "2 proton";
			h2 -> cnames [2] = "3 proton";
			h2 -> cnames [3] = "4 proton";
			h2 -> cnames [4] = "more protons";

			h3 = new mhist1D ("D14", "", "{/Symbol n}_{/Symbol m}, E = 3 GeV, CC, on argon, no meson in FS", "cos{/Symbol q} of most energetic proton", "d{/Symbol s}/d{/Symbol q}", 20, -1, 1, 0, 5);
			h3 -> cnames [0] = "1 proton";
			h3 -> cnames [1] = "2 proton";
			h3 -> cnames [2] = "3 proton";
			h3 -> cnames [3] = "4 proton";
			h3 -> cnames [4] = "more protons";
		}
		
		~D()
		{
			delete h1;
			delete h2;
			delete h3;
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
			
			h3 -> name = "E14";
			h3 -> title = "{/Symbol n}_{/Symbol m}, E = 1 GeV, CC, on argon, no meson in FS";
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
	
		F () : pattern ("NUINT12 F", 10000000)
		{
			h1 = new mhist1D ("Fnu", "", "{/Symbol n}_{/Symbol m}, E = 0.2 - 5 GeV, on CH_{2}", "neutrino xsec", "neutrino energy [GeV]", 50, 0, 5, 0, 4);
			h1 -> cnames [0] = "CC with MEC (F11)";
			h1 -> cnames [1] = "NC with MEC (F12)";
			h1 -> cnames [2] = "CC without MEC (F13)";
			h1 -> cnames [3] = "NC without MEC (F14)";

			h2 = new mhist1D ("Fanu", "", "anti {/Symbol n}_{/Symbol m}, E = 0.2 - 5 GeV, on CH_{2}", "antineutrino xsec", "neutrino energy [GeV]", 50, 0, 5, 0, 4);
			h2 -> cnames [0] = "CC with MEC (F11)";
			h2 -> cnames [1] = "NC with MEC (F12)";
			h2 -> cnames [2] = "CC without MEC (F13)";
			h2 -> cnames [3] = "NC without MEC (F14)";

			h3 = new mhist1D ("F", "", "E = 0.2 - 5 GeV, on CH_{2}", "neutrino xsec / antineutrino xsec", "neutrino energy [GeV]", 50, 0, 5, 0, 4);
			h3 -> cnames [0] = "CC with MEC (F11)";
			h3 -> cnames [1] = "NC with MEC (F12)";
			h3 -> cnames [2] = "CC without MEC (F13)";
			h3 -> cnames [3] = "NC without MEC (F14)";
		}
		
		~F ()
		{
			ratio ();
			
			delete h1;
			delete h2;
			delete h3;
		}
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
			h1 = new hist2D ("G11", "", "{/Symbol n}_{/Symbol m}, E = 600 MeV, CC QE, on CH_{2}", "{/Symbol m} kinetic energy [MeV]", "cos{/Symbol q}", "d{/Symbol s} / dTk / dcos{/Symbol q}", 40, 0, 1200, 20, -1, 1, 0);
			h2 = new hist2D ("G11b", "", "{/Symbol n}_{/Symbol m}, E = 600 MeV, CC QE + MEC, on CH_{2}", "{/Symbol m} kinetic energy [MeV]", "cos{/Symbol q}", "d{/Symbol s} / dTk / dcos{/Symbol q}", 40, 0, 1200, 20, -1, 1, 0);
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
			h1 = new mhist1D ("G2121b", "", "{/Symbol n}_{/Symbol m}, E = 800 MeV, CC, on CH_{2}", "total kinetic energy of all protons", "#events", 40, 0, 800, 100, 2);
			h1 -> cnames[0] = "only QE (G21)";
			h1 -> cnames[1] = "QE + MEC (G21b)";
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
			h1 -> cnames[0] = "only QE (G22)";
			h1 -> cnames[1] = "QE + MEC (G22b)";
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
			h1 -> cnames[0] = "only QE (G23)";
			h1 -> cnames[1] = "QE + MEC (G23b)";
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
			h1 -> cnames[0] = "only QE (G24)";
			h1 -> cnames[1] = "QE + MEC (G24b)";
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
		hist2D *h18;
		
		void set_params ();
		void calculate (event *e);
		
	public:
	
		H () : pattern ("NUINT H", 5000000)
		{
			h1 = new hist1D ("H11", "", "{/Symbol n}_{/Symbol m}, E = 5 GeV, CC, on carbon, no pions in FS", "leadin proton kinetic energy [GeV]", "#events", 100, 0, 5, 100);
			h2 = new hist1D ("H21", "", "{/Symbol n}_{/Symbol m}, E = 5 GeV, CC, on carbon, no pions in FS", "sum of all proton kinetic energy [GeV]", "#events", 100, 0, 5, 100);
			h3 = new hist2D ("H31", "", "{/Symbol n}_{/Symbol m}, E = 5 GeV, CC, on carbon, {/Symbol p}^{+} in FS", "{/Symbol p}^{+} energy [GeV]", "cos{/Symbol q}", "d{/Symbol s} / dTk / dcos{/Symbol q}", 100, 0, 5, 20, -1, 1, 0);
			h4 = new hist2D ("H41", "", "{/Symbol n}_{/Symbol m}, E = 5 GeV, CC, on carbon, {/Symbol p}^{-} in FS", "{/Symbol p}^{-} energy [GeV]", "cos{/Symbol q}", "d{/Symbol s} / dTk / dcos{/Symbol q}", 100, 0, 5, 20, -1, 1, 0);
			h5 = new hist2D ("H51", "", "{/Symbol n}_{/Symbol m}, E = 5 GeV, CC, on carbon, {/Symbol p}^{0} in FS", "{/Symbol p}^{0} energy [GeV]", "cos{/Symbol q}", "d{/Symbol s} / dTk / dcos{/Symbol q}", 100, 0, 5, 20, -1, 1, 0);
			h6 = new hist2D ("H61", "", "{/Symbol n}_{/Symbol m}, E = 5 GeV, CC, on carbon, at least 2 {/Symbol p} in FS", "leading {/Symbol p}^{+} energy [GeV]", "cos{/Symbol q}", "d{/Symbol s} / dTk / dcos{/Symbol q}", 100, 0, 5, 20, -1, 1, 0);
			h7 = new hist2D ("H71", "", "{/Symbol n}_{/Symbol m}, E = 5 GeV, CC, on carbon, at least 2 {/Symbol p} in FS", "leading {/Symbol p}^{-} energy [GeV]", "cos{/Symbol q}", "d{/Symbol s} / dTk / dcos{/Symbol q}", 100, 0, 5, 20, -1, 1, 0);
			h8 = new hist2D ("H81", "", "{/Symbol n}_{/Symbol m}, E = 5 GeV, CC, on carbon, at least 2 {/Symbol p} in FS", "leading {/Symbol p}^{0} energy [GeV]", "cos{/Symbol q}", "d{/Symbol s} / dTk / dcos{/Symbol q}", 100, 0, 5, 20, -1, 1, 0);
			h9 = new hist1D ("H91", "", "{/Symbol n}_{/Symbol m}, E = 5 GeV, CC, on carbon, at least 2 {/Symbol p} in FS", "sum of all {/Symbol p}^{+} energies [GeV]", "#events", 100, 0, 5, 100);
			h10 = new hist1D ("H101", "", "{/Symbol n}_{/Symbol m}, E = 5 GeV, CC, on carbon, at least 2 {/Symbol p} in FS", "sum of all {/Symbol p}^{-} energies [GeV]", "#events", 100, 0, 5, 100);
			h11 = new hist1D ("H111", "", "{/Symbol n}_{/Symbol m}, E = 5 GeV, CC, on carbon, at least 2 {/Symbol p} in FS", "sum of all {/Symbol p}^{0} energies [GeV]", "#events", 100, 0, 5, 100);
			h13 = new hist2D ("H131", "", "{/Symbol n}_{/Symbol m}, E = 5 GeV, CC, on carbon", "energy transfer [GeV]", "(sum of all {/Symbol p}^{0} energies) / (energy transfer)", "#events", 100, 0, 5, 10, 0, 1, 100);
			h19 = new hist2D ("H191", "", "{/Symbol n}_{/Symbol m}, E = 5 GeV, CC, on carbon", "energy transfer [GeV]", "(sum of all {/Symbol p}^{+} energies) / (energy transfer)", "#events", 100, 0, 5, 10, 0, 1, 100);
			h20 = new hist2D ("H201", "", "{/Symbol n}_{/Symbol m}, E = 5 GeV, CC, on carbon", "energy transfer [GeV]", "(sum of all {/Symbol p}^{-} energies) / (energy transfer)", "#events", 100, 0, 5, 10, 0, 1, 100);
			h14 = new hist2D ("H141", "", "{/Symbol n}_{/Symbol m}, E = 5 GeV, CC, on carbon", "energy transfer [GeV]", "{Symbol p}^{0} multiplicity", "#events", 100, 0, 5, 10, 0, 10, 100);
			h15 = new hist2D ("H151", "", "{/Symbol n}_{/Symbol m}, E = 5 GeV, CC, on carbon", "energy transfer [GeV]", "{Symbol p}^{+} multiplicity", "#events", 100, 0, 5, 10, 0, 10, 100);
			h16 = new hist2D ("H161", "", "{/Symbol n}_{/Symbol m}, E = 5 GeV, CC, on carbon", "energy transfer [GeV]", "{Symbol p}^{-} multiplicity", "#events", 100, 0, 5, 10, 0, 10, 100);
			h17 = new hist2D ("H171", "", "{/Symbol n}_{/Symbol m}, E = 5 GeV, CC, on carbon", "energy transfer [GeV]", "(total kinetic energy in all neutrons)/(energy transfer)", "#events", 100, 0, 5, 10, 0, 1, 100);
			h18 = new hist2D ("H181", "", "{/Symbol n}_{/Symbol m}, E = 5 GeV, CC, on carbon", "energy transfer [GeV]", "(total kinetic energy in all protons)/(energy transfer)", "#events", 100, 0, 5, 10, 0, 1, 100);
			
			
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
			delete h18;
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
			h18 -> title = "anti {/Symbol n}_{/Symbol m}, E = 5 GeV, CC, on carbon";
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
			h18 -> name = "H182";
			
		}
		
		~antiH () {}
};

pattern * choose (int x);
void run_command (string com);
