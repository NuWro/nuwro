#include "nuwro.h"
#include "dis/singlepion.h"
#include "beam.h"
#include <string>

const int nof_dyn = 10; //number of dynamics in NuWro
const int nof_class = 3;

const bool active_class[nof_class] =
{
	false,	//mb_ncpi0
	true,	//NUINT12 A
	true	//NUINT12 antiA
};

class histogram
{
	public:
	
		string name;
		string data;
		string xlabel;
		string ylabel;
		
		int bins_x;
		double width_x;
		double normalization; // 0 - cross section		
		
		virtual void save () {};
		virtual void plot () {};
		virtual void finalize (double *norm) {};
		
		histogram (string name_of_hist, string data_file, string x_name, string y_name, int nof_bins, double begin, double end, double norma)
			: name (name_of_hist), data (data_file), xlabel (x_name), ylabel (y_name), bins_x (nof_bins), normalization (norma)
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
		
		hist1D (string name_of_hist, string data_file, string x_name, string y_name, int nof_bins, double begin, double end, double norma);
		~hist1D ();
		
		void save ();
		void plot ();
		void finalize (double *norm);
};

class hist2D : public histogram
{
	public:
	
		string zlabel;
		
		int bins_y;
		double width_y;
		
		double **part_result [nof_dyn];
		double **result;
		
		hist2D (string name_of_hist, string data_file, string x_name, string y_name, string z_name, int nof_bins_x, double begin_x, double end_x, int nof_bins_y, double begin_y, double end_y, double norma);
		~hist2D ();
		
		void save ();
		void plot ();
		void finalize (double *norm);
};

class pattern
{
	protected:
	
		string name;
		NuWro *N;
		params P;
		int events;
		double norm [nof_dyn];

		void run ();
		virtual void set_params () {}
		virtual void calculate(event *e) {}
		virtual void finish () {}
		
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
		void finish ();		
		
	public:
	
		mb_ncpi0 () : pattern ("MB NC pi0 production", 1000000)
		{
			h1 = new hist1D ("mb_ncpi0", "data/data/mb_ncpi0.txt", "{/Symbol p}^{0} momentum [MeV]", "d{/Symbol s} / dp_{{/Symbol p}^{0}} [cm^{2} / MeV]", 30, 0, 1500, 0);
			h2 = new hist2D ("mb_ncpi0_3d", "", "{/Symbol p}^{0} momentum [MeV]", "cos{/Symbol q}", "d{/Symbol s} / dp_{{/Symbol p}^{0}} / dcos{/Symbol q} [cm^{2} / MeV]", 30, 0, 1500, 20, -1, 1, 0);
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
		void finish ();
		
	public:
	
		A () : pattern ("NUINT12 A", 5000000)
		{
			h1 = new hist1D ("A21", "", "#proton", "#events", 10, 0, 10, 100);
			h2 = new hist1D ("A23", "", "{/Symbol n}_{/Symbol m} reconstructed energy [MeV]", "#events", 40, 0, 2000, 100);
			h3 = new hist2D ("A31", "", "{/Symbol p}^{0} energy [MeV]", "cos{/Symbol q}", "#events", 40, 0, 2000, 20, -1, 1, 100);
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
			h1 = new hist1D ("A22", "", "#proton", "#events", 10, 0, 10, 100);
			h2 = new hist1D ("A24", "", "anti {/Symbol n}_{/Symbol m} reconstructed energy [MeV]", "#events", 40, 0, 2000, 100);
			h3 = new hist2D ("A32", "", "{/Symbol p}^{0} energy [MeV]", "cos{/Symbol q}", "#events", 40, 0, 2000, 20, -1, 1, 100);
		}
		
		~antiA ()
		{
			delete h1;
			delete h2;
			delete h3;
		}	
};

pattern * choose (int x);
void run_command (string com);
