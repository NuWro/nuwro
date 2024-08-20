#ifndef _nuwro_h_
#define _nuwro_h_
#include "event1.h"
#include "params.h"
#include "args.h"
#include "geomy.h"
#include "target_mixer.h"
#include "beam.h"
#include <ostream>
#include "nucleus.h"
#include "chooser.h"
#include "input_data.h"
#include "metropolis.h"

class NuWro
{
	public:
		geomy* make_detector(params &p);
		void makeevent(event* e, params &p);
		void finishevent(event* e, params &p);
		void raport(double i, double n, const char* text, int precision=1000, int k=-1, string label="", bool toFile=false);
		void init  (int argc, char **argv);
		void test_events(params &p);
		void user_events(params &p);
		void UserAction(params& p);
		void real_events(params &p);
		void kaskada_redo(string input, string output);
		void main (int argc, char **argv);
		inline int proces() {return _procesy.choose();}
		void set (params &p);
		void refresh_target (params &p);
		void refresh_dyn (params &p);
		void pot_report(ostream&, bool format);
		NuWro ();
		~NuWro();
		// Metropolis-Hastings algorithm:
		void initialize_dynamics_list();
		void real_events_mh(params &p);
		event get_event();

	private:
		params p;
		args a;
		chooser _procesy;
		ofstream _progress;
		geomy *_detector;
		beam *_beam;
		nucleus *_nucleus;
		target_mixer *_mixer;
		bool dismode;
		input_data input;
		// for metropolis-hastings algorithm
		//  - list of enabled dynamical models
		std::vector<int> enabled_dyns{};
		//  - Sampler
		Metropolis<event, size_t> sampler;
		//  - state used to calculate xsec and acceptance rate
		std::vector<double> channel_sampleing_weight{};
		std::vector<double> channel_weight_sum{}, channel_weight_sum_fraction{};
		std::unordered_map<int, double> channel_count_final{};
		size_t accepted_count{};
		bool accept{false};
};

extern NuWro nuwro;
#endif
