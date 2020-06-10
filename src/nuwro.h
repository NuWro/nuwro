#ifndef _nuwro_h_
#define _nuwro_h_
#include "event1.h"
#include "params.h"
#include "args.h"
#include <ostream>
#include "chooser.h"
#include "input_data.h"

// Forward declarations
class geomy;
class beam;
class nucleus;
class target_mixer;

namespace NUWRO {
class NuWro
{
	public:
		geomy* make_detector();
		void makeevent(event* e);
		void finishevent(event* e);
		void raport(double i, double n, const char* text, int precision=1000, int k=-1, bool toFile=false);
		void init  ();
		void init  (int argc, char **argv);
		void test_events();
		void user_events();
		void UserAction();
		void real_events();
		bool simulate_event(event*e, int k=-1);
		void kaskada_redo(string input, string output);
		void set_param_file(const char* fn);
		// Need to include function definition in header b/c c++ weirdness
		template <typename T>
		void set_param(string key, T value) {
		    stringstream param_line;
		    param_line << key <<"="<< value << endl;
		    a.add_param(param_line.str());
		}
		void set_param(string key, double x, double y, double z);
		void set_path(char* path);
		void main (int argc, char **argv);
		inline int proces() {return _procesy.choose();}
		void set (params &par);
		void refresh_target ();
		void refresh_dyn ();
		void pot_report(ostream&);
		NuWro ();
		~NuWro();
		params p;

	private:
		args a;
		chooser _procesy;
		ofstream _progress;
		geomy *_detector;
		beam *_beam;
		nucleus *_nucleus;
		target_mixer *_mixer;
		bool dismode;
		input_data input;
};
}
#endif
