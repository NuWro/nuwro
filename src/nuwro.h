#ifndef _nuwro_h_
#define _nuwro_h_
#include "../data/F1F209/F1F209Wrapper.hh"
#include "args.h"
#include "beam.h"
#include "chooser.h"
#include "event1.h"
#include "geomy.h"
#include "input_data.h"
#include "nucleus.h"
#include "params.h"
#include "target_mixer.h"
#include <ostream>

class NuWro {
public:
  geomy *make_detector(params &p);
  void makeevent(event *e, params &p);
  void finishevent(event *e, params &p);
  void raport(double i, double n, const char *text, int precision = 1000,
              int k = -1, string label = "", bool toFile = false);
  void init(int argc, char **argv);
  void test_events(params &p);
  void user_events(params &p);
  void UserAction(params &p);
  void real_events(params &p);
  void kaskada_redo(string input, string output);
  void main(int argc, char **argv);
  inline int proces() { return _procesy.choose(); }
  void set(params &p);
  void refresh_target(params &p);
  void refresh_dyn(params &p);
  void pot_report(ostream &, bool format);
  NuWro();
  ~NuWro();

private:
  params p;
  args a;
  chooser _procesy;
  ofstream _progress;
  geomy *_detector;
  beam *_beam;
  nucleus *_nucleus;
  F1F209Wrapper *BostedFit;
  target_mixer *_mixer;
  bool dismode;
  input_data input;
};

extern NuWro nuwro;
#endif
