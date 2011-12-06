#ifndef __dirs__
#define __dirs__

#include <fstream>

using namespace std;

string get_data_dir ();
string get_bin_dir ();

void set_dirs (char* exename);

bool open_data_file(ifstream & file, string name);

#endif
