//#include <iostream>
//#include "stdlib.h"
#include "params.h"
#include "args.h"
#include "dirs.h"

int main(int argc, char** argv){
	set_dirs(argv[0]);
	params p;
	args a;
	p.read (a.input);
	p.read (a.params, "command line");
	p.list ();
	p.makeFormParams();
	return 0;
}
