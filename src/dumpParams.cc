//#include <iostream>
//#include "stdlib.h"
#include "params.h"
#include "args.h"
#include "dirs.h"

int main(int argc, char** argv){
	set_dirs(argv[0]);
	params p;
	args a;
        a.read (argc, argv);
	p.read (a.input);
	p.read (a.params, "command line");
	p.list ();
	return 0;
}
