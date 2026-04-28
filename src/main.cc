#include "nuwro.h"
#include "output.h"

NuWro nuwro;

#include "tests.h"
int main(int argc, char** argv)
{
  print_nuwro_banner(VERSION);

  nuwro.main(argc,argv);

  return 0;
}
