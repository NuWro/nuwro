#include "beam.h"
#include "params.h"
#include "beam_uniform.h"

params p;


int main()
{ p.read("params.txt");
  p.list();
  beam_uniform b(p);
  for(int i=0;i<15000;i++)
     cout<< b.shoot()<<endl;

}
