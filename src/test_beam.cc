#include "beam_uniform.h"
#include "params.h"

params p;


int main()
{ //p.read("params.txt");
  p.beam_type=0;
  p.beam_particle=14;
  p.beam_energy="0 3 1 2 3";
  p.list();
  beam_uniform b(p);
  particle c;
  int ile[]={0,0,0};
  for(int i=0;i<22000000;i++)
   {  c=b.shoot(true);
//      cout<< c <<endl;
      ile[int(c.t)]++;
   }
   cout<<ile[0]<<' '<<ile[1]<<' '<<ile[2]<<endl;
}
