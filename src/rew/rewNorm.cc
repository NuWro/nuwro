#include "../params.h"
#include "rewparams.h"
#include "../event1.h"
#include "../nucleus.h"

double calcNorm(event & e, params &p, nucleus &t)
{
	double w=1;

	if(e.flag.qel) w*=rew.qelNorm();
	if(e.flag.res) w*=rew.resNorm();
	if(e.flag.dis) w*=rew.disNorm();
	if(e.flag.coh) w*=rew.cohNorm();
	if(e.flag.mec) w*=rew.mecNorm();
	if(e.flag.cc)  w*=rew.ccNorm ();
	if(e.flag.nc)  w*=rew.ncNorm ();
	if(e.flag.anty)w*=rew.antyNorm();
	
    if(0<=e.dyn && e.dyn<10)
	  w*=(&rew.dynNorm0)[e.dyn].val; 

	return w;	
}
