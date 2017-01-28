#ifndef REWEIGHTERS_H
#define REWEIGHTERS_H

#include "../event1.h"
#include "rewparams.h"
#include "../params.h"
#include "../ff.h"
#include "../nucleus.h"

typedef double (*RewCalc)(event & e, params &p, nucleus &t);

double calcDummy(event & e, params &p, nucleus &t);
double calcNorm(event & e, params &p, nucleus &t); 
double calcQEL(event & e, params &p, nucleus &t);
double calcRES(event & e, params &p, nucleus &t);


struct Reweighter
{
	bool active;
	string name;
	RewCalc calc; 

public:
	Reweighter(string name0="",RewCalc calc0=calcDummy):active(false),name(name0),calc(calc0){}
};


struct Reweighters
{	
	Reweighter rEnd;

	Reweighter rewNorm;
	Reweighter rewQEL;
	Reweighter rewRES;

	Reweighter End;

public:

	Reweighters():
		rEnd(),
		rewNorm("rewNorm",calcNorm),
		rewQEL("rewQEL",calcQEL),
		rewRES("rewRES",calcRES),
		End()
		{}

	Reweighter * begin(){return &rEnd+1;}
	Reweighter * end(){return &End;}

	double weight(event &e ,params &p, nucleus &t)
	{
		double w=1;
		for(Reweighter* r=begin();r!=end();r++)
		{	
			if(r->active)
				w*=r->calc(e,p,t);
		}
		return w;
	}

	Reweighter & operator() (string name)
	{
		for(Reweighter* r=begin();r!=end();r++)
		{	
			if(r->name==name)
				return *r;
		}
		return End;
	}

	Reweighters & init(params &p)
	{
		rew.init(p);
		ff_configure(p);
        return *this;
	}


};


extern Reweighters REW; //global reweighters object 



#endif