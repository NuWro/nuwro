#ifndef _elements_h_
#define _elements_h_

struct element
{
	int Z;
	const char* symbol;
	const char* name;
	int group;
	int period;
	double weight;
	double density;
	double melt;
	double boil;
	double heat;
	double neg;
	double abundance;
};

extern element el[];

#endif
