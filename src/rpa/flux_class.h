#ifndef _flux_class_h_
#define _flux_class_h_

#include <iostream>
#include <fstream>
#include <vector>

namespace rpa
{

using namespace std;

struct punkt{

	double data[2];

};


class data_set{

public:

	int N_bin;
	double step;

	vector <punkt> set;

	double calka;   // suma po binach

	data_set(char* nazwa_pliku, int K);

	double weight(double E);

	double N(double EMIN, double EMAX);
};


//vector <punkt> ANL;
//vector <punkt> BNL;

data_set::data_set(char* nazwa_pliku, int K)
{

	N_bin =K;
	ifstream in(nazwa_pliku);

	punkt prob;

	if(!in.good()) cout<< "brak pliku!!!"<<endl;

	double E, w;

	while(true){
		if(!in.good()) break;

		in>>E>>w;

		prob.data[0] =  E;
		prob.data[1] =  w;

		set.push_back(prob);

	}

	in.close();

	///      --------        ///
	///      liczymy calke   ///

	step =   (set[set.size()-1].data[0] - set[0].data[0])/N_bin;


	double Estart = set[0].data[0];

	double sum =0;
	for(int i = 0; i < N_bin; i ++)
		sum+=weight(Estart + i*step);

	calka = sum*step;

}


double data_set::weight(double E)
{

	double E1;
	double E2;


	if(E > set[set.size()-1].data[0] )
		return 0;


	if(E <= set[0].data[0] )
		return 0;


	int i =0;
	E2 = set[0].data[0];

	//cout<<E1<<endl;

	while( E > E2 )
	{
		E2 = set[i].data[0];


		i++;
	}
	E1 = set[i-2].data[0];


	double a  = (set[i-1].data[1] - set[i-2].data[1])/(E2 - E1);
	double b  =  set[i-1].data[1] - set[i-1].data[0]*a;

	return a*E + b;

}


double data_set::N(double EMIN, double EMAX)
{
	int k=200;
	double sum = 0;

	double step = (EMAX - EMIN)/k;

	for(int i = 0; i < k; i++ )
	sum += weight(EMIN + i*step );

	return sum*step;
}

}
/*

int main(){


data_set ANL("ANL_flux.dat", 50 );

data_set BNL("BNL_flux.dat", 50 );


for(double E=0; E<6; E+= 6./ANL.N_bin)
cout << E << "\t" << ANL.weight(E)<<endl;


//cout << BNL.calka<<endl;
cout << ANL.calka<<endl;


return 0;

}
*/
#endif
