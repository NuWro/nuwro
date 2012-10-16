#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <fstream>
using namespace std;

double mm=105.658;

int events=50000;				 // for a moment
const int kinbins=18;
const int katbins=20;
double kinmin=200;
double kinmax=2000;
double katmin=-1;
double katmax=1;

								 //0 data 1 errors 2 mec; 3-21 varying axial mass
double FINAL_RESULTS [22] [katbins] [kinbins];

								 //read cross section from the file .root.txt
void crosssections (double cross[], ifstream& Input)
{
	string y;
	double x;
	if (Input)
	{
		do
		{
			getline (Input, y);
			for (int j=0; j<8; j++)
			{
				for (int k=0; k<3; k++)
					{Input>>y;}
					Input>>cross[j];
			}
		} while (Input);
	}
}


void to_table_3D (double table[][katbins][kinbins], ifstream& Input, int model, int katbins, int kinbins, double weight)
{
	double x;
	if (Input)
	{
		do
		{
			for (int s=0; s<katbins; s++)
			{
				for (int ss=0; ss<kinbins; ss++)
				{
					Input>>x;
					cout<<s<<"  "<<ss<<"  "<<x<<endl;
					table[model][s][ss]=weight*x;
					//cout<<x<<"  "<<table [model] [s] [ss]<<endl;
				}
			}
		}
		while (Input);
	}
	Input.close();
}


void to_table_3D_bis (double table[][katbins][kinbins], ifstream& Input, int model, int katbins, int kinbins, double weight)
{
	double x;
	while (!Input.eof())
	{

		{
			for (int s=0; s<katbins; s++)
			{
				for (int ss=0; ss<kinbins; ss++)
				{
					Input>>x;
					cout<<"tutaj  "<<x<<"  "<<endl;
					table[model][s][ss]=weight*x;
					//cout<<x<<"  "<<table [model] [s] [ss]<<endl;
				}
			}
		}
	}
}


ifstream input_data ("MB_data1.dat");
//ifstream input_error ("MB_data2.dat");
//ifstream input_mec ("mec.dat");

char first[130];
char second[130];

int main ()
{
	to_table_3D_bis (FINAL_RESULTS, input_data, 0, katbins, kinbins, 1e-44);
	//to_table_3D (FINAL_RESULTS, input_error, 1, katbins, kinbins, 1e-45);
	//to_table_3D (FINAL_RESULTS, input_mec, 2, katbins, kinbins, 1e-3);

	for (int s=0; s<20; s++)
	{
		for (int ss=0; ss<18; ss++)
		{
			cout<<FINAL_RESULTS [0][s][ss]<<"  ";
		}
		cout<<endl;
	}

}
