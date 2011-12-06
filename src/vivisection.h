#ifndef _vivisection_h
#define _vivisection_h

#include <string>

using namespace std;

const int nofb = 10;
const int nofa = 40;
const int nofp = 15;
const int events = 100000;

void multi(string filename, double vivi[][nofa][nofp], double norm[][nofa], double counter[][5]);
void maketex (string file, double vivi[][nofa][nofp], double norm[][nofa], double counter[][5], string name, int fz);
int viviNomad (int fz);
int viviMultiplicity (int fz);

#endif
