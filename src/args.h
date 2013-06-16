#ifndef _args_h_
#define _args_h_
#include <iostream>
#include <sstream>
#include <cstdlib>

using namespace std;


class args
{
public:
const char * program;
const char * input;
const char * output;
const char * progress;
stringstream params;

args(const char* p="nuwro",const char* i="params.txt",const char* o="eventsout.root"):
  program(p),input(i),output(o)
{}

int read(int argc, char** argv)
{ 
  if(argc%2==0) 
     usage();
  for(int i=1;i<argc-1;i++)
     {
	if(string(argv[i])=="-i" and argv[i+1][0]!='-')
           input=argv[++i];
        else 
	if(string(argv[i])=="-o" and argv[i+1][0]!='-')
           output=argv[++i];
        else 
	if(string(argv[i])=="-progress" and argv[i+1][0]!='-')
           progress=argv[++i];
	else
	if(string(argv[i])=="-p" and argv[i+1][0]!='-')
           params<<argv[++i]<<endl;
        else 
	   usage();
     } 
}

void usage()
{
  cout <<"usage: "<< program <<" [-i input_paramaters_file] [-o output_root_file] [-p \"param1=value1\"] [-p \"param2=value2\"]...\n";
  exit(22);  
}

};
#endif
