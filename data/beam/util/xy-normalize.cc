#include <iostream>
#include <fstream>
#include <vector>

using namespace std;


int main()
{
    /* 
     * normalize the plot read from stdin 
     * such that maximum be 1
     * output to stdout
     */
     
    vector<double> X,Y;
    double key,val;
    double max_val=0;
    while(cin>>key>>val)
    {
        X.push_back(key);
        Y.push_back(val);
	if(val>max_val)
	    max_val=val;
    }
    int n=Y.size();
    for(int i=0;i<n;i++)
    {
        cout<<X[i]<<"\t"<<Y[i]/max_val<<endl;
    }
}
