#include <iostream>
#include <fstream>
#include <vector>

using namespace std;
int main()
{
    /*
     * convert one line  nuwro energy profile 
     * in form E_min E_max bin1 bin2 ... binN
     * to file of x y values (one pair per line) 
     * such that: 
     * x-coordinates are centers of bins
     * y-coordinates are heights of bins
     */
    vector<double> x,y;
    double x0,x1,val;
    cin>>x0>>x1;
    while(cin>>val)
	y.push_back(val);
    int n=y.size();
    double dx=(x1-x0)/n;
    double sx=1,sy=1;
    for(int i=0;i<n;i++)
	cout<< (x0+dx*(0.5+i))*sx<<"\t"<<y[i]*sy<<endl;
}
