#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

double interpolate(double x,int n,double X[],double Y[])
{
    /*
     * Given n points (X[i],Y[i]) and real x
     * find y such that (x,y) linearly interpolates one of the 
     * sections or y=Y[0] if x<=X[0] or y=Y[n-1] if x>=X[n-1]
     */
     
     
    if(x<=X[0]) 
        return Y[0];
    if(x>=X[n-1])
        return Y[n-1];
    int i=0;
    while(x>X[i+1]) 
        i++;
    double a=x-X[i];
    double b=X[i+1]-x;
    
    return (b*Y[i]+a*Y[i+1])/(a+b);

}

int main()
{
    /*
     * Convert any (x,y) data file to one 
     * with equally spaced x coordinates.
     */
     
    vector<double> X,Y;
    double key,val;
    while(cin>>key>>val)
    {
        X.push_back(key);
        Y.push_back(val);
    }
    int n=Y.size();
    int N=2*n-1;
    double x0=X[0];
    double x1=X[n-1];
    double dx=(x1-x0)/(N-1);
    for(int i=0;i<N;i++)
    {
        double x=x0+dx*i;
        double y=interpolate(x, n, &X[0], &Y[0]);
        cout<<x<<"\t"<<y<<endl;
    }
}
