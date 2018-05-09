/////////////////////////////////////////////////////////////////////////////////////
//     Tworzenie histogramów i zapisywanie ich w plikach dla xmgrace
/////////////////////////////////////////////////////////////////////////////////////

#ifndef _hist_h_
#define _hist_h_
#include <iostream>
#include <fstream>

/////////////////////////////////////////////////////////////////////////////////////
/// hist class provides plots of differential cross sections
/////////////////////////////////////////////////////////////////////////////////////
class hist
{
private:
  string name;                  // parameter name
  double min, max, width;       // x_min , x_max  and bin width
  int n_bins;                   // number of bins 
  double *sum;			// bins array
  int count;			// number of insertions
  double _total;                // total of insertions
  inline double middle(int i, double half=0.5); // middle x for i-th bin

public:

   inline hist (string n="", double minv = 0, double maxv = 1,  int ile = 1);
    
   inline ~hist () { delete [] sum;}
   
   inline  void insert_value (double val, double weight);
   
   inline  void plot (ostream & out, double xunit, double yunit, double half=0.5, char separator='\t');
   
   inline  void plot (string filename, double xunit, double yunit, double half=0.5, char separator='\t');
   
   inline  double total (){    return _total;  }   

};


//---------------------------------------------------------------------------------------
double hist::middle (int i, double pol)   // middle of i-th bin
{
    return (i + pol) * width + min;
}


//---------------------------------------------------------------------------------------
hist::hist (string n, double minv, double maxv, int bins)
           :name(n), 
	    min(minv), 
	    max(maxv), 
	    width((maxv-minv)/bins),
	    n_bins(bins),
	    count(0),
	    _total(0)
{ 
    sum = new double[n_bins];
    for(int i=0;i< n_bins;i++)
         sum[i] = 0;
}

//---------------------------------------------------------------------------------------
void hist::insert_value (double val, double weight)
{
    int i = int ((val - min) / width);
    count++;
    if (i >= 0 && i<n_bins && weight!=0)
      {sum[i] += weight;
       _total += weight;
      } 
}
//---------------------------------------------------------------------------------------
void hist::plot (ostream & out, double xunit, double yunit, double pol, char separator)
{
    out << "#  " << name << " dependance "<<endl;
    out << "#  mean = " << _total/count/(yunit*xunit)<<endl;
    
    for (int i = 0; i < n_bins; i++)
      {
		out << middle (i, pol) / xunit << separator;
		if (count > 0)
		  out << (sum[i] / count) /width / yunit << endl;
		else
		  out << 0 << endl;
      }
}
  
//---------------------------------------------------------------------------------------
void hist::plot (string filename, double xunit, double yunit, double pol, char separator)
{
    std::ofstream wyk (filename.c_str());
    plot (wyk, xunit, yunit, pol, separator);
}

#endif
