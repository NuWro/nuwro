/////////////////////////////////////////////////////////////////////////////////////
//     Tworzenie histogramów i zapisywanie ich w plikach dla xmgrace
/////////////////////////////////////////////////////////////////////////////////////

#ifndef _hist_h_
#define _hist_h_
#include <iostream>
#include <fstream>

/////////////////////////////////////////////////////////////////////////////////////
///   KLASA HIST - ró¿niczkowe przekroje czynne i ich wykresy
/////////////////////////////////////////////////////////////////////////////////////
class hist
{
private:
  char * name;                  // nazwa parametru
  double min, max, width;       // minimum , maximim i szeroko¶æ kube³ka
  int nkub;                     // ilo¶æ kube³ków 
  double *sum;			// tablica kube³ków
  int count;			// ilo¶æ wstawieñ
  double _total;                // suma wstawieñ
  inline int ind (double arg);   // nr kube³ka je¶li parametr wynosi 'arg'
  inline double srodek (int i, double pol=0.5);  // ¶rodek i-tego kube³ka

public:

   inline hist (char *n = NULL, double minv = 0, double maxv = 1,  int ile = 1);
    
   inline ~hist () { delete [] sum;}
   
   inline  void insert_value (double val, double weight);
   // Wstawienie warto¶ci do histogramu 
   
   inline  void wykres (ostream & out, 
                        double jednostkax, 
			double jednostkay, 
			double pol=0.5, 
			char separator='\t');
   // zrzuca do strumienia out dane potrzebne do zrobienia 
   // wykresu zale¿no¶ci warto¶ci przekroju czynnego od parametru

   
   inline  void wykres (const char* filename, 
                        double jednostkax, 
			double jednostkay, 
			double pol=0.5);
   // filename=nazwa pliku który zostanie utworzony.
   // parametr=1, 2 lub 3 numer parametru na osi x   

   
   inline  double total (){    return _total;  }   
   //ca³kowity przekrój czynny ? NIE 
   // Suma wrzuconych liczb

};

/////////////////////////////////////////////////////////////////////////////////////
///                    I M P L E M E N T A C J A                                  ///
/////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////
double hist::srodek (int i, double pol)   // ¶rodek itego kube³ka
/////////////////////////////////////////////////////////////////////////////////////
  {
    return (i + pol) * width + min;
  }

/////////////////////////////////////////////////////////////////////////////////////
int hist::ind (double arg)   // nr kube³ka je¶li parametr wynosi 'arg'
/////////////////////////////////////////////////////////////////////////////////////
  {
    if (arg < min || arg>max)
      return -1;
    int i = int ((arg - min) / width);
    if (i >= nkub)
      return -1;
    return i;
  }



/////////////////////////////////////////////////////////////////////////////////////
hist::hist (char *n, double minv, double maxv, int ile)
/////////////////////////////////////////////////////////////////////////////////////
           :name(n), 
	    min(minv), 
	    max(maxv), 
	    width((maxv-minv)/ile),
	    nkub(ile),
	    count(0),
	    _total(0)
  { 
    sum = new double[nkub];
    for(int i=0;i< nkub;i++)
         sum[i] = 0;
  }

/////////////////////////////////////////////////////////////////////////////////////
  void hist::insert_value (double val, double weight)
/////////////////////////////////////////////////////////////////////////////////////
  {
    int i = ind (val);
    count++;
    if (i >= 0)
      {sum[i] += weight;
       _total += weight;
      } 
  }

/////////////////////////////////////////////////////////////////////////////////////
  void hist::wykres (ostream & out, double jednostkax, double jednostkay, double pol, char separator)
/////////////////////////////////////////////////////////////////////////////////////
{
    out << "#  Zale¿no¶æ od " << name << endl;
    
    for (int i = 0; i < nkub; i++)
      {
	out << srodek (i, pol) / jednostkax << separator;
	if (count > 0)
	  out << (sum[i] / count) /width / jednostkay << endl;
	else
	  out << 0 << endl;
      }
  }
  
/////////////////////////////////////////////////////////////////////////////////////
  void hist::wykres (const char* filename, double jednostkax, double jednostkay, double pol)
/////////////////////////////////////////////////////////////////////////////////////
  {
    std::ofstream wyk (filename);
    wykres (wyk, jednostkax, jednostkay, pol);
  }

/////////////////////////////////////////////////////////////////////////////////////
///        K O N I E C   I M P L E M E N T A C J I                                ///
/////////////////////////////////////////////////////////////////////////////////////
#endif
