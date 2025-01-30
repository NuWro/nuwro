#include <iostream>


inline void frame_top(string text)
{
  string line (78,' ');
  string line2(76,'_');
  string line3(text.length(),' ');
  string line4(71-text.length(),' ');
  string line5(text.length(),'_');

  cout << "  "   << line2 << endl;
  cout << " |  " << line3 << "  |" << line4 << "|" << endl;
  cout << " |  " << text  << "  |" << endl;
  cout << " |__" << line5 << "__|" << endl;
  cout << " |"   << endl  << " |"  << endl;
}

inline void frame_bottom()
{
  string line (78,' ');
  string line2(76,'_');

  cout << line << "|"   << endl
       << line << "|"   << endl
       << " |" << line2 << "|" << endl << endl;
}