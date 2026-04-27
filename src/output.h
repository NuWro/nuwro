#include <iostream>
#include <string>


static const int BOX_WIDTH = 76;

inline string pad(int n, char c = ' ') { return string(n > 0 ? n : 0, c); }

inline void box_blank()
{
    cout << " |" << pad(BOX_WIDTH) << "|" << endl;
}

inline void box_line(const string& text, int left_pad = 3)
{
    int right_pad = BOX_WIDTH - left_pad - (int)text.length();
    cout << " |" << pad(left_pad) << text << pad(right_pad) << "|" << endl;
}

inline void box_bottom()
{
    cout << " |" << pad(BOX_WIDTH, '_') << "|" << endl << endl;
}

inline void box_top()
{
    cout << "  " << pad(BOX_WIDTH, '_') << endl;
    cout << " |" << pad(BOX_WIDTH)      << "|" << endl;
}


inline void print_nuwro_banner()
{
    box_top();
    box_blank();
    box_line(R"(`.``      `.       .-.        `      )", 38);
    box_line(R"(`-/+os+s  ./ohmN:.    sNNNy:`   .---.+`  )", 34);
    box_line(R"(|\ |     |  |  _  _           :oooyysy: +oodMMd-`  .MMMMM/-   -----s/   )", 4);
    box_line(R"(| \| |_| |/\| |  (_)           `.`oyy+d`   `mMMo-   yMMMyo`   .----d.   )", 4);
    box_line(R"(__   __                  .yyyoo    :MMN:.  :MMho.    `---h:    )", 13);
    box_line(R"(_) |__    /|  /|         :yyoh-    sMMh-`.mMho-     ---h:     )", 14);
    box_line(R"(/__  __| . _|_ _|_         oyy+h    `mMM+-mMho-     ---h:      )", 13);
    box_line(R"(.yyys+    -MMNmMh+-     ---h:       )", 40);
    box_line(R"(:yyod.    sMMMh+-     ---h:        )", 41);
    box_line(R"(Wrocław Neutrino Event Generator       oyy+h    `mMh-+.    ---h:         )", 3);
    box_line(R"(https://github.com/NuWro/nuwro         .yyss/  .s/y--oh   .--y-          )", 3);
    box_line(R"(:yy+d`.sy+----d/ .--y-           )", 43);
    box_line(R"(J. T. Sobczyk et al.                     oyy++syoy`.--:s.--y-            )", 3);
    box_line(R"(Institute of Theoretical Physics         .yyssyoy.  ------y-             )", 3);
    box_line(R"(University of Wrocław                     :yyyoh.   `----y-              )", 3);
    box_line(R"(Poland                                     osoh.     .-:y-               )", 3);
    box_line(R"(`-:.       .:-                )", 46);
    box_blank();
    box_bottom();
}

inline void print_cascade_mode_info()
{
    box_top();
    box_line(R"(__                                                            )", 14);
    box_line(R"(/    _   _  _  _   _|  _   |\/|  _   _|  _                     )", 13);
    box_line(R"(\__ (_| _) (_ (_| (_| (-   |  | (_) (_| (-                     )", 13);
    box_blank();
    box_blank();
    box_line("Hadrons are introduced directly to the cascade model. Only selected");
    box_line("parameters are active. The incident particle starting point follows");
    box_line("the beam_placement parameter as");
    box_line("(0) nucleus center");
    box_line("(1) random nucleon's position:           transparency mode");
    box_line("(2) just under the surface of the nucleus: scattering mode");
    box_blank();
    box_bottom();
}


inline void frame_top(string text)
{
    string line2(BOX_WIDTH, '_');
    string line3(text.length(), ' ');
    string line4(71 - text.length(), ' ');
    string line5(text.length(), '_');
    cout << "  "   << line2 << endl;
    cout << " |  " << line3 << "  |" << line4 << "|" << endl;
    cout << " |  " << text  << "  |" << endl;
    cout << " |__" << line5 << "__|" << endl;
    cout << " |"   << endl  << " |"  << endl;
}

inline void frame_bottom()
{
    string line (78, ' ');
    string line2(BOX_WIDTH, '_');
    cout << line << "|"   << endl
         << line << "|"   << endl
         << " |" << line2 << "|" << endl << endl;
}