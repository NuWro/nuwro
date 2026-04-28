#include <iostream>
#include <string>
#include <regex>


static const int BOX_WIDTH = 76;

inline string pad(int n, char c = ' ') { return string(n > 0 ? n : 0, c); }

inline int visible_len(const string& s)
{
  int len = 0;
  for (size_t i = 0; i < s.size(); )
  {
    unsigned char c = (unsigned char)s[i];
    if      (c < 0x80) { i += 1; }
    else if (c < 0xE0) { i += 2; }
    else if (c < 0xF0) { i += 3; }
    else               { i += 4; }
    len++;
  }
  return len;
}

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

inline string parse_nuwro_version(const string& tag)
{
  static const std::regex re(R"(nuwro_(\d+)\.(\d+).*)");
  std::smatch m;
  // Use regex_match so the string must start with the expected prefix
  if (!std::regex_match(tag, m, re))
    return "";
  return m[1].str() + "." + m[2].str();
}

static const string DIGIT_GLYPHS[11][3] = {
  // row0     row1      row2
  { " __ ",  "/  \\",  "\\__/" }, // 0
  { "  ",    "/|",     " |"   },  // 1
  { "__ ",   " _)",    "/__"  },  // 2
  { "__ ",   " _)",    "__)"  },  // 3
  { "    ",  "|__|",   "   |" },  // 4
  { " __",   "|_ ",    "__)"  },  // 5
  { " __ ",  "/__ ",   "\\__)" }, // 6
  { "___",   "  /",    " / "  },  // 7
  { " __ ",  "(__)",   "(__)" },  // 8
  { " __ ",  "(__\\",  " __/" },  // 9
  { " ",     " ",      "."   },   // . (index 10)
};

inline int digit_index(char c)
{
  if (c >= '0' && c <= '9') return c - '0';
  if (c == '.')             return 10;
  return -1; // unknown
}

inline void render_version_art(const string& version, string out[3])
{
  static const string GAP = " "; // single space between glyphs
  out[0] = out[1] = out[2] = "";
  bool first = true;
  for (char c : version)
  {
    int idx = digit_index(c);
    if (idx < 0) continue; // skip non-renderable chars (e.g. leading 'v')
    if (!first)
      for (int r = 0; r < 3; r++) out[r] += GAP;
      for (int r = 0; r < 3; r++) out[r] += DIGIT_GLYPHS[idx][r];
        first = false;
    }
}

static const int  VERSION_LEFT_PAD  = 13;
static const int  VERSION_WIDTHS[3] = { 25, 26, 27 };
static const string LOGO_SUFFIX[3]  = {
  ".yyyoo    :MMN:.  :MMho.    `---h:    ",
  ":yyoh-    sMMh-`.mMho-     ---h:     ",
  "oyy+h    `mMM+-mMho-     ---h:      ",
};

inline void print_nuwro_banner(const string& raw_tag)
{
  const string version = parse_nuwro_version(raw_tag);

  string vart[3];
  render_version_art(version, vart);

  box_top();
  box_blank();
  box_line(R"(`.``      `.       .-.        `      )", 38);
  box_line(R"(`-/+os+s  ./ohmN:.    sNNNy:`   .---.+`  )", 34);
  box_line(R"(|\ |     |  |  _  _           :oooyysy: +oodMMd-`  .MMMMM/-   -----s/   )", 4);
  box_line(R"(| \| |_| |/\| |  (_)           `.`oyy+d`   `mMMo-   yMMMyo`   .----d.   )", 4);
  for (int r = 0; r < 3; r++)
  {
    string text = vart[r];
    while ((int)text.size() < VERSION_WIDTHS[r]) text += ' ';
    text += LOGO_SUFFIX[r];
    box_line(text, VERSION_LEFT_PAD);
  }
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