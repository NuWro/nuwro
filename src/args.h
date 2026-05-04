#ifndef _args_h_
#define _args_h_
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <sstream>

using namespace std;

struct FileInfo {
  string _filename;
  string _path;

  // Returns filename and path
  void ExtractFilename(const std::string &fullPath) {
    filesystem::path p(fullPath);
    this->_filename = p.filename().string();

    this->_path = p.parent_path().string();

    if (!this->_path.empty())
      this->_path += "/";
  }
};

class args {
public:
  const char *program;
  const char *input;
  const char *output;
  FileInfo
      outputFileInfo; // It was found that if one wants to save test events in
                      // NuWro then one has to forcefully place it in main/
                      // directory (project directory). As the output
                      // filename would be "weighted.<name_of_the_file>.root",
                      // passing the relative/absolute path of the output nuwro
                      // file along with its name would result in production of
                      // the file in the main/ directory. This structure was
                      // created to store the name of the file and its path (if
                      // provided) in separate string such that
                      // "weighted.<name_of_the_file>.root" can be safely
                      // produced at the provided path of the file
                      // along with its name
  const char *progress;
  stringstream params;

  args(const char *p = "nuwro", const char *i = "params.txt",
       const char *o = "eventsout.root")
      : program(p), input(i), output(o) {
    outputFileInfo._filename = string(output);
  }

  int read(int argc, char **argv) {
    if (argc == 2 and
        (string(argv[1]) == "-v" or string(argv[1]) == "--version")) {
      //  cout << VERSION << endl;
      exit(0);
    }
    if (argc % 2 == 0)
      usage();
    for (int i = 1; i < argc - 1; i++) {
      if (string(argv[i]) == "-i" and argv[i + 1][0] != '-')
        input = argv[++i];
      else if (string(argv[i]) == "-o" and argv[i + 1][0] != '-')
        output = argv[++i];
      else if (string(argv[i]) == "-progress" and argv[i + 1][0] != '-')
        progress = argv[++i];
      else if (string(argv[i]) == "-p" and argv[i + 1][0] != '-')
        params << argv[++i] << endl;
      else
        usage();
    }
    return 0;
  }

  void usage() {
    cout << "usage: " << program
         << " [-i input_paramaters_file] [-o output_root_file] [-p "
            "\"param1=value1\"] [-p \"param2=value2\"]...\n";
    exit(22);
  }
};
#endif
