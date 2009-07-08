#include <fstream>

#include "pathwaytab.h"

void usage(int exit_code) {
  cout << "Usage: " << endl
       << "  pathwaytab2daifg pathway_file"<< endl;
  exit(exit_code);
}

int main(int argc, char** argv) {
  if (argc != 2) {
    usage(2);
  }

  ifstream path_stream(argv[1]);
  PathwayTab path = PathwayTab::create(path_stream);

  path.printNodeMap();
  path.printDaiFactorSection();
  return 0;
}
