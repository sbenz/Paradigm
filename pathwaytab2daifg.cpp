#include <fstream>

#include "pathwaytab.h"

void usage(int exit_code) {
  cout << "Usage: " << endl
       << "  pathwaytab2daifg pathway_file [config_file]"<< endl;
  exit(exit_code);
}

int main(int argc, char** argv) {
  if (argc != 2 && argc != 3) {
    usage(2);
  }

  ifstream path_stream(argv[1]);
  PropertySet p;
  if (argc == 3) {
    RunConfiguration c(argv[2]);
    p = c.pathwayProps();
  }
  PathwayTab path = PathwayTab::create(path_stream, p);

  path.printNodeMap();
  path.printDaiFactorSection();
  return 0;
}
