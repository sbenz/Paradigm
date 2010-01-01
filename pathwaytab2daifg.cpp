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
  RunConfiguration c(argc == 3 ? argv[2] : "/dev/null");
  PathwayTab path = PathwayTab::create(path_stream, c.pathwayProps());
  
  vector< Factor > factors;
  vector< MaximizationStep > msteps;
  
  path.constructFactors(c.emSteps(), factors, msteps);

  FactorGraph fg(factors);

  path.printNodeMap();
  cout << fg;
  return 0;
}
