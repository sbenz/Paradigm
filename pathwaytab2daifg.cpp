/********************************************************************************/
/* Copyright 2009-2011 -- The Regents of the University of California           */
/* This code is provided for research purposes to scientists at non-profit		*/
/*  organizations.  All other use is strictly prohibited.  For further			*/
/*  details please contact University of California, Santa Cruz or				*/
/*	Five3 Genomics, LLC (http://five3genomics.com).								*/
/********************************************************************************/

#include <fstream>

#include "pathwaytab.h"

void usage(int exit_code) {
  cout << "pathwaytab2daifg"
#ifdef VERSION
       << " -- " << VERSION
#endif
       << endl
       << "Usage: " << endl
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
