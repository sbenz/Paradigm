#include <getopt.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <dai/alldai.h>

#include "common.h"
#include "configuration.h"
#include "evidencesource.h"

using namespace std;
using namespace dai;

#define GENE_COL 0
#define MOLECULE_COL 1
#define EVIDENCE_COL 2

#define VAR_DIM 3

void print_usage(int signal)
{
  cerr << "hgFactorGraph" << endl
       << "Usage:" << endl
       <<"  hgFactorGraph [options] -p pathway.tab -c config.txt -b batchPrefix [-e emOutput]" << endl
       << "C++ program for taking bioInt data and performing inference using libdai" << endl
       << "Note this can't be linked in to kent src as libDAI is GPL" << endl
       << "Valid options:" << endl
       << "\t-v,--verbose : verbose mode" << endl;
  exit(signal);
}

void die(string msg) 
{
  cerr << msg << endl;
  exit(-1);
}

inline double log10odds(double post,double prior) 
{ 
  return log( ( post / (1.0 - post) ) 
	      / ( prior / (1.0 - prior) ) ) 
    / log(10); 
}

void outputFastaPerturbations(string sampleName, InfAlg* prior, InfAlg* sample,
			      FactorGraph& fg, map<long,string>& activeNodes,
			      ostream& out)
{
  out << "> " << sampleName;
  out << " loglikelihood=" << (sample->logZ() - prior->logZ()) 
       << endl;
  for (size_t i = 0; i < fg.nrVars(); ++i)
    {
      const Var& v = fg.var(i);
      if(activeNodes.count(v.label()) == 0)
	continue;
      out << activeNodes[v.label()];
      Factor priorBelief = prior->belief(v);
      Factor belief = sample->belief(v);
      vector<double> priors;
      vector<double> posteriors;
      bool beliefEqualOne = false;
      
      for (size_t j = 0; j < belief.states(); ++j)
	{
	  if(belief[j] == 1 || priorBelief[j] == 1)
	    {
	      beliefEqualOne = true;
	      break;
	    }
	  priors.push_back(priorBelief[j]);
	  posteriors.push_back(belief[j]);
	}
      out << "\t";
      if(beliefEqualOne)
	out << "NA";
      else
	{
	  double down = log10odds(posteriors[0],priors[0]); 
	  double nc = log10odds(posteriors[1],priors[1]); 
	  double up = log10odds(posteriors[2],priors[2]); 
	  
	  if (nc > down && nc > up) 
	    out << "0";
	  else if (down > up)
	    out << (-1.0*down);
	  else
	    out << up;				
	}
      out << endl;
    }
}

void outputEmInferredParams(ostream& out, EMAlg& em, PathwayTab& pathway) {
  out << "> em_iters=" << em.getCurrentIters() 
      << " logZ=" << em.getLogZ() << endl;
  for (EMAlg::s_iterator m = em.s_begin(); m != em.s_end(); ++m) {
    for (MaximizationStep::iterator pit = m->begin(); pit != m->end(); ++pit) {
      vector< Var > vars;
      vector< Real > vals;
      vector< size_t > dims;

      pit->collectParameters(em.eStep().fg(), vals, vars);
      for (size_t vi = 0; vi < vars.size(); ++vi) {
	dims.push_back(vars[vi].states());
	if (vi == 0) {
	  out << "> child='"<< pathway.getNode(vars[vi].label()).second << "'";
	} else {
	  out << " edge" << vi << "='" 
	      << pathway.getInteraction(vars[0].label(), vars[vi].label())
	      << '\'';
	}
      }
      out << endl;
      for(MultiFor s(dims); s.valid(); ++s) {
	for (size_t state = 0; state < dims.size(); ++state) {
	  out << s[state] << '\t';
	}
	out << vals[(size_t)s] << endl;
      }
    }
  }
}

int main(int argc, char *argv[])
{
  const char* const short_options = "hp:b:c:e:o:v";
  const struct option long_options[] = {
    { "batch", 0, NULL, 'b' },
    { "config", 0, NULL, 'c' },
    { "pathway", 0, NULL, 'p' },
    { "em", 0, NULL, 'e' },
    { "output", 0, NULL, 'o' },
    { "help", 0, NULL, 'h' },
	{ "verbose", 0, NULL, 'v' },
    { NULL, 0, NULL, 0 }
  };
  int next_options;
  
  string pathwayFilename;
  string batchPrefix;
  string configFile;
  string paramsOutputFile;
  string actOutFile;

  // /////////////////////////////////////////////////
  // Read in command line options
  //
  do {
    next_options=getopt_long(argc,argv,short_options,long_options,NULL);
    switch (next_options){
    case 'h': print_usage(EXIT_SUCCESS); break;
      //case 't': cout << "Found t=" << optarg << endl; break;
    case 'b': batchPrefix = optarg; break;
    case 'c': configFile = optarg; break;
    case 'p': pathwayFilename = optarg; break;
    case 'e': paramsOutputFile = optarg; break;
    case 'o': actOutFile = optarg; break;
	case 'v': VERBOSE = true; break;
    }
  } while (next_options != -1);

  ostream* outstream;
  ofstream outFileStream;
  if (actOutFile == "") {
    outstream = &cout;
  } else {
    outFileStream.open(actOutFile.c_str());
    if (!outFileStream.is_open()) {
      die("couldn't open output file");
    }
    outstream = &outFileStream;
  }

  // /////////////////////////////////////////////////
  // Verify that command line options are valid
  //
  if(pathwayFilename.length() == 0)
    {
      cerr << "Missing required arguments" << endl;
      print_usage(EXIT_FAILURE);
    }
  if(batchPrefix.length() == 0)
    {
      cerr << "Missing batch prefix" << endl;
      print_usage(EXIT_FAILURE);
    }
  if(configFile.length() == 0)
    {
      cerr << "Missing configuration file" << endl;
      print_usage(EXIT_FAILURE);
    }

  // /////////////////////////////////////////////////
  // Load configuration
  RunConfiguration conf(configFile);
  if (conf.evidenceSize() == 0) {
    die("Must have at least one evidence file in configuration.");
  }

  // /////////////////////////////////////////////////
  // Load pathway
  ifstream pathwayStream;
  pathwayStream.open(pathwayFilename.c_str());
  if (!pathwayStream.is_open()) {
    die("Could not open pathway stream");
  }
  PathwayTab pathway = PathwayTab::create(pathwayStream);

  // /////////////////////////////////////////////////
  // Read in evidence
  Evidence evidence;

  vector<EvidenceSource> evid;
  for(size_t i = 0; i < conf.evidenceSize(); i++) {
    EvidenceSource e(conf.evidence(i), batchPrefix);
    if(VERBOSE)
      cerr << "Parsing evidence file: " << e.evidenceFile() << endl;
    e.loadFromFile(pathway, evidence);
    evid.push_back(e);
    if (i > 0 && e.sampleNames() != evid[0].sampleNames()) 
      {
	die("Sample names differ in files " + e.evidenceFile() + " and " 
	    + evid[0].evidenceFile());
      }
  }
  if(VERBOSE)
    cerr << "Added evidence for " << evid[0].sampleNames().size() 
	 << " samples" << endl;

  // /////////////////////////////////////////////////
  // Construct the factor graph
  vector< Factor > factors;
  vector< MaximizationStep > msteps;
  pathway.constructFactors(conf.emSteps(), factors, msteps);
  
  FactorGraph priorFG(factors);
  PropertySet inferenceOptions = conf.getInferenceProperties(pathwayFilename);
  std::string method = inferenceOptions.GetAs<std::string>("method");
  map< long, string > outNodes = pathway.getOutputNodeMap();
  
  InfAlg* prior = newInfAlg(method, priorFG, inferenceOptions);
  prior->init();

  // /////////////////////////////////////////////////
  // Run EM
  const PropertySet& em_conf = conf.emProps();
  EMAlg em(evidence, *prior, msteps, &em_conf);
  while(!em.hasSatisfiedTermConditions()) {
    em.iterate();
    if (VERBOSE) {
      outputEmInferredParams(cerr, em, pathway);
    }
  }
  em.run();

  ofstream paramsOutputStream;
  paramsOutputStream.open(paramsOutputFile.c_str());
  if (paramsOutputStream.is_open()) {
    outputEmInferredParams(paramsOutputStream, em, pathway);
  }
  prior->run();
  
  // /////////////////////////////////////////////////
  // Run inferenec on each of the samples
  Evidence::iterator sample_iter = evidence.begin();
  for ( ; sample_iter != evidence.end(); ++sample_iter) {
    InfAlg* clamped = prior->clone();
    
    sample_iter->second.applyEvidence(*clamped);
    clamped->run();
    
    outputFastaPerturbations(sample_iter->first, prior, clamped, priorFG, 
			     outNodes, *outstream);

    delete clamped;
  }

  delete prior;

  return 0;
}
