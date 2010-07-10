#include <getopt.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <dai/alldai.h>

#include <sys/time.h>
#include <sys/resource.h>

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
  cerr << "paradigm" 
#ifdef VERSION
       << " -- " << VERSION 
#endif
       << endl
       << "Usage:" << endl
       <<"  paradigm [options] -p path.tab -c cfg.txt -b prefix" << endl
       << "C++ program for taking bioInt data and performing inference using libdai" << endl
       << "Note this can't be linked in to kent src as libDAI is GPL" << endl
       << "Valid options:" << endl
       << "\t-e emOutputFile" << endl
       << "\t-m max_mem      : The maximum number of GB that can be allocated" << endl
       << "\t-v,--verbose    : verbose mode" << endl;
  exit(signal);
}

void die(string msg) 
{
  cerr << msg << endl;
  exit(-1);
}

inline double log10odds(double post,double prior) 
{ 
  return std::log( ( post / (1.0 - post) ) 
		   / ( prior / (1.0 - prior) ) ) 
    / log(10); 
}

// returns a map of all the nodes in a single subnet, for testing against
map< long, bool > nodesInSingleNet(vector< Factor > factors)
{
  // start with a single (random) node, and follow it out, marking down all the nodes we see
  map< long, bool > nodesSeen;
  vector< long > nodesToCheck;
  vector< Factor >::iterator factorIter = factors.begin();
  const VarSet tmpVars = factorIter->vars();
  vector< Var >::const_iterator tmpVarsIter;
  for (tmpVarsIter = tmpVars.begin(); tmpVarsIter != tmpVars.end(); ++tmpVarsIter) {
	const long varLabel = tmpVarsIter->label();
	nodesToCheck.push_back(varLabel);
	nodesSeen[varLabel] = true;
  }
  while(nodesToCheck.size() > 0)
  {
	long currNode = nodesToCheck.back();
	nodesToCheck.pop_back();
	for ( factorIter = factors.begin(); factorIter != factors.end(); ++factorIter) {
	  const VarSet tmpVars = factorIter->vars();
	  //vector< Var >::const_iterator tmpVarsIter = tmpVars.begin();
	  for (tmpVarsIter = tmpVars.begin(); tmpVarsIter != tmpVars.end(); ++tmpVarsIter) {
		const long varLabel = tmpVarsIter->label();
		if(varLabel == currNode)
		{
		  vector< Var >::const_iterator tmpVarsIter2 = tmpVars.begin();
		  for ( ; tmpVarsIter2 != tmpVars.end(); ++tmpVarsIter2) {
			const long varLabel2 = tmpVarsIter2->label();
			if(!nodesSeen[varLabel2])
			{
			  nodesToCheck.push_back(varLabel2);
			  nodesSeen[varLabel2] = true;  
			}
		  }		  
		}
	  }	
	}
	
  }
  return nodesSeen;
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
      
      for (size_t j = 0; j < belief.nrStates(); ++j)
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

void outputEmInferredParams(ostream& out, EMAlg& em, PathwayTab& pathway,
			    const vector< vector < SharedParameters::FactorOrientations > > &var_orders) {
  out << "> em_iters=" << em.Iterations() 
      << " logZ=" << em.logZ() << endl;
  size_t i = 0;
  for (EMAlg::s_iterator m = em.s_begin(); m != em.s_end(); ++m, ++i) {
    size_t j = 0;
    for (MaximizationStep::iterator pit = m->begin(); pit != m->end(); 
	 ++pit, ++j) {
      SharedParameters::FactorOrientations::const_iterator fo 
	= var_orders[i][j].begin();
      const vector< Var >& vars  = fo->second;
      Permute perm(fo->second);
      const Factor f = em.eStep().fg().factor(fo->first);
      vector< size_t > dims;

      // Output column headers
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
      // Output actual parameters
      out << endl;
      for(multifor s(dims); s.valid(); ++s) {
	for (size_t state = 0; state < dims.size(); ++state) {
	  out << s[state] << '\t';
	}
	out << f[perm.convertLinearIndex((size_t)s)] << endl;
      }
    }
  }
}

void setMaxMem(unsigned long maxmem) {
  struct rlimit memlimits;
  //  memlimits.rlim_cur = memlimits.rlim_max = 0x0C0000000; // 6 * 1024 * 1024 * 1024;
  // memlimits.rlim_cur = memlimits.rlim_max = 3 * 1024 * 1024 * 1024UL;
  memlimits.rlim_cur = memlimits.rlim_max = maxmem;
  printf("rlim_t size is: %ld\n", sizeof(memlimits.rlim_cur));
  printf("rlim_max size is: %ld\n", (memlimits.rlim_max));
  setrlimit(RLIMIT_AS, &memlimits);
}

void setMaxMemGigs(char * maxmem) {
  double md = strtod(maxmem, NULL);
  unsigned long mul = (unsigned long)(md * 1024 * 1024 * 1024);
  setMaxMem(mul);
}

int main(int argc, char *argv[])
{
  const char* const short_options = "hp:b:c:e:m:o:v";
  const struct option long_options[] = {
    { "batch", 0, NULL, 'b' },
    { "config", 0, NULL, 'c' },
    { "pathway", 0, NULL, 'p' },
    { "em", 0, NULL, 'e' },
    { "maxmem", 0, NULL, 'm' },
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
    case 'b': batchPrefix = optarg; break;
    case 'c': configFile = optarg; break;
    case 'p': pathwayFilename = optarg; break;
    case 'e': paramsOutputFile = optarg; break;
    case 'm': setMaxMemGigs(optarg); break;
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
  PathwayTab pathway = PathwayTab::create(pathwayStream, conf.pathwayProps());

  // /////////////////////////////////////////////////
  // Read in evidence
  vector<EvidenceSource> evid;
  map<string,size_t> sampleMap;
  vector<Evidence::Observation> sampleData;
  for(size_t i = 0; i < conf.evidenceSize(); i++) {
    EvidenceSource e(conf.evidence(i), batchPrefix);
    if(VERBOSE)
      cerr << "Parsing evidence file: " << e.evidenceFile() << endl;
    e.loadFromFile(pathway, sampleMap, sampleData);
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
  vector< vector < SharedParameters::FactorOrientations > > var_orders;
  var_orders = pathway.constructFactors(conf.emSteps(), factors, msteps);
  map< long, string > outNodes = pathway.getOutputNodeMap();
  
  // add in additional factors to link disconnected pieces of the pathway
  FactorGraph *testGraphConnected = new FactorGraph(factors);
  while(!testGraphConnected->isConnected())
  {
    map< long, bool > nodesInSubNet = nodesInSingleNet(factors);
    long backNode = nodesInSubNet.rbegin()->first;
    VarSet I_vars;
    I_vars |= Var(backNode, PathwayTab::VARIABLE_DIMENSION);
    bool addedLink = false;
    // now go find ones that aren't connected to this
    vector< Factor >::iterator factorIter;
    for ( factorIter = factors.begin(); !addedLink && factorIter != factors.end(); ++factorIter) {
      const VarSet tmpVars = factorIter->vars();
      vector< Var >::const_iterator tmpVarsIter;
      for (tmpVarsIter = tmpVars.begin(); !addedLink && tmpVarsIter != tmpVars.end(); ++tmpVarsIter) {
	const long varLabel = tmpVarsIter->label();
	if(!nodesInSubNet[varLabel])
	  {
	    I_vars |= Var(varLabel, PathwayTab::VARIABLE_DIMENSION);
	    
	    factors.push_back( Factor( I_vars, 1.0 ) );
	    addedLink = true;
	  }
      }	
    }
    break;
    delete testGraphConnected;
    testGraphConnected = new FactorGraph(factors);
  }
  delete testGraphConnected;
  
  FactorGraph priorFG(factors);

  PropertySet inferenceOptions = conf.getInferenceProperties(pathwayFilename);
  std::string method = inferenceOptions.getAs<std::string>("method");
  
  InfAlg* prior = newInfAlg(method, priorFG, inferenceOptions);
  prior->init();

  // /////////////////////////////////////////////////
  // Run EM
  const PropertySet& em_conf = conf.emProps();
  Evidence evidence(sampleData);
  EMAlg em(evidence, *prior, msteps, em_conf);
  while(!em.hasSatisfiedTermConditions()) {
    em.iterate();
    if (VERBOSE) {
      outputEmInferredParams(cerr, em, pathway, var_orders);
    }
  }
  em.run();

  ofstream paramsOutputStream;
  paramsOutputStream.open(paramsOutputFile.c_str());
  if (paramsOutputStream.is_open()) {
    outputEmInferredParams(paramsOutputStream, em, pathway, var_orders);
  }
  prior->run();
  
  // /////////////////////////////////////////////////
  // Run inference on each of the samples
  map<string, size_t>::iterator sample_iter = sampleMap.begin();
  for ( ; sample_iter != sampleMap.end(); ++sample_iter) {
    InfAlg* clamped = prior->clone();
    Evidence::Observation *e = &sampleData[sample_iter->second];
    for (Evidence::Observation::const_iterator i = e->begin(); i != e->end(); ++i) {
      clamped->clamp( clamped->fg().findVar(i->first), i->second);
    }
    clamped->init();
    clamped->run();
    
    outputFastaPerturbations(sample_iter->first, prior, clamped, priorFG, 
			     outNodes, *outstream);

    delete clamped;
  }

  delete prior;

  return 0;
}
