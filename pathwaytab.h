#ifndef HEADER_PATHWAYTAB_H
#define HEADER_PATHWAYTAB_H

#include <iostream>
#include <set>
#include <map>
#include <vector>
#include <sstream>

#include <dai/util.h>
#include <dai/var.h>
#include <dai/factor.h>
#include <dai/emalg.h>

#include "configuration.h"

using namespace std;
using namespace dai;

class FactorGenerator {
public:
  virtual void generateValues(const vector< string >& edge_types, 
			      vector< Real >& outVals) const = 0;
  virtual ~FactorGenerator() {}
};

class RepressorDominatesVoteFactorGenerator : public FactorGenerator{
private:
  Real _epsilon;
public:
  RepressorDominatesVoteFactorGenerator(double epsilon=0.001) : 
    _epsilon(epsilon) {}
  ~RepressorDominatesVoteFactorGenerator() {}
  void generateValues(const vector< string >& edge_types, 
		      vector< Real >& outVals) const;
};

class SingleMemberNeededFactorGenerator : public FactorGenerator{
private:
  Real _epsilon;
public:
  SingleMemberNeededFactorGenerator(double epsilon=0.001) : 
    _epsilon(epsilon) {}
  ~SingleMemberNeededFactorGenerator() {}
  void generateValues(const vector< string >& edge_types, 
		      vector< Real >& outVals) const;
};

class PathwayTab;

void readInteractionMap(istream& is, map< string, vector< string > >& out_imap);

class GeneProteinExpressionModel {
private:
  set< string > _states;
  set< string > _steps;
public:
  GeneProteinExpressionModel(istream& is);
  void addGeneDogma(const string& genename, PathwayTab& pathway_graph);
};

class PathwayTab {
public:
  typedef pair< string, string > Node;

  static const size_t VARIABLE_DIMENSION; // = 3
  static const std::string DEFAULT_INTERACTION_MAP;
  static const std::string CENTRAL_DOGMA;

private:
  static const string OBSERVATION_INTERACTION;

  map< Node, size_t > _nodemap;
  vector< Node > _nodevector;
  map< Node, map< Node, string > > _parents;
  map< string, string > _entities;				
  GeneProteinExpressionModel _dogma;
  map< string, vector< string > > _imap;
  
  PropertySet _props;
  map< pair<string, string>, FactorGenerator* > _factorGenLookup;
  FactorGenerator* _defaultFactorGen;

  PathwayTab(istream& pathway_stream, 
	     istream& imap_stream, 
	     istream& dogma_stream,
	     const PropertySet& props);
public:
  static PathwayTab create(istream& pathway_stream,
			   const PropertySet& props,
			   istream* imap_stream=NULL,
			   istream* dogma_stream=NULL) {
    istringstream is(DEFAULT_INTERACTION_MAP);
    istringstream ds(CENTRAL_DOGMA);
    if (imap_stream == NULL) {
      imap_stream = &is;
    }
    if (dogma_stream == NULL) {
      dogma_stream = &ds;
    }
    return PathwayTab(pathway_stream, *imap_stream, *dogma_stream, props);
  }
		      
  ~PathwayTab() { 
    delete _defaultFactorGen; 
    map< pair<string, string>, FactorGenerator* >::iterator i;
    i= _factorGenLookup.begin();
    for ( ; i != _factorGenLookup.end(); ++i) {
      delete i->second;
    }
  }

  void addEntity(const string& entity, const string& type="protein");
  Var addObservationNode(const string& entity, const string& on_type,
			 const string& obstype);
  void addInteraction(const string& entity_from, const string& entity_to, 
		      const string& interaction);

  void addNode(Node nodename);
  void addNode(const string& entity, const string& subtype);
  void addEdge(const Node& from, const Node& to, const string& label);
  void getAppropriateEntityNode(const string& entity, const string& species,
				Node& out_node);
  Node getNode(size_t i) const { return _nodevector[i]; }

  size_t getNodeIndex(const Node& n) const {
    map< Node, size_t >::const_iterator i = _nodemap.find(n);
    return i->second;
  }

  void dumpNodeIndexMap() const;

  string getInteraction(size_t child_i, size_t parent_i) const {
    Node c = getNode(child_i);
    Node p = getNode(parent_i);
    string result = "";
    map< Node, map< Node, string > >::const_iterator ci = _parents.find(c);
    if (ci != _parents.end()) {
      map< Node, string >::const_iterator pi = ci->second.find(p);
      if (pi != ci->second.end()) {
	result = pi->second;
      }
    }
    return result;
  }

  void addFactorGenerator(const string& entity_type, 
			  const string& node_type,
			  FactorGenerator* factor_gen);

  int debugPrintParents(size_t node_i);
  void printNodeMap(ostream& to=cout, const string& prefix="# ");
  void printDaiFactorSection(ostream& to=cout);

  void splitNodeParents(const Node& n, const size_t maxParents);
  void splitHighInDegree(const int maxParents);

  void generateFactorValues(const Node& child, 
			    const vector< string >& edge_types,
			    vector< Real >& outValues) const;
  void constructFactors(const RunConfiguration::EMSteps& sp,
			vector< Factor >& outFactors, 
			vector< MaximizationStep >& outMsteps) const;
  map< long, string > getOutputNodeMap();  
};

#endif
