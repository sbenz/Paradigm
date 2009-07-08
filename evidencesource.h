#ifndef HEADER_EVIDENCESOURCE_H
#define HEADER_EVIDENCESOURCE_H

#include <vector>
#include <dai/alldai.h>
#include <dai/evidence.h>

#include "pathwaytab.h"

using namespace std;
using namespace dai;

class EvidenceFactorGen : public FactorGenerator {
public:
  Real _epsilon;
  Real _epsilon0;

  EvidenceFactorGen(Real epsilon=0.01, Real epsilon0=0.2) : 
    _epsilon(epsilon),
    _epsilon0(epsilon0) {}
  ~EvidenceFactorGen() {}
  
  void generateValues(const vector< string >& edge_types, 
		      vector< Real >& outVals) const;
};

typedef map<string, map<string,int> > SampleEvidMap;

class EvidenceSource 
{
private:
  vector<double> cutoffs;
  PropertySet options;
  double epsilon;
  double epsilon0;
  string attachPoint;
  string _suffix;
  string _evidenceFile;

  // add these from main.cpp
  vector<SampleEvidMap> _disc;
  vector<string> _sampleNames;
  vector<string> _sampleFactors;
  vector<int> _sampleFactorNum;

public:
  /// Default constructor
  EvidenceSource() : cutoffs(), 
		     epsilon(0.01), 
		     epsilon0(0.2), 
		     attachPoint(), 
		     _suffix(),
		     _evidenceFile(),
		     _disc(),
		     _sampleNames(),
		     _sampleFactors(),
		     _sampleFactorNum()

  {
    setCutoffs("-1.3;1.3");
  }
  /// Copy constructor
  EvidenceSource(const EvidenceSource &x) 
  : cutoffs(x.cutoffs), 
    epsilon(x.epsilon), 
    epsilon0(x.epsilon0), 
    attachPoint(x.attachPoint), 
    _suffix(x._suffix),
    _evidenceFile(x._evidenceFile),
    _disc(x._disc),
    _sampleNames(x._sampleNames),
    _sampleFactors(x._sampleFactors),
    _sampleFactorNum(x._sampleFactorNum) {}
  
  /// Assignment operator
  EvidenceSource& operator=(const EvidenceSource &x) {
    if (this != &x) {
      cutoffs = x.cutoffs;
      options = x.options;
      epsilon = x.epsilon;
      epsilon0 = x.epsilon0;
      attachPoint = x.attachPoint;
      _suffix = x._suffix;
      _evidenceFile = x._evidenceFile;
      _disc = x._disc;
      _sampleNames = x._sampleNames;
      _sampleFactors = x._sampleFactors;
      _sampleFactorNum = x._sampleFactorNum;
    }
    return *this;
  }

  /// Useful constructor
  EvidenceSource(PropertySet &p, string base);

  void setCutoffs(string discLimits);
  int discCutoffs (float x);  
  string genFactorString(int nodeID, int state);

  void loadFromFile(PathwayTab& p, Evidence& e);
  void addNodeNumber(int node_number, string name, string type);

  const string& evidenceFile() {return _evidenceFile;}
  const vector<string>& sampleNames() {return _sampleNames;}
  const int factorCount(size_t sample) {return _sampleFactorNum.at(sample);}
  const string& factorString(size_t sample) {return _sampleFactors.at(sample);}
};

void Tokenize(const string& str,
	      vector<string>& tokens,
	      const string& delimiters = " ");

#endif
