#include <fstream>

#include "common.h"
#include "evidencesource.h"

#define THROW(msg) throw std::runtime_error(msg)

EvidenceFactorGen::EvidenceFactorGen(const PropertySet& p) : _params() 
{
  _params.reserve(9);
  if (p.hasKey("factorParams")) {
    vector<string> paramsStr;
    Tokenize(p.getStringAs<string>("factorParams"),paramsStr,";");
    if (paramsStr.size() != 9) {
      THROW("must have 9 elements in factorParams");
    }
    for(size_t i = 0; i < paramsStr.size(); i++)
    {
      _params.push_back(atof(paramsStr[i].c_str()));
    }
    return;
  }

  if (!p.hasKey("epsilon") || !p.hasKey("epsilon0")) {
    cout << "woops" << endl;
    cout << p << endl;
    THROW("must specify either factorParams, or epsilon and epsilon0");
  }

  double epsilon = p.getStringAs<double>("epsilon");
  double epsilon0 = p.getStringAs<double>("epsilon0");
  Real major = 1 - epsilon;
  Real minor = epsilon / 2;
  Real major0 = 1 - epsilon0;
  Real minor0 = epsilon0 / 2;
  
  _params.push_back(major);
  _params.push_back(minor0);
  _params.push_back(minor);
  
  _params.push_back(minor);
  _params.push_back(major0);
  _params.push_back(minor);
  
  _params.push_back(minor);
  _params.push_back(minor0);
  _params.push_back(major);
}

void EvidenceFactorGen::generateValues(const vector< string >& edge_types, 
		      vector< Real >& outVals) const
{
  assert(edge_types.size() == 1);
  for (size_t i = 0; i < _params.size(); i++) {
    outVals.push_back(_params[i]);
  }
}

EvidenceSource::EvidenceSource(PropertySet &p, string base) : 
  cutoffs(),
  options(p),
  attachPoint(), 
  _evidenceFile()
{
  if (p.hasKey("disc"))
    setCutoffs(p.GetAs<string>("disc"));
  else
    setCutoffs("-1.3,1.3");

  if (p.hasKey("node"))
    attachPoint = p.GetAs<string>("node");
  else
    THROW("EvidenceSource conf. is missing the required property \"node\"");

  if (p.hasKey("suffix")) {
    _suffix = p.GetAs<string>("suffix");
    _evidenceFile = base + _suffix;
  } else {
    THROW("EvidenceSource conf. is missing the required property \"suffix\"");
  }
}

void EvidenceSource::setCutoffs(string discLimits)
{
  vector<string> cutoffsStr;
  Tokenize(discLimits,cutoffsStr,";");
  for(size_t i = 0; i < cutoffsStr.size(); i++)
    {
      if(VERBOSE)
		  cerr << "Added cutoff " << atof(cutoffsStr[i].c_str()) << endl;
      cutoffs.push_back(atof(cutoffsStr[i].c_str()));
    }
}

int EvidenceSource::discCutoffs  (float x)
{
  if(cutoffs.size() == 0)
    return 0;
  if(x < cutoffs[0])
    return 0;
  size_t i = 1;
  while (i < cutoffs.size())
    {
      if (x < cutoffs[i])
	return i;
      i++;
    }
  return i;
}

double stringToDouble(const string& s) {
  stringstream ss(s);
  double result;
  if (!(ss >> result).eof()) {
    THROW("String " + s + " can not be converted to double.");
  }
  return result;
}

void EvidenceSource::loadFromFile(PathwayTab& p, 
				  map<string, size_t>& sampleMap, 
				  vector<Observation>& sampleData) 
{
  ifstream infile;
  infile.open( _evidenceFile.c_str() );
  if( infile.is_open() ) {
    string line;

    if(!getline(infile,line))
      return;
    vector<string> header;
    Tokenize(line,header,"\t");
    header.erase(header.begin());

    vector<Var> vars;
    vars.reserve(header.size());

    for (size_t h = 0; h < header.size(); h++) {
      vars.push_back(p.addObservationNode(header[h], attachPoint, _suffix));
    }

    FactorGenerator* fgen = new EvidenceFactorGen(options);
    p.addFactorGenerator("protein", _suffix, fgen);
    
    while(getline(infile,line)) {
      vector<string> vals;
      Tokenize(line,vals,"\t");
      if (vals.size() > vars.size() + 1) {
	THROW("Entries in evidence line does not match header length");
      }
      string sample = vals[0];
      _sampleNames.push_back(sample);
      vals.erase(vals.begin());
      
      for(size_t i = 0; i < vals.size(); i++) {
	if(strcmp(vals[i].c_str(),"NA")==0)
	  continue;
	double evidence = stringToDouble(vals[i]);
	if (sampleMap.count(sample) == 0) {
	  sampleMap[sample] = sampleData.size();
	  sampleData.push_back(Observation());
	}
	size_t sample_idx = sampleMap[sample];
	sampleData[sample_idx].addObservation(vars[i], discCutoffs(evidence));
      }
    }
    infile.close();
  }
  return;
}


void Tokenize(const string& str,
	      vector<string>& tokens,
	      const string& delimiters)
{
  // Skip delimiters at beginning.
  string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  // Find first "non-delimiter".
  string::size_type pos     = str.find_first_of(delimiters, lastPos);

  while (string::npos != pos || string::npos != lastPos)
    {
      // Found a token, add it to the vector.
      tokens.push_back(str.substr(lastPos, pos - lastPos));
      // Skip delimiters.  Note the "not_of"
      lastPos = str.find_first_not_of(delimiters, pos);
      // Find next "non-delimiter"
      pos = str.find_first_of(delimiters, lastPos);
    }
}

