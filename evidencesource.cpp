#include <fstream>

#include "evidencesource.h"

#define THROW(msg) throw std::runtime_error(msg)

void EvidenceFactorGen::generateValues(const vector< string >& edge_types, 
		      vector< Real >& outVals) const
{
  assert(edge_types.size() == 1);
  Real major = 1 - _epsilon;
  Real minor = _epsilon / 2;
  Real major0 = 1 - _epsilon0;
  Real minor0 = _epsilon0 / 2;

  outVals.push_back(major);
  outVals.push_back(minor0);
  outVals.push_back(minor);

  outVals.push_back(minor);
  outVals.push_back(major0);
  outVals.push_back(minor);

  outVals.push_back(minor);
  outVals.push_back(minor0);
  outVals.push_back(major);
}


EvidenceSource::EvidenceSource(PropertySet &p, string base) : 
  cutoffs(),
  epsilon(0.01), 
  epsilon0(0.2), 
  attachPoint(), 
  _evidenceFile()
{
  if (p.hasKey("disc"))
    setCutoffs(p.GetAs<string>("disc"));
  else
    setCutoffs("-1.3,1.3");

  if (p.hasKey("epsilon")) {
    epsilon = p.getStringAs<double>("epsilon");
  }

  if (p.hasKey("epsilon0")) {
    epsilon0 = p.getStringAs<double>("epsilon0");
  }

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

string EvidenceSource::genFactorString(int nodeID, int state)
{
  int varDimension = 3;
  ostringstream result;
	
  double eps = (state == 1) ? epsilon0 : epsilon;
  double e = eps / (double)(varDimension - 1);
  result << "\n" << "1" << "\n" << (nodeID) << "\n" << (varDimension) << "\n" << (varDimension);
  for(int i = 0; i < varDimension; i++)
    {
      if(i == state)
	result << "\n" << (i) << "\t" << (1 - eps);
      else
	result << "\n" << (i) << "\t" << (e);
    }
  result << endl;
  return result.str();
}

double stringToDouble(const string& s) {
  stringstream ss(s);
  double result;
  if (!(ss >> result).eof()) {
    throw "String " + s + " can not be converted to double.";
  }
  return result;
}

void EvidenceSource::loadFromFile(PathwayTab& p, Evidence& e)
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

    FactorGenerator* fgen = new EvidenceFactorGen(epsilon, epsilon0);
    p.addFactorGenerator("protein", _suffix, fgen);
    
    while(getline(infile,line)) {
      vector<string> vals;
      Tokenize(line,vals,"\t");
      if (vals.size() > vars.size() + 1) {
	throw "Entries in evidence line does not match header length";
      }
      string sample = vals[0];
      e.ensureSampleExists(sample);
      _sampleNames.push_back(sample);
      vals.erase(vals.begin());
      
      for(size_t i = 0; i < vals.size(); i++) {
	if(strcmp(vals[i].c_str(),"NA")==0)
	  continue;
	double evidence = stringToDouble(vals[i]);
	e.addObservation(sample, vars[i], discCutoffs(evidence));
      }
    }
    infile.close();
  }
  return;
}

void EvidenceSource::addNodeNumber(int node_number, string name, string type)
{
  for (size_t sample = 0; sample < _disc.size(); sample++)
    {
      if (_disc[sample].count(name) > 0 && _disc[sample][name].count(type) > 0)
	{
	  int state = _disc[sample][name][type];
	  _sampleFactors[sample] += genFactorString(node_number, state);
	  _sampleFactorNum[sample]++;

	  _disc[sample][name].erase(type);
	  if (_disc[sample][name].size() == 0)
	    _disc[sample].erase(name);
	}
    }
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

