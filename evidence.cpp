#include <fstream>

#include "evidence.h"

#define THROW(msg) throw std::runtime_error(msg)

Evidence::Evidence(PropertySet &p, string base) : cutoffs(),
						  epsilon(0.01), 
						  epsilon0(0.2), 
						  attachPoint(), 
						  _evidenceFile() 
{
  if (p.hasKey("disc"))
    setCutoffs(p.GetAs<string>("disc"));
  else
    setCutoffs("-1.3,1.3");

  if (p.hasKey("epsilon"))
    epsilon = p.GetAs<double>("epsilon");

  if (p.hasKey("epsilon0"))
    epsilon0 = p.GetAs<double>("epsilon0");

  if (p.hasKey("node"))
    attachPoint = p.GetAs<string>("node");
  else
    THROW("Evidence conf. is missing the required property \"node\"");

  if (p.hasKey("suffix"))
    _evidenceFile = base + p.GetAs<string>("suffix");
  else
    THROW("Evidence conf. is missing the required property \"suffix\"");
}

void Evidence::setCutoffs(string discLimits)
{
  vector<string> cutoffsStr;
  Tokenize(discLimits,cutoffsStr,";");
  for(size_t i = 0; i < cutoffsStr.size(); i++)
    {
      cerr << "Added cutoff " << atof(cutoffsStr[i].c_str()) << endl;
      cutoffs.push_back(atof(cutoffsStr[i].c_str()));
    }
}

int Evidence::discCutoffs  (float x)
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

string Evidence::genFactorString(int nodeID, int state)
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

void Evidence::loadFile()
{
  ifstream infile;
  infile.open( _evidenceFile.c_str() );
  if( infile.is_open() ) 
    {
      string line;
      if(!getline(infile,line))
	return;
      vector<string> header;
      Tokenize(line,header,"\t");
      header.erase(header.begin());
      //cout << line << endl;
      size_t lineCount = 0;
      while(getline(infile,line))
	{
	  map<string, map<string,int> > evidenceMap;
	  if(_disc.size() > lineCount)
	    evidenceMap = _disc[lineCount];
	  //cout << line << endl;
	  vector<string> vals;
	  Tokenize(line,vals,"\t");
	  string patient = vals[0];
	  vals.erase(vals.begin());
			
	  for(size_t i = 0; i < vals.size(); i++)
	    {
	      if(strcmp(vals[i].c_str(),"NA")==0)
		continue;
	      double evidence = atof(vals[i].c_str());
	      if(evidence == 0.0)
		{
		  cerr << "Evidence column must be numeric" << endl;
		  return;
		}
	      //cerr << header[i] << endl;
	      evidenceMap[header[i]][attachPoint] = discCutoffs(evidence);
	    }
	  _disc.push_back(evidenceMap);
	  _sampleNames.push_back(patient);
	  _sampleFactors.push_back("");
	  _sampleFactorNum.push_back(0);
	  lineCount++;
	}
      infile.close();
    }
  return;
}

void Evidence::addNodeNumber(int node_number, string name, string type)
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

