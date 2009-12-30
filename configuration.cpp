#include "configuration.h"

const std::string RunConfiguration::INFERENCE_CONF_TOKEN("inference");
const std::string RunConfiguration::EVIDENCE_CONF_TOKEN("evidence");
const std::string RunConfiguration::PATHWAY_CONF_TOKEN("pathway");
const std::string RunConfiguration::EM_STEP_CONF_TOKEN("em_step");
const std::string RunConfiguration::EM_CONF_TOKEN("em");

const std::string RunConfiguration::INFERENCE_MATCH_TOKEN("pathway_match");

RunConfiguration::RunConfiguration(const std::string& configure_filename)
{
  std::ifstream is(configure_filename.c_str());
  if (!is.is_open())
    THROW("couldn't open configuration file");
  addConfigurations(is);
}

PropertySet& 
RunConfiguration::getInferenceProperties(const std::string& pathway_filename) 
{
  for (size_t i = 0; i < _inferences.size(); i++)
    {
      if (!_inferences[i].hasKey(INFERENCE_MATCH_TOKEN))
	{
	  return _inferences[i];
	}

      std::string m = _inferences[i].GetAs<std::string>(INFERENCE_MATCH_TOKEN);
      if (pathway_filename.find(m) != std::string::npos)
	{
	  return _inferences[i];
	}
    }
  THROW("No matching inference configurations");
}

void
RunConfiguration::addConfigurations(std::istream& is)
{
  std::string type;
  PropertySet conf;
  while(is >> type)
    {
      is >> conf;
      if (type == INFERENCE_CONF_TOKEN) {
	_inferences.push_back(conf);
      } else if (type == EVIDENCE_CONF_TOKEN) {
	_evidences.push_back(conf);
      } else if (type == EM_STEP_CONF_TOKEN) {
	std::set<PropertyKey> keys = conf.keys();
	std::set<PropertyKey>::iterator i = keys.begin();
	EMStep e;
	for ( ; i != keys.end(); ++i) {
	  std::vector< std::string > edges;
	  tokenizeString(conf.GetAs<std::string>(*i), edges, ";");
	  SmallSet< std::string > s(edges.begin(), edges.end(), edges.size());
	  e[*i] = s;
	}
	_emsteps.push_back(e);
      } else if (type == EM_CONF_TOKEN) {
	_em = conf;
      } else if (type == PATHWAY_CONF_TOKEN) {
	_path = conf;
      } else {
	THROW("Expecting an inference or evidence token in conf file");
      }
    }
}

size_t RunConfiguration::evidenceSize() {return _evidences.size();}
PropertySet& RunConfiguration::evidence(size_t i) {return _evidences.at(i);}
