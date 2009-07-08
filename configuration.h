#ifndef HEADER_CONFIGURATION_H
#define HEADER_CONFIGURATION_H

#include <fstream>
#include <vector>
#include <dai/properties.h>
#include <dai/smallset.h>
#include <dai/util.h>

#define THROW(msg) throw std::runtime_error(msg)

using namespace dai;

class RunConfiguration 
{
public:
  static const std::string INFERENCE_CONF_TOKEN;
  static const std::string EVIDENCE_CONF_TOKEN;
  static const std::string PATHWAY_STRUCTURE_CONF_TOKEN;
  static const std::string EM_STEP_CONF_TOKEN;
  static const std::string EM_CONF_TOKEN;

  static const std::string INFERENCE_MATCH_TOKEN;

  typedef std::map< std::string, SmallSet< std::string > > EMStep;
  typedef std::vector< EMStep > EMSteps;

private:
  std::vector<PropertySet> _inferences;
  std::vector<PropertySet> _evidences;
  EMSteps _emsteps;
  PropertySet _em;
  
public:

  /// Default constructor
  RunConfiguration() : _inferences(), _evidences(), _emsteps(), _em() {
    _em.Set("max_iters", 0);
  }

  /// Copy constructor
  RunConfiguration(const RunConfiguration &x) : 
    _inferences(x._inferences), 
    _evidences(x._evidences),
    _emsteps(x._emsteps),
    _em(x._em)
  {}

  /// Assignment operator
  RunConfiguration& operator=(const RunConfiguration &x) {
    if (this != &x) {
      _inferences = x._inferences;
      _evidences = x._evidences;
      _emsteps = x._emsteps;
      _em = x._em;
    }
    return *this;
  }
  
  /// Useful constructor
  RunConfiguration(const std::string& configure_filename);

  PropertySet& getInferenceProperties(const std::string& pathway_filename);

  void addConfigurations(std::istream& is);

  PropertySet& evidence(size_t i);

  size_t evidenceSize();

  const PropertySet& emProps() { return _em; }
  
  const EMSteps emSteps() const {return _emsteps;}
};

#endif
