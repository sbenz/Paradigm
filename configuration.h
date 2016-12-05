/********************************************************************************/
/* Copyright 2009-2011 -- The Regents of the University of California           */
/* This code is provided for research purposes to scientists at non-profit		*/
/*  organizations.  All other use is strictly prohibited.  For further			*/
/*  details please contact University of California, Santa Cruz or				*/
/*	Five3 Genomics, LLC (http://five3genomics.com).								*/
/********************************************************************************/

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
  static const std::string PATHWAY_CONF_TOKEN;
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
  PropertySet _path;

public:

  /// Default constructor
  RunConfiguration() : _inferences(), _evidences(), _emsteps(), _em(), _path(){
    _em.set("max_iters", 0);
  }

  /// Copy constructor
  RunConfiguration(const RunConfiguration &x) :
    _inferences(x._inferences),
    _evidences(x._evidences),
    _emsteps(x._emsteps),
    _em(x._em),
    _path(x._path)
  {}

  /// Assignment operator
  RunConfiguration& operator=(const RunConfiguration &x) {
    if (this != &x) {
      _inferences = x._inferences;
      _evidences = x._evidences;
      _emsteps = x._emsteps;
      _em = x._em;
      _path = x._path;
    }
    return *this;
  }

  /// Useful constructor
  RunConfiguration(const std::string& configure_filename);

  PropertySet& getInferenceProperties(const std::string& pathway_filename);

  void addConfigurations(std::istream& is);

  PropertySet& evidence(size_t i);

  PropertySet& pathwayProps() {return _path;}

  size_t evidenceSize();

  const PropertySet& emProps() { return _em; }

  const EMSteps emSteps() const {return _emsteps;}
};

#endif
