#include <map>
#include <vector>
#include <iostream>
#include <sstream>

#include <dai/index.h>
#include "pathwaytab.h"

using namespace std;

const size_t PathwayTab::VARIABLE_DIMENSION = 3;

const std::string PathwayTab::DEFAULT_INTERACTION_MAP = 
  "-dt>	genome	mRNA	positive\n"
  "-dr>	mRNA	protein	positive\n"
  "-dp>	protein	active	positive\n"
  "-t>	active	mRNA	positive\n"
  "-t|	active	mRNA	negative\n"
  "-a>	active	active	positive\n"
  "-a|	active	active	negative\n"
  "-ap>	active	active	positive\n"
  "-ap|	active	active	negative\n"
  "->	active	active	positive\n"
  "-|	active	active	negative\n"
  "<->	active	active	positive\n"
  "component>	active	active	positive\n"
  "member>	active	active	positive\n"
  ;

const std::string PathwayTab::CENTRAL_DOGMA = 
  "genome	mRNA	-dt>\n"
  "mRNA	protein	-dr>\n"
  "protein	active	-dp>\n"
  ;

const std::string PathwayTab::OBSERVATION_INTERACTION = "-obs>";

size_t countVotesRepressorDominates(size_t down, size_t up) {
  if (up > 0 & up > down) {
    return 2;
  } else if (down > 0 & down >= up) {
    return 0;
  } else {
    return 1;
  }
}

void RepressorDominatesVoteFactorGenerator::generateValues(const vector< string >& edge_types, 
							   vector< Real >& v) const {
  Real minor = _epsilon / 2;
  Real major = 1 - _epsilon;
  vector< size_t > dims;
  dims.reserve(edge_types.size());
  for (size_t i = 0; i < edge_types.size(); ++i) {
    dims.push_back(PathwayTab::VARIABLE_DIMENSION);
  }
  dai::MultiFor s(dims);
  for ( ; s.valid(); ++s) {
    std::vector< size_t > votes(PathwayTab::VARIABLE_DIMENSION,0);
    for (size_t i = 0; i < dims.size(); ++i) {
      if (edge_types[i] == "negative") {
	++votes[ votes.size() - 1 - s[i] ];
      } else {
	++votes[s[i]];
      }
    }
    size_t expected_state = countVotesRepressorDominates(votes[0], votes[2]);
    for (size_t i = 0; i < PathwayTab::VARIABLE_DIMENSION; ++i) {
      v.push_back(i == expected_state ? major : minor);
    }
  }
}

void SingleMemberNeededFactorGenerator::generateValues(const vector< string >& edge_types, 
							   vector< Real >& v) const {
  Real minor = _epsilon / 2;
  Real major = 1 - _epsilon;
  vector< size_t > dims;
  dims.reserve(edge_types.size());
  for (size_t i = 0; i < edge_types.size(); ++i) {
    dims.push_back(PathwayTab::VARIABLE_DIMENSION);
  }
  dai::MultiFor s(dims);

  for ( ; s.valid(); ++s) {
	size_t highestIndex = s[0];
    for (size_t i = 0; i < edge_types.size(); ++i) {
	  if(s[i] > highestIndex)
		highestIndex = s[i];
    }

    for (size_t i = 0; i < PathwayTab::VARIABLE_DIMENSION; ++i) {
	  v.push_back(i == highestIndex ? major : minor);
    }
  }
}
void readInteractionMap(istream& is, 
			map< string, vector< string > >& out_imap) {
  string line;
  while(getline(is, line)) {
    vector< string > vals;
    dai::tokenizeString(line, vals);
    if (vals.size() != 4) {
      THROW("Interaction map lines must have 4 entries");
    }
    string interaction = *vals.begin();
    vals.erase(vals.begin());
    out_imap[interaction] = vals;
  }
}

GeneProteinExpressionModel::GeneProteinExpressionModel(istream& is)
  : _states(), _steps() {
  string line;
  while(getline(is, line)) {
    vector<string> vals;
    dai::tokenizeString(line, vals);
    if (vals.size() != 3) {
      THROW("Must have three values per line in central dogma");
    }
    _steps.insert(vals[2]);
    _states.insert(vals[0]);
    _states.insert(vals[1]);
  }
}

void GeneProteinExpressionModel::addGeneDogma(const string& genename, 
					      PathwayTab& pathway_graph) {
  set< string >::iterator state_iterator = _states.begin();
  for ( ; state_iterator != _states.end(); ++state_iterator) {
    PathwayTab::Node node(genename, *state_iterator);
    pathway_graph.addNode(node);
  }
  set< string >::iterator step_iterator = _steps.begin();
  for ( ; step_iterator != _steps.end(); ++step_iterator) {
    pathway_graph.addInteraction(genename, genename, *step_iterator);
  }
}


PathwayTab::PathwayTab(istream& pathway_stream, 
		       istream& imap_stream, 
		       istream& dogma_stream,
		       const PropertySet& props)
  : _nodemap(),
    _nodevector(),
    _parents(),
    _entities(),
    _dogma(dogma_stream),
    _imap(),
    _props(props),
    _factorGenLookup(),
    _defaultFactorGen(new RepressorDominatesVoteFactorGenerator()) {
  vector< vector< string > > entity_lines;
  vector< vector< string > > interaction_lines;
  string line;

  readInteractionMap(imap_stream, _imap);

  while(getline(pathway_stream, line)) {
    vector< string > vals;
    dai::tokenizeString(line, vals);
    if (vals.size() == 2) {
      entity_lines.push_back(vals);
    } else if (vals.size() == 3) {
      interaction_lines.push_back(vals);
    } else {
      THROW("Must have either two or three entries per line");
    }
  }
  vector< vector< string > >::iterator v = entity_lines.begin();
  for ( ; v != entity_lines.end(); ++v ) {
    addEntity(v->at(1), v->at(0));
  }
  for (v = interaction_lines.begin(); v != interaction_lines.end(); ++v) {
    addInteraction(v->at(0), v->at(1), v->at(2));
  }

  if (_props.hasKey("max_in_degree")) {
    int max_degree = props.getStringAs<int>("max_in_degree");
    splitHighInDegree(max_degree);
  }
}

void PathwayTab::addEntity(const string& entity, const string& type) {
  if (_entities.count(entity) == 0) {
    _entities[entity] = type;
    if (type == "protein") {
      _dogma.addGeneDogma(entity, *this);
    } else {
      Node n(entity, "active");
      addNode(n);
    }
  }
}

Var PathwayTab::addObservationNode(const string& entity, 
				   const string& on_type, 
				   const string& obs_type) {
  Node obs_node(entity, obs_type);
  Node hidden_node;
  addNode(obs_node);
  getAppropriateEntityNode(entity, on_type, hidden_node);
  addEdge(hidden_node, obs_node, OBSERVATION_INTERACTION);

  return Var(_nodemap[obs_node], VARIABLE_DIMENSION);
}

void PathwayTab::addInteraction(const string& entity_from, 
				const string& entity_to, 
				const string& interaction) {
  if (_imap.count(interaction) == 0) {
    THROW("Unrecognized interaction type: " + interaction);
  }
  vector< string > i = _imap[interaction];
  assert(i.size() == 3);
  addEntity(entity_from);
  addEntity(entity_to);
  Node node_from;
  Node node_to;
  getAppropriateEntityNode(entity_from, i[0], node_from);
  getAppropriateEntityNode(entity_to, i[1], node_to);
  if (node_from == node_to) {
    return;
  }
  addNode(node_from);
  addNode(node_to);
  addEdge(node_from, node_to, i[2]);
}

void PathwayTab::addNode(Node nodename) {
  if (_nodemap.count(nodename) == 0) {
    map< Node, string > empty;
    _nodemap[nodename] = _nodevector.size();
    _nodevector.push_back(nodename);
    _parents[nodename] = empty;
  }
}

void PathwayTab::addEdge(const Node& from, const Node& to, const string& lbl) {
  _parents[to][from] = lbl;
}

void PathwayTab::getAppropriateEntityNode(const string& entity, 
					  const string& species,
					  Node& out_node) {
  out_node.first = entity;
  if (_entities[entity] == "protein") {
    out_node.second = species;
  } else {
    out_node.second = "active";
  }
}

void PathwayTab::addFactorGenerator(const string& entity_type, 
				    const string& node_type,
				    FactorGenerator* factor_gen) {
  pair< string, string > entry(entity_type, node_type);
  _factorGenLookup[entry] = factor_gen;
}			  

void PathwayTab::printNodeMap(ostream& to, const string& prefix) {
  for (size_t i = 0; i < _nodevector.size(); ++i) {
    to << prefix << i 
       << '\t' << _nodevector[i].first 
       << '\t' << _nodevector[i].second << endl;
  }
}

void PathwayTab::printDaiFactorSection(ostream& to) {
  size_t factor_count = 0;
  map< Node, map< Node, string > >::iterator c_iter = _parents.begin();
  vector< Node > corder;
  corder.reserve(_parents.size());
  for ( ; c_iter != _parents.end(); ++c_iter) {
    factor_count += (c_iter->second.size() > 0);
    corder.push_back(c_iter->first);
  }
  to.precision(6);
  to << fixed;
  to << factor_count << endl;

  sort(corder.begin(), corder.end());
  vector< Node >::iterator corder_iter = corder.begin();
  for ( ; corder_iter != corder.end(); ++corder_iter) {
    map< Node, string >& pmap = _parents[*corder_iter];

    if (pmap.size() == 0) {
      continue;
    }

    vector< Node > parents;
    vector< string > edge_types;


    if (pmap.count(*corder_iter) > 0) {
      // is this causing a seg fault?
      // pmap.erase(c_iter->first);
    }
    map< Node, string >::iterator p_iter =  pmap.begin();
    for ( ; p_iter != pmap.end(); ++p_iter) {
      parents.push_back(p_iter->first);
    }
    sort(parents.begin(), parents.end());

    /// Factor line: number of variables in factor
    to << endl << (parents.size() + 1) << endl;
    /// Output variable ids
    to << _nodemap[*corder_iter];
    vector< Node >::iterator pnode_iter = parents.begin();
    for ( ; pnode_iter != parents.end(); ++pnode_iter) {
      to << ' ' << _nodemap[*pnode_iter];
      edge_types.push_back(pmap[*pnode_iter]);
    }
    to << endl;
    /// Output variable dimensions
    to << VARIABLE_DIMENSION;
    size_t total_dimension = VARIABLE_DIMENSION;
    for (size_t i = 0; i < parents.size(); ++i) {
      to << ' ' << VARIABLE_DIMENSION;
      total_dimension *= VARIABLE_DIMENSION;
    }
    to << endl;
    
    vector< Real > factor_vals;
    factor_vals.reserve(total_dimension);
    _defaultFactorGen->generateValues(edge_types, factor_vals);
    to << factor_vals.size() << endl;
    for (size_t i = 0; i < factor_vals.size(); ++i) {
      to << i << '\t' << factor_vals[i] << endl;
    }
  }
}

void PathwayTab::splitNodeParents(const Node& n, const size_t maxParents) {
  map< Node, map< Node, string > >::iterator pmap_i = _parents.find(n);
  assert(pmap_i != _parents.end());
  unsigned int numNodes = (pmap_i->second.size() / maxParents) 
    + (pmap_i->second.size() % maxParents > 0);
  if (numNodes > maxParents) {numNodes = maxParents;}
  if (numNodes > 1) {
    vector<Node> newNodes;

    // Step 1 of 4: create intermediate nodes
    for (size_t i = 0; i < numNodes; ++i) {
      stringstream s;
      s << n.first << "__" << i;
      Node newNode(s.str(), n.second);

      newNodes.push_back(newNode);
      addEntity(newNode.first, _entities[n.second]);
    }
    
    // Step 2 of 4: connect parents to intermediate nodes
    map< Node, string >::iterator parent_i = pmap_i->second.begin();
    for (int nodeIndex = 0; parent_i != pmap_i->second.end(); ++parent_i) {
      addEdge(parent_i->first, newNodes[nodeIndex], parent_i->second);
      nodeIndex = (nodeIndex + 1) % numNodes;
    }

    // Step 3 of 4: remove edges to child and replace with intermediate edges
    _parents[n].clear();
    for (size_t i = 0; i < numNodes; ++i) {
      addEdge(newNodes[i], n, _parents[newNodes[i]].begin()->second);
    }

    // Step 4 of 4: recursively break up nodes, if necessary
    for (size_t i = 0; i < numNodes; ++i) {
      splitNodeParents(newNodes[i], maxParents);
    }
  }
}

void PathwayTab::splitHighInDegree(const int maxParents) {
  vector< Node > origN = _nodevector;
  for (vector<Node>::iterator i = origN.begin(); i != origN.end(); ++i) {
    splitNodeParents(*i, maxParents);
  }
}
void PathwayTab::generateFactorValues(const Node& child, 
				      const vector< string >& edge_types,
				      vector< Real >& outValues) const {
  /// \todo make this polymorphic based on the entity type and node sub-type
  map< string, string >::const_iterator entity_type_i;
  entity_type_i = _entities.find(child.first);
  if (entity_type_i == _entities.end()) {
    THROW("Could not find entity in generateFactorValues");
  }
  const string& entity_type = entity_type_i->second;
  pair<string, string> lookup(entity_type, child.second);
  map< pair<string, string>, FactorGenerator* >::const_iterator fgen_lookup;
  fgen_lookup = _factorGenLookup.find(lookup);
  if (fgen_lookup != _factorGenLookup.end()) {
    FactorGenerator* f = fgen_lookup->second;
    f->generateValues(edge_types, outValues);
  } else {
    _defaultFactorGen->generateValues(edge_types, outValues);
  }
}

void PathwayTab::constructFactors(const RunConfiguration::EMSteps& sp,
				  vector< Factor >& outFactors,
				  vector< MaximizationStep >& outMsteps) const{
  vector< vector < SharedParameters::FactorOrientations > > var_orders;
  vector< vector < size_t > > sp_total_dim;
  var_orders.resize(sp.size());
  sp_total_dim.resize(sp.size());
  for (size_t i = 0; i < var_orders.size(); ++i) {
    var_orders[i].resize(sp[i].size());
    sp_total_dim[i] = vector< size_t >(sp[i].size());
  }

  map< Node,map< Node,string > >::const_iterator child_iter = _parents.begin();
  for ( ; child_iter != _parents.end(); ++child_iter) {
    const Node& child_node = child_iter->first;
    const map< Node, string >& pmap = child_iter->second;
    
    if (pmap.size() == 0) {
      continue;
    }

    vector< Var > factor_vars;
    factor_vars.reserve(pmap.size() + 1);
    Var child_var(Var(getNodeIndex(child_node), VARIABLE_DIMENSION));
    factor_vars.push_back(child_var);

    vector< string > edge_types;
    edge_types.reserve(pmap.size());
    size_t total_dimension = VARIABLE_DIMENSION;
    map< Node, string >::const_iterator p_iter = pmap.begin();
    for ( ; p_iter != pmap.end(); ++p_iter) {
      const Node& parent_node = p_iter->first;
      const string& edge_type = p_iter->second;
      Var parent_var(getNodeIndex(parent_node), VARIABLE_DIMENSION);
      factor_vars.push_back(parent_var);
      edge_types.push_back(edge_type);
      total_dimension *= VARIABLE_DIMENSION;
    }
    
    vector< Real > factor_vals;
    factor_vals.reserve(total_dimension);
    generateFactorValues(child_node, edge_types, factor_vals);
    assert(factor_vals.size() == total_dimension);
    
    Factor f(factor_vars, factor_vals);
    outFactors.push_back(f);

    for (size_t i = 0; i < sp.size(); ++i) {
      map< string, SmallSet<string> >::const_iterator jit = sp[i].begin();
      for(size_t j = 0; j < sp[i].size(); ++j, ++jit) {
	SmallSet< string > eset(edge_types.begin(), edge_types.end(), 
				edge_types.size());
	const string& spec_subtype = jit->first;
	const string& node_subtype = child_node.second;
	if (spec_subtype == node_subtype
	    && edge_types.size() == eset.size()
	    && jit->second == eset) {
	  vector< Var > o;
	  o.reserve(factor_vars.size());
	  o.push_back(factor_vars[0]);
	  SmallSet< string >::const_iterator sp_iter = jit->second.begin();
	  for ( ; sp_iter != jit->second.end(); ++sp_iter) {
	    for (size_t k = 0; k < edge_types.size(); ++k) {
	      if (edge_types[k] == *sp_iter) {
		o.push_back(factor_vars[k+1]);
		break;
	      }
	    }
	  }
	  assert (o.size() == eset.size() + 1);
	  var_orders[i][j][outFactors.size()-1] = o;
	  sp_total_dim[i][j] = total_dimension;
	}
      }
    }
  }

  // Construct all Msteps
  for (size_t i= 0; i < var_orders.size(); ++i) {
    vector< SharedParameters> spvec;
    spvec.reserve(var_orders[i].size());
    RunConfiguration::EMStep::const_iterator spec_iterator = sp[i].begin();
    for (size_t j = 0; j < var_orders[i].size(); ++j, ++spec_iterator) {
      if (var_orders[i][j].size() == 0) {
	cerr << "!! Did not find any variables of sub-type '" 
	     << spec_iterator->first << "' with incoming edges matching: "
	     << endl;
	SmallSet< string >::const_iterator it = spec_iterator->second.begin();
	for ( ; it != spec_iterator->second.end(); ++it) {
	  cerr << "!!  " << *it  << endl;
	}
	continue;
      }
      PropertySet props;
      props.Set("total_dim", sp_total_dim[i][j]);
      props.Set("target_dim", VARIABLE_DIMENSION);
      ParameterEstimation* pe;
      pe = ParameterEstimation::construct("ConditionalProbEstimation", props);
      spvec.push_back(SharedParameters(var_orders[i][j], pe, 1));
    }
    if (spvec.size() > 0) {
      outMsteps.push_back(MaximizationStep(spvec));
    } else {
      cerr << "!! em_step number " << i 
	   << " had no matching nodes in the pathway" << endl;
    }
  }
}

map< long, string > PathwayTab::getOutputNodeMap() {
  map< long, string > result;
  for (size_t i = 0; i < _nodevector.size(); ++i) {
    if (_nodevector[i].second == "active") {
      result[i] = _nodevector[i].first;
    }
  }
  return result;
}

void PathwayTab::dumpNodeIndexMap() const {
  map< Node, size_t>::const_iterator i = _nodemap.begin();
  for( ; i != _nodemap.end(); ++i) {
    const Node& n = i->first;
    size_t idx = i->second;
    const Node& nv = _nodevector[idx];
    cerr << idx << '\t' 
	 << n.first << '\t' << n.second << '\t' 
	 << nv.first << '\t' << nv.second << endl;
  }
}


int PathwayTab::debugPrintParents(size_t node_i) {
  Node c = _nodevector[node_i];
  if (_parents.count(c) > 0) {
    map<Node, string>::iterator i = _parents[c].begin();
    for (; i != _parents[c].end(); ++i) {
      cout << _nodemap[i->first] << '\t' << i->first.first << "\t" << i->first.second << '\t' << i->second <<endl;
    }
  }
  return -1;
}
