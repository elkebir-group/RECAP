/*
 * inputinstance.cpp
 *
 *  Created on: 10-jun-2019
 *      Author: M. El-Kebir
 */

#include "inputinstance.h"

InputInstance::InputInstance()
  : _scriptT()
{
}

bool InputInstance::hasMutationClusters(const std::string& sep) const
{
  const int mutSize = getMutationSetSize(sep);
  
  const int nrPatients = getNrPatients();
  for (int i = 0; i < nrPatients; ++i)
  {
    for (const CloneTree& T : _scriptT[i])
    {
      if (lemon::countNodes(T.tree()) < mutSize)
      {
        return true;
      }
      assert(lemon::countNodes(T.tree()) == mutSize);
    }
  }
  
  return false;
}

void InputInstance::expand(const std::string& sep)
{
  const int n = _scriptT.size();
  for (int i = 0; i < n; ++i)
  {
    CloneTreeVector ctv_i;
    
    for (const CloneTree& T : _scriptT[i])
    {
      CloneTree::expandTree(T, sep, ctv_i);
    }
    
    _scriptT[i] = ctv_i;
  }
}

void InputInstance::randPermutePerTree()
{
  const int n = _scriptT.size();
  for (int i = 0; i < n; ++i)
  {
   
    CloneTreeVector ctv_i;
    
    for (const CloneTree& T : _scriptT[i])
    {
      CloneTree::permuteTree(T, ctv_i);
    }
    
    _scriptT[i] = ctv_i;
  }
}

void InputInstance::randPermutePerPatient()
{
  const int n = _scriptT.size();
  for (int i = 0; i < n; ++i)
  {
   
    CloneTreeVector ctv_i;
      
    // Get one random permutation for this patient
    // Use this permutation for all trees
    std::map<std::string, std::string> labelMap;
    
    if (_scriptT[i].size() > 0){
      CloneTree::permuteMap(_scriptT[i][0], labelMap);
    
      for (const CloneTree& T : _scriptT[i])
      {
        CloneTree::permuteTree(T, ctv_i, labelMap);
      }
    
      _scriptT[i] = ctv_i;
    }
  }
}

int InputInstance::getMutationSetSize(const std::string& sep) const
{
  return getMutationSet(sep).size();
}

const StringSet& InputInstance::getMutationSet(const std::string& sep) const
{
  if (_mutationSet.empty())
  {
    updateMutationSet(sep);
  }
  return _mutationSet;
}

const StringSet& InputInstance::updateMutationSet(const std::string& sep) const
{
  _patient2mutations.clear();
  _patient2missingMutations.clear();
  _mutationCounts.clear();
  _mutationSet.clear();
  
    const int nrPatients = getNrPatients();
    _patient2mutations = std::vector<StringSet>(nrPatients);
    _patient2missingMutations = std::vector<StringSet>(nrPatients);
    for (int i = 0; i < nrPatients; ++i)
    {
      const CloneTreeVector& ctv_i = _scriptT[i];
      for (const CloneTree& T : ctv_i)
      {
        StringSet muts = T.getIdSet();
        for (const std::string& mut : muts)
        {
          StringVector s;
          boost::split(s, mut, boost::is_any_of(sep));
          
          _mutationSet.insert(s.begin(), s.end());
          _patient2mutations[i].insert(s.begin(), s.end());
        }
      }
      
      // Add mutations to count dictionary for this patient only once
      for (auto x: _patient2mutations[i]){
        if (_mutationCounts.find(x) == _mutationCounts.end()){
          _mutationCounts[x] = 1;
        }
        else{
          _mutationCounts[x] = _mutationCounts[x]+1;
        }
      }
    }
    
    
    for (int i = 0; i < nrPatients; ++i)
    {
      std::set_difference(_mutationSet.begin(), _mutationSet.end(),
                          _patient2mutations[i].begin(), _patient2mutations[i].end(),
                          std::inserter(_patient2missingMutations[i], _patient2missingMutations[i].begin()));
    }

  return _mutationSet;
}

void InputInstance::completePatientTrees(const std::string& sep, const std::string& dummyNode){
   
  // Find total number of patients
  const int nrPatients = getNrPatients();
  
  // Ensure mutation set has been identified
  updateMutationSet(sep);
  int numMut(_mutationSet.size());
  _mutationSet.insert(dummyNode);
  _mutationCounts["dummy"] = nrPatients;
    
  // Iterate over each patient
  for (int i = 0; i < nrPatients; ++i){
        
      // Iterate over each tree for each patient
      for (CloneTree& T : _scriptT[i]){
            
          // Create edge dictionary
          std::vector<StringPair> edgeMap;
            
          // Get root from T
          Node root = T.root();
          std::string rLabel = T.label(root);
            
          // Add edge from root to "missing" dummy node
          edgeMap.push_back(make_pair(rLabel, dummyNode));
            
          // Iterate over the missing mutations for each patient
          for (const std::string& mut : _patient2missingMutations[i]){
                edgeMap.push_back(make_pair(dummyNode, mut));
          }
            
          // Add edge to tree
          T.addEdges(edgeMap, sep);
      }
  }
  
  for (int i = 0; i < nrPatients; ++i){
    for (const CloneTree& T : _scriptT[i]){
            
      int c(0);
      StringSet muts = T.getIdSet();
      for (const std::string& mut : muts)
      {
        StringVector s;
        boost::split(s, mut, boost::is_any_of(sep));
              
        c += s.size();
      }
            
      if (c != numMut+1){
        throw std::runtime_error("Error: completing patient trees with missing mutations");
      }
    }
  }
}

void InputInstance::removeRareMutations(const std::string& sep, const int t){
  
  // Find total number of patients
  const int nrPatients = getNrPatients();
  
  // Make sure the minimum threshold is less than the number of patients
  if (t > nrPatients)
    throw std::runtime_error("Error: The minimum count threshold exceeds the number of patients.");
  
  // Ensure mutation set has been identified
  getMutationSet(sep);
  
  // Identify rare mutations
  StringSet rareMuts;
  for (auto x: _mutationCounts){
    if (x.second < t){
      rareMuts.insert(x.first);
    }
  }
  
  // Iterate over each patient
  for (int i = 0; i < nrPatients; ++i){
    
    // Find rare mutations for this patient
    StringSet patientRareMuts;
    std::set_intersection(rareMuts.begin(), rareMuts.end(),
                          _patient2mutations[i].begin(), _patient2mutations[i].end(),
                          std::inserter(patientRareMuts, patientRareMuts.begin()));
    
    if (patientRareMuts.size() > 0 ){
    
    // Create mapping from old node labels to new labels, removing rare mutations
    std::unordered_map<std::string, std::string> newLabelMap;
    
    // Iterate over nodes of the tree
    for (NodeIt v(_scriptT[i][0].tree()); v != lemon::INVALID; ++v){
      std::string label = _scriptT[i][0].getIdMap()[v];
      StringVector s;
      boost::split(s, label, boost::is_any_of(sep));
      
      // Find rare mutations at node
      StringVector keep;
      for (auto x: s){
        if (patientRareMuts.find(x) == patientRareMuts.end()){
          keep.push_back(x);
        }
      }
      
      // Only add if at least one mutation still present
      if (keep.size() > 0){
      
        // Get new label 
        std::string joinedString = boost::algorithm::join(keep, sep);
      
        // Add to newLabelMap Dictionary
        newLabelMap[label] = joinedString;
      }
    }
    
    // Create new clone tree vector
    CloneTreeVector newTrees;
    
    // Iterate over each tree for each patient
    for (const CloneTree& T : _scriptT[i]){
        
      // Create new clone tree
      CloneTree newT(T, newLabelMap);
      newTrees.push_back(newT);
    }
    _scriptT[i] = newTrees;
  }
  }
  
  // Update mutation sets
  updateMutationSet(sep);
}

std::ostream& operator<<(std::ostream& out, const InputInstance& scriptT)
{
  const int n = scriptT.getNrPatients();
  out << n << " #patients" << std::endl;
  
  for (int p = 0; p < n; ++p)
  {
    const CloneTreeVector& scriptT_p = scriptT.getPatientTrees(p);
    out << scriptT_p.size() << " #trees for patient " << p << std::endl;
    
    int idx = 0;
    for (const CloneTree& T : scriptT_p)
    {
      const Digraph& TT = T.tree();
      out << lemon::countArcs(TT) << " #edges, tree " << idx++ << std::endl;
      
      for (ArcIt a(TT); a != lemon::INVALID; ++a)
      {
        out << T.label(TT.source(a)) << " " << T.label(TT.target(a)) << std::endl;
      }
    }
  }
  
  return out;
}

std::istream& operator>>(std::istream& in, InputInstance& scriptT)
{
  scriptT._scriptT.clear();
  
  std::string line;
  getline(in, line);
  
  std::stringstream ss(line);
  int n = -1;
  ss >> n;
  
  if (n <= 0)
  {
    throw std::runtime_error("Error: invalid number of patients");
  }
  
  for (int p = 0; p < n; ++p)
  {
    scriptT._scriptT.push_back(CloneTreeVector());
    in >> scriptT._scriptT.back();
  }
  
  return in;
}
