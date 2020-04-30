/*
 * solution.cpp
 *
 *  Created on: 20-oct-2019
 *      Author: M. El-Kebir
 */

#include "solution.h"

Solution::Solution()
  : _consensusTrees()
  , _selectedTrees()
  , _clustering()
  , _selection()
{
}

Solution::Solution(const CloneTreeVector& consensusTrees,
                   const CloneTreeVector& selectedTrees,
                   const IntVector& clustering,
                   const IntVector& selection)
  : _consensusTrees(consensusTrees)
  , _selectedTrees(selectedTrees)
  , _clustering(clustering)
  , _selection(selection)
{
}

std::ostream& operator<<(std::ostream& out, const Solution& sol)
{
  const int k = sol._consensusTrees.size();
  out << k << " # clusters" << std::endl;
  for (int i = 0; i < k; ++i)
  {
    const auto& edges = sol._consensusTrees[i].getEdgeSet();
    out << edges.size() << " # edges of consensus tree " << i << std::endl;
    for (const auto& edge : edges)
    {
      out << edge.first << " " << edge.second << std::endl;
    }
  }
  
  const int n = sol._selectedTrees.size();
  
  out << n << " # patients" << std::endl;
  
  for (int i = 0; i < n; ++i)
  {
    out << sol._clustering[i] << " # cluster index for patient " << i << std::endl;
    out << sol._selection[i] << " # selected input tree index for patient " << i << std::endl;
    const auto& edges = sol._selectedTrees[i].getEdgeSet();
    out << edges.size() << " # edges in selected input tree for patient " << i << std::endl;
    for (const auto& edge : edges)
    {
      out << edge.first << " " << edge.second << std::endl;
    }
  }
  
  return out;
}

std::istream& operator>>(std::istream& in, Solution& sol)
{
  std::stringstream ss;
  std::string line;
  getline(in, line);
  
  StringVector s;
  boost::split(s, line, boost::is_any_of("\t "));
  int k = boost::lexical_cast<int>(s[0]);
  if (k <= 0)
  {
    throw std::runtime_error("Error: invalid number of clusters");
  }
  
  sol._consensusTrees.clear();
  for (int i = 0; i < k; ++i)
  {
    BaseTree T;
    
    getline(in, line);
    boost::split(s, line, boost::is_any_of("\t "));
    int nrEdges = boost::lexical_cast<int>(s[0]);
    if (nrEdges < 0)
    {
      throw std::runtime_error("Error: invalid number of edges");
    }
    
    if (nrEdges > 0)
    {
      ss.clear();
      for (int ii = 0; ii < nrEdges; ++ii)
      {
        getline(in, line);
        ss << line << std::endl;
      }
      
      T.read(ss);
    }
    sol._consensusTrees.push_back(T);
  }
  
  getline(in, line);
  boost::split(s, line, boost::is_any_of("\t "));
  int n = boost::lexical_cast<int>(s[0]);
  if (n <= 0)
  {
    throw std::runtime_error("Error: invalid number of patients");
  }
  
  sol._clustering = IntVector(n, -1);
  sol._selection = IntVector(n, -1);
  
  sol._selectedTrees.clear();
  for (int i = 0; i < n; ++i)
  {
    getline(in, line);
    boost::split(s, line, boost::is_any_of("\t "));
    int j = boost::lexical_cast<int>(s[0]);
    if (!(0 <= j && j < k))
    {
      throw std::runtime_error("Error: invalid cluster index");
    }
    
    sol._clustering[i] = j;
    
    getline(in, line);
    boost::split(s, line, boost::is_any_of("\t "));
    int sel = boost::lexical_cast<int>(s[0]);
//    if (!(0 <= sel))
//    {
//      throw std::runtime_error("Error: tree selection");
//    }

    sol._selection[i] = sel;
    
    BaseTree T;
    getline(in, line);
    boost::split(s, line, boost::is_any_of("\t "));
    int nrEdges = boost::lexical_cast<int>(s[0]);
    if (nrEdges <= 0)
    {
      throw std::runtime_error("Error: invalid number of edges");
    }
    ss.clear();
    for (int ii = 0; ii < nrEdges; ++ii)
    {
      getline(in, line);
      ss << line << std::endl;
    }
    
    T.read(ss);
    sol._selectedTrees.push_back(T);
  }
  
  return in;
}

DoublePair Solution::computeBIC(int nrMutations,
                                bool ancestorDescendantDistance,
                                bool mixtures,
                                bool trim) const
{
  const int n = getNrSelectedTrees();
  const int k = getNrClusters();
    
  if (!mixtures)
  {
    double numerator = 0;
    
    for (int i = 0; i < n; ++i)
    {
      if (ancestorDescendantDistance)
      {
        if (trim)
          numerator += getAncestorDescendantDistance(i, true);
        else
          numerator += getAncestorDescendantDistance(i);
      }
      else
      {
        if (trim)
          numerator += getParentChildDistance(i, true);
        else
          numerator += getParentChildDistance(i);
      }
    }
  
    if (numerator == 0){
      return DoublePair((double)k * log(n), 0);
    }
    double sim = 1 - numerator / (n);
  
    return DoublePair((double)k * log(n), - 2*n*log(sim));
  }
  else
  {
    double sum = 0;
    for (int j = 0; j < k; ++j)
    {
      double numerator = 0;
      int nrTrees_j = 0;
      
      for (int i = 0; i < n; ++i)
      {
        if (_clustering[i] == j)
        {
          ++nrTrees_j;
          if (ancestorDescendantDistance)
          {
            if (trim)
              numerator += getAncestorDescendantDistance(i, true);
            else
              numerator += getAncestorDescendantDistance(i);
          }
          else
          {
            if (trim)
              numerator += getParentChildDistance(i, true);
            else
              numerator += getParentChildDistance(i);
          }
        }
      }
           
      if (nrTrees_j > 0)
      {
        sum += nrTrees_j * log(1 - numerator / (nrTrees_j));
      }
    }
  
    if (sum == 0){
      return DoublePair((double)k * log(n), 0);
    }
    return DoublePair((double)k * log(n), - 2*sum);
  }
}

  
StringPairSet Solution::removeDummyNodes(const StringPairSet& edges){
  
  StringPairSet newEdges;
  for (auto x: edges){
    if (x.first != "MISSING" and x.second != "MISSING" and x.first != "dummy" and x.second != "dummy"){
      newEdges.insert(x);
    }
  }
  
  return newEdges;
}
