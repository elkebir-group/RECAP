/*
 * solution.h
 *
 *  Created on: 20-oct-2019
 *      Author: M. El-Kebir
 */

#ifndef SOLUTION_H
#define SOLUTION_H

#include "clonetree.h"
#include "utils.h"

class Solution
{
public:
  Solution();
  
  Solution(const CloneTreeVector& consensusTrees,
           const CloneTreeVector& selectedTrees,
           const IntVector& clustering,
           const IntVector& selection);
  
  const CloneTreeVector& getConsensusTrees() const
  {
    return _consensusTrees;
  }
  
  const CloneTreeVector& getSelectedTrees() const
  {
    return _selectedTrees;
  }
  
  const IntVector& getClustering() const
  {
    return _clustering;
  }
  
  const IntVector& getSelection() const
  {
    return _selection;
  }
  
  int getNrSelectedTrees() const
  {
    return _selection.size();
  }
  
  int getNrClusters() const
  {
    return _consensusTrees.size();
  }
  
  const CloneTree& getConsensusTree(int i) const
  {
    assert(0 <= i && i < _consensusTrees.size());
    return _consensusTrees[i];
  }
  
  const CloneTree& getSelectedTree(int i) const
  {
    assert(0 <= i && i < _selectedTrees.size());
    return _selectedTrees[i];
  }
  
  double getParentChildDistance(int i, const bool trim=false) const
  {
    assert(0 <= i && i < getNrSelectedTrees());
    
    int res;
    
    StringPairSet consensusEdges = _consensusTrees[_clustering[i]].getParentalPairs();
    StringPairSet selectedEdges = _selectedTrees[i].getParentalPairs();
    
    // REMOVING DUMMY NODES
    if (trim){
      consensusEdges = Solution::removeDummyNodes(consensusEdges);
      selectedEdges = Solution::removeDummyNodes(selectedEdges);
      if (consensusEdges.size() + selectedEdges.size() == 0){
        return 0;
      }
    }
    
    StringPairSet resSet;
    std::set_symmetric_difference(consensusEdges.begin(), consensusEdges.end(),
                                  selectedEdges.begin(), selectedEdges.end(),
                                  std::inserter(resSet, resSet.begin()));
    
    res = resSet.size();
    
    return res/((double)(consensusEdges.size() + selectedEdges.size()));
  }
  
  int getAncestorDescendantDistance(int i, const bool trim=false) const
  {
    assert(0 <= i && i < getNrSelectedTrees());
    
    int res;
    
    StringPairSet consensusEdges = _consensusTrees[_clustering[i]].getAncestralPairs();
    StringPairSet selectedEdges = _selectedTrees[i].getAncestralPairs();
    
    
    // REMOVING DUMMY NODES
    if (trim){
      consensusEdges = Solution::removeDummyNodes(consensusEdges);
      selectedEdges = Solution::removeDummyNodes(selectedEdges);
      if (consensusEdges.size() + selectedEdges.size() == 0){
        return 0;
      }
    }
    
    StringPairSet resSet;
    std::set_symmetric_difference(consensusEdges.begin(), consensusEdges.end(),
                                  selectedEdges.begin(), selectedEdges.end(),
                                  std::inserter(resSet, resSet.begin()));
    
    res = resSet.size();
    
    return res/((double)(consensusEdges.size() + selectedEdges.size()));
  }
  
  double getParentChildDistance(const bool trim=false) const
  {
    double res = 0;
    
    const int n = getNrSelectedTrees();
    for (int i = 0; i < n; ++i)
    {
      res += getParentChildDistance(i, trim);
    }
    
    return res;
  }
  
  int getAncestorDescendantDistance(const bool trim=false) const
  {
    int res = 0;
    
    const int n = getNrSelectedTrees();
    for (int i = 0; i < n; ++i)
    {
      res += getAncestorDescendantDistance(i, trim);;
    }
    
    return res;
  }
  
  DoublePair computeBIC(int nrMutations,
                        bool ancestorDescendantDistance,
                        bool mixtures,
                        bool trim = false) const;
    
  static StringPairSet removeDummyNodes(const StringPairSet& edges);
  
private:
  CloneTreeVector _consensusTrees;
  
  CloneTreeVector _selectedTrees;
  
  IntVector _clustering;
  
  IntVector _selection;
  
  friend std::ostream& operator<<(std::ostream& out, const Solution& sol);
  friend std::istream& operator>>(std::istream& in, Solution& sol);
};

std::ostream& operator<<(std::ostream& out, const Solution& sol);

std::istream& operator>>(std::istream& in, Solution& sol);

#endif // SOLUTION_H
