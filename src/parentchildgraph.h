/*
 * parentchildgraph.h
 *
 *  Created on: 6-dec-2018
 *      Author: N. Aguse
 */

#ifndef PARENTCHILDGRAPH_H
#define PARENTCHILDGRAPH_H

#include "clonetree.h"
#include <lemon/min_cost_arborescence.h>
#include <lemon/maps.h>
#include <algorithm>

class ParentChildGraph
{
public:
  /// Constructor
  ///
  /// @param ctv Clone tree vector
  ParentChildGraph(const CloneTreeVector& ctv,
                   const StringSet& mutationSet,
                   const std::string& rootMutation);
  
  /// Clear selected edges and transitive edges
  void clear()
  {
    lemon::mapFill(_G, _arborescence, false);
    _selectedEdgeSet.clear();
    _selectedTransEdgeSet.clear();
  }
  
  /// Write parent-child graph
  ///
  /// @param out Output stream
  void writeDOT(std::ostream & out) const;
  
  /// Write parent-child graph
  ///
  /// @param numTrees Number of trees
  /// @param cost Total distance
  void writeDOT(std::ostream & out,
                int numTrees,
                double cost) const;
  
  /// Compute maximum weight spanning arborescence
  void maxWeightSpanningArborescence();
  
  /// Compute parent-child distance to consensus tree of given clone tree
  ///
  /// @param T Clone tree
  double parentChildDistance(const CloneTree& T);
  
  /// Compute ancestor-descendant distance to consensus tree of given clone tree
  ///
  /// @param T Clone tree
  double ancestorDescendantDistance(const CloneTree& T);
  
  /// Compute sum of parent-child distances to consensus tree of given clone trees
  ///
  /// @param cluster Clone trees
  double parentChildDistance(const CloneTreeVector& cluster);
  
  /// Compute sum of ancestor-descendant distances to consensus tree of given clone trees
  ///
  /// @param cluster Clone trees
  double ancestorDescendantDistance(const CloneTreeVector& cluster);
  
  /// Return the selected edge set
  const StringPairSet& getSelectedEdgeList() const
  {
    return _selectedEdgeSet;
  }
  
  /// Return the selected trans edge set
  const StringPairSet& getSelectedTransEdgeList() const
  {
    return _selectedTransEdgeSet;
  }
  
  /// Set selected edge set
  ///
  /// @param selectedEdgeSet Edge set
  void setSelectedEdgeList(const StringPairSet& selectedEdgeSet)
  {
    _selectedEdgeSet = selectedEdgeSet;
  }
  
  /// Set selected edge set (_arborescence and _selectEdgeSet are not updated)
  ///
  /// @param selectedTransEdgeSet Transitive edge set
  void setSelectedTransEdgeList(const StringPairSet& selectedTransEdgeSet)
  {
    _selectedTransEdgeSet = selectedTransEdgeSet;
  }
  
  bool isActiveNode(const std::string& mutation) const
  {
    assert(_mutationSet.count(mutation) == 1);
    Node v = _labelToNode.find(mutation)->second;
    return _activeNode[v];
  }
  
  bool isActiveArc(const std::string& mutSource,
                   const std::string& mutTarget) const
  {
    assert(_mutationSet.count(mutSource) == 1);
    assert(_mutationSet.count(mutTarget) == 1);
    
    Node v = _labelToNode.find(mutSource)->second;
    Node w = _labelToNode.find(mutTarget)->second;

    return _activeNode[v] && _activeNode[w];
  }
  
  double getArcOccurrence(const std::string& mutSource,
                          const std::string& mutTarget) const
  {
//    assert(_mutationSet.count(mutSource) == 1);
//    assert(_mutationSet.count(mutTarget) == 1);
    
    Node v = _labelToNode.find(mutSource)->second;
    Node w = _labelToNode.find(mutTarget)->second;
    
    Arc a = _arcLookUp(v, w);
    if (a == lemon::INVALID)
    {
      return 0;
    }
    else
    {
      assert(a != lemon::INVALID);
      return -_arcWeight[a];
    }
  }

private:
  /// Number of input clone trees
  const int _nrCloneTrees;
  /// Mutation set
  const StringSet& _mutationSet;
  /// Root mutation
  const std::string& _rootMutation;
  /// Parent-child graph
  Digraph _G;
  /// Arc look up
  DynArcLookUp _arcLookUp;
  /// String to node map
  StringToNodeMap _labelToNode;
  /// Node to string map
  StringNodeMap _nodeToLabel;
  /// Node to index map
  IntNodeMap _nodeToIndex;
  /// Index to node map
  NodeVector _indexToNode;
  /// Arc weights
  DoubleArcMap _arcWeight;
  /// Trans arc weights
  DoubleMatrix _transArcWeights;
  /// Indicator of real/fake node
  BoolNodeMap _activeNode;
  /// Arborescence indicator
  Digraph::ArcMap<bool> _arborescence;
  /// Edges of arborescence
  StringPairSet _selectedEdgeSet;
  /// Transitive closure of edges of arborescence
  StringPairSet _selectedTransEdgeSet;
  
  void init(const CloneTreeVector& ctv);
};

#endif // PARENTCHILDGRAPH_H
