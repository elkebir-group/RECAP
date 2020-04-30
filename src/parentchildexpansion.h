/*
 * parentchildexpansion.h
 *
 *  Created on: 28-oct-2019
 *      Author: M. El-Kebir
 */

#ifndef PARENTCHILDEXPANSION_H
#define PARENTCHILDEXPANSION_H

#include "utils.h"
#include "clonetree.h"
#include <unordered_map>

class ParentChildExpansion
{
public:
  typedef std::unordered_map<std::string, std::pair<StringVector, std::vector<StringPair>>> LookUpTableType;
  
  /// Constructor
  ///
  /// @param clusteredCloneTree Clustered clone tree
  /// @param expandedCloneTree Expanded clone tree
  /// @param sep Mutation cluster separator
  /// @param lookUpTable Look-up table to identify expansion and matchedEdges corresponding to <start,stop> pairs (against a fixed consensus)
  ParentChildExpansion(const CloneTree& clusteredCloneTree,
                       const StringPairSet& expandedParentChildPairs,
                       const std::string& sep,
                       LookUpTableType& lookUpTable);
  
  CloneTree getExpandedCloneTree();
  
  void backTrace(Node v,
                 const std::string& stop,
                 std::unordered_map<std::string, StringVector>& expansionMap) const;
  
private:
  typedef std::map<StringPair, int> StringPairToIntMap;
  typedef std::map<StringPair, StringVector> StringPairToStringVectorMap;
  typedef std::map<StringPair, std::vector<StringPair>> StringPairToStringPairVectorMap;
  typedef Digraph::NodeMap<StringPairToIntMap> DpTableType;
  typedef Digraph::NodeMap<StringPairToStringVectorMap> DpTableExpType;
  typedef Digraph::NodeMap<StringPairToStringPairVectorMap> DpTableEdgesType;
  
  /// Clustered clone tree
  const CloneTree& _clusteredCloneTree;
  /// Expanded clone tree
  const StringPairSet& _expandedParentChildPairs;
  /// Mutation cluster separator
  const std::string& _sep;
  /// Look-up table to identify expansion and matchedEdges corresponding to <start,stop> pairs (against a fixed consensus)
  LookUpTableType& _lookUpTable;
  /// Dynamic programming table: maximum number of matching parent-child pairs in subtree rooted at node v when assigning s as start and t as end
  DpTableType _table;
  /// Dynamic programming table: expansion with maximum number of matching parent-child pairs in subtree rooted at node v when assigning s as start and t as end
  DpTableExpType _expTable;
  /// Matched edges
  DpTableEdgesType _matchedEdgesTable;
  
  void fillTable(Node v);
  
  int getMaxSimExpansion(const StringSet& mutationsInCluster,
                         const std::string& start,
                         const std::string& stop,
                         StringVector& expansion,
                         std::vector<StringPair>& matchesEdges) const;
  
  void fillMatchingDescendants(const std::string& rootMut,
                               std::unordered_map<std::string, std::map<std::string, StringVector>>& rootToMatchingEdges,
                               const std::string& mut,
                               std::unordered_map<std::string, std::map<std::string, StringSet>>& rootToMatchingDescendants) const
  {
    assert(rootToMatchingEdges.count(rootMut) == 1);
    assert(rootToMatchingEdges[rootMut].count(mut) == 1);
    const StringVector& children = rootToMatchingEdges[rootMut][mut];
    
    rootToMatchingDescendants[rootMut][mut].clear();
    rootToMatchingDescendants[rootMut][mut].insert(mut);
    for (const std::string& childMut : children)
    {
      fillMatchingDescendants(rootMut, rootToMatchingEdges, childMut, rootToMatchingDescendants);
      rootToMatchingDescendants[rootMut][mut].insert(rootToMatchingDescendants[rootMut][childMut].begin(),
                                                     rootToMatchingDescendants[rootMut][childMut].end());
    }
  }
  
  void getComponents(const StringSet& mutationsInCluster,
                     const std::vector<StringPair>& matchedEdgesInConsensus,
                     std::unordered_map<std::string, std::map<std::string, StringVector>>& rootToMatchingEdges,
                     std::unordered_map<std::string, std::map<std::string, StringSet>>& rootToMatchingDescendants,
                     std::unordered_map<std::string, StringSet>& rootToMatchingMutations,
                     StringSet& rootMutations) const
  {
    StringSet matchingMutations;
    for (const StringPair& matchingEdge : matchedEdgesInConsensus)
    {
      matchingMutations.insert(matchingEdge.first);
      matchingMutations.insert(matchingEdge.second);
    }
    
    for (const std::string& mut : matchingMutations)
    {
      bool root = true;
      for (const StringPair& matchingEdge : matchedEdgesInConsensus)
      {
        if (matchingEdge.second == mut)
        {
          root = false;
          break;
        }
      }
      if (root)
      {
        rootToMatchingEdges[mut][mut];
        rootToMatchingMutations[mut].insert(mut);
        rootMutations.insert(mut);
      }
    }
    
    // partition matching edges into components
    StringPairSet matchedEdgesInConsensusSet(matchedEdgesInConsensus.begin(), matchedEdgesInConsensus.end());
    while (!matchedEdgesInConsensusSet.empty())
    {
      for (auto it = matchedEdgesInConsensusSet.begin(); it != matchedEdgesInConsensusSet.end();)
      {
        bool reassigned = false;
        for (const std::string& rootMut : rootMutations)
        {
          if (rootToMatchingEdges[rootMut].count(it->first) == 1)
          {
            rootToMatchingEdges[rootMut][it->first].push_back(it->second);
            rootToMatchingEdges[rootMut][it->second];
            rootToMatchingMutations[rootMut].insert(it->second);
            it = matchedEdgesInConsensusSet.erase(it);
            reassigned = true;
            break;
          }
        }
        if (!reassigned)
        {
          ++it;
        }
      }
    }
    
    // put unmatched nodes as separate components
    for (const std::string& mut : mutationsInCluster)
    {
      if (rootToMatchingEdges.count(mut) == 0)
      {
        bool found = false;
        for (const std::string& rootMut : rootMutations)
        {
          if (rootToMatchingMutations[rootMut].count(mut) == 1)
          {
            found = true;
            break;
          }
        }
        if (!found)
        {
          rootToMatchingMutations[mut].insert(mut);
          rootToMatchingEdges[mut][mut];
          rootMutations.insert(mut);
        }
      }
    }
    for (const std::string& rootMut : rootMutations)
    {
      fillMatchingDescendants(rootMut, rootToMatchingEdges, rootMut, rootToMatchingDescendants);
    }
  }
};

#endif // PARENTCHILDEXPANSION_H
