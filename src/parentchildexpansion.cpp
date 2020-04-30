/*
 * parentchildexpansion.cpp
 *
 *  Created on: 28-oct-2019
 *      Author: M. El-Kebir
 */

#include "parentchildexpansion.h"

ParentChildExpansion::ParentChildExpansion(const CloneTree& clusteredCloneTree,
                                           const StringPairSet& expandedParentChildPairs,
                                           const std::string& sep,
                                           LookUpTableType& lookUpTable)
  : _clusteredCloneTree(clusteredCloneTree)
  , _expandedParentChildPairs(expandedParentChildPairs)
  , _sep(sep)
  , _lookUpTable(lookUpTable)
  , _table(_clusteredCloneTree.tree())
  , _expTable(_clusteredCloneTree.tree())
  , _matchedEdgesTable(_clusteredCloneTree.tree())
{
}

CloneTree ParentChildExpansion::getExpandedCloneTree()
{
  const Node root = _clusteredCloneTree.root();
  fillTable(root);
  
//  Node w = _clusteredCloneTree.tree().target(OutArcIt(_clusteredCloneTree.tree(), root));
  
  // construct expansionMap
  std::unordered_map<std::string, StringVector> expansionMap;
  
  int max_sim = std::numeric_limits<int>::min();
  std::vector<StringPair> startStopPairList;
  for (const auto& kv : _table[root])
  {
    if (kv.second > max_sim)
    {
      startStopPairList.clear();
      max_sim = kv.second;
      startStopPairList.push_back(kv.first);
    }
    else if (kv.second == max_sim)
    {
      startStopPairList.push_back(kv.first);
    }
  }
  
  std::shuffle(startStopPairList.begin(), startStopPairList.end(), g_rng);
  expansionMap[_clusteredCloneTree.label(root)] = _expTable[root][startStopPairList.front()];
  
  backTrace(root, startStopPairList.front().second, expansionMap);
  CloneTree res(_clusteredCloneTree, expansionMap);
  
#ifdef DEBUG
  StringPairSet GEdges = _expandedParentChildPairs;
  StringPairSet treeEdges = res.getEdgeSet();
  
  StringPairSet resSet;
  std::set_symmetric_difference(GEdges.begin(), GEdges.end(),
                                treeEdges.begin(), treeEdges.end(),
                                std::inserter(resSet, resSet.begin()));
  
  assert(resSet.size() == (GEdges.size() + treeEdges.size() - 2*max_sim));
//  {
//    std::cout << std::endl;
//  }
#endif
  
  return res;
}

void ParentChildExpansion::backTrace(Node v,
                                     const std::string& stop,
                                     std::unordered_map<std::string, StringVector>& expansionMap) const
{
  const Digraph& T = _clusteredCloneTree.tree();
  
  for (OutArcIt a(T, v); a != lemon::INVALID; ++a)
  {
    Node w = T.target(a);
    
    int max_sim = std::numeric_limits<int>::min();
    std::vector<StringPair> startStopPairList;
    for (const auto& kv : _table[w])
    {
      int sim = kv.second;
      if (_expandedParentChildPairs.count(StringPair(stop, kv.first.first)) == 1)
      {
        sim += 1;
      }
          
      if (sim > max_sim)
      {
        startStopPairList.clear();
        max_sim = sim;
        startStopPairList.push_back(kv.first);
      }
      else if (sim == max_sim)
      {
        startStopPairList.push_back(kv.first);
      }
    }
    
    std::shuffle(startStopPairList.begin(), startStopPairList.end(), g_rng);
    expansionMap[_clusteredCloneTree.label(w)] = _expTable[w].find(startStopPairList.front())->second;
    backTrace(w, startStopPairList.front().second, expansionMap);
  }
}

void ParentChildExpansion::fillTable(Node v)
{
  const Digraph& T = _clusteredCloneTree.tree();
  
  // Post-order traversal
  for (OutArcIt a(T, v); a != lemon::INVALID; ++a)
  {
    Node w = T.target(a);
    fillTable(w);
  }
  
  StringVector s;
  boost::split(s, _clusteredCloneTree.label(v), boost::is_any_of(_sep));
  StringSet mutationsInCluster(s.begin(), s.end());
  assert(s.size() >= 1);
  
  if (s.size() == 1)
  {
    std::string start = s.front();
    std::string stop = s.front();
    
    int total_sim = 0;
    _expTable[v][StringPair(start, stop)] = StringVector(1, start);
    for (OutArcIt a(T, v); a != lemon::INVALID; ++a)
    {
      Node w = T.target(a);
      
      // find minimum distance for w
      int max_sim_w = std::numeric_limits<int>::min();
      StringPair startStop_w;
      for (const auto& kv : _table[w])
      {
        // set sim to maximum number of common PC pairs of subtree rooted w with expansion (kv)
        int sim = kv.second;
        // need to do a +1 if (stop, kv.first.first) occurs in consensus
        if (_expandedParentChildPairs.count(StringPair(stop, kv.first.first)) == 1)
        {
          sim += 1;
        }
        // find maximum
        if (sim > max_sim_w)
        {
          max_sim_w = sim;
        }
      }
      total_sim += max_sim_w;
    }
    _table[v][StringPair(start, stop)] = total_sim;
  }
  else
  {
    // make all pairs
    for (const std::string& start : s)
    {
      for (const std::string& stop : s)
      {
        if (start == stop) continue;
        
        // next do +whatever we have in best expansion
        StringVector expansion;
        int total_sim = getMaxSimExpansion(mutationsInCluster, start, stop, expansion, _matchedEdgesTable[v][StringPair(start, stop)]);
        _expTable[v][StringPair(start, stop)] = expansion;
        
        for (OutArcIt a(T, v); a != lemon::INVALID; ++a)
        {
          Node w = T.target(a);
          
          // find minimum distance for w
          int max_sim_w = std::numeric_limits<int>::min();
          for (const auto& kv : _table[w])
          {
            // set sim to maximum number of common PC pairs of subtree rooted w with expansion (kv)
            int sim = kv.second;
            // need to do a +1 if (stop, kv.first.first) occurs in consensus
            if (_expandedParentChildPairs.count(StringPair(stop, kv.first.first)) == 1)
            {
              sim += 1;
            }
            
            // find maximum
            if (sim > max_sim_w)
            {
              max_sim_w = sim;
            }
          }
          total_sim += max_sim_w;
        }
        
        _table[v][StringPair(start, stop)] = total_sim;
      }
    }
  }
}

int ParentChildExpansion::getMaxSimExpansion(const StringSet& mutationsInCluster,
                                             const std::string& start,
                                             const std::string& stop,
                                             StringVector& expansion,
                                             std::vector<StringPair>& matchedEdgesInConsensus) const
{
//  StringPair startStop(start, stop);
  std::string startStop = start;
  startStop += stop;
  if (_lookUpTable.count(startStop) == 1)
  {
    expansion = _lookUpTable[startStop].first;
    matchedEdgesInConsensus = _lookUpTable[startStop].second;
    return matchedEdgesInConsensus.size();
  }

  matchedEdgesInConsensus.clear();
  for (const StringPair& consensusEdge : _expandedParentChildPairs)
  {
    if ((mutationsInCluster.count(consensusEdge.first) == 1) && (mutationsInCluster.count(consensusEdge.second) == 1))
    {
      if (consensusEdge.first != stop && consensusEdge.second != start)
      {
        if (!(mutationsInCluster.size() > 2 && consensusEdge.first == start && consensusEdge.second == stop))
        {
          matchedEdgesInConsensus.push_back(consensusEdge);
        }
      }
    }
  }
  
  // generate subgraph using selected edges (those in matchedEdgesInConsensus)
  // if there are >= 2 components and one of them contains a path from start to stop then remove one edge on this path
  // while there is a node with more than two children, remove one edge
  // first identify roots
  StringSet matchingMutations;
  for (const StringPair& matchingEdge : matchedEdgesInConsensus)
  {
    matchingMutations.insert(matchingEdge.first);
    matchingMutations.insert(matchingEdge.second);
  }
  
  std::unordered_map<std::string, std::map<std::string, StringVector>> rootToMatchingEdges;
  std::unordered_map<std::string, std::map<std::string, StringSet>> rootToMatchingDescendants;
  std::unordered_map<std::string, StringSet> rootToMatchingMutations;
  StringSet rootMutations;
  
  getComponents(mutationsInCluster,
                matchedEdgesInConsensus,
                rootToMatchingEdges,
                rootToMatchingDescendants,
                rootToMatchingMutations,
                rootMutations);
  
  // if there are more than two components, remove an edge on the path from start to stop if present in one component
  for (const std::string& rootMut : rootMutations)
  {
    if (rootToMatchingMutations[rootMut].count(start) == 1
        && rootToMatchingMutations[rootMut].count(stop) == 1)
    {
      // we have such a component, reconstruct path
      std::vector<StringPair> path;
      std::string startMut = start;
      while (startMut != stop)
      {
        for (const StringPair& edge : matchedEdgesInConsensus)
        {
          if (edge.first == startMut && rootToMatchingDescendants[rootMut][edge.second].count(stop) == 1)
          {
            path.push_back(edge);
            startMut = edge.second;
            break;
          }
        }
      }
      if (path.size() != mutationsInCluster.size() - 1)
      {
        // score edges in path by number of chains they induce upon removal
        // since our goal is to maximize number of matched edges, we are only interested in edges
        // whose removal induces largest number of chains
        StringPair edgeToRemove;
        
        IntVector desiredRemovalIndices;
        for (int i = 0; i < path.size(); ++i)
        {
          if (rootToMatchingEdges[rootMut][path[i].first].size() >= 2)
          {
            desiredRemovalIndices.push_back(i);
          }
        }
        if (!desiredRemovalIndices.empty())
        {
          std::shuffle(desiredRemovalIndices.begin(), desiredRemovalIndices.end(), g_rng);
          edgeToRemove = path[desiredRemovalIndices.front()];
        }
        else
        {
          // shuffle path
          std::shuffle(path.begin(), path.end(), g_rng);
          edgeToRemove = path[0];
        }
          
        // remove edge
        auto it = std::find(matchedEdgesInConsensus.begin(), matchedEdgesInConsensus.end(), edgeToRemove);
        matchedEdgesInConsensus.erase(it);
        
        auto it2 = std::find(rootToMatchingEdges[rootMut][edgeToRemove.first].begin(),
                             rootToMatchingEdges[rootMut][edgeToRemove.first].end(),
                             edgeToRemove.second);
        rootToMatchingEdges[rootMut][edgeToRemove.first].erase(it2);
        
        if (rootToMatchingEdges[rootMut][edgeToRemove.second].size() == 0)
        {
          rootToMatchingMutations[rootMut].erase(edgeToRemove.second);
        }
        break;
      }
      else
      {
        path.clear();
      }
    }
  }
  
  // update our components again
  rootToMatchingEdges.clear();
  rootToMatchingMutations.clear();
  rootToMatchingDescendants.clear();
  rootMutations.clear();
  getComponents(mutationsInCluster,
                matchedEdgesInConsensus,
                rootToMatchingEdges,
                rootToMatchingDescendants,
                rootToMatchingMutations,
                rootMutations);
  
  // now we need to greedily match edges, a parent can only be used once!
  // that is why we keep track of available parents
  // initially all parents are available except stop
  std::shuffle(matchedEdgesInConsensus.begin(), matchedEdgesInConsensus.end(), g_rng);
  
  StringSet availableParents(mutationsInCluster);
  availableParents.erase(stop);
  
  // at the same time we identify mutations in v that are left unmatched in consensus tree
  StringSet missingMutationSet(mutationsInCluster);
  missingMutationSet.erase(start);
  missingMutationSet.erase(stop);
  
  for (auto it = matchedEdgesInConsensus.begin(); it != matchedEdgesInConsensus.end();)
  {
    // match edge if we can (parent not matched previously)
    if (availableParents.count(it->first) == 1)
    {
      missingMutationSet.erase(it->first);
      missingMutationSet.erase(it->second);
      availableParents.erase(it->first);
      ++it;
    }
    else
    {
      // if we can't, erase edge
      it = matchedEdgesInConsensus.erase(it);
    }
  }
  
  // update our components again, because we just deleted some stuff
  rootToMatchingEdges.clear();
  rootToMatchingMutations.clear();
  rootToMatchingDescendants.clear();
  rootMutations.clear();
  getComponents(mutationsInCluster,
                matchedEdgesInConsensus,
                rootToMatchingEdges,
                rootToMatchingDescendants,
                rootToMatchingMutations,
                rootMutations);
  
  // each component is now a chain of matched parent child edges
  std::vector<StringList> matchedChains;
  for (const std::string& rootMut : rootMutations)
  {
    matchedChains.push_back(StringList());
    std::string startMut = rootMut;
    while (true)
    {
      matchedChains.back().push_back(startMut);
      if (!rootToMatchingEdges[rootMut][startMut].empty())
      {
        assert(rootToMatchingEdges[rootMut][startMut].size() == 1);
        startMut = rootToMatchingEdges[rootMut][startMut][0];
      }
      else
      {
        break;
      }
    }
    if (matchedChains.back().size() == 1)
    {
      matchedChains.pop_back();
    }
  }
  
  // extract the two chains containing start and stop (if any)
  StringList matchedChainStart;
  StringList matchedChainStop;
  
  for (auto it = matchedChains.begin(); it != matchedChains.end();)
  {
    if (it->front() == start)
    {
      matchedChainStart = *it;
      it = matchedChains.erase(it);
    }
    else if (it->back() == stop)
    {
      matchedChainStop = *it;
      it = matchedChains.erase(it);
    }
    else
    {
      ++it;
    }
  }
  
  if (matchedChainStart.empty() &&
      std::find(matchedChainStop.begin(), matchedChainStop.end(), start) == matchedChainStop.end())
  {
    matchedChainStart.push_back(start);
  }
  if (matchedChainStop.empty() &&
      std::find(matchedChainStart.begin(), matchedChainStart.end(), stop) == matchedChainStart.end())
  {
    matchedChainStop.push_back(stop);
  }
  
  // now shuffle matchedChains
  std::shuffle(matchedChains.begin(), matchedChains.end(), g_rng);
  
  // shuffle missing mutations
  StringVector missingMutations(missingMutationSet.begin(),
                                missingMutationSet.end());
  std::shuffle(missingMutations.begin(), missingMutations.end(), g_rng);
  
  // determine where to place missing mutations:
  // they can be:   (i) in front of matchedChains but after matchedChainStart (pos = 0)
  //               (ii) at the end of matchedChains but prior to matchedChainStop (pos = |matchedChains|)
  //              (iii) in between matchedChains (0 < pos < |matchedChains|)
  std::uniform_int_distribution<> u_int_dist(0, matchedChains.size());
  IntVector insertionPositions(missingMutations.size(), -1);
  for (int idx = 0; idx < missingMutations.size(); ++idx)
  {
    insertionPositions[idx] = u_int_dist(g_rng);
  }
  
  for (const std::string& mut : matchedChainStart)
  {
    expansion.push_back(mut);
  }
  
  for (int pos = 0; pos <= matchedChains.size(); ++pos)
  {
    for (int idx = 0; idx < missingMutations.size(); ++idx)
    {
      if (insertionPositions[idx] == pos)
      {
        expansion.push_back(missingMutations[idx]);
      }
    }
    // insert matchedChains[pos]
    if (pos < matchedChains.size())
    {
      for (const std::string& mut : matchedChains[pos])
      {
        expansion.push_back(mut);
      }
    }
    else
    {
      for (const std::string& mut : matchedChainStop)
      {
        expansion.push_back(mut);
      }
    }
  }
  
#ifdef DEBUG
  StringSet testSet(expansion.begin(), expansion.end());
  assert(testSet.size() == expansion.size());
#endif
  assert(expansion.front() == start);
  assert(expansion.back() == stop);
  assert(expansion.size() == mutationsInCluster.size());
  
  _lookUpTable[startStop] = std::make_pair(expansion, matchedEdgesInConsensus);
  return matchedEdgesInConsensus.size();
}
