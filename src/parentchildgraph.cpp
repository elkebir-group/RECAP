/*
 * inputgraph.cpp
 *
 *  Created on: 6-dec-2018
 *      Author: N. Aguse
 */

#include "parentchildgraph.h"

ParentChildGraph::ParentChildGraph(const CloneTreeVector& ctv,
                                   const StringSet& mutationSet,
                                   const std::string& rootMutation)
  : _nrCloneTrees(ctv.size())
  , _mutationSet(mutationSet)
  , _rootMutation(rootMutation)
  , _G()
  , _arcLookUp(_G)
  , _labelToNode()
  , _nodeToLabel(_G)
  , _nodeToIndex(_G)
  , _indexToNode()
  , _arcWeight(_G)
  , _activeNode(_G)
  , _arborescence(_G)
{
  init(ctv);
}

void ParentChildGraph::init(const CloneTreeVector& ctv)
{
  StringSet presentMutationSet;
  for (const CloneTree& T : ctv)
  {
    for (NodeIt v(T.tree()); v != lemon::INVALID; ++v)
    {
      presentMutationSet.insert(T.label(v));
    }
  }
  
//  StringSet absentMutationSet;
//  std::set_difference(_mutationSet.begin(), _mutationSet.end(),
//                      presentMutationSet.begin(), presentMutationSet.end(),
//                      std::inserter(absentMutationSet, absentMutationSet.begin()));
  
  int nrNodes = 0;
  
  // Construct complete graph _G
  _indexToNode = NodeVector(_mutationSet.size(), lemon::INVALID);
  for (const std::string& mutation : _mutationSet)
//  _indexToNode = NodeVector(presentMutationSet.size(), lemon::INVALID);
//  for (const std::string& mutation : presentMutationSet)
  {
    Node n = _G.addNode();
    _nodeToLabel[n] = mutation;
    _activeNode[n] = presentMutationSet.count(mutation) == 1 ? true : false;
    _labelToNode[mutation] = n;
    _nodeToIndex[n] = nrNodes;
    _indexToNode[nrNodes] = n;
    ++nrNodes;
  }
  assert(nrNodes == _mutationSet.size());
  
  // add missing mutation node
  // Node missingNode = _G.addNode();
  // _nodeToLabel[missingNode] = "MISSING";
  // _activeNode[missingNode] = false;
  // _labelToNode["MISSING"] = missingNode;
  // _nodeToIndex[missingNode] = -1;
    
 
  const Node root = _labelToNode[_rootMutation];
  // _arcWeight[_G.addArc(root, missingNode)] = 0;
  
  for (int i = 0; i < nrNodes; ++i)
  {
    Node v_i = _indexToNode[i];
    for (int j = 0; j < nrNodes; ++j)
    {
      Node v_j = _indexToNode[j];
      if (v_j == root)
        continue;
      if (i == j)
        continue;
      
      if (_activeNode[v_j])
      {
        Arc newarc = _G.addArc(v_i, v_j);
        _arcWeight[newarc] = 0;
      }
      else
      {
        if (ctv.size() != 0)
          throw std::runtime_error("Error: There is a missing node in the PC graph.");
        // Arc newarc = _G.addArc(missingNode, v_j);
        // _arcWeight[newarc] = 0;
      }
    }
  }

  for (const CloneTree& T : ctv)
  {
    const double nrEdgesT = lemon::countArcs(T.tree());
    for (ArcIt a(T.tree()); a != lemon::INVALID; ++a)
    {
      Node src = T.tree().source(a);
      Node tgt = T.tree().target(a);
      const std::string& label_src = T.label(src);
      const std::string& label_tgt = T.label(tgt);
      
      assert(_labelToNode.count(label_src) == 1);
      assert(_labelToNode.count(label_tgt) == 1);
      
      Arc a_exist = _arcLookUp(_labelToNode[label_src], _labelToNode[label_tgt]);
      assert(a_exist != lemon::INVALID);
      _arcWeight[a_exist] -= 1. / (nrNodes + 1 + nrEdgesT);
      // assymmetric difference
//      _arcWeight[a_exist] -= 1. / (nrEdgesT);
    }
  }
  
  double def = 0.;
  for (const CloneTree& T : ctv)
  {
    const double denom = T.getAncestralPairs().size();
    def += (-.5 / denom);
  }
  
  _transArcWeights = DoubleMatrix(nrNodes, DoubleVector(nrNodes, def));
  for (const CloneTree& T : ctv)
  {
    const double denom = T.getAncestralPairs().size();
    for (const StringPair& pair : T.getAncestralPairs())
    {
      assert(_labelToNode.count(pair.first) == 1);
      assert(_labelToNode.count(pair.second) == 1);
      Node u = _labelToNode[pair.first];
      Node v = _labelToNode[pair.second];
      
      int i = _nodeToIndex[u];
      int j = _nodeToIndex[v];
      
      _transArcWeights[i][j] += (1) / denom;
    }
  }
}

void ParentChildGraph::writeDOT(std::ostream &out,
                                int numTrees, double cost) const
{
  out << "digraph T {" << std::endl;
  out << "\tlabel=\"Number of trees: " << numTrees << "\\nNormalized cost: " << cost << "\"" << std::endl;
  for (NodeIt u(_G); u != lemon::INVALID; ++u)
  {
    std::string nl(_nodeToLabel[u]);
    
    /// Only include nodes that have an incoming edge other than dummy (i.e. aren't missing in all trees)
    bool include(false);
    if (nl == "GL"){
      include = true;
    }
    for (ArcIt i(_G); i!=lemon::INVALID; ++i){
      if (_nodeToLabel[_G.source(i)] != "dummy" and _G.target(i)==u and _arcWeight[i]!=0)
        include = true;
    }

    if (nl != "dummy" and include){
      out << "\t" << _G.id(u) << " [label=\"" << nl << "\"]" << std::endl;
    }
    else if (nl != "dummy" and !include){
      out << "\t" << _G.id(u) << " [style=dashed, label=\"" << nl << "\"]" << std::endl;
    }
  }
  
  for (ArcIt a(_G); a != lemon::INVALID; ++a)
  {
    std::pair<std::string, std::string> edgePair = std::make_pair(_nodeToLabel[_G.source(a)], _nodeToLabel[_G.target(a)]);
    double edgeSupport(-2.0*_arcWeight[a]*_mutationSet.size());
    if (edgeSupport > 0 and edgePair.first != "dummy" and edgePair.second != "dummy"){
      out << "\t" << _G.id(_G.source(a)) << " -> " << _G.id(_G.target(a));
      out << " [label=\"" << edgeSupport << "\"";
      //if (_arborescence[a])
      if (_selectedEdgeSet.find(edgePair) != _selectedEdgeSet.end())
      {
        out << ",penwidth=3,color=red";
      }
      out << "]" << std::endl;
    }
    
  }
  
  out << "}" << std::endl;
}

void ParentChildGraph::writeDOT(std::ostream &out) const
{
  out << "digraph T {" << std::endl;
  
  for (NodeIt u(_G); u != lemon::INVALID; ++u)
  {
    out << "\t" << _G.id(u) << " [label=\"" << _nodeToLabel[u] << "\"";
    if (_activeNode[u])
    {
      out << ",color=blue";
    }
    out << "]" << std::endl;
  }
  
  for (ArcIt a(_G); a != lemon::INVALID; ++a)
  {
    out << "\t" << _G.id(_G.source(a)) << " -> " << _G.id(_G.target(a));
    out << " [label=\"" << _arcWeight[a] << "\"";
    if (_arborescence[a])
    {
      out << ",penwidth=3,color=red";
    }
    out << "]" << std::endl;
    
  }
  
  out << "}" << std::endl;
}

double ParentChildGraph::parentChildDistance(const CloneTree & T)
{
  StringPairSet GEdges = getSelectedEdgeList();
  StringPairSet treeEdges = T.getEdgeSet();
  
  StringPairSet resSet;
  std::set_symmetric_difference(GEdges.begin(), GEdges.end(),
                                treeEdges.begin(), treeEdges.end(),
                                std::inserter(resSet, resSet.begin()));
  
//  std::set_difference(treeEdges.begin(), treeEdges.end(),
//                      GEdges.begin(), GEdges.end(),
//                      std::inserter(resSet, resSet.begin()));
  
  
  return resSet.size() / double(GEdges.size() + treeEdges.size());
  
//  return resSet.size() / double(treeEdges.size());
}

double ParentChildGraph::parentChildDistance(const CloneTreeVector& cluster)
{
  double total = 0;
  for (const CloneTree& T : cluster)
  {
    total += parentChildDistance(T);
  }
  
  return total;
}

double ParentChildGraph::ancestorDescendantDistance(const CloneTree& T)
{
  StringPairSet GEdges = getSelectedTransEdgeList();
  StringPairSet treeEdges = T.getAncestralPairs();
  
  StringPairSet resSet;
  std::set_symmetric_difference(GEdges.begin(), GEdges.end(),
                                treeEdges.begin(), treeEdges.end(),
                                std::inserter(resSet, resSet.begin()));
  
  return resSet.size() / double(treeEdges.size());
}

double ParentChildGraph::ancestorDescendantDistance(const CloneTreeVector& cluster)
{
  int total = 0;
  for (const CloneTree& T : cluster)
  {
    total += ancestorDescendantDistance(T);
  }
  
  return total;
}

void ParentChildGraph::maxWeightSpanningArborescence()
{
//  if (lemon::countNodes(_G) == 0)
//  {
//    _selectedEdgeSet.clear();
//    _selectedTransEdgeSet.clear();
//  }
//  else
  {
    assert(_labelToNode.count(_rootMutation) == 1);
    Node root = _labelToNode[_rootMutation];
    
    lemon::mapFill(_G, _arborescence, false);
//    {
//      static int idx = 0;
//      char buf[1024];
//      snprintf(buf, 1024, "/tmp/G%d.dot", idx);
//      std::ofstream outF(buf);
//      writeDOT(outF);
//      outF.close();
//      ++idx;
//    }
    double arborescence_cost = lemon::minCostArborescence(_G, _arcWeight, root, _arborescence);
    
    _selectedEdgeSet.clear();
    _selectedTransEdgeSet.clear();
    for (ArcIt a_ij(_G); a_ij != lemon::INVALID; ++a_ij)
    {
      Node v_i = _G.source(a_ij);
      Node v_j = _G.target(a_ij);

      if (_arborescence[a_ij])
      {
        _selectedEdgeSet.insert(StringPair(_nodeToLabel[v_i], _nodeToLabel[v_j]));
      }
    }
  }
}
