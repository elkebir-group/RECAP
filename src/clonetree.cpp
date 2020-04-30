/*
 * clonetree.cpp
 *
 *  Created on: 19-oct-2016
 *      Author: M. El-Kebir
 */

#include "clonetree.h"
#include <iomanip>

CloneTree::CloneTree()
  : BaseTree()
  , _l(_tree)
{
}

CloneTree::CloneTree(const CloneTree& other)
  : BaseTree(other)
  , _l(_tree)
{
  for (NodeIt u(_tree); u != lemon::INVALID; ++u)
  {
    const std::string& lbl = label(u);
    Node other_u = other.getNodeByLabel(lbl);
    
    _l[u] = other._l[other_u];
  }
}

CloneTree& CloneTree::operator =(const CloneTree& other)
{
  if (this != &other)
  {
    lemon::digraphCopy(other._tree, _tree)
      .node(other._root, _root)
      .nodeMap(other._nodeToId, _nodeToId)
      .run();
    
    for (NodeIt u(_tree); u != lemon::INVALID; ++u)
    {
      const std::string& str = _nodeToId[u];
      _idToNode[str] = u;
    }
    
    _edgeSet = other._edgeSet;
    _transEdgeSet = other._transEdgeSet;
    
    init();
    
    for (NodeIt u(_tree); u != lemon::INVALID; ++u)
    {
      const std::string& lbl = label(u);
      Node other_u = other.getNodeByLabel(lbl);
      
      _l[u] = other._l[other_u];
    }
  }
  
  return *this;
}

CloneTree::CloneTree(const Digraph& T,
                     Node root,
                     const StringNodeMap& label,
                     const StringNodeMap& l)
  : BaseTree(T, root, label)
  , _l(_tree)
{
  for (NodeIt v(T); v != lemon::INVALID; ++v)
  {
    if (OutArcIt(T, v) == lemon::INVALID)
    {
      const std::string& label_v = label[v];
      _l[getNodeByLabel(label_v)] = l[v];
    }
  }
}

CloneTree::CloneTree(const CloneTree& clusterTree,
                     const std::unordered_map<std::string, StringVector>& expansionMap)
  : BaseTree()
  , _l(_tree)
{
    // Create a set to keep track of nodes we have added
    std::unordered_set<std::string> addedLabels;
    
    // Iterate over nodes in clusterTree
    for (NodeIt n(clusterTree.tree()); n != lemon::INVALID; ++n){
        
        // Get label 
        std::string clusterLabel = clusterTree.label(n);
        
        // First create a new node for each item in cluster node expansion
        StringVector expansionV = expansionMap.at(clusterLabel);
        for (auto & newLabel : expansionV){
            if (addedLabels.count(newLabel)==0){
                Digraph::Node x = _tree.addNode();
                _nodeToId[x] = newLabel;
                _idToNode[newLabel] = x;
                addedLabels.insert(newLabel);
            }
        }
        
        // Now add an edge between nodes in this expansion
        Digraph::Node prevNode(_idToNode[expansionV[0]]);
        for(unsigned int i = 1; i != expansionV.size(); i++) {
            Digraph::Node nextNode(_idToNode[expansionV[i]]);
            _tree.addArc(prevNode, nextNode);
            prevNode = nextNode;
        }
        
        // Now add an edge to first element of children clusters
        for (Digraph::OutArcIt a(clusterTree.tree(), n); a != lemon::INVALID; ++a){
            
            // Get this child cluster label
            std::string clusterLabel = clusterTree.label(clusterTree.tree().target(a));
            
            // Get first new label for this cluster
            std::string firstLabel = expansionMap.at(clusterLabel)[0];
            
            // If this has not yet been added as a node, add it
            if (addedLabels.count(firstLabel)==0){
                Digraph::Node x = _tree.addNode();
                _nodeToId[x] = firstLabel;
                _idToNode[firstLabel] = x;
                addedLabels.insert(firstLabel);
            }
            
            // Now add new directed edge
            _tree.addArc(prevNode, _idToNode[firstLabel]);
        }
        
    }
    
    // Set root node
    _root = lemon::INVALID;
    for (NodeIt v(_tree); v != lemon::INVALID; ++v){
        if (InArcIt(_tree, v) == lemon::INVALID){
            assert(_root == lemon::INVALID);
            _root = v;
        }
    }

    init();
}

CloneTree::CloneTree(const CloneTree& oldTree,
                     const std::unordered_map<std::string, std::string>& relabelMap)
  : BaseTree()
  , _l(_tree)
{
    // Create a set to keep track of nodes we have added
    std::unordered_set<std::string> addedLabels;
    
    // Iterate over nodes in oldTree
    for (NodeIt n(oldTree.tree()); n != lemon::INVALID; ++n){
        
      // Get old label
      std::string oldLabel = oldTree.label(n);
      
      // Get new label and add node if present in new tree
      if (relabelMap.find(oldLabel) != relabelMap.end()){
        std::string newLabel = relabelMap.at(oldLabel);
        Digraph::Node x = _tree.addNode();
        _nodeToId[x] = newLabel;
        _idToNode[newLabel] = x;
        addedLabels.insert(newLabel);
      }
    }
  
    // Iterate over edges in oldTree
    for (ArcIt a(oldTree.tree()); a != lemon::INVALID; ++a)
    {
      Node u = oldTree.tree().source(a);
      Node v = oldTree.tree().target(a);
      std::string oldSource = oldTree.label(u);
      std::string oldTarget = oldTree.label(v);
      
      // Only add edges from sources that are still present
      if (relabelMap.find(oldSource) != relabelMap.end()){
        std::string newSource = relabelMap.at(oldSource);
        renameEdgesRecursive(_idToNode[newSource], v, oldTree, relabelMap);
      }
    }
    
    // Set root node
    _root = lemon::INVALID;
    for (NodeIt v(_tree); v != lemon::INVALID; ++v){
        if (InArcIt(_tree, v) == lemon::INVALID){
            assert(_root == lemon::INVALID);
            _root = v;
        }
    }

    init();
}

void CloneTree::renameEdgesRecursive(const Node& newSourceNode, const Node& currentOldNode, const CloneTree& oldTree, const std::unordered_map<std::string, std::string>& relabelMap){
  
  // Get child target label
  std::string oldTargetLabel = oldTree.label(currentOldNode);
  
  // If present, add edge
  if (relabelMap.find(oldTargetLabel) != relabelMap.end()){
    std::string newTargetLabel = relabelMap.at(oldTargetLabel);
    _tree.addArc(newSourceNode, _idToNode[newTargetLabel]);
  }
  
  // Otherwise recurse over children
  else{
    for (Digraph::ListDigraph::OutArcIt a(oldTree.tree(), currentOldNode); a != lemon::INVALID; ++a){
      renameEdgesRecursive(newSourceNode, oldTree.tree().target(a), oldTree, relabelMap);
    }
  }
}

CloneTree::CloneTree(const CloneTree& clusteredCloneTree,
                     const std::string& sep)
  : BaseTree()
  , _l(_tree)
{
  Digraph::NodeMap<Digraph::Node> oldNodeToClusterStart(clusteredCloneTree._tree,
                                                        lemon::INVALID);
  Digraph::NodeMap<Digraph::Node> oldNodeToClusterEnd(clusteredCloneTree._tree,
                                                      lemon::INVALID);
  Digraph::NodeMap<StringVector> expansion(_tree);
  
  for (NodeIt v(clusteredCloneTree._tree); v != lemon::INVALID; ++v)
  {
    //const std::string& lbl_v = clusteredCloneTree._nodeToId[v];
    const std::string& lbl_v = clusteredCloneTree.label(v);
    StringVector exp_v;
    randomNodeExpansion(lbl_v, sep, exp_v);
    assert(!exp_v.empty());
    
    Node prev_node = lemon::INVALID;
    for (const std::string& mut : exp_v)
    {
      Node new_node = _tree.addNode();
      _nodeToId[new_node] = mut;
      _idToNode[mut] = new_node;
      _l[new_node] = clusteredCloneTree._l[v];
      
      if (prev_node != lemon::INVALID)
      {
        _tree.addArc(prev_node, new_node);
      }
      
      prev_node = new_node;
    }
    assert(_idToNode.count(exp_v.front()) == 1);
    assert(_idToNode.count(exp_v.back()) == 1);
    
    oldNodeToClusterStart[v] = _idToNode[exp_v.front()];
    oldNodeToClusterEnd[v] = _idToNode[exp_v.back()];
  }
  
  // Set root node
  _root = oldNodeToClusterStart[clusteredCloneTree._root];
  
  // Insert edges between start and end nodes of path expansions
  for (ArcIt a(clusteredCloneTree._tree); a != lemon::INVALID; ++a)
  {
    Node u = clusteredCloneTree._tree.source(a);
    Node v = clusteredCloneTree._tree.target(a);
    
    Node end_u = oldNodeToClusterEnd[u];
    Node start_v = oldNodeToClusterStart[v];
    
    _tree.addArc(end_u, start_v);
  }
  
  // Initialize the other stuff now that we have a valid expanded tree
  init();
}

bool CloneTree::readLeafLabeling(std::istream& in)
{
  for (NodeIt u(_tree); u != lemon::INVALID; ++u)
  {
    _l[u] = "";
  }
  
  while (in.good())
  {
    std::string line;
    getline(in, line);
    
    if (line.empty())
      continue;
    
    StringVector s;
    boost::split(s, line, boost::is_any_of("\t "));
    
    std::string label_u = s[0];
    std::string label_s = s[1];
    
    if (_idToNode.count(label_u) == 0)
    {
      std::cerr << "Error: clone-tree vertex with label '" << label_u << "' does not exist" << std::endl;
      return false;
    }
    
    Node u = _idToNode[label_u];
//    if (!_isLeaf[u])
//    {
//      std::cerr << "Error: clone-tree vertex with label '" << label_u << "' is not a leaf" << std::endl;
//      return false;
//    }
    
    _l[u] = label_s;
  }
  
  for (Node u : _leafSet)
  {
    if (_l[u].empty())
    {
      std::cerr << "Error: leaf '" << label(u) << "' left unlabeled" << std::endl;
      return false;
    }
  }
  
  return true;
}

void CloneTree::writeLeafLabeling(std::ostream& out) const
{
  for (Node u : _leafSet)
  {
    out << _nodeToId[u] << " " << _l[u] << std::endl;
  }
}

void CloneTree::writeVertexLabeling(std::ostream& out,
                                    const StringNodeMap& lPlus) const
{
  for (NodeIt u(_tree); u != lemon::INVALID; ++u)
  {
    out << _nodeToId[u] << " " << lPlus[u] << std::endl;
  }
}

void CloneTree::writeDOT(std::ostream& out,
                         const StringNodeMap& lPlus) const
{
  StringToIntMap colorMap = generateColorMap();
  writeDOT(out, lPlus, colorMap);
}

void CloneTree::writeDOT(std::ostream& out,
                         const StringNodeMap& lPlus,
                         const StringToIntMap& colorMap) const
{
  out << "digraph T {" << std::endl;
  out << "\t{" << std::endl;
  out << "\t\trank=same" << std::endl;
  for (NodeIt u(_tree); u != lemon::INVALID; ++u)
  {
    if (_isLeaf[u])
    {
      out << "\t\t" << _tree.id(u) << " [penwidth=3,colorscheme=set19,color="
          << colorMap.find(lPlus[u])->second
          << ",label=\"" << _nodeToId[u] << "\\n" << lPlus[u] << "\"]" << std::endl;
    }
  }
  out << "\t}" << std::endl;
  
  for (NodeIt u(_tree); u != lemon::INVALID; ++u)
  {
    if (!_isLeaf[u])
    {
      out << "\t" << _tree.id(u) << " [penwidth=3,colorscheme=set19,color="
          << colorMap.find(lPlus[u])->second << ",label=\"" << _nodeToId[u] << "\"]" << std::endl;
    }
  }
  
  for (ArcIt a(_tree); a != lemon::INVALID; ++a)
  {
    Node u = _tree.source(a);
    Node v = _tree.target(a);
    
    const std::string& s_u = lPlus[u];
    const std::string& s_v = lPlus[v];
    
    out << "\t" << _tree.id(u) << " -> " << _tree.id(v);
    if (s_u == s_v)
    {
      out << " [penwidth=3,colorscheme=set19,color=" << colorMap.find(s_u)->second << "]";
    }
    else
    {
      out << " [penwidth=3,colorscheme=set19,color=\"" << colorMap.find(s_u)->second << ";0.5:" << colorMap.find(s_v)->second << "\"]";
    }
    
    out << std::endl;
  }
  
  out << "}" << std::endl;
}

void CloneTree::writeDOT(std::ostream& out,
                         const StringNodeMap& lPlus,
                         const StringToIntMap& colorMap,
                         const DoubleNodeMap& U,
                         const IntNodeMap& characterLabel) const
{
  out << "digraph T {" << std::endl;
  out << "\t{" << std::endl;
  out << "\t\trank=same" << std::endl;
  for (NodeIt u(_tree); u != lemon::INVALID; ++u)
  {
    if (_isLeaf[u])
    {
      out << "\t\t" << _tree.id(u) << " [width=0.2,height=0.2,penwidth=3,colorscheme=set19,color="
          << colorMap.find(lPlus[u])->second
          << ",label=\"";
      out << std::setprecision(3) << 100 * U[u] << "%";
      out << "\"]" << std::endl;
    }
  }
  out << "\t}" << std::endl;
  
  for (NodeIt u(_tree); u != lemon::INVALID; ++u)
  {
    if (!_isLeaf[u])
    {
      out << "\t" << _tree.id(u) << " [width=0.2,height=0.2,penwidth=3,colorscheme=set19,color="
      << colorMap.find(lPlus[u])->second << ",label=\"\"]" << std::endl;
    }
  }
  out << "\tinv [style=\"invis\"]" << std::endl;
  out << "\tinv -> "
  << _tree.id(_root) <<"[penwidth=3,colorscheme=set19,color="
  << colorMap.find(lPlus[_root])->second << ",label=\"" << characterLabel[_root]
  << "\",fontsize=18]" << std::endl;
  
  for (ArcIt a(_tree); a != lemon::INVALID; ++a)
  {
    Node u = _tree.source(a);
    Node v = _tree.target(a);
    
    const std::string& s_u = lPlus[u];
    const std::string& s_v = lPlus[v];
    
    out << "\t" << _tree.id(u) << " -> " << _tree.id(v);
    if (s_u == s_v)
    {
      out << " [fontsize=18,penwidth=3,colorscheme=set19,color=" << colorMap.find(s_u)->second;
    }
    else
    {
      out << " [fontsize=18,penwidth=3,colorscheme=set19,color=\"" << colorMap.find(s_u)->second << ";0.5:" << colorMap.find(s_v)->second << "\"";
    }
    if (characterLabel[_tree.target(a)] != -1)
    {
      out << ",label=\"" << characterLabel[_tree.target(a)] << "\"";
    }
    out << "]";
    out << std::endl;
  }
  
  out << "}" << std::endl;
}


void CloneTree::writeDOT(std::ostream& out,
                         const StringNodeMap& lPlus,
                         const StringToIntMap& colorMap,
                         const DoubleVectorNodeMap& U,
                         const StringNodeMap& characterLabel) const
{
  out << "digraph T {" << std::endl;
  out << "\t{" << std::endl;
  out << "\t\trank=same" << std::endl;
  for (NodeIt u(_tree); u != lemon::INVALID; ++u)
  {
    if (_isLeaf[u])
    {
      out << "\t\t" << _tree.id(u) << " [width=0.2,height=0.2,penwidth=3,colorscheme=set19,color="
          << colorMap.find(lPlus[u])->second
          << ",label=\"";
      
      const DoubleVector& vecU = U[u];
      bool first = true;
      for (int p = 0; p < vecU.size(); ++p)
      {
//        if (g_tol.nonZero(vecU[p]))
        {
          if (first)
            first = false;
          else
            out << "\\n";
          
          out << std::setprecision(3) << 100 * vecU[p] << "%";
        }
      }
      out << "\"]" << std::endl;
    }
  }
  out << "\t}" << std::endl;
  
  for (NodeIt u(_tree); u != lemon::INVALID; ++u)
  {
    if (!_isLeaf[u])
    {
      out << "\t" << _tree.id(u) << " [width=0.2,height=0.2,penwidth=3,colorscheme=set19,color="
          << colorMap.find(lPlus[u])->second << ",label=\"\"]" << std::endl;
    }
  }
  out << "\tinv [style=\"invis\"]" << std::endl;
  out << "\tinv -> "
      << _tree.id(_root) <<"[penwidth=3,colorscheme=set19,color="
      << colorMap.find(lPlus[_root])->second << ",label=\"" << characterLabel[_root]
      << "\",fontsize=18]" << std::endl;
  
  for (ArcIt a(_tree); a != lemon::INVALID; ++a)
  {
    Node u = _tree.source(a);
    Node v = _tree.target(a);
    
    const std::string& s_u = lPlus[u];
    const std::string& s_v = lPlus[v];
    
    out << "\t" << _tree.id(u) << " -> " << _tree.id(v);
    if (s_u == s_v)
    {
      out << " [fontsize=18,penwidth=3,colorscheme=set19,color=" << colorMap.find(s_u)->second;
    }
    else
    {
      out << " [fontsize=18,penwidth=3,colorscheme=set19,color=\"" << colorMap.find(s_u)->second << ";0.5:" << colorMap.find(s_v)->second << "\"";
    }
    if (!characterLabel[_tree.target(a)].empty())
    {
      out << ",label=\"" << characterLabel[_tree.target(a)] << "\"";
    }
    out << "]";
    out << std::endl;
  }
  
  out << "}" << std::endl;
}

void CloneTree::writeDOT(std::ostream& out,
                         const StringNodeMap& lPlus,
                         const StringToIntMap& colorMap,
                         const DoubleVectorNodeMap& F,
                         const DoubleNodeMap& U) const
{
  out << "digraph T {" << std::endl;
  out << "\t{" << std::endl;
  out << "\t\trank=same" << std::endl;
  for (NodeIt u(_tree); u != lemon::INVALID; ++u)
  {
    if (_isLeaf[u])
    {
      out << "\t\t" << _tree.id(u) << " [penwidth=3,colorscheme=set19,color="
          << colorMap.find(lPlus[u])->second
          << ",label=\"" << _nodeToId[u];
      out << "\\n" << U[u] << "\\n" << lPlus[u] << "\"]" << std::endl;
    }
  }
  out << "\t}" << std::endl;
  
  for (NodeIt u(_tree); u != lemon::INVALID; ++u)
  {
    if (!_isLeaf[u])
    {
      out << "\t" << _tree.id(u) << " [penwidth=3,colorscheme=set19,color="
          << colorMap.find(lPlus[u])->second << ",label=\"" << _nodeToId[u];
      for (int s = 0; s < F[u].size(); ++s)
      {
        out << "\\n" << F[u][s];
      }
      out << "\"]" << std::endl;
    }
  }
  
  for (ArcIt a(_tree); a != lemon::INVALID; ++a)
  {
    Node u = _tree.source(a);
    Node v = _tree.target(a);
    
    const std::string& s_u = lPlus[u];
    const std::string& s_v = lPlus[v];
    
    out << "\t" << _tree.id(u) << " -> " << _tree.id(v);
    if (s_u == s_v)
    {
      out << " [penwidth=3,colorscheme=set19,color=" << colorMap.find(s_u)->second << "]";
    }
    else
    {
      out << " [penwidth=3,colorscheme=set19,color=\"" << colorMap.find(s_u)->second << ";0.5:" << colorMap.find(s_v)->second << "\"]";
    }
    
    out << std::endl;
  }
  
  out << "}" << std::endl;
}

void CloneTree::writeDOT(std::ostream& out) const
{
  StringToIntMap colorMap = generateColorMap();
  writeDOT(out, colorMap);
}

void CloneTree::writeDOT(std::ostream& out,
                         const StringToIntMap& colorMap) const
{
  out << "digraph T {" << std::endl;
  out << "\t{" << std::endl;
  out << "\t\trank=same" << std::endl;
  for (NodeIt u(_tree); u != lemon::INVALID; ++u)
  {
    if (_isLeaf[u])
    {
      out << "\t\t" << _tree.id(u) << " [penwidth=3,colorscheme=set19,color="
          << colorMap.find(l(u))->second << ",label=\"" << _nodeToId[u] << "\\n" << _l[u] << "\"]" << std::endl;
    }
  }
  out << "\t}" << std::endl;
  
  for (NodeIt u(_tree); u != lemon::INVALID; ++u)
  {
    if (!_isLeaf[u])
    {
      out << "\t" << _tree.id(u) << " [label=\"" << _nodeToId[u] << "\"]" << std::endl;
    }
  }
  
  for (ArcIt a(_tree); a != lemon::INVALID; ++a)
  {
    out << "\t" << _tree.id(_tree.source(a)) << " -> " << _tree.id(_tree.target(a)) << std::endl;
  }
  
  out << "}" << std::endl;
}

void CloneTree::renameNode(const Node& node, const std::string& label){
  _nodeToId[node] = label;
  _idToNode[label] = node;
}

void CloneTree::addNodes(const std::vector<std::string>& labels){
  
  for (auto x: labels) {
    Node new_node = _tree.addNode();
    _nodeToId[new_node] = x;
    _idToNode[x] = new_node;
    _l[new_node] = "";
  }
  
  // Update dictionaries associated with tree
  init();
}

void CloneTree::addEdges(const std::vector<StringPair>& edgeMap, std::string sep){
    
    // Split mutation clusters and get mutations for this tree
    std::unordered_set<std::string> muts;
    for (NodeIt v(_tree); v != lemon::INVALID; ++v){
      std::string label = _nodeToId[v];
      StringVector s;
      boost::split(s, label, boost::is_any_of(sep));
      muts.insert(s.begin(), s.end());
    }
    
    for (const StringPair& element : edgeMap){
        
        // Add starting node if it does not exist
        if (muts.find(element.first) == muts.end()) {
          Node new_node = _tree.addNode();
          _nodeToId[new_node] = element.first;
          _idToNode[element.first] = new_node;
          _l[new_node] = "";
          muts.insert(element.first);
        }
        
        // Add ending node if it does not exist
        if (muts.find(element.second) == muts.end()) {
          Node new_node = _tree.addNode();
          _nodeToId[new_node] = element.second;
          _idToNode[element.second] = new_node;
          _l[new_node] = "";
          muts.insert(element.second);
        }
        
        // Add edge
        _tree.addArc(_idToNode[element.first], _idToNode[element.second]);
    }
    
    // Update dictionaries associated with tree
    init();
}


void CloneTree::removeNodes(const std::vector<std::string>& nodeVect){
  
  for (auto x: nodeVect){
    
    if (_idToNode.find(x) != _idToNode.end()){
      _tree.erase(_idToNode[x]);
    
      // Other NodeMaps should automatically be updated according to documentation
      _idToNode.erase(x);
    }
  }
  
  // Update dictionaries associated with tree
  init();
}

void CloneTree::mergeSameSiblingLeaves()
{
  typedef std::map<std::string, NodeList> LocationToNodeList;
  char buf[1024];
  
  for (NodeIt v(_tree); v != lemon::INVALID; ++v)
  {
    if (_isLeaf[v])
    {
      continue;
    }
    
    LocationToNodeList map;
    for (OutArcIt a(_tree, v); a != lemon::INVALID; ++a)
    {
      Node w = _tree.target(a);
      if (_isLeaf[w])
      {
        map[l(w)].push_back(w);
      }
    }
    
    int j = 0;
    for (const auto& sV : map)
    {
      const NodeList& V = sV.second;
      Node grandParent = _tree.source(InArcIt(_tree, V.front()));
      
      if (V.size() > 1)
      {
        // make a backbone
        NodeVector backbone;
        for (int i = 0; i < V.size() - 1; ++i)
        {
          backbone.push_back(_tree.addNode());
          
          snprintf(buf, 1024, "%s^%d^%d", _nodeToId[grandParent].c_str(), j, i);
          std::string lbl_parent =  buf;
          _nodeToId[backbone.back()] = lbl_parent;
          _idToNode[lbl_parent] = backbone.back();
        }
        
        _tree.addArc(grandParent, backbone.front());
        
        for (int i = 1; i < V.size() - 1; ++i)
        {
          _tree.addArc(backbone[i-1], backbone[i]);
        }
        
        int i = 0;
        for (Node w : V)
        {
          _tree.erase(InArcIt(_tree, w));
          _tree.addArc(backbone[i], w);
          if (i < V.size() - 2)
            ++i;
        }
        ++j;
      }
    }
  }
  
  init();
}

void CloneTree::writeDOT(std::ostream& out,
                         const StringNodeMap& lPlus,
                         const StringToIntMap& colorMap,
                         const DoubleNodeMap& U,
                         const StringNodeMap& characterLabel) const
{
  DoubleVectorNodeMap newU(_tree);
  for (NodeIt v(_tree); v != lemon::INVALID; ++v)
  {
    newU[v] = DoubleVector(1, U[v]);
  }
  
  writeDOT(out, lPlus, colorMap, newU, characterLabel);
}

void CloneTree::writeDOT(std::ostream& out,
                         const StringToIntMap& colorMap,
                         const DoubleNodeMap& U) const
{
  out << "digraph T {" << std::endl;
  out << "\t{" << std::endl;
  out << "\t\trank=same" << std::endl;
  for (NodeIt u(_tree); u != lemon::INVALID; ++u)
  {
    if (_isLeaf[u])
    {
      out << "\t\t" << _tree.id(u) << " [penwidth=3,colorscheme=set19,color="
          << colorMap.find(l(u))->second << ",label=\""
          << _nodeToId[u] << "\\n"
          << _l[u] << "\\n"
          << std::setprecision(2) << 100 * U[u]
          << "%"<< "\"]" << std::endl;
    }
  }
  out << "\t}" << std::endl;
  
  for (NodeIt u(_tree); u != lemon::INVALID; ++u)
  {
    if (!_isLeaf[u])
    {
      out << "\t" << _tree.id(u) << " [label=\"" << _nodeToId[u] << "\"]" << std::endl;
    }
  }
  
  for (ArcIt a(_tree); a != lemon::INVALID; ++a)
  {
    out << "\t" << _tree.id(_tree.source(a)) << " -> " << _tree.id(_tree.target(a)) << std::endl;
  }
  
  out << "}" << std::endl;
}

void CloneTree::writeDOT(std::ostream& out,
                         const StringNodeMap& lPlus,
                         const StringToIntMap& colorMap,
                         const DoubleVectorNodeMap& U) const
{
  out << "digraph T {" << std::endl;
  out << "\t{" << std::endl;
  out << "\t\trank=same" << std::endl;
  for (NodeIt u(_tree); u != lemon::INVALID; ++u)
  {
    if (_isLeaf[u])
    {
      out << "\t\t" << _tree.id(u) << " [penwidth=3,colorscheme=set19,color="
          << colorMap.find(l(u))->second << ",label=\""
          << _nodeToId[u] << "\\n"
          << _l[u] << "\\n";
      
      bool first = true;
      for (double usage : U[u])
      {
        if (first)
          first = false;
        else
          out << " ";
        out << std::setprecision(2) << 100 * usage << "%";
      }
      out << "\"]" << std::endl;
    }
  }
  out << "\t}" << std::endl;
  
  for (NodeIt u(_tree); u != lemon::INVALID; ++u)
  {
    if (!_isLeaf[u])
    {
      out << "\t" << _tree.id(u) << " [penwidth=3,colorscheme=set19,color="
          << colorMap.find(lPlus[u])->second << ",label=\""
          << _nodeToId[u] << "\"]" << std::endl;
    }
  }
  
  for (ArcIt a(_tree); a != lemon::INVALID; ++a)
  {
    Node u = _tree.source(a);
    Node v = _tree.target(a);
    
    const std::string& s_u = lPlus[u];
    const std::string& s_v = lPlus[v];
    
    out << "\t" << _tree.id(u) << " -> " << _tree.id(v);
    if (s_u == s_v)
    {
      out << " [penwidth=3,colorscheme=set19,color=" << colorMap.find(s_u)->second << "]";
    }
    else
    {
      out << " [penwidth=3,colorscheme=set19,color=\"" << colorMap.find(s_u)->second << ";0.5:" << colorMap.find(s_v)->second << "\"]";
    }
    
    out << std::endl;
  }
  
  out << "}" << std::endl;
}

void CloneTree::writeDOT(std::ostream& out,
                         const StringToIntMap& colorMap,
                         const DoubleVectorNodeMap& U) const
{
  out << "digraph T {" << std::endl;
  out << "\t{" << std::endl;
  out << "\t\trank=same" << std::endl;
  for (NodeIt u(_tree); u != lemon::INVALID; ++u)
  {
    if (_isLeaf[u])
    {
      out << "\t\t" << _tree.id(u) << " [penwidth=3,colorscheme=set19,color="
      << colorMap.find(l(u))->second << ",label=\""
      << _nodeToId[u] << "\\n"
      << _l[u] << "\\n";
      
      bool first = true;
      for (double usage : U[u])
      {
        if (first)
          first = false;
        else
          out << " ";
        out << std::setprecision(2) << 100 * usage << "%";
      }
      out << "\"]" << std::endl;
    }
  }
  out << "\t}" << std::endl;
  
  for (NodeIt u(_tree); u != lemon::INVALID; ++u)
  {
    if (!_isLeaf[u])
    {
      out << "\t" << _tree.id(u) << " [label=\"" << _nodeToId[u] << "\"]" << std::endl;
    }
  }
  
  for (ArcIt a(_tree); a != lemon::INVALID; ++a)
  {
    out << "\t" << _tree.id(_tree.source(a)) << " -> " << _tree.id(_tree.target(a)) << std::endl;
  }
  
  out << "}" << std::endl;
}

SplitSet CloneTree::getSplits() const
{
  SplitSet splits;
  const NodeSet allLeaves = leafSet();
  for (NodeIt v(_tree); v != lemon::INVALID; ++v)
  {
    const NodeSet& L_v = leafSubset(v);
    NodeSet L_complement_v;
    std::set_difference(allLeaves.begin(), allLeaves.end(),
                        L_v.begin(), L_v.end(),
                        std::inserter(L_complement_v,
                                      L_complement_v.begin()));
    
    StringSet split_v;
    StringSet complement_v;
    for (Node u : L_v)
    {
      const std::string& sStr = l(u);
      split_v.insert(sStr);
    }
    for (Node u : L_complement_v)
    {
      const std::string& sStr = l(u);
      complement_v.insert(sStr);
    }
    
    Split split;
    split.insert(split_v);
    split.insert(complement_v);
    splits.insert(split);
  }
  
  return splits;
}

void CloneTree::randomNodeExpansion(const std::string& label,
                                    const std::string& sep,
                                    StringVector& expansion)
{
  // Split string
  expansion.clear();
  boost::split(expansion, label, boost::is_any_of(sep));
  
  // Shuffle the string vector
  if (expansion.size() > 1)
  {
    std::shuffle(expansion.begin(), expansion.end(), g_rng);
  }
}

void CloneTree::levelNodeExpansion(const std::string& label, const std::string& sep, StringVector& expansion, const std::map<std::string, int> mutationToLevelMap){
    
    // Split string
    StringVector s;
    boost::split(s, label, boost::is_any_of(sep));
    
    // Get number of mutations
    int numMut(s.size());
    
    // Initialize map from depth to vector of mutations in cluster at this depth
    std::map<int, std::vector<std::string>> leveltoMutationMap;
    
    // Initialize vector of mutations not in consensus
    std::vector<std::string> extraMutations;
    
    // Loop over mutations in this node
    for (auto& mutation : s){
        
        // Check to see if mutation is in consensus depth map
        auto iter = mutationToLevelMap.find(mutation);
        
        // If so, we get the depth
        if (iter != mutationToLevelMap.end() ){
            int d = iter->second;
            
            // We then add this mutation to our level map
            auto iter2 = leveltoMutationMap.find(d);
            
            if (iter2 != leveltoMutationMap.end()){
                std::vector<std::string>& strvect = iter2->second;
                strvect.push_back(mutation);
            }
            else{
                std::vector<std::string> vect;
                vect.push_back(mutation);
                leveltoMutationMap[d]=vect;
            }
        }
        // Otherwise this mutation is missing from the consensus
        else{
            extraMutations.push_back(mutation);
        }
    }
    
    // GET LEVEL ORDER MUTATIONS INTO ONE QUEUE
    // Loop over map and insert into expansion vector
    int check(-1);
    std::queue<std::string> tmp;
    for (auto const& x : leveltoMutationMap){
        
        // Make sure map is iterating in order of depth
        assert(x.first > check);
        check = x.first;
        
        // Shuffle mutations at this depth
        std::vector<std::string> v(x.second);
        std::shuffle(v.begin(), v.end(), g_rng);
        for (auto const m: v){
            tmp.push(m);
        }
    }
    
    // SELECT PLACES FOR EXRA MUTATIONS
    // Now select random subset of indices for where to place missing mutations in final ordering
    std::vector<int> index(numMut);
    std::iota(index.begin(), index.end(), 0);
    std::shuffle(index.begin(), index.end(), g_rng);
    std::unordered_set<int> missingIndices(index.begin(), index.begin () + extraMutations.size());
    assert(missingIndices.size() == extraMutations.size());
    
    std::shuffle(extraMutations.begin(), extraMutations.end(), g_rng);
    
    // FILL IN EXPANSION VECTOR WITH LEVEL ORDER AND EXTRA MUTATIONS
    expansion.resize(numMut);
    for(int i=0; i < expansion.size(); i++){
        // If this is not an index to place a missing mutation
        // Then we pop off the next level-order mutation
        if (missingIndices.find(i) == missingIndices.end()){
            expansion.at(i) = tmp.front();
            tmp.pop();
        }
        // Otherwise, we add the next missing mutation
        else{
            expansion.at(i) = extraMutations.back();
            extraMutations.pop_back();
        }
    }
}


void CloneTree::expandNode(const std::string& label,
                           const std::string& sep,
                           std::vector<StringVector>& expansion)
{
    // Split string
    StringVector s;
    boost::split(s, label, boost::is_any_of(sep));
    
    // Sorting the string vector
    sort(s.begin(), s.end());
    
    // Add all permutations of the string vector to expansion
    do {
        expansion.push_back(s);
    } while ( std::next_permutation(s.begin(),s.end()) );
    
}

void CloneTree::permuteTree(const CloneTree& C,
                           std::vector<CloneTree>& permutation)
{
  // Get tree from C
  const Digraph& t = C.tree();
  
  // Get root from C
  Node root = C.root();
  std::string rLabel = C.label(root);
  
  // Get node labels excluding the root
  StringVector labels;
  for (NodeIt u(t); u != lemon::INVALID; ++u)
  {
    std::string x = C.label(u);
    if (x != rLabel){
      labels.push_back(x);
    }
  }
  
  // Randomly shuffle node labels
  std::shuffle(labels.begin(), labels.end(), g_rng);
  
  // Create new label map that will house permutation
  StringNodeMap new_l(t);
  
  // Iterate over nodes and relabel nodemap
  for (NodeIt n(t); n != lemon::INVALID; ++n){
    
    // If node is root we keep same label
    if (C.label(n) == rLabel){
      new_l[n] = rLabel;
    }
    
    // If node is not root we will relabel it
    else {
      
      // Get new random label
      std::string x = labels.back();
      labels.pop_back();
      
      // Assign label
      new_l[n] = x;
    }
  }
  
  if (labels.size() > 0)
  {
    throw std::runtime_error("Error: number of labels remaining after permutation should be 0");
  }
  
  // Create permutation tree
  CloneTree P(t, root, new_l, new_l);
  
  // Append new clone tree to vector
  permutation.push_back(P);
    
}

void CloneTree::permuteTree(const CloneTree& C,
                           std::vector<CloneTree>& permutation, std::map<std::string, std::string>& map)
{
  // Get tree from C
  const Digraph& t = C.tree();
    
  // Get root from C
  Node root = C.root();
  std::string rLabel = C.label(root);
  
  // Create new label map that will house permutation
  StringNodeMap new_l(t);
  
  // Iterate over nodes and relabel nodemap
  for (NodeIt n(t); n != lemon::INVALID; ++n){
      
      // Assign label
      new_l[n] = map[C.label(n)];
  }
  
  if (new_l[root] != rLabel)
  {
    throw std::runtime_error("Error: cannot permute root label.");
  }
  
  // Create permutation tree
  CloneTree P(t, root, new_l, new_l);
  
  // Append new clone tree to vector
  permutation.push_back(P);
    
}

void CloneTree::permuteMap(const CloneTree& C, std::map<std::string, std::string>& map){
    
    // Get tree from C
    const Digraph& t = C.tree();
    
    // Get root from C
    Node root = C.root();
    std::string rLabel = C.label(root);
    
    // Get node labels excluding the root
    StringVector labels;
    for (NodeIt u(t); u != lemon::INVALID; ++u)
    {
      std::string x = C.label(u);
      if (x != rLabel){
        labels.push_back(x);
      }
    }
    
    // Randomly shuffle node labels
    std::shuffle(labels.begin(), labels.end(), g_rng);
    
    // Iterate over nodes and relabel nodemap
    for (NodeIt n(t); n != lemon::INVALID; ++n){
      
      // If node is root we keep same label
      if (C.label(n) == rLabel){
        map[rLabel] = rLabel;
      }
      
      // If node is not root we will relabel it
      else {
        
        // Get new random label
        std::string x = labels.back();
        labels.pop_back();
        
        // Assign label
        map[C.label(n)] = x;
      }
    }
    
    if (labels.size() > 0)
    {
      throw std::runtime_error("Error: number of labels remaining after permutation should be 0");
    }
    
}

void CloneTree::expandTree(const CloneTree& C,
                           const std::string& sep,
                           std::vector<CloneTree>& expansion)
{
    // Get tree from C
    const Digraph& t = C.tree();
    
    // Initialize vector of nodes and expansions
    std::vector<std::string> nodeV;
    std::vector<std::vector<StringVector>> expV;
    std::vector<std::vector<int>> indV;
    
    // Iterate over nodes of clone tree and expand each node
    for (NodeIt u(t); u != lemon::INVALID; ++u)
    {
        // Create vector to pass in
        std::vector<StringVector> tmp;
        
        // Expand node
        CloneTree::expandNode(C.label(u), sep, tmp);
        
        // Update vectors
        nodeV.push_back(C.label(u));
        expV.push_back(tmp);
        std::vector<int> v(tmp.size());
        std::iota(v.begin(), v.end(), 0);
        indV.push_back(v);
    }
    
    // Get cartesian product of vector of indices
    std::vector<std::vector<int>> prodV(CloneTree::product(indV));
    
    // For each product, construct a dictionary of node labels to expansion
    std::vector<std::unordered_map<std::string, StringVector>> expDicts;
    for (const auto& p : prodV)
    {
        std::unordered_map<std::string, StringVector> tmp;
        
        // For each node, add expansion based on product indices
        for (unsigned int i=0; i<p.size(); i++){
            
            // i.e., for this node we will use the p[i]th expansion in this map
            tmp[nodeV[i]] = expV[i][p[i]];
        }
        
        expDicts.push_back(tmp);
    }
    
    // For each map, we now construct a clone tree
    for (auto& x: expDicts)
    {
        expansion.push_back(CloneTree(C, x));
    }
}

CloneTree CloneTree::levelExpandTree(const std::string sep, const std::map<std::string, int>& mutationToLevelMap, const CloneTree& clusterTree){
  
  // Get tree from C
  const Digraph& t = clusterTree.tree();
  
  // Initialize vector of nodes and expansions
  std::unordered_map<std::string, StringVector> expansionMap;
  
  // Iterate over nodes of clone tree and expand each node
  for (NodeIt u(t); u != lemon::INVALID; ++u)
  {
    // Create expansion vector to pass in
    StringVector tmp;
      
    // Expand node
    CloneTree::levelNodeExpansion(clusterTree.label(u), sep, tmp, mutationToLevelMap);
      
    // Update dictionary
    expansionMap[clusterTree.label(u)]=tmp;
  }
  
  // Build and return clone tree
  CloneTree exp(clusterTree, expansionMap);
  return exp;
}

std::vector<std::vector<int>> CloneTree::product(const std::vector<std::vector<int>>& lists)
{
    std::vector<std::vector<int>> result;
    
    if (std::find_if(std::begin(lists), std::end(lists),
                     [](std::vector<int> e) -> bool { return e.size() == 0; }) != std::end(lists)) {
        return result;
    }
    for (auto& e : lists[0]) {
        result.push_back({ e });
    }
    for (size_t i = 1; i < lists.size(); ++i) {
        std::vector<std::vector<int>> temp;
        for (auto& e : result) {
            for (auto f : lists[i]) {
                auto e_tmp = e;
                e_tmp.push_back(f);
                temp.push_back(e_tmp);
            }
        }
        result = temp;
    }
    return result;
}

std::ostream& operator<<(std::ostream& out, const CloneTree& T)
{
  out << lemon::countArcs(T.tree()) << " #edges" << std::endl;
  T.write(out);
  out << T.leafSet().size() << " #leaves" << std::endl;
  T.writeLeafLabeling(out);
  
  return out;
}

std::istream& operator>>(std::istream& in, CloneTree& T)
{
  std::string line;
  getline(in, line);
  
  StringVector s;
  boost::split(s, line, boost::is_any_of("\t "));
  int nrEdges = boost::lexical_cast<int>(s[0]);
  if (nrEdges < 0)
  {
    throw std::runtime_error("Error: invalid number of edges");
  }
  
  std::stringstream ss;
  for (int i = 0; i < nrEdges; ++i)
  {
    getline(in, line);
    ss << line << std::endl;
  }
  
  if (!T.read(ss))
  {
    throw std::runtime_error("Error: incorrect edge list");
  }
  
  getline(in, line);
  
  boost::split(s, line, boost::is_any_of("\t "));
  int nrLeaves = boost::lexical_cast<int>(s[0]);
  if (nrLeaves < 0)
  {
    throw std::runtime_error("Error: invalid number of leaves");
  }
  
  ss.clear();
  for (int i = 0; i < nrLeaves; ++i)
  {
    getline(in, line);
    ss << line << std::endl;
  }
  
  if (!T.readLeafLabeling(ss))
  {
    throw std::runtime_error("Error: incorrect leaf labeling");
  }
  
  return in;
}

std::ostream& operator<<(std::ostream& out,
                         const CloneTreeVector& trees)
{
  out << trees.size() << " #trees" << std::endl;
  int idx = 0;
  for (const CloneTree& T : trees)
  {
    out << lemon::countArcs(T.tree()) << " #edges, tree " << idx + 1 << std::endl;
    T.write(out);
    ++idx;
  }
  return out;
}

std::istream& operator>>(std::istream& in,
                         CloneTreeVector& trees)
{
  int nrTrees = -1;
  
  std::string line;
  getline(in, line);
  
  while (line.empty() || line[0] == '#')
  {
    getline(in, line);
  }
  
  std::stringstream ss(line);
  ss >> nrTrees;
  
  if (nrTrees < 0)
  {
    throw std::runtime_error("Error: number of trees should be nonnegative");
  }
  
  for (int i = 0; i < nrTrees; ++i)
  {
    int nrEdges = -1;
    
    getline(in, line);
    ss.clear();
    ss.str(line);
    ss >> nrEdges;
    
    if (nrEdges < 0)
    {
      throw std::runtime_error("Error: number of edges should be nonnegative");
    }
    
    ss.clear();
    ss.str("");
    for (int j = 0; j < nrEdges; ++j)
    {
      getline(in, line);
      ss << line << std::endl;
    }
    
    CloneTree T;
    if (!T.read(ss))
      throw std::runtime_error("");
    
//    // Expand Tree
////    std::cout << "Expanding Tree Now" << std::endl;
//    std::vector<CloneTree> exp;
//    std::string sep(";");
//    CloneTree::expandTree(T, sep, exp);
//    int k(0);
//    //std::cout << T << std::endl;
//    for (auto x: exp){
//          trees.push_back(x);
//          //std::cout << x << std::endl;
//          k=k+1;
//    }
////    std::cout << k << std::endl;
    
    trees.push_back(T);
  }
  
  return in;
}

void expandCloneTreeVector(const CloneTreeVector& ctv,
                           const std::string& sep,
                           CloneTreeVector& expCtv)
{
  for (const CloneTree& cloneTree : ctv)
  {
    CloneTree::expandTree(cloneTree, sep, expCtv);
  }
}

