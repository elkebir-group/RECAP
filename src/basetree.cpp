/*
 * basetree.cpp
 *
 *  Created on: 24-oct-2016
 *      Author: M. El-Kebir
 */

#include "basetree.h"
#include <lemon/bfs.h>
#include <lemon/connectivity.h>

BaseTree::BaseTree()
  : _tree()
  , _root(lemon::INVALID)
  , _nodeToId(_tree)
  , _idToNode()
  , _isLeaf(_tree, false)
  , _leafSet()
  , _leafSubset(_tree)
  , _level(_tree)
  , _edgeSet()
  , _transEdgeSet()
{
}

BaseTree::BaseTree(const BaseTree& other)
  : _tree()
  , _root(lemon::INVALID)
  , _nodeToId(_tree)
  , _idToNode()
  , _isLeaf(_tree, false)
  , _leafSet()
  , _leafSubset(_tree)
  , _level(_tree)
  , _edgeSet(other._edgeSet)
  , _transEdgeSet(other._transEdgeSet)
{
  if (other._idToNode.size() > 0)
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
    
    init();
  }
}

BaseTree::BaseTree(const StringSet& idSet,
                   const std::string& rootId,
                   const StringVector& pruferSeq)
  : _tree()
  , _nodeToId(_tree)
  , _idToNode()
  , _isLeaf(_tree, false)
  , _leafSet()
  , _leafSubset(_tree)
  , _level(_tree)
  , _edgeSet()
  , _transEdgeSet()
{
  // 1. Construct undirected tree from pruferSeq
  Graph T;
  std::map<std::string, Graph::Node> idToT;
  Graph::NodeMap<std::string> TtoId(T);
  
  // 1a. Construct nodes and their degrees of T
  Graph::NodeMap<int> degT(T, 0);
  for (const std::string& id : idSet)
  {
    Graph::Node v = T.addNode();
    idToT[id] = v;
    TtoId[v] = id;
    degT[v] = 1;
    
    for (const std::string& id2 : pruferSeq)
    {
      if (id2 == id)
      {
        ++degT[v];
      }
    }
  }
  
  for (const std::string& id_i : pruferSeq)
  {
    assert(idToT.count(id_i));
    Graph::Node v_i = idToT[id_i];
    
    // find the first node, j, with degree equal to 1
    std::string min_id_j;
    for (const std::string& id_j : idSet)
    {
      assert(idToT.count(id_j));
      Graph::Node v_j = idToT[id_j];
      if (degT[v_j] == 1 && (min_id_j.empty() || min_id_j > id_j))
      {
        min_id_j = id_j;
      }
    }
    
    assert(!min_id_j.empty());
    assert(idToT.count(min_id_j));
    Graph::Node v_j = idToT[min_id_j];
    
    T.addEdge(v_j, v_i);
    --degT[v_i];
    --degT[v_j];
  }
  
  // 1b. Add remaining edge
  Graph::Node v_i = lemon::INVALID;
  Graph::Node v_j = lemon::INVALID;
  for (const std::string& id : idSet)
  {
    assert(idToT.count(id));
    Graph::Node v = idToT[id];
    if (degT[v] == 1 && v_i == lemon::INVALID)
    {
      v_i = v;
    }
    else if (degT[v] == 1 && v_j == lemon::INVALID)
    {
      v_j = v;
    }
    else if (degT[v] == 1)
    {
      assert(false);
    }
  }
  T.addEdge(v_i, v_j);
  
  assert(lemon::acyclic(T));
  
  // 2. Root T
  _root = _tree.addNode();
  _nodeToId[_root] = rootId;
  _idToNode[rootId] = _root;
  rootUndirectedTree(T, TtoId, idToT[rootId], _root);
  
  // 3. Initialize data structures
  init();
  
  // 4. Initialize edges and transEdge
  _edgeSet = getParentalPairs();
  _transEdgeSet = getAncestralPairs();
}

BaseTree::BaseTree(const Digraph& T,
                   Node root,
                   const StringNodeMap& label)
  : _tree()
  , _nodeToId(_tree)
  , _idToNode()
  , _isLeaf(_tree, false)
  , _leafSet()
  , _leafSubset(_tree)
  , _level(_tree)
{
  lemon::digraphCopy(T, _tree)
    .node(root, _root)
    .nodeMap(label, _nodeToId)
    .run();
  
  for (NodeIt u(_tree); u != lemon::INVALID; ++u)
  {
    const std::string& str = _nodeToId[u];
    assert(_idToNode.count(str) == 0);
    _idToNode[str] = u;
  }
  
  init();
  
  // Initialize edges and transEdge
  _edgeSet = getParentalPairs();
  _transEdgeSet = getAncestralPairs();
}

void BaseTree::rootUndirectedTree(const Graph& T,
                                  const Graph::NodeMap<std::string>& TtoId,
                                  Graph::Node& v,
                                  Node vv)
{
  for (Graph::IncEdgeIt e(T, v); e != lemon::INVALID; ++e)
  {
    Graph::Node w = T.oppositeNode(v, e);
    const std::string& id_w = TtoId[w];
    
    if (_idToNode.count(id_w) == 0)
    {
      Node ww = _tree.addNode();
      _idToNode[id_w] = ww;
      _nodeToId[ww] = id_w;
      
      _tree.addArc(vv, ww);
      rootUndirectedTree(T, TtoId, w, ww);
    }
  }
}

StringVector BaseTree::getPruferSequence() const
{
  const int n = lemon::countNodes(_tree);
 
  // construct Prufer sequence
  NodeSet leaves = _leafSet;
  BoolNodeMap visitedMap(_tree, false);
  StringVector pruferSeq;
  while (pruferSeq.size() < n - 2)
  {
    std::string min_id;
    Node min_leaf = lemon::INVALID;
    
    // Find leaf with smallest label
    for (Node u : leaves)
    {
      if (min_leaf == lemon::INVALID || min_id > _nodeToId[u])
      {
        min_id = _nodeToId[u];
        min_leaf = u;
      }
    }
    
    visitedMap[min_leaf] = true;
    Node pi = parent(min_leaf);
    
    pruferSeq.push_back(_nodeToId[pi]);
    
    // remove leaf
    leaves.erase(min_leaf);
    
    // check whether to add pi to leaves
    bool add = true;
    for (OutArcIt a(_tree, pi); a != lemon::INVALID; ++a)
    {
      Node v = _tree.target(a);
      add &= visitedMap[v];
    }
    
    if (add)
    {
      leaves.insert(pi);
    }
  }
  
  return pruferSeq;
}

void BaseTree::write(std::ostream& out) const
{
  for (ArcIt a(_tree); a != lemon::INVALID; ++a)
  {
    Node u = _tree.source(a);
    Node v = _tree.target(a);
    
    out << _nodeToId[u] << " " << _nodeToId[v] << std::endl;
  }
}

bool BaseTree::read(std::istream& in)
{
  _idToNode.clear();
  
  while (in.good())
  {
    std::string line;
    getline(in, line);
    
    if (line.empty())
      break;
    
    StringVector s;
    boost::split(s, line, boost::is_any_of("\t "));
    
    if (s.size() < 2)
    {
      std::cerr << "Error: line '" << line << "' incorrect number of vertices" << std::endl;
      return false;
    }
    
    std::string label_u = s[0];
    std::string label_v = s[1];
    
    if (_idToNode.count(label_u) == 0)
    {
      Node u = _tree.addNode();
      _idToNode[label_u] = u;
      _nodeToId[u] = label_u;
    }
    if (_idToNode.count(label_v) == 0)
    {
      Node v = _tree.addNode();
      _idToNode[label_v] = v;
      _nodeToId[v] = label_v;
    }
    
    Node u = _idToNode[label_u];
    Node v = _idToNode[label_v];
    
    _tree.addArc(u, v);
  }
  
  _root = lemon::INVALID;
  for (NodeIt node(_tree); node != lemon::INVALID; ++node)
  {
    if (InArcIt(_tree, node) == lemon::INVALID)
    {
      if (_root != lemon::INVALID)
      {
        std::cerr << "Error: multiple root node '" << _nodeToId[node]
                  << "' and '" << _nodeToId[_root] << "'" << std::endl;
        return false;
      }
      _root = node;
    }
  }
  
  if (_idToNode.empty())
  {
    std::cerr << "Error: empty tree" << std::endl;
    return false;
  }
  
  init();
  
  if (!isValid())
  {
    std::cerr << "Error: tree is not a valid tree" << std::endl;
    return false;
  }
  
  return true;
}

NodeList BaseTree::pathFromRoot(const Digraph& T, Node u)
{
  NodeList result;
  while ((InArcIt(T, u) != lemon::INVALID))
  {
    result.push_front(u);
    u = T.source(InArcIt(T, u));
  }
  result.push_front(u);
  return result;
}

NodeList BaseTree::path(Node u, Node v) const
{
  NodeList result;
  
  if (!isAncestor(u, v))
  {
    return result;
  }
  
  while (u != v)
  {
    result.push_front(v);
    v = _tree.source(InArcIt(_tree, v));
  }
  result.push_front(u);
  
  return result;
}

void BaseTree::init()
{
  _leafSet.clear();
  for (NodeIt u(_tree); u != lemon::INVALID; ++u)
  {
    _leafSubset[u].clear();
    
    OutArcIt a(_tree, u);
    if (a == lemon::INVALID)
    {
      _leafSet.insert(u);
      _isLeaf[u] = true;
    }
  }
  
  lemon::bfs(_tree).distMap(_level).run(_root);
  
  initLeafSubset(_root);
    
  // Update edge set
  _edgeSet.clear();
  getParentalPairs();
  
  // Update ancestor set
  _transEdgeSet.clear();
  getAncestralPairs();
}

void BaseTree::initLeafSubset(Node u)
{
  NodeSet merged;
  for (OutArcIt a(_tree, u); a != lemon::INVALID; ++a)
  {
    Node v = _tree.target(a);
    initLeafSubset(v);
    merged.insert(_leafSubset[v].begin(), _leafSubset[v].end());
  }
  
  if (merged.empty())
  {
    // node is a leaf
    merged.insert(u);
  }
  _leafSubset[u] = merged;
}

bool BaseTree::readVertexLabeling(std::istream& in,
                                  const BaseTree& T,
                                  StringNodeMap& label)
{
  const Digraph& tree = T.tree();
  for (NodeIt u(tree); u != lemon::INVALID; ++u)
  {
    label[u] = "";
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
    
    Node u = T.getNodeByLabel(label_u);
    if (u == lemon::INVALID)
    {
      std::cerr << "Error: clone-tree vertex with label '" << label_u << "' does not exist" << std::endl;
      return false;
    }
    
    label[u] = label_s;
  }
  
  for (NodeIt u(tree); u != lemon::INVALID; ++u)
  {
    if (label[u].empty())
    {
      std::cerr << "Error: vertex '" << T.label(u) << "' left unlabeled" << std::endl;
      return false;
    }
  }
  
  return true;
}

bool BaseTree::readColorMap(std::istream& in,
                            StringToIntMap& colorMap)
{
  colorMap.clear();
  
  while (in.good())
  {
    std::string line;
    getline(in, line);
    
    if (line.empty())
      continue;
    
    StringVector s;
    boost::split(s, line, boost::is_any_of("\t "));
    
    std::string location = s[0];
    int color = boost::lexical_cast<int>(s[1]);
    
    if (colorMap.count(location) != 0)
    {
      std::cerr << "Error: location '" << location << "' is already assigned a color" << std::endl;
      return false;
    }

    colorMap[location] = color;
  }
  
  return true;
}

Node BaseTree::getLCA(const Digraph& T, const NodeSet& nodes)
{
  if (nodes.size() == 1)
  {
    return *nodes.begin();
  }
  
  NodeListVector allPaths;
  NodeListItVector iteratorList;
  for (Node node : nodes)
  {
    allPaths.push_back(pathFromRoot(T, node));
    iteratorList.push_back(allPaths.back().begin());
  }
  
  // all first elements in iteratorList are set to the root and hence are the same
  Node lca;
  bool same = true;
  while (same)
  {
    lca = *iteratorList.front();
    
    for (int i = 0; i < iteratorList.size(); ++i)
    {
      if (++iteratorList[i] == allPaths[i].end())
        same = false;
    }
    
    if (same)
    {
      Node first = *iteratorList.front();
      for (const NodeListIt& it : iteratorList)
      {
        if (*it != first)
        {
          same = false;
        }
      }
    }
  }
  
  return lca;
}

bool BaseTree::isValid() const
{
  if (_root == lemon::INVALID)
    return false;
  
  // check connectivity
  lemon::Bfs<Digraph> bfs(_tree);
  bfs.run(_root);
  
  for (NodeIt v(_tree); v != lemon::INVALID; ++v)
  {
    if (!bfs.reached(v))
      return false;
  }
  
  return true;
}

bool BaseTree::isConnected(const NodeSet& nodes) const
{
  Node lca = getLCA(nodes);
  
  if (nodes.count(lca) == 0)
  {
    return false;
  }
  
  for (Node u : nodes)
  {
    NodeList P = path(lca, u);
    for (Node v : P)
    {
      if (nodes.count(v) == 0)
      {
        return false;
      }
    }
  }
  
  return true;
}

void BaseTree::writeDOT(std::ostream& out) const
{
  out << "digraph T {" << std::endl;
  out << "\t{" << std::endl;
  out << "\t\trank=same" << std::endl;
  for (NodeIt u(_tree); u != lemon::INVALID; ++u)
  {
    if (_isLeaf[u])
    {
      out << "\t\t" << _tree.id(u) << " [label=\"" << _nodeToId[u] << "\"]" << std::endl;
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

