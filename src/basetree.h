/*
 * basetree.h
 *
 *  Created on: 13-oct-2016
 *      Author: M. El-Kebir
 */

#ifndef BASETREE_H
#define BASETREE_H

#include "utils.h"

/// This class models a tree whose nodes are labeled by unique identifiers
class BaseTree
{
public:
  /// Default constructor
  BaseTree();
  
  /// Constructor
  ///
  /// @param idSet Set of node identifiers
  /// @param rootId Identifier of root node
  /// @param pruferSeq Prufer sequence
  BaseTree(const StringSet& idSet,
           const std::string& rootId,
           const StringVector& pruferSeq);
  
  /// Constructor
  ///
  /// @param T Directed graph
  /// @param root Root node
  /// @param id Node identifier
  BaseTree(const Digraph& T,
           Node root,
           const StringNodeMap& id);
  
  /// Copy constructor
  BaseTree(const BaseTree& other);
  
  /// Destructor
  virtual ~BaseTree()
  {
  }
  
  /// Deserialize
  ///
  /// @param in Input stream
  virtual bool read(std::istream& in);
  
  /// Serialize
  ///
  /// @param out Output stream
  void write(std::ostream& out) const;
  
  /// Print tree in DOT format
  ///
  /// @param out Output stream
  virtual void writeDOT(std::ostream& out) const;
  
  /// Return number of locations
  virtual int getNrLocations() const
  {
    return 0;
  }
  
  /// Return location labels
  virtual StringSet getLocations() const
  {
    return StringSet();
  }
  
  /// Return underlying LEMON tree
  const Digraph& tree() const
  {
    return _tree;
  }
  
  /// Decide whether two nodes occur in distinct branches
  ///
  /// @param u Node
  /// @param v Node
  bool areIncomparable(Node u, Node v) const
  {
    return !isAncestor(u, v) && !isAncestor(v, u);
  }
  
  /// Decide whether one node is an ancestor (or self) of another
  ///
  /// @param u Node
  /// @param v Node
  bool isAncestor(Node u, Node v) const
  {
    while (v != _root)
    {
      if (u == v)
        return true;
      else
        v = _tree.source(InArcIt(_tree, v));
    }
    
    return u == v;
  }
  
  static bool isAncestor(const Digraph& tree,
                         Node u,
                         Node v)
  {
    while (v != lemon::INVALID)
    {
      if (u == v)
      {
        return true;
      }
      else
      {
        InArcIt a(tree, v);
        if (a == lemon::INVALID)
        {
          v = lemon::INVALID;
        }
        else
        {
          v = tree.source(a);
        }
      }
    }
    
    return u == v;
  }

  /// Return parent node. If u is the root lemon::INVALID is returned.
  ///
  /// @param u Node
  Node parent(Node u) const
  {
    if (u == _root)
    {
      return lemon::INVALID;
    }
    else
    {
      return _tree.source(InArcIt(_tree, u));
    }
  }
  
  /// Return the root node
  Node root() const
  {
    return _root;
  }
  
  /// Return the node identifier map
  const StringNodeMap& getLabelMap() const
  {
    return _nodeToId;
  }
  
  /// Return the identifier of a node
  ///
  /// @param u Node
  const std::string& label(Node u) const
  {
    return _nodeToId[u];
  }
  
  /// Decide whether a node is a leaf
  ///
  /// @param u Node
  bool isLeaf(Node u) const
  {
    return _isLeaf[u];
  }
  
  /// Return the BFS level of a node. The root node has BFS level 0.
  ///
  /// @param u Node
  int level(Node u) const
  {
    return _level[u];
  }
  
  /// Return the set of leaves that occur in the subtree rooted at the provided node
  ///
  /// @param u Node
  const NodeSet& leafSubset(Node u) const
  {
    return _leafSubset[u];
  }
  
  /// Return the leaf set
  const NodeSet& leafSet() const
  {
    return _leafSet;
  }
  
  /// Return the LCA of the provided set of nodes
  ///
  /// @param T Tree
  /// @param nodes Node set
  static Node getLCA(const Digraph& T,
                     const NodeSet& nodes);
  
  /// Return the LCA of the provided set of nodes
  ///
  /// @param nodes Node set
  Node getLCA(const NodeSet& nodes) const
  {
    return getLCA(_tree, nodes);
  }
  
  /// Decides whether the given set of nodes is connected
  ///
  /// @param nodes Node set
  bool isConnected(const NodeSet& nodes) const;
  
  /// Initialize auxilliary data structures
  /// (_isLeaf, _leafSet, _leafSubset, _level)
  void init();
  
  /// Return a node corresponding to the given identifier.
  /// If there is no node with the given identifier, lemon::INVALID is returned.
  ///
  /// @param lbl Identifier
  Node getNodeByLabel(const std::string& lbl) const
  {
    auto it = _idToNode.find(lbl);
    if (it == _idToNode.end())
    {
      return lemon::INVALID;
    }
    else
    {
      return it->second;
    }
  }
  
  /// Return the unique path from the root to u
  ///
  /// @param u Node
  NodeList pathFromRoot(Node u) const
  {
    return pathFromRoot(_tree, u);
  }
  
  /// Return the unique path from the root to u
  ///
  /// @param T Tree
  /// @param u Node
  static NodeList pathFromRoot(const Digraph& T,
                               Node u);
  
  /// Return the unique path from u to v
  NodeList path(Node u, Node v) const;
  
  /// Read vertex labeling
  ///
  /// @param in Input stream
  /// @param T Tree
  /// @param label Output labeling
  static bool readVertexLabeling(std::istream& in,
                                 const BaseTree& T,
                                 StringNodeMap& label);
  
  /// Read color map
  ///
  /// @param in Input stream
  /// @param colorMap Output color map
  static bool readColorMap(std::istream& in,
                           StringToIntMap& colorMap);
  
  /// Return a generated color map based on the locations present in the current tree
  StringToIntMap generateColorMap() const
  {
    StringSet Sigma = getLocations();
    
    StringToIntMap colorMap;
    int idx = 0;
    for (const std::string& s : Sigma)
    {
      colorMap[s] = ++idx;
    }
    
    return colorMap;
  }
  
  /// Return the edge set
  const StringPairSet& getEdgeSet() const
  {
    return getParentalPairs();
  }
  
  /// Return the edge list
  StringPairList getEdgeList() const
  {
    StringPairList res;
    for (ArcIt a(_tree); a != lemon::INVALID; ++a)
    {
      Node u = _tree.source(a);
      Node v = _tree.target(a);
      
      res.push_back(StringPair(_nodeToId[u], _nodeToId[v]));
    }
    return res;
  }
  
  /// Return the edge list
  std::string getEdgeListString() const
  {
    std::string res;
    for (ArcIt a(_tree); a != lemon::INVALID; ++a)
    {
      Node u = _tree.source(a);
      Node v = _tree.target(a);
      
      if (!res.empty())
        res += " ; ";
      
      res += "(" + _nodeToId[u] + "," + _nodeToId[v] + ")";
    }
    return res;
  }
  
  /// Return node id map
  const StringNodeMap& getIdMap() const
  {
    return _nodeToId;
  }
  
  /// Return Prufer sequence
  StringVector getPruferSequence() const;
  
  StringSet getIdSet() const
  {
    StringSet IDs;
    
    for (NodeIt v(_tree); v != lemon::INVALID; ++v)
    {
      IDs.insert(_nodeToId[v]);
    }
    
    assert(IDs.size() == lemon::countNodes(_tree));
    return IDs;
  }
  
  /// Return set of ordered pairs of ancestral vertices
  const StringPairSet& getAncestralPairs() const
  {
    if (_transEdgeSet.empty())
    {
      for (NodeIt u(_tree); u != lemon::INVALID; ++u)
      {
        for (NodeIt v(_tree); v != lemon::INVALID; ++v)
        {
          if (u == v) continue;
          if (isAncestor(u, v))
          {
            assert(_transEdgeSet.count(StringPair(_nodeToId[u], _nodeToId[v])) == 0);
            _transEdgeSet.insert(StringPair(_nodeToId[u], _nodeToId[v]));
          }
        }
      }
    }
    return _transEdgeSet;
  }
  
  /// Return set of ordered pairs of ancestral vertices
  StringPairSet getAncestralPairsNoRoot() const
  {
    StringPairSet res;
    for (NodeIt u(_tree); u != lemon::INVALID; ++u)
    {
      if (u == _root || parent(u) == _root) continue;
      for (NodeIt v(_tree); v != lemon::INVALID; ++v)
      {
        if (u == v) continue;
        if (OutArcIt(_tree, v) == lemon::INVALID) continue;
        if (isAncestor(u, v))
        {
          assert(res.count(StringPair(_nodeToId[u], _nodeToId[v])) == 0);
          res.insert(StringPair(_nodeToId[u], _nodeToId[v]));
        }
      }
    }
    return res;
  }
  
  /// Return set of directed edges
  const StringPairSet& getParentalPairs() const
  {
    if (_edgeSet.empty())
    {
      for (ArcIt a(_tree); a != lemon::INVALID; ++a)
      {
        Node u = _tree.source(a);
        Node v = _tree.target(a);
        
        assert(_edgeSet.count(StringPair(_nodeToId[u], _nodeToId[v])) == 0);
        _edgeSet.insert(StringPair(_nodeToId[u], _nodeToId[v]));
      }
    }
    return _edgeSet;
  }

  /// Return set of unordered pairs of incomparable vertices
  StringPairSet getIncomparablePairs() const
  {
    StringPairSet res;
    for (NodeIt u(_tree); u != lemon::INVALID; ++u)
    {
      for (NodeIt v(_tree); v != lemon::INVALID; ++v)
      {
        if (_nodeToId[u]  >= _nodeToId[v]) continue;
        if (areIncomparable(u, v))
        {
          assert(res.count(StringPair(_nodeToId[u], _nodeToId[v])) == 0);
          res.insert(StringPair(_nodeToId[u], _nodeToId[v]));
        }
      }
    }
    return res;
  }
  
  /// Return recall
  ///
  /// @param inferredSet Set of inferred pairs
  /// @param trueSet Set of true pairs
  static double recall(const StringPairSet& inferredSet,
                       const StringPairSet& trueSet)
  {
    StringPairSet intersection;
    std::set_intersection(inferredSet.begin(), inferredSet.end(),
                          trueSet.begin(), trueSet.end(),
                          std::inserter(intersection, intersection.begin()));
    
    return (double) intersection.size() / (double) trueSet.size();
  }
  
protected:
  /// Tree
  Digraph _tree;
  /// Root
  Node _root;
  /// Node label
  StringNodeMap _nodeToId;
  /// Get node by label
  StringToNodeMap _idToNode;
  /// Indicates whether a node is a leaf
  BoolNodeMap _isLeaf;
  /// Leaf set
  NodeSet _leafSet;
  /// Leaf set of subtree rooted at a node
  NodeNodeSetMap _leafSubset;
  /// BFS level of a node
  IntNodeMap _level;
  /// Edge list
  mutable StringPairSet _edgeSet;
  /// Transitive closure of edge list
  mutable StringPairSet _transEdgeSet;
  
  /// Initialize _leafSubset
  void initLeafSubset(Node node);
  
  /// Decide whether the tree is valid
  virtual bool isValid() const;
  
  /// Root provided undirected tree
  ///
  /// @param T Undirected tree
  /// @param TtoId Identifier node map
  /// @param v Node of undirected tree
  /// @param vv Corresponding node rooted/directed tree
  void rootUndirectedTree(const Graph& T,
                          const Graph::NodeMap<std::string>& TtoId,
                          Graph::Node& v,
                          Node vv);
};

#endif // BASETREE_H
