/*
 * mcctsolver.cpp
 *
 *  Created on: 12-dec-2018
 *      Author: N. Aguse
 */

#include "mcctsolver.h"

McctSolver::McctSolver(const InputInstance& scriptT,
                       const std::string& rootMutation,
                       const std::string& sep,
                       int k,
                       int timeLimit,
                       int minClusterSize)
  : _scriptT(scriptT)
  , _rootMutation(rootMutation)
  , _sep(sep)
  , _k(k)
  , _timeLimit(timeLimit)
  , _minClusterSize(minClusterSize)
  , _hasMutationClusters(scriptT.hasMutationClusters(sep))
  , _startTimePoint(std::chrono::high_resolution_clock::now())
  , _clusteringCost(-1)
  , _cluster2patients(k)
  , _cluster2consensus(k)
  , _cluster2cost(k)
  , _cluster2nrpatients(k)
  , _patient2cost()
  , _patient2cluster()
  , _selection()
  , _selectedExpandedTrees()
  , _consensusEdges()
  , _consensusTransEdges()
  , _id()
  , _mutationToLevelMap(k)
{
}

McctSolver::~McctSolver()
{
  clearConsensusTrees();
}

void McctSolver::run(McctSolver& solver,
                     const std::string& outputPrefix)
{
  solver.solve();
  
  // print output (need to do)
  {
    std::ofstream outFile;
    outFile.open(outputPrefix + "_summary.tsv");
    solver.writeSummaryHeader(outFile, true);
    solver.writeSummary(outFile, true);
  }
  
  {
    std::ofstream outFile;
    outFile.open(outputPrefix + ".solution.txt");
    solver.writeSolution(outFile);
  }
  
  solver.displayConsensusTrees();

  for (int clusterIdx = 0; clusterIdx < solver._k; ++clusterIdx)
  {
    char buf[1024];
    snprintf(buf, 1024, "%s_cluster%d.dot", outputPrefix.c_str(), clusterIdx);
    std::ofstream outFile(buf);
    solver._cluster2consensus[clusterIdx]->writeDOT(outFile,
                                                    solver._cluster2nrpatients[clusterIdx],
                                                    solver._cluster2cost[clusterIdx]);
    outFile.close();
  }
}

void McctSolver::displayConsensusTrees() const
{
  for (int i = 0; i < _k; ++i)
  {
    auto edges = _cluster2consensus[i]->getSelectedEdgeList();
    for (const auto& edge : edges)
    {
      std::cerr << edge.first << " " << edge.second << " " << _cluster2consensus[i]->getArcOccurrence(edge.first, edge.second) << std::endl;
    }
    std::cerr << std::endl << std::endl;
  }
}

void McctSolver::resetClusteringAndSelection()
{
  _selection = IntVector(_scriptT.getNrPatients(), -1);
  _patient2cluster = IntVector(_scriptT.getNrPatients(), -1);
  _patient2cost = DoubleVector(_scriptT.getNrPatients(), std::numeric_limits<double>::max());
  _cluster2patients = IntSetVector(_k);
}

void McctSolver::clearConsensusTrees()
{
  _mutationToLevelMap = std::vector<std::map<std::string, int>>(_k);
  for (auto it = _cluster2consensus.begin();
       it != _cluster2consensus.end();
       it++)
  {
    delete *it;
    *it = NULL;
  }
  _consensusTransEdges.clear();
  _consensusEdges.clear();
}

void McctSolver::generateParentChildGraphs()
{
  clearConsensusTrees();
  assert(_cluster2patients.size() == _k);
  assert(_selection.size() == _scriptT.getNrPatients());
  
  int clusterIdx = 0;
  for (const IntSet& patients : _cluster2patients)
  {
    std::vector<CloneTree> ctvcopy;
    
    for (int p : patients)
    {
      assert(0 <= p && p < _scriptT.getNrPatients());
      ctvcopy.push_back(_selectedExpandedTrees[p]);
    }
    _cluster2consensus[clusterIdx] = new ParentChildGraph(ctvcopy,
                                                          _scriptT.getMutationSet(_sep),
                                                          _rootMutation);
    
    ++clusterIdx;
  }
}

void McctSolver::computeConsensusTrees()
{
  assert(_cluster2patients.size() == _k);
  
  _consensusEdges = std::vector<StringPairSet>(_k);
  _consensusTransEdges = std::vector<StringPairSet>(_k);
  
  for (int s = 0; s < _k; ++s)
  {
    computeConsensusTree(s);
  }
}

void McctSolver::computeConsensusTree(const int s){
  
  _cluster2consensus[s]->maxWeightSpanningArborescence();
  
  _consensusEdges[s] = _cluster2consensus[s]->getSelectedEdgeList();
  _consensusTransEdges[s] = _cluster2consensus[s]->getSelectedTransEdgeList();
  const StringPairSet& edges = _cluster2consensus[s]->getSelectedEdgeList();
  
  
  // FINDING ROOT OF CONSENSUS TREE
  // Identify set of mutations
  std::unordered_set<std::string> mutations;
  for (auto elem : edges){
    mutations.insert(std::get<0>(elem));
    mutations.insert(std::get<1>(elem));
  }
  
  // Remove mutations that have a parent, leaving root
  for (auto elem : edges){
    mutations.erase(std::get<1>(elem));
  }
  
  // Clear consensus map
  _mutationToLevelMap[s].clear();
  
  if (mutations.size() != 0){
    assert(mutations.size()==1);
    std::string root(*mutations.begin());
  
    // UPDATING mutationToLevelMap RECURSIVELY
    _mutationToLevelMap[s][root] = 0;
    fillMutationToLevelMap(edges, s, root);
  }
}

void McctSolver::fillMutationToLevelMap(const StringPairSet& edges, const int cluster, const std::string& vertex)
{
  int currentDepth(_mutationToLevelMap[cluster][vertex]);
    
  // Loop over edges to identify children of current node
  for (const auto elem : edges){
      
    // If this edge involves our current vertex as a parent
    if (std::get<0>(elem) == vertex){
    
      // Identify child of current vertex
      std::string child(std::get<1>(elem));
            
      // Update depth of child in dictionary
      _mutationToLevelMap[cluster][child] = currentDepth+1;
            
      // Call recursively on child
      fillMutationToLevelMap(edges, cluster, child);
    }
  }
}

void McctSolver::updateClusteringCost()
{
  const int nrPatients = _scriptT.getNrPatients();
  _clusteringCost = 0;
  _patient2cost = DoubleVector(nrPatients, 0);
  
  assert(_cluster2patients.size() == _k);
  assert(_cluster2consensus.size() != 0);
  
  int clusterIdx = 0;
  for (const IntSet& patients : _cluster2patients)
  {
    _cluster2cost[clusterIdx] = 0;
    
    for (int p : patients)
    {
      assert(0 <= p && p < _scriptT.getNrPatients());
      
      double cost = _cluster2consensus[clusterIdx]->parentChildDistance(_selectedExpandedTrees[p]);
      _patient2cost[p] = cost;
      _cluster2cost[clusterIdx] += cost;
    }
    _cluster2nrpatients[clusterIdx] = patients.size();
    _clusteringCost += _cluster2cost[clusterIdx];
  
    ++clusterIdx;
  }
}

void McctSolver::setSolution(const IntVector& clustering,
                             const IntVector& selection,
                             const CloneTreeVector& selectedExpandedTrees,
                             const std::vector<StringPairSet>& consensusEdges,
                             const std::vector<StringPairSet>& consensusTransEdges)
{
  assert(clustering.size() == _scriptT.getNrPatients());
  assert(selection.size() == _scriptT.getNrPatients());
  _patient2cluster = clustering;
  _selection = selection;
 
  // update _cluster2trees
  const int n = _scriptT.getNrPatients();
  _cluster2patients = IntSetVector(_k);
  for (int p = 0; p < n; ++p)
  {
    _cluster2patients[clustering[p]].insert(p);
  }
  
  _selectedExpandedTrees = selectedExpandedTrees;
 
  generateParentChildGraphs();
  for (int j = 0; j < _k; ++j)
  {
    _cluster2consensus[j]->setSelectedEdgeList(consensusEdges[j]);
    _cluster2consensus[j]->setSelectedTransEdgeList(consensusTransEdges[j]);
  }
  
  _consensusEdges = consensusEdges;
  _consensusTransEdges = consensusTransEdges;
  
  updateClusteringCost();
}

void McctSolver::writeSelection(std::ostream& out) const
{
  for (const auto &e : _selection) out << e << "\n";
}

void McctSolver::writeClustering(std::ostream& out) const
{
  for (const auto &e : _patient2cluster) out << e << "\n";
}

void McctSolver::writeSolution(std::ostream& out) const
{
  out << _k << " # clusters" << std::endl;
  for (int i = 0; i < _k; ++i)
  {
    auto edges = _cluster2consensus[i]->getSelectedEdgeList();
    out << edges.size() << " # edges in consensus tree " << i << ", distance " << _cluster2cost[i] << std::endl;
    for (const auto& edge : edges)
    {
        out << edge.first << " " << edge.second << " " << _cluster2consensus[i]->getArcOccurrence(edge.first, edge.second) << std::endl;
    }
  }

  const int n = _scriptT.getNrPatients();
  assert(n == _selection.size() && n == _patient2cluster.size());
  out << n << " # patients" << std::endl;
  for (int i = 0; i < n; ++i)
  {
    out << _patient2cluster[i] << " # cluster index for patient " << i << ", distance " << _patient2cost[i] << std::endl;
    out << _selection[i] << " # selected input tree index for patient " << i << std::endl;
    const CloneTree& S_i = _selectedExpandedTrees[i];
    out << S_i.getEdgeList().size() << " # edges in selected input tree for patient " << i << std::endl;
    for (const auto& edge : S_i.getEdgeList())
    {
        out << edge.first << " " << edge.second << std::endl;
    }
  }
}

void McctSolver::writeSummaryHeader(std::ostream& out, bool newLine) const
{
  out << "method\ttime\tk\tcost";
  for (int s = 0; s < _k; ++s)
  {
    out << "\tcost-" << s;
  }
  out << "\ttrees";
  for (int s = 0; s < _k; ++s)
  {
    out << "\ttrees-" << s;
  }
  
  if (newLine)
  {
    out << std::endl;
  }
}

void McctSolver::writeSummary(std::ostream& out, bool newLine) const
{
  auto finish = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> diff = finish - _startTimePoint;
  double secondsElapsed = diff.count();
  
  out << getMethodName() << "\t" << secondsElapsed << "\t" << _k << "\t " << getClusteringCost();
  
  for (int s = 0; s < _k; ++s)
  {
    out << "\t" << _cluster2cost[s];
  }
  
  out << "\t" << _scriptT.getNrPatients();
  for (int s = 0; s < _k; ++s)
  {
    out << "\t" << _cluster2nrpatients[s];
  }
  
  if (newLine)
  {
    out << std::endl;
  }
}
