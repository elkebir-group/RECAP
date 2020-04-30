/*
 * mcctsolverca.cpp
 *
 *  Created on: 22-sep-2019
 *      Author: J. Kim
 */

#include "utils.h"
#include "mcctsolverca.h"
#include "inputinstance.h"
#include "parentchildexpansion.h"

McctSolverCA::McctSolverCA(const InputInstance& scriptT,
                           const std::string& rootMutation,
                           const std::string& sep,
                           int k,
                           int timeLimit,
                           int minClusterSize,
                           int nrMaxRestarts)
  : McctSolver(scriptT, rootMutation, sep, k, timeLimit, minClusterSize)
  , _nrMaxRestarts(nrMaxRestarts)
  , _nrRestarts(0)
{
  //assert(k <= ctv.size());
  init();
}

void McctSolverCA::writeSummaryHeader(std::ostream& out, bool newLine) const
{
  McctSolver::writeSummaryHeader(out, false);
  out << "\t" << "restarts" << std::endl;
}

void McctSolverCA::writeSummary(std::ostream& out, bool newLine) const
{
  McctSolver::writeSummary(out, false);
  out << "\t" << _nrRestarts << std::endl;
}

void McctSolverCA::init()
{

}

void McctSolverCA::updateClusteringAndSelection()
{
  assert(_cluster2patients.size() != 0);  
  resetClusteringAndSelection();
  
  const int nrPatients = _scriptT.getNrPatients();
  for (int s = 0; s < _k; ++s)
  {
    for (int i = 0; i < nrPatients; i++)
    {
      ParentChildExpansion::LookUpTableType lookUpTable;
      for (int j = 0; j < _scriptT.getPatientTrees(i).size(); j++)
      {
        if (_hasMutationClusters)
        {
          CloneTree exp = ParentChildExpansion(_scriptT.getPatientTrees(i)[j],
                                               _consensusEdges[s],
                                               _sep, lookUpTable).getExpandedCloneTree();
          double dist = _cluster2consensus[s]->parentChildDistance(exp);
          if (dist < _patient2cost[i])
          {
            _patient2cost[i] = dist;
            _patient2cluster[i] = s;
            _selection[i] = j;
            _selectedExpandedTrees[i] = exp;
          }
        }
        else
        {
          double dist = _cluster2consensus[s]->parentChildDistance(_scriptT.getPatientTrees(i)[j]);
          if (dist < _patient2cost[i])
          {
            _patient2cost[i] = dist;
            _patient2cluster[i] = s;
            _selection[i] = j;
            _selectedExpandedTrees[i] = _scriptT.getPatientTrees(i)[j];
          }
        }
      }
    }
  }
  
  // assign patients to clusters
  for (int i = 0; i < nrPatients; i++)
  {
    _cluster2patients[_patient2cluster[i]].insert(i);
  }
///
}

void McctSolverCA::fillEmptyClustersSingleton(){
    /// The following code ensures that all clusters achieve the minimum cluster size
    /// For each cluster not achieving the minimum size, we "steal" the closest patient
    /// The patient must be from a cluster that is above the minimum cluster size
        
  const int nrPatients = _scriptT.getNrPatients();
    
  /// For each cluster
  for (int s = 0; s < _k; ++s){
        
    /// While the cluster is still below the minimum size
    while (_cluster2patients[s].size() < _minClusterSize){
          
      /// We maintain the closest patient, the patient's current cluster assignment, and the patient's distance
      int bestPatient = -1;
      int bestSelect = -1;
      //int bestDist = maxDist;
      double bestDist = 0;
          
      /// Now we loop over all of the patients
      for (int i = 0; i < nrPatients; i++){
            
        /// We check if this patient is elligible to be stolen by making sure the patient isn't already in this cluster
        /// And checking the current size of the patient's cluster
        if ((_patient2cluster[i] != s) and (_cluster2patients[_patient2cluster[i]].size() > _minClusterSize)){
                
          int j = _selection[i];
          /// If so, we loop over ALL of the patient's trees, not just the current selected tree
          //for (int j = 0; j < _scriptT.getPatientTrees(i).size(); j++)
          //{
            /// Calculate the distance for this tree to consensus tree for current cluster
            //int dist = _cluster2consensus[s]->parentChildDistance(_scriptT.getPatientTrees(i)[j]);
            /// Calculate the distance for this tree to current consensus tree
            double dist = _patient2cost[i];
                  
            /// If the distance is better than current distance, we take this patient
            //if (dist < bestDist)
            if (dist >= bestDist){
              bestDist = dist;
              bestPatient = i;
              bestSelect = j;
            }
              //}
        }
      }

      assert(bestPatient != -1);
      assert(bestSelect != -1);

      // Now we need to update member variables w/r/t this selection
      int oldCluster(_patient2cluster[bestPatient]);
            
      // Change cluster sets of patients to reflect change
      _cluster2patients[oldCluster].erase(bestPatient);
      _cluster2patients[s].insert(bestPatient);
      _cluster2nrpatients[s] = 1;
            
      // Change patient cluster
      _patient2cluster[bestPatient] = s;
      _selection[bestPatient] = bestSelect;
            
      // Update consensus tree
      std::vector<CloneTree> ctvcopy;
      ctvcopy.push_back(_selectedExpandedTrees[bestPatient]);
      _cluster2consensus[s] = new ParentChildGraph(ctvcopy, _scriptT.getMutationSet(_sep), _rootMutation);
      computeConsensusTree(s);
      double dist = _cluster2consensus[s]->parentChildDistance(_selectedExpandedTrees[bestPatient]);
      assert(dist == 0);
      _patient2cost[bestPatient] = dist;
    }
  }
}

void McctSolverCA::fillEmptyClustersSplit(){
  
  /// For each cluster
  for (int s = 0; s < _k; ++s){
        
    /// If the cluster is empty
    if (_cluster2patients[s].size() == 0){
  
      /// The goal of this function is to split the worst cluster in half to seed the empty cluster
      /// To identify the worst cluster we need to update the current cost
      updateClusteringCost();
  
      /// We now iterate over all clusters to find the worst score
      int nrClusters(_cluster2cost.size());
      double worstCost = 0;
      int worstCluster = -1;
      for (int i = 0; i < nrClusters; i++){
    
        // if this cluster has the highest cost so far
        // and if there is more than one person we choose it
        if (_cluster2cost[i]>=worstCost and _cluster2patients[i].size()>1){
          worstCost = _cluster2cost[i];
          worstCluster = i;
        }
      }
      
      assert(worstCluster != -1);
      int worstSize = _cluster2patients[worstCluster].size();
      

      /// We will pick who to reassign based on sharing a common mutation
      StringSet presentMutationSet;
      std::unordered_map<std::string, int> presentMutationMap;
      for (auto p: _cluster2patients[worstCluster])
      {
        auto T(_selectedExpandedTrees[p]);
        for (NodeIt v(T.tree()); v != lemon::INVALID; ++v)
        {
          presentMutationSet.insert(T.label(v));
          if (presentMutationMap.count(T.label(v)) == 0){
            presentMutationMap[T.label(v)] = 1;
          }
          else{
            presentMutationMap[T.label(v)] = presentMutationMap[T.label(v)]+1;
          }
        }
      }
      
      /// Pick a random mutation
      std::vector<std::string> mv(presentMutationSet.begin(), presentMutationSet.end());
      std::shuffle(mv.begin(), mv.end(), g_rng);
      std::string m(mv.back());
      for (auto candidateM: presentMutationSet){
        if (presentMutationMap[candidateM] < worstSize and presentMutationMap[candidateM] >1){
          m = candidateM;
        }
      }
      
      /// Keep track of selected trees and patients
      std::vector<CloneTree> ctv;
      std::vector<int> splitPatients;
      int fillCount(0);
      
      /// Move all patients with this mutation
      /// Make sure to leave at least one patient
      for (auto p: _cluster2patients[worstCluster])
      {
        auto T(_selectedExpandedTrees[p]);
        for (NodeIt v(T.tree()); v != lemon::INVALID; ++v)
        {
          if (T.label(v) == m and fillCount < worstSize-1){
                        
            /// Add selected tree and patient
            ctv.push_back(_selectedExpandedTrees[p]);
            splitPatients.push_back(p);
            
            fillCount++;
          }
        }
      }
      
      /// Update dictionaries for the selected patients
      for (auto p: splitPatients){
        /// Change cluster sets of patients to reflect change
        _cluster2patients[worstCluster].erase(p);
        _cluster2patients[s].insert(p);
      
        /// Change patient cluster
        _patient2cluster[p] = s;
      }
      
      // Update consensus tree for this cluster
      // At worst it could be the same as before and the score would not change
      // At best it could improve
      _cluster2consensus[s] = new ParentChildGraph(ctv, _scriptT.getMutationSet(_sep), _rootMutation);
      computeConsensusTree(s);
      
      _cluster2nrpatients[s] = splitPatients.size();
      _cluster2nrpatients[worstCluster] = worstSize - splitPatients.size();
      
      /// Finally update distance to the consensus tree for this cluster
      for (auto p: splitPatients){
      double dist = _cluster2consensus[s]->parentChildDistance(_selectedExpandedTrees[p]);
      _patient2cost[p] = dist;
      }
    }
  }
}

void McctSolverCA::generateInitialClusteringAndSelection()
{
  assert(_k >= 1);
  resetClusteringAndSelection();
  
  /// Our goal is to populate the following two vectors
  const int nrPatients = _scriptT.getNrPatients();
  IntVector selection(nrPatients, -1);
  IntVector clustering(nrPatients,-1);
    
  /// SELECTION
  /// First, we randomly select a tree for each patient
  for (int i = 0; i < nrPatients; i++)
  {
    const CloneTreeVector& scriptT_i = _scriptT.getPatientTrees(i);
    std::uniform_int_distribution<> unif(0, scriptT_i.size()-1);
    selection[i] = unif(g_rng);
  }

  /// CLUSTERING
  if (_minClusterSize > 0)
  {
    /// Now, for each cluster we randomly select the minimum number of patients to be assigned to the cluster
    /// We start by initializing a vector of patient indices and randomly shuffling the order
    IntVector v(nrPatients);
    std::iota(v.begin(), v.end(), 0);
    std::shuffle(v.begin(), v.end(), g_rng);

    /// For every cluster
    for (int i = 0; i < _k; i++){
        int fillCount(0);
        /// While the cluster does not meet minimum size requirement
        while (fillCount < _minClusterSize){
            /// Take last patient and assign to this cluster
            int p(v.back());
            clustering[p] = i;
            v.pop_back();
            fillCount++;
        }
    }

    /// Now that the minimum number of patients is satisfied, we assign the remaing patients randomly
    std::uniform_int_distribution<> unif(0, _k - 1);
    /// While some patients have not been assigned
    while(v.size() > 0){
        /// Randomly select a cluster
        int c(unif(g_rng));
        /// Take last patient and assign to this cluster
        int p(v.back());
        clustering[p] = c;
        v.pop_back();
    }
  }
  else
  {
    if (_k == 1)
    {
      clustering = IntVector(nrPatients, 0);
    }
    else
    {
      IntVector v(nrPatients);
      std::iota(v.begin(), v.end(), 0);
      std::shuffle(v.begin(), v.end(), g_rng);
      
      // pick k breakpoints
      IntVector breakPoints(nrPatients-1);
      std::iota(breakPoints.begin(), breakPoints.end(), 1);
      std::shuffle(breakPoints.begin(), breakPoints.end(), g_rng);
      breakPoints.erase(breakPoints.begin() + _k - 1, breakPoints.end());

      std::sort(breakPoints.begin(), breakPoints.end());
      
      for (int clusterIdx = 0; clusterIdx < _k; ++clusterIdx)
      {
        if (0 < clusterIdx && clusterIdx < _k - 1)
        {
          for (int i = breakPoints[clusterIdx-1]; i < breakPoints[clusterIdx]; ++i)
          {
            clustering[v[i]] = clusterIdx;
          }
        }
        else if (clusterIdx == 0)
        {
          for (int i = 0; i < breakPoints[clusterIdx]; ++i)
          {
            clustering[v[i]] = 0;
          }
        }
        else if (clusterIdx == _k - 1)
        {
          for (int i = breakPoints[clusterIdx-1]; i < nrPatients; ++i)
          {
            clustering[v[i]] = clusterIdx;
          }
        }
      }
    }
  }


  /// Check that everything has been set
  for (int i = 0; i < nrPatients; i++){
      assert(selection[i] != -1);
      assert(clustering[i] != -1);
  }
  
  // Randomly expand selected trees
  _selectedExpandedTrees.clear();
  for (int i = 0; i < nrPatients; i++)
  {
    const CloneTree& T = _scriptT.getPatientTrees(i)[selection[i]];
    CloneTree exp_T(T, _sep);
    
    _selectedExpandedTrees.push_back(exp_T);
  }

//  std::cout << _selectedExpandedTrees << std::endl;
  
  StringPairSet initialConsensusEdges;
//  for (const std::string& mut : _scriptT.getMutationSet(_sep))
//  {
//    initialConsensusEdges.insert(StringPair("MISSING", mut));
//  }
//  initialConsensusEdges.insert(StringPair(_rootMutation, "MISSING"));
  std::vector<StringPairSet> emptySet(_k, initialConsensusEdges);
  
  setSolution(clustering,
              selection,
              _selectedExpandedTrees,
              emptySet,
              emptySet);
}

void McctSolverCA::solve()
{
  
  if (_k*_minClusterSize > _scriptT.getNrPatients())
  {
    throw std::runtime_error("Error: There are not enough patients to satisfy the minimum cluster size for the specified number of clusters.");
  }
    
  double minCost = std::numeric_limits<double>::max();
  IntVector bestSelection;
  IntVector bestClustering;
  CloneTreeVector bestSelectedExpandedTrees;
  std::vector<StringPairSet> bestConsensusEdges, bestConsensusTransEdges;
  
  _nrRestarts = 1;
  SecondsType timeLimit(_timeLimit);
  
  while (true)
  {
    generateInitialClusteringAndSelection();
    checkMinClusters();
    
    //int currCost = -1;
    double currCost = std::numeric_limits<double>::max();
    while (true)
    {
      // 2a. solve the maximum weight arborescence problem for each graph
      computeConsensusTrees();
      
      // 2b. update cost
      updateClusteringCost();
      std::cerr << "Restart " << _nrRestarts << " -- identified R_1, ..., R_k --  distance " << _clusteringCost << " :";
      for (double c : _cluster2cost)
      {
        std::cerr << " " << c;
      }
      std::cerr << std::endl;
      
      // 3. update clustering and tree selection
      updateClusteringAndSelection();
      
      // 3.5. Deal with empty clusters
      //fillEmptyClustersSingleton();
      fillEmptyClustersSplit();
      checkMinClusters();
        
      // 4. update cost
      updateClusteringCost();
      std::cerr << "Restart " << _nrRestarts << " -- identified S_1, ..., S_n and sigma --  distance " << _clusteringCost << " :";
      for (double c : _cluster2cost)
      {
        std::cerr << " " << c;
      }
      std::cerr << std::endl;
 
      double newCost = getClusteringCost();
      
      assert(!g_tol.less(currCost, newCost));
      if (!g_tol.different(currCost, newCost))
        break;
      currCost = newCost;
      
      // 1. generate parent-child graphs for each cluster using tree selection
      generateParentChildGraphs();
    }
    
    if (currCost < minCost)
    {
      minCost = currCost;
      bestClustering = _patient2cluster;
      bestSelection = _selection;
      bestSelectedExpandedTrees = _selectedExpandedTrees;
      bestConsensusEdges = _consensusEdges;
      bestConsensusTransEdges = _consensusTransEdges;
    }
    
    if (_nrMaxRestarts > 0 && _nrRestarts == _nrMaxRestarts)
    {
      break;
    }
    _nrRestarts++;
    
    if (_nrMaxRestarts <= 0)
    {
      auto finish = std::chrono::high_resolution_clock::now();
      if (finish - _startTimePoint >= timeLimit)
      {
        break;
      }
    }
  }

  setSolution(bestClustering,
              bestSelection,
              bestSelectedExpandedTrees,
              bestConsensusEdges,
              bestConsensusTransEdges);
}
