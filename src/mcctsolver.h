/*
 * mcctsolver.h
 *
 *  Created on: 12-dec-2018
 *      Author: N. Aguse
 */

#ifndef MCCTSOLVER_H
#define MCCTSOLVER_H

#include "inputinstance.h"
#include "clonetree.h"
#include "parentchildgraph.h"

class McctSolver
{
public:
  /// Constructor
  ///
  /// @param scriptT Input instance
  /// @param rootMutation Root mutation
  /// @param sep Mutation cluster separator
  /// @param k Number of clusters
  /// @param timeLimit Time limit in seconds
  /// @param minClusterSize Minimum cluster size
  McctSolver(const InputInstance& scriptT,
             const std::string& rootMutation,
             const std::string& sep,
             int k,
             int timeLimit,
             int minClusterSize);
  
  /// Destructor
  virtual ~McctSolver();
  
  /// Construct parent-child graphs given patient clustering and tree selection
  void generateParentChildGraphs();
  
  /// Generate consensus tree in each parent-child graph
  ///
  /// @param ancestorDescendantDistance Use ancestor-descendant distance
  void computeConsensusTrees();
    
  /// Generate consensus tree in each parent-child graph
  ///
  /// @param ancestorDescendantDistance Use ancestor-descendant distance
  void computeConsensusTree(const int s);
  
  /// Print consensus trees to standard out
  void displayConsensusTrees() const;
 
  /// Write patient clustering
  ///
  /// @param out Output stream
  void writeClustering(std::ostream& out) const;
  
  /// Write tree selection
  ///
  /// @param out Output stream
  void writeSelection(std::ostream& out)  const;
  
  /// Cost of consensus trees, patient clustering and tree selection
  double getClusteringCost() const
  {
    return _clusteringCost;
  }
  
  /// Patient clustering
  const IntVector& getClustering() const
  {
    return _patient2cluster;
  }

  /// Tree selection
  const IntVector& getSelection() const
  {
    return _selection;
  }
  
  /// Set patient clustering and tree selection and compute optimal consensus trees
  ///
  /// @param clustering Patient clustering
  /// @param selection Tree selection
  /// @param selectedExpandedTrees Selected expanded trees
  /// @param consensusEdges Consensus tree edges
  /// @param consensusTransEdges Consensus tree transitive edges
  /// @param ancestorDescendantDistance Use ancestor-descendant distance
  void setSolution(const IntVector& clustering,
                   const IntVector& selection,
                   const CloneTreeVector& selectedExpandedTrees,
                   const std::vector<StringPairSet>& consensusEdges,
                   const std::vector<StringPairSet>& consensusTransEdges);

  /// Clear current consensus trees
  void clearConsensusTrees();

  /// Get cluster-specific distances
  const DoubleVector& getCluster2Cost() const
  {
    return _cluster2cost;
  }
  
  /// Get number of patients per cluster
  const IntVector& getCluster2NrPatients() const
  {
    return _cluster2nrpatients;
  }

  /// Get number of patients
  int getNrPatients() const
  {
    return _scriptT.getNrPatients();
  }
  
  /// Get parent-child graph and consensus tree
  ///
  /// @param j Cluster index
  const ParentChildGraph* getConsensus(int j) const
  {
    assert(0 <= j && j < _k);
    return _cluster2consensus[j];
  }
  
  /// Get number of patients
  ///
  /// @param j Cluster index
  int getClusterNrPatients(int j) const
  {
    assert(0 <= j && j < _k);
    return _cluster2nrpatients[j];
  }
  
  /// Get cluster distance
  ///
  /// @param j Cluster index
  int getClusterCost(int j) const
  {
    assert(0 <= j && j < _k);
    return _cluster2cost[j];
  }
  
  /// Get total cost
  double getTotalCost() const
  {
    double total=0;
    
    for(int i=0; i < _cluster2cost.size(); i++)
    {
      total += _cluster2cost[i];
    }
    return total;
  }
  
  /// Check if minimum cluster size constraint is satisfied
  bool checkMinClusters()
  {
    bool res = true;
    for(int i=0; i < _cluster2nrpatients.size(); i++)
    {
      res &= _cluster2nrpatients[i] >= _minClusterSize;
      assert(_cluster2nrpatients[i] >= _minClusterSize);
    }
    return res;
  }
  
  /// Recursively fill out mutationToLevelMap
  void fillMutationToLevelMap(const StringPairSet& edges,
                              const int cluster,
                              const std::string& vertex);
    
  /// Write solution (consensus trees, patient clustering and tree selection)
  ///
  /// @param out Output stream
  virtual void writeSolution(std::ostream& out) const;
  
  /// Write header for summary table
  ///
  /// @param out Output stream
  /// @param newLine Terminate line with new line character
  virtual void writeSummaryHeader(std::ostream& out, bool newLine) const;
  
  /// Write row in summary table
  ///
  /// @param out Output stream
  /// @param newLine Terminate line with new line character
  virtual void writeSummary(std::ostream& out, bool newLine) const;

  /// Get method name
  virtual std::string getMethodName() const = 0;
  
  /// Solve
  ///
  /// @param minClusterSize Minimum number of patients per cluster
  /// @param ancestorDescendantDistance Use ancestor-descendant distance
  virtual void solve() = 0;
  
  /// Run solver
  ///
  /// @param solver Solver
  /// @param outputPrefix Prefix used in filename for output files
  static void run(McctSolver& solver,
                  const std::string& outputPrefix);
  
protected:
  /// Update clustering cost using current patient clustering, tree selection and consensus trees
  void updateClusteringCost();
  
  /// Clear patient clustering and tree selection
  void resetClusteringAndSelection();
  
protected:
  /// Input trees
  const InputInstance& _scriptT;
  /// Root mutation
  const std::string& _rootMutation;
  /// Mutation cluster separator
  const std::string& _sep;
  /// Number of clusters
  const int _k;
  /// Timelimit for the algorithm to run
  const int _timeLimit;
  /// Minimum number of patients per cluster
  const int _minClusterSize;
  ///  Indicates presence of mutations clutsers
  const bool _hasMutationClusters;
  /// Starting time point
  const TimePointType _startTimePoint;
  /// Sum of distances between each input tree and its corresponding consensus tree
  double _clusteringCost;
  /// Set of patient indices for each cluster
  IntSetVector _cluster2patients;
  /// Consensus tree of each cluster
  std::vector<ParentChildGraph*> _cluster2consensus;
  /// Total distance of the trees to the consensus tree in each cluster
  DoubleVector _cluster2cost;
  /// Number of trees in each cluster
  IntVector _cluster2nrpatients;
  /// Distance of selected patient tree to its consensus tree
  DoubleVector _patient2cost;
  /// Patient index to cluster index
  IntVector _patient2cluster;
  /// Selected tree index for each patient
  IntVector _selection;
  /// Selected expanded trees for each patient
  CloneTreeVector _selectedExpandedTrees;
  /// Consensus tree edges
  std::vector<StringPairSet> _consensusEdges;
  /// Consensus tree transitive edges
  std::vector<StringPairSet> _consensusTransEdges;
  /// Output prefix
  std::string _id;
  /// Assigns level to every mutation in consensus tree
  std::vector<std::map<std::string, int>> _mutationToLevelMap;
};

#endif // MCCTSOLVER_H
