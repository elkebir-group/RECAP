/*
 * mcctsolverca.h
 *
 *  Created on: 12-dec-2018
 *      Author: N. Aguse
 */

#ifndef MCCTSOLVERCA_H
#define MCCTSOLVERCA_H

#include "mcctsolver.h"

class McctSolverCA: public McctSolver
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
  /// @param nrMaxRestarts Maximum number of restarts
  McctSolverCA(const InputInstance& scriptT,
               const std::string& rootMutation,
               const std::string& sep,
               int k,
               int timeLimit,
               int minClusterSize,
               int nrMaxRestarts);

  virtual void solve();
  
  virtual std::string getMethodName() const
  {
    char buf[1024];
    snprintf(buf, 1024, "CA-r%d_t%d", _nrMaxRestarts, _timeLimit);
    return buf;
  }
  
  virtual void writeSummaryHeader(std::ostream& out, bool newLine) const;
  
  virtual void writeSummary(std::ostream& out, bool newLine) const;

private:
  void init();
  
  void generateInitialClusteringAndSelection();
  
  /// Update clustering when consensus trees are fixed
  void updateClusteringAndSelection();
    
  /// Fill empty clusters with patient with worst score
  void fillEmptyClustersSingleton();
  
  /// Fill empty clusters by splitting worst cluster
  void fillEmptyClustersSplit();
  
private:
  const int _nrMaxRestarts;
  int _nrRestarts;
};

#endif // MCSCTSOLVERCA_H
