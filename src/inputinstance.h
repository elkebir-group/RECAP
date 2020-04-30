/*
 * inputinstance.h
 *
 *  Created on: 10-jun-2019
 *      Author: M. El-Kebir
 */

#ifndef INPUTINSTANCE_H
#define INPUTINSTANCE_H

#include "utils.h"
#include "clonetree.h"

class InputInstance
{
public:
  /// Constructor
  InputInstance();
  
  /// Return clone trees of the specified patient
  ///
  /// @param p Patient
  const CloneTreeVector& getPatientTrees(int p) const
  {
    assert(0 <= p && p < _scriptT.size());
    return _scriptT[p];
  }
  
  /// Return number of patients
  int getNrPatients() const
  {
    return _scriptT.size();
  }
  
  /// Return minimum number of trees per patient
  int getMinNrTrees() const
  {
    int res = std::numeric_limits<int>::max();

    for (const CloneTreeVector& ctv : _scriptT)
    {
      if (res > ctv.size())
      {
        res = ctv.size();
      }
    }

    return res;
  }
  
  /// Return maximum number of trees per patient
  int getMaxNrTrees() const
  {
    int res = 0;
    
    for (const CloneTreeVector& ctv : _scriptT)
    {
      if (res < ctv.size())
      {
        res = ctv.size();
      }
    }
    
    return res;
  }
  
  /// Return median number of trees per patient
  int getMedianNrTrees() const
  {
    IntVector nrTrees;
    
    for (const CloneTreeVector& ctv : _scriptT)
    {
      nrTrees.push_back(ctv.size());
    }
    
    std::sort(nrTrees.begin(), nrTrees.end());
    
    return nrTrees[getNrPatients() / 2];
  }
  
  /// Return total number of trees across all patients
  int getTotalNrTrees() const
  {
    int res = 0;
    
    for (const CloneTreeVector& ctv : _scriptT)
    {
      res += ctv.size();
    }
    
    return res;
  }
  
  /// Returns whether there are mutation clusters
  ///
  /// @param sep Separator
  bool hasMutationClusters(const std::string& sep) const;
  
  /// Return number of mutations
  ///
  /// @param sep Separator
  int getMutationSetSize(const std::string& sep) const;
  
  /// Return mutations
  ///
  /// @param sep Separator
  const StringSet& getMutationSet(const std::string& sep) const;
  
  /// Update mutations
  ///
  /// @param sep Separator
  const StringSet& updateMutationSet(const std::string& sep) const;
    
  /// Adds dummy mutations to trees
  ///
  /// @param sep Separator
  /// @param dummyNode Dummy node label
  void completePatientTrees(const std::string& sep, const std::string& dummyNode="dummy");
  
  /// Removes nodes from patient trees below cohort count threshold
  ///
  /// @param sep Separator
  /// @param t threshold
  void removeRareMutations(const std::string& sep, const int t);
  
  /// Exhausively enumerate all possible tree expansions
  ///
  /// @param sep Separator
  void expand(const std::string& sep);
  
  /// Randomly permute labels on clone trees (excluding root)
  void randPermutePerTree();
    
  /// Randomly permute labels on clone trees (excluding root)
  void randPermutePerPatient();

private:
  typedef std::vector<CloneTreeVector> CloneTreeMatrix;
  
private:
  /// Patient trees
  CloneTreeMatrix _scriptT;
  /// Patient to mutations
  mutable std::vector<StringSet> _patient2mutations;
  /// Patient to missing mutations
  mutable std::vector<StringSet> _patient2missingMutations;
  /// Set of mutations present among patients
  mutable StringSet _mutationSet;
  /// Dictionary of counts of mutations
  mutable std::unordered_map<std::string, int> _mutationCounts;
  
  friend std::ostream& operator<<(std::ostream& out, const InputInstance& scriptT);
  friend std::istream& operator>>(std::istream& in, InputInstance& scriptT);
};

std::ostream& operator<<(std::ostream& out, const InputInstance& scriptT);

std::istream& operator>>(std::istream& in, InputInstance& scriptT);

#endif // INPUTINSTANCE_H
