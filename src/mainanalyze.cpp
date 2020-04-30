/*
 * mainanalyze.cpp
 *
 *  Created on: 20-oct-2019
 *      Author: M. El-Kebir
 */

#include <lemon/arg_parser.h>
#include <lemon/matching.h>
#include "utils.h"
#include "solution.h"
#include "inputinstance.h"

double symmDiffDist(const StringPairSet& set1, const StringPairSet& set2)
{
  StringPairSet resSet;
  std::set_symmetric_difference(set1.begin(), set1.end(),
                                set2.begin(), set2.end(),
                                std::inserter(resSet, resSet.begin()));
  
  return double(resSet.size()) / double(set1.size() + set2.size());
}

std::pair<double, double> computeSelectedDistance(const Solution& solInferred,
                                                  const Solution& solTrue)
{
  if (solInferred.getNrSelectedTrees() != solTrue.getNrSelectedTrees())
  {
    throw std::runtime_error("Error: varying number of trees");
  }
  
  const int n = solInferred.getNrSelectedTrees();
  
  double distSelectedPC = 0;
  double distSelectedAD = 0;
  for (int i = 0; i < n; ++i)
  {
    {
      // REMOVING DUMMY NODES
      StringPairSet transEdgesInferred = Solution::removeDummyNodes(solInferred.getSelectedTree(i).getAncestralPairs());
      StringPairSet transEdgesSolution = Solution::removeDummyNodes(solTrue.getSelectedTree(i).getAncestralPairs());
      distSelectedAD += symmDiffDist(transEdgesInferred, transEdgesSolution);
    }
    {
      // REMOVING DUMMY NODES
      StringPairSet edgesInferred = Solution::removeDummyNodes(solInferred.getSelectedTree(i).getParentalPairs());
      StringPairSet edgesSolution = Solution::removeDummyNodes(solTrue.getSelectedTree(i).getParentalPairs());
      
      distSelectedPC += symmDiffDist(edgesInferred, edgesSolution);
    }
  }
  
  return std::make_pair(distSelectedPC / n, distSelectedAD / n);
}

std::pair<double, double> computeConsensusDistance(const Solution& solInferred,
                                                   const Solution& solTrue)
{
  const int inferred_k = solInferred.getNrClusters();
  const int true_k = solTrue.getNrClusters();
  
  // 1. Compute distances of consensus trees
  Graph G_consensus;
  Graph::EdgeMap<double> distMapConsensusPC(G_consensus);
  Graph::EdgeMap<double> distMapConsensusAD(G_consensus);
  std::vector<Graph::Node> nodesConsensusInferred;
  std::vector<Graph::Node> nodesConsensusSolution;
  for (int i = 0; i < inferred_k; ++i)
  {
    nodesConsensusInferred.push_back(G_consensus.addNode());
  }
  for (int i = inferred_k; i < true_k; ++i)
  {
    nodesConsensusInferred.push_back(G_consensus.addNode());
  }
  
  for (int i = 0; i < true_k; ++i)
  {
    nodesConsensusSolution.push_back(G_consensus.addNode());
  }
  for (int i = true_k; i < inferred_k; ++i)
  {
    nodesConsensusSolution.push_back(G_consensus.addNode());
  }
  
  for (int i1 = 0; i1 < inferred_k; ++i1)
  {
    StringPairSet tmp = solInferred.getConsensusTree(i1).getAncestralPairs();
    StringPairSet tmp2 = solInferred.getConsensusTree(i1).getParentalPairs();
    
    // REMOVING DUMMY NODES
    StringPairSet transEdgesInferred = Solution::removeDummyNodes(tmp);
    StringPairSet edgesInferred = Solution::removeDummyNodes(tmp2);

    for (int i2 = 0; i2 < true_k; ++i2)
    {
      Graph::Edge edge = G_consensus.addEdge(nodesConsensusInferred[i1], nodesConsensusSolution[i2]);
      {
        const StringPairSet& transEdgesSolution = solTrue.getConsensusTree(i2).getAncestralPairs();
        
        distMapConsensusAD[edge] = 1 - symmDiffDist(transEdgesInferred, transEdgesSolution);
      }
      {
        const StringPairSet& edgesSolution = solTrue.getConsensusTree(i2).getParentalPairs();
        
        distMapConsensusPC[edge] = 1 - symmDiffDist(edgesInferred, edgesSolution);
      }
    }
  }
  
  for (int i1 = inferred_k; i1 < true_k; ++i1)
  {
    for (int i2 = 0; i2 < true_k; ++i2)
    {
      Graph::Edge edge = G_consensus.addEdge(nodesConsensusInferred[i1], nodesConsensusSolution[i2]);
      distMapConsensusAD[edge] = 0;
      distMapConsensusPC[edge] = 0;
    }
  }
  
  for (int i2 = true_k; i2 < inferred_k; ++i2)
  {
    for (int i1 = 0; i1 < inferred_k; ++i1)
    {
      Graph::Edge edge = G_consensus.addEdge(nodesConsensusInferred[i1], nodesConsensusSolution[i2]);
      distMapConsensusAD[edge] = 0;
      distMapConsensusPC[edge] = 0;
    }
  }
  
  double distConsensusPC;
  double distConsensusAD;
  
  {
    lemon::MaxWeightedPerfectMatching<Graph, Graph::EdgeMap<double>> mwmPC(G_consensus, distMapConsensusPC);
    mwmPC.run();
    distConsensusPC = std::max(inferred_k, true_k) - mwmPC.matchingWeight();
  }
  {
    lemon::MaxWeightedPerfectMatching<Graph, Graph::EdgeMap<double>> mwmAD(G_consensus, distMapConsensusAD);
    mwmAD.run();
    distConsensusAD = std::max(inferred_k, true_k) - mwmAD.matchingWeight();
  }
  
  return std::make_pair(distConsensusPC / std::max(inferred_k, true_k), distConsensusAD / std::max(inferred_k, true_k));
}

void getClusteringPairs(const IntVector& clustering,
                        IntPairSet& positives,
                        IntPairSet& negatives)
{
  const int n = clustering.size();
  for (int i = 0; i < n; ++i)
  {
    for (int j = i + 1; j < n; ++j)
    {
      if (clustering[i] == clustering[j])
      {
        positives.insert(IntPair(i, j));
      }
      else
      {
        negatives.insert(IntPair(i, j));
      }
    }
  }
}

void getClusteringStats(const IntVector& inferredClustering,
                        const IntVector& trueClustering,
                        double& recall,
                        double& precision,
                        double& rand)
{
  IntPairSet inferredPositives, inferredNegatives;
  IntPairSet solutionPositives, solutionNegatives;
  
  getClusteringPairs(inferredClustering, inferredPositives, inferredNegatives);
  getClusteringPairs(trueClustering, solutionPositives, solutionNegatives);
  
  IntPairSet TP, FP, TN, FN;

  // true positives
  std::set_intersection(inferredPositives.begin(), inferredPositives.end(),
                        solutionPositives.begin(), solutionPositives.end(),
                        std::inserter(TP, TP.begin()));
  
  // false positives
  std::set_intersection(inferredPositives.begin(), inferredPositives.end(),
                        solutionNegatives.begin(), solutionNegatives.end(),
                        std::inserter(FP, FP.begin()));
  
  // true negatives
  std::set_intersection(inferredNegatives.begin(), inferredNegatives.end(),
                        solutionNegatives.begin(), solutionNegatives.end(),
                        std::inserter(TN, TN.begin()));
  
  // false negatives
  std::set_intersection(inferredNegatives.begin(), inferredNegatives.end(),
                        solutionPositives.begin(), solutionPositives.end(),
                        std::inserter(FN, FN.begin()));
  
  recall = solutionPositives.size() == 0 ? 1. : (double)TP.size() / solutionPositives.size();
  precision = inferredPositives.size() == 0 ? 1. : (double)TP.size() / inferredPositives.size();
  rand = (double)(TP.size() + TN.size()) / (TP.size() + TN.size() + FP.size() + FN.size());
}

double getSelectionAccuracy(const IntVector& inferredSelection,
                            const IntVector& trueSelection)
{
  const int n = inferredSelection.size();
  
  double accuracy = 0;
  
  for (int i = 0; i < n; ++i)
  {
    if (inferredSelection[i] == trueSelection[i])
    {
      accuracy += 1;
    }
  }
  
  return accuracy / n;
}

double getSelectionAccuracy(const Solution& solInferred,
                            const Solution& solTrue)
{
  const int n = solInferred.getNrSelectedTrees();
  
  double accuracy = 0;
  
  for (int i = 0; i < n; ++i)
  {
    
    // REMOVING DUMMY NODES
    StringPairSet edgesInferred = Solution::removeDummyNodes(solInferred.getSelectedTree(i).getParentalPairs());
    
    if (edgesInferred == solTrue.getSelectedTree(i).getParentalPairs())
    {
      accuracy += 1;      
    }
  }
  
  return accuracy / n;
}

int main(int argc, char** argv)
{
  bool header = false;
  std::string inputFilename;
  std::string groundTruthFilename;
  
  lemon::ArgParser ap (argc, argv);
  ap.other("sol", "Solutions");
  ap.refOption("i", "Input file name", inputFilename, true);
  ap.refOption("s", "Ground truth file name", groundTruthFilename);
  ap.refOption("H", "Output header", header);
  ap.parse();
  
  std::vector<Solution> inferredSolutions;
  for (const std::string& filename : ap.files())
  {
    Solution sol;
    try
    {
      std::ifstream inSol(filename.c_str());
      if (!inSol.good())
      {
        std::cerr << "Error: could not open '" << filename << "' for reading" << std::endl;
        return 1;
      }
      inSol >> sol;
      inferredSolutions.push_back(sol);
    }
    catch (std::runtime_error& e)
    {
      std::cerr << "Error: " << e.what() << std::endl;
      return 1;
    }
  }
  
  int minTrees = -1, maxTrees = -1, medianTrees = -1, totalTrees = -1;
  int nrMutations = -1;
  if (!inputFilename.empty())
  {
    InputInstance input;
    std::ifstream inInput(inputFilename.c_str());
    if (!inInput.good())
    {
      std::cerr << "Error: could not open '" << inputFilename << "' for reading" << std::endl;
      return 1;
    }
    inInput >> input;
    minTrees = input.getMinNrTrees();
    maxTrees = input.getMaxNrTrees();
    medianTrees = input.getMedianNrTrees();
    totalTrees = input.getTotalNrTrees();
    nrMutations = input.getMutationSetSize(";");
    
    if (inferredSolutions.empty())
    {
      for (int i = 0; i < input.getNrPatients(); ++i)
      {
        for (int j = 0; j < input.getPatientTrees(i).size(); ++j)
        {
          CloneTree exp(input.getPatientTrees(i)[j], ";");
          std::cout << "Patient " << i << ", tree " << j << ", AD pairs: " << exp.getAncestralPairs().size() << std::endl;
        }
      }
    }
  }
  else
  {
    std::cerr << "Error: input file name missing" << std::endl;
    return 1;
  }
  
  if (!groundTruthFilename.empty())
  {
    Solution groundTruthSolution;
    try
    {
      std::ifstream inSol(groundTruthFilename.c_str());
      if (!inSol.good())
      {
        std::cerr << "Error: could not open '" << groundTruthFilename << "' for reading" << std::endl;
        return 1;
      }
      inSol >> groundTruthSolution;
    }
    catch (std::runtime_error& e)
    {
      std::cerr << "Error: " << e.what() << std::endl;
      return 1;
    }
    
    if (header)
    {
      std::cout << "filename" << "\t";
      
      if (!inputFilename.empty())
      {
        std::cout << "min_trees" << "\t" << "median_trees" << "\t" << "max_trees" << "\t" << "total_trees" << "\t";
      }
      
      std::cout << "PC_dist" << "\t"
                << "patients" << "\t"
                << "inferred_k" << "\t" << "true_k" << "\t"
                << "consensus_PC" << "\t"
                << "selected_PC" << "\t"
                << "clustering_recall" << "\t" << "clustering_precision" << "\t"
                << "clustering_Rand" << "\t" << "selection_accuracy" << "\t"
                << "sol_PC_dist" << std::endl;
    }
    
    for (int solIdx = 0; solIdx < inferredSolutions.size(); ++solIdx)
    {
      std::cout << ap.files()[solIdx] << "\t";
      if (!inputFilename.empty())
      {
        std::cout << minTrees << "\t" << medianTrees << "\t" << maxTrees << "\t" << totalTrees << "\t";
      }
      
      // These functions use original trees output by RECAP
      std::cout << inferredSolutions[solIdx].getParentChildDistance()
                << "\t" << inferredSolutions[solIdx].getNrSelectedTrees()
                << "\t" << inferredSolutions[solIdx].getNrClusters()
                << "\t" << groundTruthSolution.getNrClusters();
      
      
      // REMOVE DUMMY NODES FOR REMAINING FUNCTIONS
      auto consensusDist = computeConsensusDistance(inferredSolutions[solIdx], groundTruthSolution);
      std::cout << "\t" << consensusDist.first;
      
      // 2. Compute distances of selected trees
      
      auto selectionDist = computeSelectedDistance(inferredSolutions[solIdx], groundTruthSolution);
      std::cout << "\t" << selectionDist.first;
      
      double recall;
      double precision;
      double rand;
      
      getClusteringStats(inferredSolutions[solIdx].getClustering(), groundTruthSolution.getClustering(),
                         recall, precision, rand);
      
      std::cout << "\t" << recall << "\t" << precision << "\t" << rand;
      
      std::cout << "\t" << getSelectionAccuracy(inferredSolutions[solIdx], groundTruthSolution);
      
      std::cout << "\t" << groundTruthSolution.getParentChildDistance(true);
      
      std::cout << std::endl;
    }
  }
  else
  {
    if (header)
    {
      std::cout << "filename" << "\t";
    
      if (!inputFilename.empty())
      {
        std::cout << "min_trees" << "\t" << "median_trees" << "\t" << "max_trees" << "\t" << "total_trees" << "\t";
      }
    
      std::cout << "PC_dist" << "\t"
                << "patients" << "\t"
                << "inferred_k" << std::endl;
    }
    
    for (int solIdx = 0; solIdx < inferredSolutions.size(); ++solIdx)
    {
      std::cout << ap.files()[solIdx] << "\t";
      if (!inputFilename.empty())
      {
        std::cout << minTrees << "\t" << medianTrees << "\t" << maxTrees << "\t" << totalTrees << "\t";
      }
      
      std::cout << inferredSolutions[solIdx].getParentChildDistance()
                << "\t" << inferredSolutions[solIdx].getNrSelectedTrees()
                << "\t" << inferredSolutions[solIdx].getNrClusters()
                << std::endl;
    }
  }
  
  return 0;
}
