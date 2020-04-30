
#include "utils.h"
#include <lemon/arg_parser.h>
#include "clonetree.h"
#include "parentchildgraph.h"
#include "mcctsolverca.h"

int main(int argc, char** argv)
{ 
  int k = 0;
  int restarts = 50;
  int seed = 0;
  int timeLimit = -1;
  int minClusterSize = 1;
  int minCount = 1;
  std::string resultspath;
  std::string sep = ";";
  std::string rootMutation;
  
  lemon::ArgParser ap (argc, argv);
  ap.other("trees", "Input trees");
  ap.refOption("k", "Number of clusters", k, true);
  ap.refOption("r", "Number of restarts (default: 50)", restarts);
  ap.refOption("t", "Time limit in seconds (default: -1, no time limit)", timeLimit);
  ap.refOption("s", "Seed for random number generator (default: 0)", seed);
  ap.refOption("p", "Path to results (make it unique)", resultspath, true);
  ap.refOption("S", "Mutation cluster separator (default: ;)", sep);
  ap.refOption("R", "Root mutation", rootMutation, true);
  ap.refOption("c", "Minimum mutation count in cohort (default: 1)", minCount);
  ap.parse();
  
  g_rng = std::mt19937(seed);
  
  if (ap.files().size() != 1)
  {
    std::cerr << "Error: <trees> must be specified" << std::endl;
    return 1;
  }

  
  std::string filenameTrees = ap.files()[0];
    
  InputInstance scriptT;
 
  try{
    if (filenameTrees != "-"){
      std::ifstream inScriptT(filenameTrees.c_str());
      if (!inScriptT.good()){
        std::cerr << "Error: could not open '" << filenameTrees << "' for reading" << std::endl;
        return 1;
      }
      inScriptT >> scriptT; 
    }
    else {
      std::cin >> scriptT;
    }
  }
  catch (std::runtime_error& e){
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }
  
  if (timeLimit > 0)
  {
    restarts = -1;
  }
  
  if (minCount > 1){
    scriptT.removeRareMutations(sep, minCount);
  }
  scriptT.completePatientTrees(sep);
    
  McctSolverCA solver(scriptT, rootMutation, sep,
                      k, timeLimit,
                      minClusterSize, restarts);
  McctSolver::run(solver, resultspath);
}
