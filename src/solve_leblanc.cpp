#include "selfish_routing.h"
#include "fr_solver.h"

void usage()
{
  std::cout << std::endl;
  std::cout << "USAGE: leblanc_solver nodeFileName edgeFileName SR/FR" << std::endl;
}

int main(int argc, char **argv)
{
  srand(time(NULL));

  if (argc == 4)
  {
    std::string nodeFileName = argv[1];
    std::string edgeFileName = argv[2];
    NetworkDesignProblem problem(nodeFileName, edgeFileName);
    SelfishRoutingSolver solver(problem, true);
 
    std::string solveMethod = argv[3];
    if (solveMethod == "SR")
    {
      std::vector<std::vector<double>> emptyInitialSrFlows;
      std::cout << "solve user optimal" << std::endl;
      solver.solve(problem.initialFrEdges, FunctionType::USER_OPTIMAL, emptyInitialSrFlows);
      std::cout << "solve system optimal" << std::endl;
      solver.solve(problem.initialFrEdges, FunctionType::SYSTEM_OPTIMAL, emptyInitialSrFlows);
      std::cout << "done solving" << std::endl;
    }
    else
    {
      FRSolver frSolver(problem);
      frSolver.solveBranchAndBound();
      std::cout << "done solving" << std::endl;
    }
  }
  else
  {
    usage();
    return -1;
  }

  return 0;
};
