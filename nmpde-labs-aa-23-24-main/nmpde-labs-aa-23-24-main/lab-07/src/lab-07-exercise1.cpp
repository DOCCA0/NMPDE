#include "NonLinearDiffusion.hpp"
#include <deal.II/base/convergence_table.h>


// Main function.
int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);
  ConvergenceTable table;


  const std::string  mesh_file_name = "../mesh/mesh-square-h0.012500.msh";
  const unsigned int degree         = 1;

  NonLinearDiffusion problem(mesh_file_name, degree);

  problem.setup();
  problem.solve_newton();
  problem.output();


  return 0;
}