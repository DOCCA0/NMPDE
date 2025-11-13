#include <fstream>

#include "Parabolic.hpp"
#include <deal.II/base/convergence_table.h>


#ifdef CONVERGENCE
// Main function.
int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);
  const unsigned int               mpi_rank =
    Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  const unsigned int degree = 2;

  double T     = 2.0;
  double theta = 0.5;

  const std::vector<double> deltat_vector = {
    0.25, 0.125, 0.0625, 0.03125, 0.015625};
  std::vector<double> errors_L2;
  std::vector<double> errors_H1;

  for (const auto &deltat : deltat_vector)
    {
      Parabolic problem("../mesh/mesh-square-h0.012500.msh", degree, T, deltat, theta);

      problem.setup();
      problem.solve();

      errors_L2.push_back(problem.compute_error(VectorTools::L2_norm));
      errors_H1.push_back(problem.compute_error(VectorTools::H1_norm));
    }

  // Print the errors and estimate the convergence order.
  if (mpi_rank == 0)
    {
      std::cout << "==============================================="
                << std::endl;

      std::ofstream convergence_file("convergence.csv");
      convergence_file << "dt,eL2,eH1" << std::endl;

      for (unsigned int i = 0; i < deltat_vector.size(); ++i)
        {
          convergence_file << deltat_vector[i] << "," << errors_L2[i] << ","
                           << errors_H1[i] << std::endl;

          std::cout << std::scientific << "dt = " << std::setw(4)
                    << std::setprecision(2) << deltat_vector[i];

          std::cout << std::scientific << " | eL2 = " << errors_L2[i];

          // Estimate the convergence order.
          if (i > 0)
            {
              const double p =
                std::log(errors_L2[i] / errors_L2[i - 1]) /
                std::log(deltat_vector[i] / deltat_vector[i - 1]);

              std::cout << " (" << std::fixed << std::setprecision(2)
                        << std::setw(4) << p << ")";
            }
          else
            std::cout << " (  - )";

          std::cout << std::scientific << " | eH1 = " << errors_H1[i];

          // Estimate the convergence order.
          if (i > 0)
            {
              const double p =
                std::log(errors_H1[i] / errors_H1[i - 1]) /
                std::log(deltat_vector[i] / deltat_vector[i - 1]);

              std::cout << " (" << std::fixed << std::setprecision(2)
                        << std::setw(4) << p << ")";
            }
          else
            std::cout << " (  - )";

          std::cout << "\n";
        }
    }


  
  #ifdef SPATIAL_CONVERGENCE
    
  ConvergenceTable table;

  const std::vector<std::string> meshes = {"../mesh/mesh-square-h0.100000.msh",
                                           "../mesh/mesh-square-h0.050000.msh",
                                           "../mesh/mesh-square-h0.025000.msh",
                                           "../mesh/mesh-square-h0.012500.msh"};
  const std::vector<double>      h_vals = {0.1,
                                           0.05,
                                           0.025,
                                           0.0125};

  // Only for Exercise 1:
  std::ofstream convergence_file("convergence.csv");
  convergence_file << "h,eL2,eH1" << std::endl;

  T =      5e-5;
  double deltat = 1e-6;

  for (unsigned int i = 0; i < meshes.size(); ++i)
    {
  
    Parabolic problem(meshes[i], degree, T, deltat,theta);

      problem.setup();
      problem.solve();

      const double error_L2 = problem.compute_error(VectorTools::L2_norm);
      const double error_H1 = problem.compute_error(VectorTools::H1_norm);

      table.add_value("h", h_vals[i]);
      table.add_value("L2", error_L2);
      table.add_value("H1", error_H1);

      convergence_file << h_vals[i] << "," << error_L2 << "," << error_H1
                       << std::endl;
    }

  table.evaluate_all_convergence_rates(ConvergenceTable::reduction_rate_log2);
  table.set_scientific("L2", true);
  table.set_scientific("H1", true);
  table.write_text(std::cout);

  #endif //SPATIAL_CONVERGENCE

  return 0;
}

#endif

#ifndef CONVERGENCE
// Main function.
int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

  const std::string  mesh_file_name = "../mesh/mesh-square-h0.012500.msh";
  
  //const std::string  mesh_file_name = "../mesh/square_mesh.msh";
  const unsigned int degree         = 2;
  const double T      = 2.0;
  const double deltat = 0.25;
  const double theta  = 0.5;

  Parabolic problem(mesh_file_name, degree, T, deltat, theta);

  problem.setup();
  problem.solve();


  return 0;
}
#endif //CONVERGENCE