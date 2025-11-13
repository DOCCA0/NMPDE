#include "HeatNonLinear.hpp"
#include <deal.II/base/convergence_table.h>




#ifdef CONVERGENCE
// Main function.
int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);
  const unsigned int               mpi_rank =
    Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  const unsigned int degree = 1;

  double T     = 0.5;
  double theta = 1.0;

  const std::vector<double> deltat_vector = {
     0.125, 0.0625, 0.03125, 0.015625, 0.0078125};
  std::vector<double> errors_L2;
  std::vector<double> errors_H1;

  for (const auto &deltat : deltat_vector)
    {
      HeatNonLinear problem("../mesh/mesh-cube-20.msh", degree, T, deltat);

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

  #endif //CONVERGENCE

  
  #ifdef SPATIAL_CONVERGENCE
    
  ConvergenceTable table;

  const std::vector<std::string> meshes = {"../mesh/mesh-cube-5.msh",
                                           "../mesh/mesh-cube-10.msh",
                                           "../mesh/mesh-cube-20.msh",
                                           "../mesh/mesh-cube-40.msh"};
  const std::vector<double>      h_vals = {0.1,
                                           0.05,
                                           0.025,
                                           0.0125};

  // Only for Exercise 1:
  std::ofstream convergence_file("convergence.csv");
  convergence_file << "h,eL2,eH1" << std::endl;

  T =      5e-4;
  double deltat = 1e-4;

  for (unsigned int i = 0; i < meshes.size(); ++i)
    {
  
  HeatNonLinear problem(meshes[i], degree, T, deltat);

      problem.setup();
      //problem.assemble();
      problem.solve();
      //problem.output();

      // Only for Exercise 1:
      const double error_L2 = problem.compute_error(VectorTools::L2_norm);
      const double error_H1 = problem.compute_error(VectorTools::H1_norm);

      table.add_value("h", h_vals[i]);
      table.add_value("L2", error_L2);
      table.add_value("H1", error_H1);

      convergence_file << h_vals[i] << "," << error_L2 << "," << error_H1
                       << std::endl;
    }

  // // Only for Exercise 1:
  table.evaluate_all_convergence_rates(ConvergenceTable::reduction_rate_log2);
  table.set_scientific("L2", true);
  table.set_scientific("H1", true);
  table.write_text(std::cout);

  #endif //SPATIAL_CONVERGENCE


  return 0;
}

#ifndef CONVERGENCE

// Main function.
int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

  const unsigned int degree = 1;

  const double T      = 1.0;
  const double deltat = 0.05;

  HeatNonLinear problem("../mesh/mesh-cube-20.msh", degree, T, deltat);

  problem.setup();
  problem.solve();

  return 0;
}
#endif //CONVERGENCE