#include <fstream>

#include "Parabolic.hpp"
#include <deal.II/base/convergence_table.h>
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

  double T     = 0.5;
  double theta = 1.0;

  const std::vector<double> deltat_vector = {
    0.25, 0.125, 0.0625, 0.03125, 0.015625};
  std::vector<double> errors_L2;
  std::vector<double> errors_H1;

  for (const auto &deltat : deltat_vector)
    {
      Parabolic problem("../mesh/m1.msh", degree, T, deltat, theta);

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

  // Default run: Exercise 1.5 settings
  std::string mesh_file_name = "AUTO_1D_10"; // auto-generate 10-element 1D mesh
  unsigned int degree = 1;                    // piecewise linear
  double T = 1.0;
  double deltat = 0.01;                       // time step
  double theta = 1.0;                         // backward Euler
  double kappa = 1.0;                        // transport coefficient
  bool use_supg = false;                      // SUPG toggle


  if (argc > 1)
    mesh_file_name = argv[1];
  if (argc > 2)
    degree = std::stoi(argv[2]);
  if (argc > 3)
    T = std::stod(argv[3]);
  if (argc > 4)
    deltat = std::stod(argv[4]);
  if (argc > 5)
    theta = std::stod(argv[5]);
  if (argc > 6)
    kappa = std::stod(argv[6]);
  if (argc > 7)
    use_supg = (std::stoi(argv[7]) != 0);

  Parabolic problem(mesh_file_name, degree, T, deltat, theta, kappa, use_supg);

  problem.setup();
  problem.solve();

  // const double error_L2 = problem.compute_error(VectorTools::L2_norm);
  // const double error_H1 = problem.compute_error(VectorTools::H1_norm);
  // printf("L2 error: %e\n", error_L2);
  // printf("H1 error: %e\n", error_H1);



  return 0;
}
#endif //CONVERGENCE