#ifndef HEAT_NON_LINEAR_HPP
#define HEAT_NON_LINEAR_HPP

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/fully_distributed_tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_values_extractors.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/grid_in.h>

#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>

using namespace dealii;

#define NEUMANN

// Class representing the non-linear diffusion problem.
class HeatNonLinear
{
public:
  // Physical dimension (1D, 2D, 3D)
  static constexpr unsigned int dim = 3;

  // Function for the mu_0 coefficient.
  class FunctionMu0 : public Function<dim>
  {
  public:
    virtual double
    value(const Point<dim> & /*p*/,
          const unsigned int /*component*/ = 0) const override
    {
      return 0.1;
    }
  };

  // Function for the mu_1 coefficient.
  class FunctionMu1 : public Function<dim>
  {
  public:
    virtual double
    value(const Point<dim> & /*p*/,
          const unsigned int /*component*/ = 0) const override
    {
      return 1.0;
    }
  };

  // Function for the forcing term.
  class ForcingTerm : public Function<dim>
  {
  public:
    virtual double
    value(const Point<dim> & /*p*/,
          const unsigned int /*component*/ = 0) const override
    {
      if (get_time() < 0.25)
        return 2.0;
      else
        return 0.0;
    }
  };

  // Function for Dirichlet boundary conditions.
  class FunctionG : public Function<dim>
  {
  public:
    virtual double
    value(const Point<dim> & /*p*/,
          const unsigned int /*component*/ = 0) const override
    {
      return 0.0;
    }
  };

       class FunctionH : public Function<dim>
{
public:
    FunctionH() : Function<dim>() {}
    virtual double value(const Point<dim> & p, const unsigned int /*component*/ = 0) const override
    {
        double t = 0.0; // replace with your time variable
        double uex = std::sin(5 * M_PI * t) * std::sin(2 * M_PI * p[0]) * std::sin(3 * M_PI * p[1]) * std::sin(4 * M_PI * p[2]);
        double scalar = 1.0 + 10.0 * uex * uex;
        if (p[0] == 0.0)
            return -1.0 * scalar * (2 * M_PI * std::cos(2 * M_PI * p[0]) * std::sin(5 * M_PI * t) * std::sin(3 * M_PI * p[1]) * std::sin(4 * M_PI * p[2]));
        else if (p[0] == 1.0)
            return 1.0 * scalar * (2 * M_PI * std::cos(2 * M_PI * p[0]) * std::sin(5 * M_PI * t) * std::sin(3 * M_PI * p[1]) * std::sin(4 * M_PI * p[2]));
        else if (p[1] == 0.0)
            return -1.0 * scalar * (3 * M_PI * std::cos(3 * M_PI * p[1]) * std::sin(5 * M_PI * t) * std::sin(2 * M_PI * p[0]) * std::sin(4 * M_PI * p[2]));
        else if (p[1] == 1.0)
            return 1.0 * scalar * (3 * M_PI * std::cos(3 * M_PI * p[1]) * std::sin(5 * M_PI * t) * std::sin(2 * M_PI * p[0]) * std::sin(4 * M_PI * p[2]));
        else if (p[2] == 0.0)
            return -1.0 * scalar * (4 * M_PI * std::cos(4 * M_PI * p[2]) * std::sin(5 * M_PI * t) * std::sin(2 * M_PI * p[0]) * std::sin(3 * M_PI * p[1]));
        else if (p[2] == 1.0)
            return 1.0 * scalar * (4 * M_PI * std::cos(4 * M_PI * p[2]) * std::sin(5 * M_PI * t) * std::sin(2 * M_PI * p[0]) * std::sin(3 * M_PI * p[1]));
        else
            return 0.0;
    }
};
  // Function for initial conditions.
  class FunctionU0 : public Function<dim>
  {
  public:
    virtual double
    value(const Point<dim> & /*p*/,
          const unsigned int /*component*/ = 0) const override
    {
      return 0.0;
    }
  };

  // Constructor. We provide the final time, time step Delta t and theta method
  // parameter as constructor arguments.
  HeatNonLinear(const std::string  &mesh_file_name_,
                const unsigned int &r_,
                const double       &T_,
                const double       &deltat_)
    : mpi_size(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD))
    , mpi_rank(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD))
    , pcout(std::cout, mpi_rank == 0)
    , T(T_)
    , mesh_file_name(mesh_file_name_)
    , r(r_)
    , deltat(deltat_)
    , mesh(MPI_COMM_WORLD)
  {}

  // Initialization.
  void
  setup();

  // Solve the problem.
  void
  solve();

protected:
  // Assemble the tangent problem.
  void
  assemble_system();

  // Solve the linear system associated to the tangent problem.
  void
  solve_linear_system();

  // Solve the problem for one time step using Newton's method.
  void
  solve_newton();

  // Output.
  void
  output(const unsigned int &time_step) const;

  // MPI parallel. /////////////////////////////////////////////////////////////

  // Number of MPI processes.
  const unsigned int mpi_size;

  // This MPI process.
  const unsigned int mpi_rank;

  // Parallel output stream.
  ConditionalOStream pcout;

  // Problem definition. ///////////////////////////////////////////////////////

  // mu_0 coefficient.
  FunctionMu0 mu_0;

  // mu_1 coefficient.
  FunctionMu1 mu_1;

  // Forcing term.
  ForcingTerm forcing_term;

  // Dirichlet boundary conditions.
  FunctionG function_g;



  #ifdef NEUMANN
  // Quadrature formula used on boundary lines.
  std::unique_ptr<Quadrature<dim - 1>> quadrature_boundary;

  // h(x).
  FunctionH function_h;
#endif //NEUMANN

  // Initial conditions.
  FunctionU0 u_0;

  // Current time.
  double time;

  // Final time.
  const double T;

  // Discretization. ///////////////////////////////////////////////////////////

  // Mesh file name.
  const std::string mesh_file_name;

  // Polynomial degree.
  const unsigned int r;

  // Time step.
  const double deltat;

  // Mesh.
  parallel::fullydistributed::Triangulation<dim> mesh;

  // Finite element space.
  std::unique_ptr<FiniteElement<dim>> fe;

  // Quadrature formula.
  std::unique_ptr<Quadrature<dim>> quadrature;

  // DoF handler.
  DoFHandler<dim> dof_handler;

  // DoFs owned by current process.
  IndexSet locally_owned_dofs;

  // DoFs relevant to the current process (including ghost DoFs).
  IndexSet locally_relevant_dofs;

  // Jacobian matrix.
  TrilinosWrappers::SparseMatrix jacobian_matrix;

  // Residual vector.
  TrilinosWrappers::MPI::Vector residual_vector;

  // Increment of the solution between Newton iterations.
  TrilinosWrappers::MPI::Vector delta_owned;

  // System solution (without ghost elements).
  TrilinosWrappers::MPI::Vector solution_owned;

  // System solution (including ghost elements).
  TrilinosWrappers::MPI::Vector solution;

  // System solution at previous time step.
  TrilinosWrappers::MPI::Vector solution_old;
};

#endif