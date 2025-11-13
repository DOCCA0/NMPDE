#ifndef NON_LINEAR_DIFFUSION_HPP
#define NON_LINEAR_DIFFUSION_HPP

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/fully_distributed_tria.h>

#include <deal.II/fe/mapping_fe.h>


#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_in.h>

#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>

//#define NEUMANN
//#define CONVERGENCE

using namespace dealii;

// Class representing the non-linear diffusion problem.
class NonLinearDiffusion
{
public:
  // Physical dimension (1D, 2D, 3D)
  static constexpr unsigned int dim = 2;

  // Function for the mu_0 coefficient.
  class FunctionMu0 : public Function<dim>
  {
  public:
    virtual double
    value(const Point<dim> & /*p*/,
          const unsigned int /*component*/ = 0) const override
    {
      return 1.0;
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

// Function for the forcing term using the Expression class.
class ForcingTerm : public Function<dim>
{
public:
  virtual double
  value(const Point<dim> & p,
        const unsigned int /*component*/ = 0) const override
  {
    double x = p[0];
    double pi = M_PI;

    return -2 * std::pow(pi, 2) * std::cos(pi * x) * std::cos(pi * x) * std::sin(pi * x)
            +std::pow(pi, 2) * std::sin(pi * x) * (1 + std::sin(pi * x) * std::sin(pi * x));
  }
};


  #ifdef NEUMANN
  // Function for Neumann boundary conditions.
class FunctionH : public Function<dim>
{
public:
    FunctionH() : Function<dim>() {}
    virtual double value(const Point<dim> & p, const unsigned int component = 0) const override
    {
        if(p[0]==1)
          return M_PI;
        if(p[0]==0)
          return -M_PI;
        else
          return 0.0;
    }
};

  #endif // NEUMANN

    // Function for Dirichlet boundary conditions.
  class FunctionG : public Function<dim>
  {
  public:
    virtual double
    value(const Point<dim> & p,
          const unsigned int component = 0) const override
    {
      return std::sin(M_PI*p[0]);
      // return 0.0;
    }
  };

     // Exact solution.
  class ExactSolution : public Function<dim>
  {
  public:
    virtual double
    value(const Point<dim> &p,
          const unsigned int component = 0) const override
    {
      return std::sin(M_PI*p[0]);
    }

    virtual Tensor<1, dim>
    gradient(const Point<dim> &p,
             const unsigned int /*component*/ = 0) const override
    {
      Tensor<1, dim> result;

      // duex / dx
      result[0] = M_PI * std::cos(M_PI * p[0]);

      // duex / dy
      result[1] = 0;

      // duex / dz
      // result[2] = 4 * M_PI * std::sin(5 * M_PI * get_time()) *
      //             std::sin(2 * M_PI * p[0]) * std::sin(3 * M_PI * p[1]) *
      //             std::cos(4 * M_PI * p[2]);

      return result;
    }
  };
  

  // Constructor.
  NonLinearDiffusion(const std::string &mesh_file_name_, const unsigned int &r_)
    : mpi_size(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD))
    , mpi_rank(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD))
    , pcout(std::cout, mpi_rank == 0)
    , mesh_file_name(mesh_file_name_)
    , r(r_)
    , mesh(MPI_COMM_WORLD)
  {}

  // Initialization.
  void
  setup();

  // Solve the problem using Newton's method.
  void
  solve_newton();

  // Output.
  void
  output() const;

  #ifdef CONVERGENCE
  // Compute the error.
  double
  compute_error(const VectorTools::NormType &norm_type);
 
 #endif //CONVERGENCE

protected:
  // Assemble the tangent problem.
  void
  assemble_system();

  // Solve the tangent problem.
  void
  solve_system();

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

  FunctionG function_g;

#ifdef NEUMANN
  // Quadrature formula used on boundary lines.
  std::unique_ptr<Quadrature<dim - 1>> quadrature_boundary;

  // h(x).
  FunctionH function_h;
#endif //NEUMANN

#ifdef CONVERGENCE
  // Exact solution.
  ExactSolution exact_solution;
#endif //CONVERGENCE

  // Discretization. ///////////////////////////////////////////////////////////

  // Mesh file name.
  const std::string &mesh_file_name;

  // Polynomial degree.
  const unsigned int r;

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

  // Solution increment (without ghost elements).
  TrilinosWrappers::MPI::Vector delta_owned;

  // System solution (without ghost elements).
  TrilinosWrappers::MPI::Vector solution_owned;

  // System solution (including ghost elements).
  TrilinosWrappers::MPI::Vector solution;
};

#endif