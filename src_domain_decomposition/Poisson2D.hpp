#ifndef POISSON_1D_HPP
#define POISSON_1D_HPP

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/solver_gmres.h>


#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>

using namespace dealii;

//#define CG
#define TRANSPORT_COEFFICIENT
#define REACTION_COEFFICIENT
#define NEUMANN
//#define ROBIN
#define CONSERVATIVE_TRANSPORT_COEFFICIENT


/**
 * Class managing the differential problem.
 */
class Poisson2D
{
public:
  // Physical dimension (1D, 2D, 3D)
  static constexpr unsigned int dim = 2;

  // transportin coefficient.
  #ifdef TRANSPORT_COEFFICIENT
  class TransportCoefficient : public Function<dim>
  {
  public:
    virtual void
    vector_value(const Point<dim> & /*p*/,
                 Vector<double> &values) const override
    {
      for (unsigned int i = 0; i < dim - 1; ++i)
      values[0] = 1.0;

      values[1] = 1.0;
    }

    virtual double
    value(const Point<dim> & /*p*/,
          const unsigned int component = 0) const override
    {
      if (component == 0)
        return 1.0;
      else
        return 1.0;
    }
  };
  #endif //TRANSPORT_COEFFICIENT

  #ifdef CONSERVATIVE_TRANSPORT_COEFFICIENT

  class ConsTransportCoefficient : public Function<dim>
  {
  public:
    virtual void
    vector_value(const Point<dim> & /*p*/,
                 Vector<double> &values) const override
    {
      for (unsigned int i = 0; i < dim - 1; ++i)
      values[0] = 1.0;

      values[1] = 1.0;
    }

    virtual double
    value(const Point<dim> & /*p*/,
          const unsigned int component = 0) const override
    {
      if (component == 0)
        return 1.0;
      else
        return 1.0;
    }
  };

  #endif //CONSERVATIVE_TRANSPORT_COEFFICIENT

  class FunctionG : public Function<dim>
  {
  public:
    // Constructor.
    FunctionG()
    {}

    // Evaluation.
    virtual double
    value(const Point<dim> &p,
          const unsigned int /*component*/ = 0) const override
    {
      return p[0];
    }
  };


  // Forcing term.
  class ForcingTerm : public Function<dim>
  {
  public:
    // Constructor.
    ForcingTerm()
    {}

    // Evaluation.
    virtual double
    value(const Point<dim> &p,
          const unsigned int /*component*/ = 0) const override
    {
      return p[0]+2.0;

    }
  };

  #ifdef ROBIN

  class FunctionGamma : public Function<dim>
  { public:
    // Constructor.
    FunctionGamma()
    {}

    // Evaluation.
    virtual double
    value(const Point<dim> &p,
          const unsigned int /*component*/ = 0) const override
    {
      if(p[0] ==0.0)
        return 1.0;
      if(p[0] ==1.0)
        return 1.0;
      else
        return 0.0;
    }
  };

  #endif //ROBIN

  #ifdef NEUMANN
   // Neumann boundary conditions.
  class FunctionH : public Function<dim>
  {
  public:
    // Constructor.
    FunctionH()
    {}

    // Evaluation:
    virtual double
    value(const Point<dim> &p, const unsigned int /*component*/ = 0) const
    {
      if(p[0] ==0.0)
        return -1.0;
      //if(p[0] ==1.0)
        //return 1.0;
      else
        return 0.0;
      
    }
  };
#endif //NEUMANN

#ifdef REACTION_COEFFICIENT
  // Reaction coefficient.
  class ReactionCoefficient : public Function<dim>
  {
  public:
    // Constructor.
    ReactionCoefficient()
    {}

    // Evaluation.
    virtual double
    value(const Point<dim> & /*p*/,
          const unsigned int /*component*/ = 0) const override
    {
      return 1.0;
    }
  };
#endif //REACTION_COEFFICIENT

  // Constructor.
  Poisson2D(const unsigned int &subdomain_id_)
    : subdomain_id(subdomain_id_)
  {}

  // Initialization.
  void
  setup();

  // System assembly.
  void
  assemble();

  // System solution.
  void
  solve();

  // Output.
  void
  output(const unsigned int &iter) const;

  // Apply Dirichlet conditions on the interface with another Poisson2D problem.
  void
  apply_interface_dirichlet(const Poisson2D &other);

  // Apply Neumann conditions on the interface with another Poisson2D problem.
  void
  apply_interface_neumann(Poisson2D &other);

  // Get the solution vector.
  const Vector<double> &
  get_solution() const
  {
    return solution;
  }

  // Apply relaxation.
  void
  apply_relaxation(const Vector<double> &old_solution, const double &lambda);

protected:
  // Build an interface map, that is construct a map that to each DoF on the
  // interface for this subproblem associates the corresponding interface DoF on
  // the other subdomain.
  std::map<types::global_dof_index, types::global_dof_index>
  compute_interface_map(const Poisson2D &other) const;

  // ID of current subdomain (0 or 1).
  const unsigned int subdomain_id;

  FunctionG function_g;

  #ifdef REACTION_COEFFICIENT

  ReactionCoefficient reaction_coefficient;

  #endif //REACTION_COEFFICIENT

#ifdef TRANSPORT_COEFFICIENT

  TransportCoefficient transport_coefficient;

#endif //TRANSPORT_COEFFICIENT

#ifdef CONSERVATIVE_TRANSPORT_COEFFICIENT

ConsTransportCoefficient cons_transport_coefficient;

#endif //CONSERVATIVE_TRANSPORT_COEFFICIENT


  ForcingTerm forcing_term;

  // Triangulation.
  Triangulation<dim> mesh;

  // Support points.
  std::map<types::global_dof_index, Point<dim>> support_points;

  // Finite element space.
  std::unique_ptr<FiniteElement<dim>> fe;

  // Quadrature formula.
  std::unique_ptr<Quadrature<dim>> quadrature;

#ifdef NEUMANN
  // Quadrature formula used on boundary lines.
  std::unique_ptr<Quadrature<dim - 1>> quadrature_boundary;

  // h(x).
  FunctionH function_h;
#endif //NEUMANN

#ifdef ROBIN

FunctionGamma function_gamma;

#endif //ROBIN


  // DoF handler.
  DoFHandler<dim> dof_handler;

  // Sparsity pattern.
  SparsityPattern sparsity_pattern;

  // System matrix.
  SparseMatrix<double> system_matrix;

  // System right-hand side.
  Vector<double> system_rhs;

  // System solution.
  Vector<double> solution;
};

#endif