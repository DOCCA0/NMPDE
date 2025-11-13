#include "Poisson2D.hpp"

void
Poisson2D::setup()
{
  // Create the mesh.
  {
    GridIn<dim> grid_in;
    grid_in.attach_triangulation(mesh);

    std::ifstream grid_in_file("../mesh/mg" +
                               std::to_string(subdomain_id) + ".msh");

    grid_in.read_msh(grid_in_file);
  }

  // Initialize the finite element space.
  {
    fe         = std::make_unique<FE_SimplexP<dim>>(1);
    quadrature = std::make_unique<QGaussSimplex<dim>>(2);

    
#ifdef NEUMANN
    quadrature_boundary = std::make_unique<QGaussSimplex<dim - 1>>(2);

    std::cout << "  Quadrature points per boundary cell = "
              << quadrature_boundary->size() << std::endl;
#endif //NEUMANN
  }

  // Initialize the DoF handler.
  {
    dof_handler.reinit(mesh);
    dof_handler.distribute_dofs(*fe);

    // Compute support points for the DoFs.
    FE_SimplexP<dim> fe_linear(1);
    MappingFE        mapping(fe_linear);
    DoFTools::map_dofs_to_support_points(mapping, dof_handler, support_points);
  }

  // Initialize the linear system.
  {
    DynamicSparsityPattern dsp(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, dsp);
    sparsity_pattern.copy_from(dsp);

    system_matrix.reinit(sparsity_pattern);
    system_rhs.reinit(dof_handler.n_dofs());
    solution.reinit(dof_handler.n_dofs());
  }
}

void
Poisson2D::assemble()
{
  const unsigned int dofs_per_cell = fe->dofs_per_cell;
  const unsigned int n_q           = quadrature->size();

  FEValues<dim> fe_values(*fe,
                          *quadrature,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);

  #ifdef NEUMANN
  // Since we need to compute integrals on the boundary for Neumann conditions,
  // we also need a FEValues object to compute quantities on boundary edges
  // (faces).
  FEFaceValues<dim> fe_values_boundary(*fe,
                                       *quadrature_boundary,
                                       update_values |
                                         update_quadrature_points |
                                         update_JxW_values);
#endif //NEUMANN

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_rhs(dofs_per_cell);

  std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

  system_matrix = 0.0;
  system_rhs    = 0.0;

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      fe_values.reinit(cell);

      cell_matrix = 0.0;
      cell_rhs    = 0.0;

      

      for (unsigned int q = 0; q < n_q; ++q)
        {
          #ifdef TRANSPORT_COEFFICIENT
          Vector<double> transport_coefficient_loc(dim);
          transport_coefficient.vector_value(fe_values.quadrature_point(q),
                                    transport_coefficient_loc);

          Tensor<1, dim> transport_coefficient_tensor;
          for (unsigned int d = 0; d < dim; ++d)
            transport_coefficient_tensor[d] = transport_coefficient_loc[d];
#endif

          #ifdef CONSERVATIVE_TRANSPORT_COEFFICIENT
          Vector<double> cons_transport_coefficient_loc(dim);
          cons_transport_coefficient.vector_value(fe_values.quadrature_point(q),
                                    cons_transport_coefficient_loc);

          Tensor<1, dim> cons_transport_coefficient_tensor;
          for (unsigned int d = 0; d < dim; ++d)
            cons_transport_coefficient_tensor[d] = cons_transport_coefficient_loc[d];
#endif

          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
              for (unsigned int j = 0; j < dofs_per_cell; ++j)
                {
                  cell_matrix(i, j) += fe_values.shape_grad(i, q) *
                                       fe_values.shape_grad(j, q) *
                                       fe_values.JxW(q);
                  #ifdef TRANSPORT_COEFFICIENT
                  cell_matrix(i, j) += scalar_product(transport_coefficient_tensor,
                                                      fe_values.shape_grad(j,q))
                                    * fe_values.shape_value(i,q)
                                    * fe_values.JxW(q);
                  #endif //TRANSPORT_COEFFICIENT

                  #ifdef CONSERVATIVE_TRANSPORT_COEFFICIENT
              cell_matrix(i, j) -= scalar_product(cons_transport_coefficient_tensor,
                                                  fe_values.shape_grad(i,q)) 
                                * fe_values.shape_value(j,q) 
                                * fe_values.JxW(q);

#endif //CONSERVATIVE_TRANSPORT_COEFFICIENT

#ifdef REACTION_COEFFICIENT
                  cell_matrix(i, j) +=
                    reaction_coefficient.value(
                      fe_values.quadrature_point(q)) * // sigma(x)
                    fe_values.shape_value(i, q) *      // phi_i
                    fe_values.shape_value(j, q) *      // phi_j
                    fe_values.JxW(q);                  // dx
#endif //REACTION_COEFFICIENT
                }

                              cell_rhs(i) += forcing_term.value(fe_values.quadrature_point(q)) *
                             fe_values.shape_value(i, q) * fe_values.JxW(q);
            }
        }

        #ifdef NEUMANN
  // If the cell is adjacent to the boundary...
  if (cell->at_boundary())
    {
      // ...we loop over its edges (referred to as faces in the deal.II
      // jargon).
      for (unsigned int face_number = 0; face_number < cell->n_faces();
           ++face_number)
        {
          // If current face lies on the boundary...
          if (cell->face(face_number)->at_boundary()&&(
              (cell->face(face_number)->boundary_id() == 0 && subdomain_id == 0) ||
              (cell->face(face_number)->boundary_id() == 1 && subdomain_id == 1) ))
            {
              // Functions::ConstantFunction<dim> one_function(1);
              // const Function<dim>* function_h = &one_function;
               fe_values_boundary.reinit(cell, face_number);

              for (unsigned int q = 0; q < quadrature_boundary->size(); ++q)
                for (unsigned int i = 0; i < dofs_per_cell; ++i)
                {
                  
                  #ifdef ROBIN
                  for (unsigned int j = 0; j < dofs_per_cell; ++j){
                    cell_matrix(i, j) +=
                     function_gamma.value(fe_values_boundary.quadrature_point(q)) * 
                     fe_values_boundary.shape_value(j, q) *
                     fe_values_boundary.shape_value(i, q) * 
                     fe_values_boundary.JxW(q);  
                  } 
                  #endif //ROBIN
                  cell_rhs(i) +=
                    function_h.value(
                    fe_values_boundary.quadrature_point(q)) * // h(xq)
                    fe_values_boundary.shape_value(i, q) *      // v(xq)
                    fe_values_boundary.JxW(q);                  // Jq wq
                }
            }
        }
    }
#endif //NEUMANN


      cell->get_dof_indices(dof_indices);

      system_matrix.add(dof_indices, cell_matrix);
      system_rhs.add(dof_indices, cell_rhs);
    }

  // Boundary conditions.
// Boundary conditions.
// Boundary conditions.
// Boundary conditions.
{
  std::map<types::global_dof_index, double>           boundary_values;
  std::map<types::boundary_id, const Function<dim> *> boundary_functions;

  // Define boundary functions for each face of each subdomain.
  Functions::ConstantFunction<dim> zero_function(0);
  Functions::ConstantFunction<dim> one_function(1);
  Functions::ConstantFunction<dim> function_bc_1(0);


  // Assign the boundary functions to the faces of the subdomain.
  if (subdomain_id == 0) {
    //boundary_functions[0] = &function_g; //     Left
    //boundary_functions[1] = &zero_function; //Right
    boundary_functions[2] = &function_g; //     Bottom
    boundary_functions[3] = &function_g; //     Top
  } else {
    //boundary_functions[0] = &zero_function; // Left
    //boundary_functions[1] = &function_g; // Right
    boundary_functions[2] = &function_g; // Bottom
    boundary_functions[3] = &function_g; // Top
  }

  // interpolate_boundary_values fills the boundary_values map.
  VectorTools::interpolate_boundary_values(dof_handler,
                                           boundary_functions,
                                           boundary_values);

  MatrixTools::apply_boundary_values(
    boundary_values, system_matrix, solution, system_rhs, false);
}



}

void
Poisson2D::solve()
{
  #ifdef CG
  std::cout << "===============================================" << std::endl;

  SolverControl solver_control(10000,  1e-6 * system_rhs.l2_norm() );

  // Since the system matrix is symmetric and positive definite, we solve the
  // system using the conjugate gradient method.
  SolverCG<Vector<double>> solver(solver_control);

  PreconditionSSOR preconditioner;
  preconditioner.initialize(
    system_matrix, PreconditionSOR<SparseMatrix<double>>::AdditionalData(1.0));

  solver.solve(system_matrix, solution, system_rhs, PreconditionIdentity());
  std::cout << "  " << solver_control.last_step() << " CG iterations"
            << std::endl;
  
  #endif //CG
#ifndef CG
  // Here we specify the maximum number of iterations of the iterative solver,
  // and its tolerance.
  SolverControl solver_control(5000, 1.0e-12*system_rhs.l2_norm());

  // Since the system matrix is symmetric and positive definite, we solve the
  // system using the conjugate gradient method.
  SolverGMRES<Vector<double>> solver(solver_control);
  PreconditionSOR preconditioner;
  preconditioner.initialize(
    system_matrix
  );

  std::cout << "  Solving the linear system" << std::endl;
  // We don't use any preconditioner for now, so we pass the identity matrix
  // as preconditioner.
  solver.solve(system_matrix, solution, system_rhs, preconditioner);
  std::cout << "  " << solver_control.last_step() << " GMRES iterations"
            << std::endl;
 #endif
}

void
Poisson2D::output(const unsigned int &iter) const
{
  DataOut<dim> data_out;

  data_out.add_data_vector(dof_handler, solution, "solution");
  data_out.build_patches();

  const std::string output_file_name = "output-" +
                                       std::to_string(subdomain_id) + "-" +
                                       std::to_string(iter) + ".vtk";
  std::ofstream output_file(output_file_name);
  data_out.write_vtk(output_file);
}

void
Poisson2D::apply_interface_dirichlet(const Poisson2D &other)
{
  const auto interface_map = compute_interface_map(other);

  // We use the interface map to build a boundary values map for interface DoFs.
  std::map<types::global_dof_index, double> boundary_values;
  for (const auto &dof : interface_map)
    boundary_values[dof.first] = other.solution[dof.second];

  // Then, we apply those boundary values.
  MatrixTools::apply_boundary_values(
    boundary_values, system_matrix, solution, system_rhs, false);
}

void
Poisson2D::apply_interface_neumann(Poisson2D &other)
{
  const auto interface_map = compute_interface_map(other);

  // We assemble the interface residual of the other subproblem. Indeed,
  // directly computing the normal derivative of the solution on the other
  // subdomain has extremely poor accuracy. This is due to the fact that the
  // trace of the derivative has very low regularity. Therefore, we compute the
  // (weak) normal derivative as the residual of the system of the other
  // subdomain, excluding interface conditions.
  Vector<double> interface_residual;
  other.assemble();
  interface_residual = other.system_rhs;
  interface_residual *= -1;
  other.system_matrix.vmult_add(interface_residual, other.solution);

  // Then, we add the elements of the residual corresponding to interface DoFs
  // to the system rhs for current subproblem.
  for (const auto &dof : interface_map)
    system_rhs[dof.first] -= interface_residual[dof.second];
}

std::map<types::global_dof_index, types::global_dof_index>
Poisson2D::compute_interface_map(const Poisson2D &other) const
{
  // Retrieve interface DoFs on the current and other subdomain.
  IndexSet current_interface_dofs;
  IndexSet other_interface_dofs;

  if (subdomain_id == 0)
    {
      current_interface_dofs =
        DoFTools::extract_boundary_dofs(dof_handler, ComponentMask(), {1});
      other_interface_dofs = DoFTools::extract_boundary_dofs(other.dof_handler,
                                                             ComponentMask(),
                                                             {0});
    }
  else
    {
      current_interface_dofs =
        DoFTools::extract_boundary_dofs(dof_handler, ComponentMask(), {0});
      other_interface_dofs = DoFTools::extract_boundary_dofs(other.dof_handler,
                                                             ComponentMask(),
                                                             {1});
    }

  // For each interface DoF on current subdomain, we find the corresponding one
  // on the other subdomain.
  std::map<types::global_dof_index, types::global_dof_index> interface_map;
  for (const auto &dof_current : current_interface_dofs)
    {
      const Point<dim> &p = support_points.at(dof_current);

      types::global_dof_index nearest = *other_interface_dofs.begin();
      for (const auto &dof_other : other_interface_dofs)
        {
          if (p.distance_square(other.support_points.at(dof_other)) <
              p.distance_square(other.support_points.at(nearest)))
            nearest = dof_other;
        }

      interface_map[dof_current] = nearest;
    }

  return interface_map;
}

void
Poisson2D::apply_relaxation(const Vector<double> &old_solution,
                            const double &        lambda)
{
  solution *= lambda;
  solution.add(1.0 - lambda, old_solution);
}