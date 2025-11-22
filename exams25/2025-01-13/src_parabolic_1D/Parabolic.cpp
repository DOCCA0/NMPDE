#include "Parabolic.hpp"

void
Parabolic::setup()
{
  // Create the mesh.
  {
    pcout << "Initializing the mesh" << std::endl;

    Triangulation<dim> mesh_serial;

    GridIn<dim> grid_in;
    grid_in.attach_triangulation(mesh_serial);

    std::ifstream grid_in_file(mesh_file_name);
    grid_in.read_msh(grid_in_file);

    GridTools::partition_triangulation(mpi_size, mesh_serial);
    const auto construction_data = TriangulationDescription::Utilities::
      create_description_from_triangulation(mesh_serial, MPI_COMM_WORLD);
    mesh.create_triangulation(construction_data);

    pcout << "  Number of elements = " << mesh.n_global_active_cells()
          << std::endl;
  }

  pcout << "-----------------------------------------------" << std::endl;

  // Initialize the finite element space.
  {
    pcout << "Initializing the finite element space" << std::endl;

    fe = std::make_unique<FE_Q<dim>>(r);

    pcout << "  Degree                     = " << fe->degree << std::endl;
    pcout << "  DoFs per cell              = " << fe->dofs_per_cell
          << std::endl;

    quadrature = std::make_unique<QGaussSimplex<dim>>(r + 1);

    pcout << "  Quadrature points per cell = " << quadrature->size()
          << std::endl;

#ifdef NEUMANN
    quadrature_boundary = std::make_unique<QGaussSimplex<dim - 1>>(r + 1);

    std::cout << "  Quadrature points per boundary cell = "
              << quadrature_boundary->size() << std::endl;
#endif //NEUMANN
  }

  pcout << "-----------------------------------------------" << std::endl;

  // Initialize the DoF handler.
  {
    pcout << "Initializing the DoF handler" << std::endl;

    dof_handler.reinit(mesh);
    dof_handler.distribute_dofs(*fe);

    locally_owned_dofs = dof_handler.locally_owned_dofs();
    DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);

    pcout << "  Number of DoFs = " << dof_handler.n_dofs() << std::endl;
  }

  pcout << "-----------------------------------------------" << std::endl;

  // Initialize the linear system.
  {
    pcout << "Initializing the linear system" << std::endl;

    pcout << "  Initializing the sparsity pattern" << std::endl;

    TrilinosWrappers::SparsityPattern sparsity(locally_owned_dofs,
                                               MPI_COMM_WORLD);
    DoFTools::make_sparsity_pattern(dof_handler, sparsity);
    sparsity.compress();

    pcout << "  Initializing the matrices" << std::endl;
    mass_matrix.reinit(sparsity);
    stiffness_matrix.reinit(sparsity);
    lhs_matrix.reinit(sparsity);
    rhs_matrix.reinit(sparsity);

    pcout << "  Initializing the system right-hand side" << std::endl;
    system_rhs.reinit(locally_owned_dofs, MPI_COMM_WORLD);
    pcout << "  Initializing the solution vector" << std::endl;
    solution_owned.reinit(locally_owned_dofs, MPI_COMM_WORLD);
    solution.reinit(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
  }
}

void
Parabolic::assemble_matrices(const double &time)
{
  pcout << "===============================================" << std::endl;
  pcout << "Assembling the system matrices" << std::endl;

  const unsigned int dofs_per_cell = fe->dofs_per_cell;
  const unsigned int n_q           = quadrature->size();

  FEValues<dim> fe_values(*fe,
                          *quadrature,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);
  #ifdef ROBIN
  // Since we need to compute integrals on the boundary for Neumann conditions,
  // we also need a FEValues object to compute quantities on boundary edges
  // (faces).
  FEFaceValues<dim> fe_values_boundary(*fe,
                                       *quadrature_boundary,
                                       update_values |
                                         update_quadrature_points |
                                         update_JxW_values);
#endif 

  FullMatrix<double> cell_mass_matrix(dofs_per_cell, dofs_per_cell);
  FullMatrix<double> cell_stiffness_matrix(dofs_per_cell, dofs_per_cell);

  std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

  mass_matrix      = 0.0;
  stiffness_matrix = 0.0;

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (!cell->is_locally_owned())
        continue;

      fe_values.reinit(cell);

      cell_mass_matrix      = 0.0;
      cell_stiffness_matrix = 0.0;

      for (unsigned int q = 0; q < n_q; ++q)
        {
#ifdef TRANSPORT_COEFFICIENT
          Vector<double> transport_coefficient_loc(dim);
          transport_coefficient.vector_value(fe_values.quadrature_point(q),
                                    transport_coefficient_loc);

          Tensor<1, dim> transport_coefficient_tensor;
          for (unsigned int d = 0; d < dim; ++d)
            transport_coefficient_tensor[d] = transport_coefficient_loc[d];
#endif //TRANSPORT_COEFFICIENT

#ifdef CONSERVATIVE_TRANSPORT_COEFFICIENT
          Vector<double> cons_transport_coefficient_loc(dim);
          cons_transport_coefficient.vector_value(fe_values.quadrature_point(q),
                                    cons_transport_coefficient_loc);

          Tensor<1, dim> cons_transport_coefficient_tensor;
          for (unsigned int d = 0; d < dim; ++d)
            cons_transport_coefficient_tensor[d] = cons_transport_coefficient_loc[d];
#endif //TRANSPORT_COEFFICIENT
          // Evaluate coefficients on this quadrature node.
#ifdef MUCOEFFICIENT
          const double mu_loc = mu.value(fe_values.quadrature_point(q));
#endif //MUCOEFFICIENT

          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
              for (unsigned int j = 0; j < dofs_per_cell; ++j)
                {
                  cell_mass_matrix(i, j) += fe_values.shape_value(i, q) *
                                            fe_values.shape_value(j, q) /
                                            deltat * fe_values.JxW(q);

#ifdef MUCOEFFICIENT
                  cell_stiffness_matrix(i, j) +=
                    mu_loc * fe_values.shape_grad(i, q) *
                    fe_values.shape_grad(j, q) * fe_values.JxW(q);
#endif //MUCOEFFICIENT

#ifdef TRANSPORT_COEFFICIENT
                  cell_stiffness_matrix(i, j) += scalar_product(transport_coefficient_tensor,
                                                      fe_values.shape_grad(j,q))
                                    * fe_values.shape_value(i,q)
                                    * fe_values.JxW(q);
#endif //TRANSPORT_COEFFICIENT

#ifdef CONSERVATIVE_TRANSPORT_COEFFICIENT
              cell_stiffness_matrix(i, j) -= scalar_product(cons_transport_coefficient_tensor,
                                                  fe_values.shape_grad(i,q)) 
                                * fe_values.shape_value(j,q) 
                                * fe_values.JxW(q);

#endif //CONSERVATIVE_TRANSPORT_COEFFICIENT

#ifdef REACTION_COEFFICIENT
                  cell_stiffness_matrix(i, j) +=
                    reaction_coefficient.value(
                    fe_values.quadrature_point(q)) * // sigma(x)
                    fe_values.shape_value(i, q) *      // phi_i
                    fe_values.shape_value(j, q) *      // phi_j
                    fe_values.JxW(q);                  // dx
#endif //REACTION_COEFFICIENT

                }
            }
        }
      #ifdef ROBIN
  // If the cell is adjacent to the boundary...
  if (cell->at_boundary())
    {
      // ...we loop over its edges (referred to as faces in the deal.II
      // jargon).
      for (unsigned int face_number = 0; face_number < cell->n_faces();
           ++face_number)
        {
          // If current face lies on the boundary, and its boundary ID (or
          // tag) is that of one of the Neumann boundaries, we assemble the
          // boundary integral.
          if (cell->face(face_number)->at_boundary() &&
              (cell->face(face_number)->boundary_id() == 1))
            {
              fe_values_boundary.reinit(cell, face_number);

              for (unsigned int q = 0; q < quadrature_boundary->size(); ++q)
                for (unsigned int i = 0; i < dofs_per_cell; ++i)
                { 
                  for (unsigned int j = 0; j < dofs_per_cell; ++j){
                    function_gamma.set_time(time);
                    cell_stiffness_matrix(i, j) +=
                     function_gamma.value(fe_values_boundary.quadrature_point(q)) * 
                     fe_values_boundary.shape_value(j, q) *
                     fe_values_boundary.shape_value(i, q) * 
                     fe_values_boundary.JxW(q);  
                  }        

                }
            }
        }
    }
#endif //ROBIN

      cell->get_dof_indices(dof_indices);

      mass_matrix.add(dof_indices, cell_mass_matrix);
      stiffness_matrix.add(dof_indices, cell_stiffness_matrix);
    }

  mass_matrix.compress(VectorOperation::add);
  stiffness_matrix.compress(VectorOperation::add);

  // We build the matrix on the left-hand side of the algebraic problem (the one
  // that we'll invert at each timestep).
  lhs_matrix.copy_from(mass_matrix);
  lhs_matrix.add(theta, stiffness_matrix);

  // We build the matrix on the right-hand side (the one that multiplies the old
  // solution un).
  rhs_matrix.copy_from(mass_matrix);
  rhs_matrix.add(-(1.0 - theta), stiffness_matrix);
}

void
Parabolic::assemble_rhs(const double &time)
{
  const unsigned int dofs_per_cell = fe->dofs_per_cell;
  const unsigned int n_q           = quadrature->size();

  FEValues<dim> fe_values(*fe,
                          *quadrature,
                          update_values | update_quadrature_points |
                            update_JxW_values);

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

  Vector<double> cell_rhs(dofs_per_cell);

  std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

  system_rhs = 0.0;

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (!cell->is_locally_owned())
        continue;

      fe_values.reinit(cell);

      cell_rhs = 0.0;

      for (unsigned int q = 0; q < n_q; ++q)
        {

          // We need to compute the forcing term at the current time (tn+1) and
          // at the old time (tn). deal.II Functions can be computed at a
          // specific time by calling their set_time method.

          // Compute f(tn+1)
          forcing_term.set_time(time);
          const double f_new_loc =
            forcing_term.value(fe_values.quadrature_point(q));

          // Compute f(tn)
          forcing_term.set_time(time - deltat);
          const double f_old_loc =
            forcing_term.value(fe_values.quadrature_point(q));

          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
              cell_rhs(i) += (theta * f_new_loc + (1.0 - theta) * f_old_loc) *
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
              // If current face lies on the boundary, and its boundary ID (or
              // tag) is that of one of the Neumann boundaries, we assemble the
              // boundary integral.
              if (cell->face(face_number)->at_boundary() &&
                  (cell->face(face_number)->boundary_id() == 1) ||
                  (cell->face(face_number)->boundary_id() == 3))
                {
                  fe_values_boundary.reinit(cell, face_number);
                  

                  for (unsigned int q = 0; q < quadrature_boundary->size(); ++q)
                  {
                      function_h.set_time(time);
                      const double h_new_loc =function_h.value(fe_values_boundary.quadrature_point(q));
                      function_h.set_time(time-deltat);
                      const double h_old_loc =function_h.value(fe_values_boundary.quadrature_point(q));
                    for (unsigned int i = 0; i < dofs_per_cell; ++i)
                      cell_rhs(i) += (theta * h_new_loc + (1.0 - theta) * h_old_loc) *
                                  fe_values_boundary.shape_value(i, q) * fe_values_boundary.JxW(q);
                    }
                }
            }
        }
#endif //NEUMANN

      cell->get_dof_indices(dof_indices);
      system_rhs.add(dof_indices, cell_rhs);
    }

  system_rhs.compress(VectorOperation::add);

  // Add the term that comes from the old solution.
  rhs_matrix.vmult_add(system_rhs, solution_owned);

  //We apply Dirichlet boundary conditions to the algebraic system.
  {
    std::map<types::global_dof_index, double> boundary_values;

    Functions::ZeroFunction<dim> zero_function(dim);

    std::map<types::boundary_id, const Function<dim> *> boundary_functions;
    //for (unsigned int i = 0; i < 4; ++i)
    //  boundary_functions[i] = &function_g;
    //boundary_functions[0] = &zero_function;
    function_g.set_time(time);
    boundary_functions[0] = &function_g;
    //boundary_functions[1] = &function_g;
    //boundary_functions[2] = &function_g;




    VectorTools::interpolate_boundary_values(dof_handler,
                                             boundary_functions,
                                             boundary_values);

    MatrixTools::apply_boundary_values(
      boundary_values, lhs_matrix, solution_owned, system_rhs, false);
  }
}

void
Parabolic::solve_time_step()
{
  SolverControl solver_control(10000, 1e-12 * system_rhs.l2_norm());

  //SolverCG<TrilinosWrappers::MPI::Vector> solver(solver_control);
  SolverGMRES<TrilinosWrappers::MPI::Vector> solver(solver_control);
  TrilinosWrappers::PreconditionSSOR      preconditioner;
  preconditioner.initialize(
    lhs_matrix, TrilinosWrappers::PreconditionSSOR::AdditionalData(1.0));

  solver.solve(lhs_matrix, solution_owned, system_rhs, preconditioner);
  pcout << "  " << solver_control.last_step() << " CG iterations" << std::endl;

  solution = solution_owned;
}

void
Parabolic::output(const unsigned int &time_step) const
{
  DataOut<dim> data_out;
  data_out.add_data_vector(dof_handler, solution, "u");

  std::vector<unsigned int> partition_int(mesh.n_active_cells());
  GridTools::get_subdomain_association(mesh, partition_int);
  const Vector<double> partitioning(partition_int.begin(), partition_int.end());
  data_out.add_data_vector(partitioning, "partitioning");

  data_out.build_patches();

  data_out.write_vtu_with_pvtu_record(
    "./", "output", time_step, MPI_COMM_WORLD, 3);
}

void
Parabolic::solve()
{
  pcout << "===============================================" << std::endl;

  time = 0.0;
  assemble_matrices(time);


  // Apply the initial condition.
  {
    pcout << "Applying the initial condition" << std::endl;

    //exact_solution.set_time(time);
    exact_solution.set_time(time);
    VectorTools::interpolate(dof_handler, u_0, solution_owned);
    solution = solution_owned;

    // Output the initial solution.
    output(0);
    pcout << "-----------------------------------------------" << std::endl;
  }

  unsigned int time_step = 0;

  while (time < T - 0.5 * deltat)
    {
      time += deltat;
      ++time_step;

      pcout << "n = " << std::setw(3) << time_step << ", t = " << std::setw(5)
            << time << ":" << std::flush;
      
      assemble_matrices(time);
      assemble_rhs(time);
      solve_time_step();
      output(time_step);
    }
}

double
Parabolic::compute_error(const VectorTools::NormType &norm_type)
{
  FE_Q<dim> fe_linear(1);
  MappingFE        mapping(fe_linear);

  const QGaussSimplex<dim> quadrature_error = QGaussSimplex<dim>(r + 2);

  exact_solution.set_time(time);

  Vector<double> error_per_cell;
  VectorTools::integrate_difference(mapping,
                                    dof_handler,
                                    solution,
                                    exact_solution,
                                    error_per_cell,
                                    quadrature_error,
                                    norm_type);

  const double error =
  VectorTools::compute_global_error(mesh, error_per_cell, norm_type);

  return error;
}

