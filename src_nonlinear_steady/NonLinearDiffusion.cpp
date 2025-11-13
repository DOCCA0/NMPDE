#include "NonLinearDiffusion.hpp"

void
NonLinearDiffusion::setup()
{
  // Create the mesh.
  {
    pcout << "Initializing the mesh" << std::endl;

    // First we read the mesh from file into a serial (i.e. not parallel)
    // triangulation.
    Triangulation<dim> mesh_serial;

    {
      GridIn<dim> grid_in;
      grid_in.attach_triangulation(mesh_serial);

      std::ifstream grid_in_file(mesh_file_name);
      grid_in.read_msh(grid_in_file);
    }

    // Then, we copy the triangulation into the parallel one.
    {
      GridTools::partition_triangulation(mpi_size, mesh_serial);
      const auto construction_data = TriangulationDescription::Utilities::
        create_description_from_triangulation(mesh_serial, MPI_COMM_WORLD);
      mesh.create_triangulation(construction_data);
    }

    // Notice that we write here the number of *global* active cells (across all
    // processes).
    pcout << "  Number of elements = " << mesh.n_global_active_cells()
          << std::endl;
  }

  pcout << "-----------------------------------------------" << std::endl;

  // Initialize the finite element space. This is the same as in serial codes.
  {
    pcout << "Initializing the finite element space" << std::endl;

    fe = std::make_unique<FE_SimplexP<dim>>(r);

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

    // We retrieve the set of locally owned DoFs, which will be useful when
    // initializing linear algebra classes.
    locally_owned_dofs = dof_handler.locally_owned_dofs();
    DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);

    pcout << "  Number of DoFs = " << dof_handler.n_dofs() << std::endl;
  }

  pcout << "-----------------------------------------------" << std::endl;

  // Initialize the linear system.
  {
    pcout << "Initializing the linear system" << std::endl;

    pcout << "  Initializing the sparsity pattern" << std::endl;

    // To initialize the sparsity pattern, we use Trilinos' class, that manages
    // some of the inter-process communication.
    TrilinosWrappers::SparsityPattern sparsity(locally_owned_dofs,
                                               MPI_COMM_WORLD);
    DoFTools::make_sparsity_pattern(dof_handler, sparsity);

    // After initialization, we need to call compress, so that all process
    // retrieve the information they need for the rows they own (i.e. the rows
    // corresponding to locally owned DoFs).
    sparsity.compress();

    // Then, we use the sparsity pattern to initialize the system matrix. Since
    // the sparsity pattern is partitioned by row, so will the matrix.
    pcout << "  Initializing the system matrix" << std::endl;
    jacobian_matrix.reinit(sparsity);

    // Finally, we initialize the right-hand side and solution vectors.
    pcout << "  Initializing the system right-hand side" << std::endl;
    residual_vector.reinit(locally_owned_dofs, MPI_COMM_WORLD);
    pcout << "  Initializing the solution vector" << std::endl;
    solution_owned.reinit(locally_owned_dofs, MPI_COMM_WORLD);
    solution.reinit(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
    delta_owned.reinit(locally_owned_dofs, MPI_COMM_WORLD);
  }
}

void
NonLinearDiffusion::assemble_system()
{
  const unsigned int dofs_per_cell = fe->dofs_per_cell;
  const unsigned int n_q           = quadrature->size();

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

  FEValues<dim> fe_values(*fe,
                          *quadrature,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_rhs(dofs_per_cell);

  std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

  jacobian_matrix = 0.0;
  residual_vector = 0.0;

  // We use these vectors to store the old solution (i.e. at previous Newton
  // iteration) and its gradient on quadrature nodes of the current cell.
  std::vector<double>         solution_loc(n_q);
  std::vector<Tensor<1, dim>> solution_gradient_loc(n_q);

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (!cell->is_locally_owned())
        continue;

      fe_values.reinit(cell);

      cell_matrix = 0.0;
      cell_rhs    = 0.0;

      // We need to compute the Jacobian matrix and the residual for current
      // cell. This requires knowing the value and the gradient of u^{(k)}
      // (stored inside solution) on the quadrature nodes of the current
      // cell. This can be accomplished through
      // FEValues::get_function_values and FEValues::get_function_gradients.
      fe_values.get_function_values(solution, solution_loc);
      fe_values.get_function_gradients(solution, solution_gradient_loc);

      for (unsigned int q = 0; q < n_q; ++q)
        {
          const double mu_0_loc = mu_0.value(fe_values.quadrature_point(q));
          const double mu_1_loc = mu_1.value(fe_values.quadrature_point(q));
          const double f_loc =
            forcing_term.value(fe_values.quadrature_point(q));

          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
              for (unsigned int j = 0; j < dofs_per_cell; ++j)
                {
                  cell_matrix(i, j) +=
                    (2.0 * mu_1_loc * solution_loc[q] *
                     fe_values.shape_value(j, q)) *
                    scalar_product(solution_gradient_loc[q],
                                   fe_values.shape_grad(i, q)) *
                    fe_values.JxW(q);

                  cell_matrix(i, j) +=
                    (mu_0_loc + mu_1_loc * solution_loc[q] * solution_loc[q]) *
                    scalar_product(fe_values.shape_grad(j, q),
                                   fe_values.shape_grad(i, q)) *
                    fe_values.JxW(q);
                }

              // -F(v)
              cell_rhs(i) +=
                f_loc * fe_values.shape_value(i, q) * fe_values.JxW(q);

              // a(u)(v)
              cell_rhs(i) -=
                (mu_0_loc + mu_1_loc * solution_loc[q] * solution_loc[q]) *
                scalar_product(solution_gradient_loc[q],
                               fe_values.shape_grad(i, q)) *
                fe_values.JxW(q);
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
              if (cell->face(face_number)->at_boundary() /* && (
                  (cell->face(face_number)->boundary_id() == 0)||
                  (cell->face(face_number)->boundary_id() == 1)) */)
                {
                  fe_values_boundary.reinit(cell, face_number);

                  for (unsigned int q = 0; q < quadrature_boundary->size(); ++q)
                    for (unsigned int i = 0; i < dofs_per_cell; ++i)
                      cell_rhs(i) +=
                        function_h.value(
                          fe_values_boundary.quadrature_point(q)) * // h(xq)
                        fe_values_boundary.shape_value(i, q) *      // v(xq)
                        fe_values_boundary.JxW(q);                  // Jq wq
                }
            }
        }
#endif //NEUMANN

      cell->get_dof_indices(dof_indices);

      jacobian_matrix.add(dof_indices, cell_matrix);
      residual_vector.add(dof_indices, cell_rhs);
    }

  jacobian_matrix.compress(VectorOperation::add);
  residual_vector.compress(VectorOperation::add);

  {
    std::map<types::global_dof_index, double> boundary_values;

    std::map<types::boundary_id, const Function<dim> *> boundary_functions;
    Functions::ZeroFunction<dim>                        zero_function;
    
    boundary_functions[0] = &zero_function;
    boundary_functions[1] = &zero_function;
    boundary_functions[2] = &zero_function;
    boundary_functions[3] = &zero_function;



    VectorTools::interpolate_boundary_values(dof_handler,
                                             boundary_functions,
                                             boundary_values);

    MatrixTools::apply_boundary_values(
      boundary_values, jacobian_matrix, delta_owned, residual_vector, true);
  }
}

void
NonLinearDiffusion::solve_system()
{
  SolverControl solver_control(10000, 1e-6 * residual_vector.l2_norm());

  SolverGMRES<TrilinosWrappers::MPI::Vector> solver(solver_control);
  TrilinosWrappers::PreconditionSSOR         preconditioner;
  preconditioner.initialize(
    jacobian_matrix, TrilinosWrappers::PreconditionSSOR::AdditionalData(1.0));

  solver.solve(jacobian_matrix, delta_owned, residual_vector, preconditioner);
  pcout << "   " << solver_control.last_step() << " GMRES iterations"
        << std::endl;
}

void
NonLinearDiffusion::solve_newton()
{
  pcout << "===============================================" << std::endl;

  const unsigned int n_max_iters        = 10000;
  const double       residual_tolerance = 1e-6;

  unsigned int n_iter        = 0;
  double       residual_norm = residual_tolerance + 1;

    {
    std::vector<bool> boundary_components(4, false); // Assuming there are 4 faces
    boundary_components[0] = true; // Face 0
    boundary_components[1] = true; // Face 1
    boundary_components[2] = true; // Face 2
    boundary_components[3] = true; // Face 3
    //boundary_components[4] = true; // Face 4
    //boundary_components[5] = true; // Face 5

    IndexSet dirichlet_dofs = DoFTools::extract_boundary_dofs(
      dof_handler,
      ComponentMask(boundary_components)
    );
    dirichlet_dofs = dirichlet_dofs & dof_handler.locally_owned_dofs();

    TrilinosWrappers::MPI::Vector vector_dirichlet(solution_owned);
    VectorTools::interpolate(dof_handler, function_g, vector_dirichlet);

    for (const auto &idx : dirichlet_dofs)
      solution_owned[idx] = vector_dirichlet[idx];

    solution_owned.compress(VectorOperation::insert);
    solution = solution_owned;
  }

  //  class InitialSolutionFunction : public Function<dim>
  // {
  // public:
  //   virtual double value(const Point<dim> &p, const unsigned int component = 0) const override
  //   {
  //     (void)component;
  //     return  0.0;
  //   }
  // };

  // InitialSolutionFunction initial_solution_function;
  // VectorTools::interpolate(dof_handler, initial_solution_function, solution_owned);
  // solution = solution_owned;


  while (n_iter < n_max_iters && residual_norm > residual_tolerance)
    {
      assemble_system();
      residual_norm = residual_vector.l2_norm();

      pcout << "Newton iteration " << n_iter << "/" << n_max_iters
            << " - ||r|| = " << std::scientific << std::setprecision(6)
            << residual_norm << std::flush;

      // We actually solve the system only if the residual is larger than the
      // tolerance.
      if (residual_norm > residual_tolerance)
        {
          solve_system();

          solution_owned += delta_owned;
          solution = solution_owned;
        }
      else
        {
          pcout << " < tolerance" << std::endl;
        }

      ++n_iter;
    }

  pcout << "===============================================" << std::endl;
}

void
NonLinearDiffusion::output() const
{
  DataOut<dim> data_out;
  data_out.add_data_vector(dof_handler, solution, "u");

  std::vector<unsigned int> partition_int(mesh.n_active_cells());
  GridTools::get_subdomain_association(mesh, partition_int);
  const Vector<double> partitioning(partition_int.begin(), partition_int.end());
  data_out.add_data_vector(partitioning, "partitioning");

  data_out.build_patches();

  const std::string output_file_name = "output-nonlineardiffusion";
  data_out.write_vtu_with_pvtu_record("./",
                                      output_file_name,
                                      0,
                                      MPI_COMM_WORLD);

  pcout << "Output written to " << output_file_name << "." << std::endl;

  pcout << "===============================================" << std::endl;
}

#ifdef CONVERGENCE
double
NonLinearDiffusion::compute_error(const VectorTools::NormType &norm_type)
{
  FE_SimplexP<dim> fe_linear(1);
  MappingFE        mapping(fe_linear);

  const QGaussSimplex<dim> quadrature_error = QGaussSimplex<dim>(r + 2);


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

#endif //CONVERGENCE