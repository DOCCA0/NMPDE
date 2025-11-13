#include <iostream>
#include <fstream>
#include <string>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>

int main()
{
    // Define the parameters
    double start = 0.0;  // Starting domain
    double end = 1.0;    // Ending domain
    double h = 0.1;      // Size of each element

    // Calculate the number of elements
    unsigned int num_elements = static_cast<unsigned int>((end - start) / h) + 1;

    // Create the mesh
    dealii::Triangulation<1> mesh;
    dealii::GridGenerator::subdivided_hyper_cube(mesh, num_elements, start, end, true);

    // Print the number of elements
    std::cout << "Number of elements = " << mesh.n_active_cells() << std::endl;

    // Write the mesh to file
    std::string mesh_file_name = "mesh-" + std::to_string(h) + ".vtk";
    dealii::GridOut grid_out;
    std::ofstream grid_out_file(mesh_file_name);
    grid_out.write_vtk(mesh, grid_out_file);
    std::cout << "Mesh saved to " << mesh_file_name << std::endl;

    return 0;
}