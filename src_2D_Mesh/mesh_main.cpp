#include <iostream>
#include <fstream>
#include <string>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>

int main()
{
    // Define the parameters
    double a = 0.75;  // Starting domain in x-direction
    double b = 0.0;  // Starting domain in y-direction
    double c = 1.0;  // Ending domain in x-direction
    double d = 1.0;  // Ending domain in y-direction
    double h = 0.05;  // Size of each element

    // Calculate the number of elements in x and y directions
    unsigned int num_elements_x = static_cast<unsigned int>((c - a) / h) + 1;
    unsigned int num_elements_y = static_cast<unsigned int>((d - b) / h) + 1;

    // Create the mesh
    dealii::Triangulation<2> mesh;
    std::vector<unsigned int> repetitions = {num_elements_x, num_elements_y};
    dealii::GridGenerator::subdivided_hyper_rectangle(mesh, repetitions, dealii::Point<2>(a, b), dealii::Point<2>(c, d), true);

    // Print the number of elements
    std::cout << "Number of elements = " << mesh.n_active_cells() << std::endl;

    // Write the mesh to file
    std::string mesh_file_name = "mesh-" + std::to_string(c) + ".vtk";
    dealii::GridOut grid_out;
    std::ofstream grid_out_file(mesh_file_name);
    grid_out.write_vtk(mesh, grid_out_file);
    std::cout << "Mesh saved to " << mesh_file_name << std::endl;

    return 0;
}
