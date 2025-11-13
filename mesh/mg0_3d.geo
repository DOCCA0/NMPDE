Lx = 1.0; // Parallelepiped length along x-axis.
Ly = 1.0; // Parallelepiped length along y-axis.
Lz = 1.0; // Parallelepiped length along z-axis.
h = 0.0125; // Mesh size.

// Create one point at (0.75, 0, 0).
Point(1) = {0, 0, 0, h};

// Extrude the point along x to create one side. The Layers option indicates the 
// number of mesh subdivisions along the extrusion.
Extrude {Lx, 0, 0} { Point{1}; Layers{Round(Lx / h)}; }

// Extrude that side along y to create a rectangle.
Extrude {0, Ly, 0} { Line{1}; Layers{Round(Ly / h)}; }

// Extrude the rectangle along z to create the parallelepiped.
Extrude {0, 0, Lz} { Surface{1}; Layers{Round(Lz / h)}; }

// Define the tags.
Physical Surface(0) = {5};
Physical Surface(1) = {6};
Physical Surface(2) = {3};
Physical Surface(3) = {4};
Physical Surface(4) = {1};
Physical Surface(5) = {2};

Physical Volume(10) = {7};

// Generate a 3D mesh.
Mesh 3;

// Save mesh to file.
str_h = Sprintf("%f", h);
Save StrCat("m0", ".msh");
