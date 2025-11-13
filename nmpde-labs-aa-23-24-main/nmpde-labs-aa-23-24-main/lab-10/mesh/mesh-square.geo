Lx = 0.75; // Rectangle length along x-axis.
Ly = 1.0; // Rectangle length along y-axis.
h = 0.0125; // Mesh size.

// Create one point in the origin.
Point(1) = {0, 0, 0, h};

// Extrude the point along x to create one side. The Layers option indicates the 
// number of mesh subdivisions along the extrusion.
Extrude {1, 0, 0} { Point{1}; Layers{Lx / h}; }

// Extrude that side along y to create the rectangle.
Extrude {0, 1, 0} { Line{1};  Layers{Ly / h}; }

// Define the tags.
Physical Line(0) = {3};
Physical Line(1) = {4};
Physical Line(2) = {1};
Physical Line(3) = {2};

Physical Surface(10) = {5};

// Generate a 2D mesh.
Mesh 2;

// Save mesh to file.
str_h = Sprintf("%f", h);
Save StrCat("mesh-h", str_h, ".msh");
