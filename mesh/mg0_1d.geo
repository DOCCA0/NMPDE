L = 0.75; // Line length.
h = 0.0125; // Mesh size.

// Create one point at (0, 0, 0).
Point(1) = {0, 0, 0, h};

// Extrude the point along x to create a line. The Layers option indicates the 
// number of mesh subdivisions along the extrusion.
Extrude {L, 0, 0} { Point{1}; Layers{Round(L / h)}; }

// Define the tags.
Physical Line(10) = {1};

// Generate a 1D mesh.
Mesh 1;

// Save mesh to file.
str_h = Sprintf("%f", h);
Save StrCat("m0", ".msh");
