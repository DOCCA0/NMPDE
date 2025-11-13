SetFactory("OpenCASCADE");
Lx = 2.0; // Length along the x-axis.
Ly = 2.0;  // Length along the y-axis.
h = 0.0125; // Mesh size.

// Create rectangle starting from (0,0) extending to (Lx, Ly).
Rectangle(1) = {-1, -1, 0, Lx, Ly};

// Set mesh size.
Mesh.CharacteristicLengthMax = h;

// Tag edges with specific identifiers.
Physical Line("BottomEdge") = {1}; // Bottom edge of the rectangle
Physical Line("RightEdge") = {2}; // Right edge of the rectangle
Physical Line("TopEdge") = {3}; // Top edge of the rectangle
Physical Line("LeftEdge") = {4}; // Left edge of the rectangle

// Tag the entire surface for element grouping.
Physical Surface("WholeSurface") = {1};

// Specify MSH file version for compatibility
Mesh.MshFileVersion = 2.2;

// Generate a 2D mesh.
Mesh 2;

// Save mesh to a file, naming it based on the mesh size.
Save "square_mesh.msh";
