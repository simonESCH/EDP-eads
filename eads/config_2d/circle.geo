// Gmsh project created on Mon Apr 16 13:41:19 2018
SetFactory("OpenCASCADE");

Point(1) = {0, 0, 0, 1.0};// centre du cercle
Point(2) = {1, 0, 0, 1.0};
Point(3) = {0, 1, 0, 1.0};
Point(4) = {-1, 0, 0, 1.0};

Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {2, 1, 4};

Line Loop(1) = {3, -2, -1};
Plane Surface(1) = {1};

Physical Line("Dirichlet") = {1};
Physical Line("Neumann") = {2};
Physical Line("Robin") = {3};
