// Gmsh project created on Wed Jul  1 14:24:32 2015
xmin = -4;  xmax = 4;
ymin = -.5;  ymax = .5; // 3
zmin = -4;  zmax = 4;

nx = 36; // 6
ny = 4; // 6
nz = 36; // 6

Point(0) = {xmin,ymin,zmin};
Point(1) = {xmax,ymin,zmin};
Point(2) = {xmax,ymax,zmin};
Point(3) = {xmin,ymax,zmin};

Point(4) = {xmin,ymin,zmax};
Point(5) = {xmax,ymin,zmax};
Point(6) = {xmax,ymax,zmax};
Point(7) = {xmin,ymax,zmax};

Line(0) = {0,1};
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,0}; 

Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,4};

Line(8) = {0,4};
Line(9) = {1,5};
Line(10) = {2,6};
Line(11) = {3,7};

Transfinite Line{0, -2, 4, -6} = nx;
Transfinite Line{1, -3, 5, -7} = ny;
Transfinite Line{8, 9, 10, 11} = nz;

Line Loop(1) = {0,1,2,3};
Line Loop(2) = {4,5,6,7};
Line Loop(3) = {3,8,-7,-11};
Line Loop(4) = {-1,9,5,-10};
Line Loop(5) = {0,9,-4,-8};
Line Loop(6) = {-2,10,6,-11};

Plane Surface(1) = {-1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Plane Surface(5) = {-5};
Plane Surface(6) = {6};

Transfinite Surface{1,2,3,4,5,6};
Recombine Surface{1,2,3,4,5,6};

Surface Loop(1) = {1, 2, 3, 4, 5, 6};
Volume(1) = {1};

Transfinite Volume{1};
Recombine Volume{1};

Physical Surface("Bottom") = {1};
Physical Surface("Top") = {2};
Physical Surface("Left") = {3};
Physical Surface("Right") = {4};
Physical Surface("Front") = {5};
Physical Surface("Back") = {6};
Physical Volume("FLUID") = {1};
