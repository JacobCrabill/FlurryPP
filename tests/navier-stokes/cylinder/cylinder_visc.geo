inSize = 6;  // farfield mesh size, inlet side
outSize = 5;  // farfield mesh size, inlet side
nfSize = .2;  // near-field mesh size

nRadial = 13;
nTagent = 20;

// Interior box of mesh
Point(1) = {-2, -2, 0, nfSize};
Point(2) = { 9, -2, 0, nfSize};
Point(3) = { 9,  2, 0, nfSize};
Point(4) = {-2,  2, 0, nfSize};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Circle & surrounding structured-quad region
Point(5) = {0,   0, 0};
Point(6) = {0,  .5, 0};
Point(7) = {0, -.5, 0};
Point(8) = {0,  1.25, 0};
Point(9) = {0, -1.25, 0};
Circle(5) = {7, 5, 6};
Circle(6) = {6, 5, 7};
Circle(7) = {8, 5, 9};
Circle(8) = {9, 5, 8};
Line(9)  = {6, 8};
Line(10) = {7, 9};
Transfinite Line {5,6,7,8} = nTagent;
Transfinite Line {9,10} = nRadial Using Progression 1.2;

// Exterior (bounding box) of mesh
Point(10) = {-30, -30, 0, inSize};
Point(11) = { 50, -30, 0, outSize};
Point(12) = { 50,  30, 0, outSize};
Point(13) = {-30,  30, 0, inSize};
Line(11) = {10, 11};
Line(12) = {11, 12};
Line(13) = {12, 13};
Line(14) = {13, 10};

Line Loop(1) = {1, 2, 3, 4, 7, 8}; // Inner Exterior
Line Loop(2) = {10, 8, -9, -5}; // RH side of quad region - note ordering
Line Loop(3) = {7, -10, -6, 9}; // LH side of quad region - note ordering
Line Loop(4) = {11,12,13,14};   // Outer Exterior

Plane Surface(1) = {1}; // Outer unstructured region
Plane Surface(2) = {2}; // RH inner structured region
Plane Surface(3) = {3}; // LH inner structured region
Plane Surface(4) = {4,-1}; // Outer outer region

Transfinite Surface{2,3};

Recombine Surface {2,3};
Recombine Surface {1,4};

Physical Line("Bottom") = {11};
Physical Line("Right")  = {12};
Physical Line("Top")    = {13};
Physical Line("Left")   = {14};
Physical Line("Circle") = {5,6};

Physical Surface("FLUID") = {1,2,3,4};
