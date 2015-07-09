// ---- 3D Circular Cylinder Gmsh Tutorial ----
// cylinder_extrude_3D.geo
// Creates a mesh of hexahedrons with an inner 'structured' region and 
// an outer unstructured region
//
// Created July 03, 2015 by Jacob Crabill
// Aerospace Computing Lab, Stanford University
// --------------------------------------------

// Gmsh allows variables; these will be used to set desired
// element sizes at various Points
cl1 = 8;     // Farthest exterior size   6
cl2 = .02;   // Near-body size    .02
cl3 = 6;     // Outlet/wake size   5
cl4 = .4;    // Nearfield-box size 4

R2 = 2;  // Outer overset-boundary radius

// Interior box of mesh
//Point(1) = {-2.5, -2.5, 0, cl4};
//Point(2) = { 3, -2.5, 0, cl4};
//Point(3) = { 3,  2.5, 0, cl4};
//Point(4) = {-2.5,  2.5, 0, cl4};
//Line(1) = {1, 2};
//Line(2) = {2, 3};
//Line(3) = {3, 4};
//Line(4) = {4, 1};

// Circle & surrounding structured-quad region
Point(5) = {0,   0, 0, cl2};
Point(6) = {0,  .5, 0, cl2};
Point(7) = {0, -.5, 0, cl2};
Point(8) = {0,  R2, 0, cl2};
Point(9) = {0, -R2, 0, cl2};
Point(10) = {.5, 0, cl2};
Point(11) = {-.5, 0, cl2};
Point(12) = {R2, 0, cl2};
Point(13) = {-R2, 0, cl2};

Circle(5) = {7, 5, 10};
Circle(6) = {10, 5, 6};
Circle(7) = {6, 5, 11};
Circle(8) = {11, 5, 7};
Circle(9) = {9, 5, 12};
Circle(10) = {12, 5, 8};
Circle(11) = {8, 5, 13};
Circle(12) = {13, 5, 9};
Line(13)  = {6, 8};
Line(14) = {7, 9};
Line(15) = {10,12};
Line(16) = {11,13};
Transfinite Line {5,6,7,8,9,10,11,12} = 10; // We want xx points along each of these lines
Transfinite Line {13,14,15,16} = 8 Using Progression 1.1;    // And xx points along each of these lines

// Exterior (bounding box) of mesh
//Point(14) = {-30, -30, 0, cl1};
//Point(15) = { 50, -30, 0, cl3};
//Point(16) = { 50,  30, 0, cl3};
//Point(17) = {-30,  30, 0, cl1};
//Line(17) = {14, 15};
//Line(18) = {15, 16};
//Line(19) = {16, 17};
//Line(20) = {17, 14};

// Each region which to be independantly meshed must have a line loop
// Regions which will be mesh with Transfinite Surface must have 4 lines
// and be labeled in CCW order, with the correct orientation of each edge
//Line Loop(1) = {1,2,3,4}; // Inner Exterior
Line Loop(2) = {-5,14,9,-15};  // RH btm side of quad region - note ordering
Line Loop(3) = {-6,15,10,-13}; // RH top side of quad region - note ordering
Line Loop(4) = {-7,13,11,-16}; // LH top side of quad region - note ordering
Line Loop(5) = {-8,16,12,-14}; // LH btm side of quad region - note ordering
//Line Loop(6) = {17,18,19,20};   // Outer Exterior

//Plane Surface(1) = {1,-2,-3,-4,-5}; // Outer unstructured region
Plane Surface(2) = {2}; // RH btm inner structured region
Plane Surface(3) = {3}; // RH top inner structured region
Plane Surface(4) = {4}; // LH top inner structured region
Plane Surface(5) = {5}; // LH btm inner structured region
//Plane Surface(6) = {6,-1}; // Outer outer region

// Mesh these surfaces in a structured manner
Transfinite Surface{2,3,4,5};

// Turn inner into quads
Recombine Surface {2,3,4,5};
// Turn outer region into unstructured quads (optional?)
//Recombine Surface {1,6};

// Extrude in Z-direction
Extrude {0,0,.5} {
  Surface{2,3,4,5}; 
  Layers{2};  
  Recombine;
}

// Apply Boundadry Conditions
// (Plane Surface #'s found by using GUI, under 'Tools->Visibility')
//Physical Surface ("Left") = {235};
//Physical Surface ("Right") = {227};
//Physical Surface ("Top") = {231};
//Physical Surface ("Bottom") = {223};
//Physical Surface ("Front") = {252,122,166,188,210,144};
//Physical Surface ("Wall") = {93,77,61,109};

Physical Surface ("Overset") = {33, 55, 77, 99};
Physical Surface ("Back") = {2,3,4,5};
Physical Surface ("Front") = {38, 60, 82, 104};
Physical Surface ("Wall") = {91, 25, 69, 47};

// IMPORTANT: "FLUID" MUST contain all fluid surfaces(2D)/volumes(3D)
Physical Volume ("FLUID") = {1,2,3,4};

