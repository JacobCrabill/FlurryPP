// Gmsh project created on Mon Jul 27 19:34:52 2015
nPtsY  = 41;
nPtsX1 = 20;
nPtsX2 = 37;

xmin = -1.25;
xmid = 0;
xmax = 1;
ymin = 0;
ymax = 2;

xprog1 = 1.49;
xprog2 = 1.22;
yprog = 1.21;

Point(1) = {xmin,ymin,0};
Point(2) = {xmid,ymin,0};
Point(3) = {xmid,ymax,0};
Point(4) = {xmin,ymax,0};
Point(5) = {xmax,ymin,0};
Point(6) = {xmax,ymax,0};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {1,4};

Line(5) = {2,5};
Line(6) = {5,6};
Line(7) = {6,3};

Line Loop(1) = {1,2,3,-4};
Line Loop(2) = {5,6,7,-2};

Transfinite Line{-1,3} = nPtsX1 Using Progression xprog1;
Transfinite Line{5,-7} = nPtsX2 Using Progression xprog2;
Transfinite Line{2,4,6} = nPtsY Using Progression yprog;

Plane Surface (1) = {1};
Plane Surface (2) = {2};

Transfinite Surface{1:2};
Recombine Surface{1:2};

Physical Surface("FLUID") = {1,2};
Physical Line("SYMMETRY") = {1};
Physical Line("WALL") = {5};
Physical Line("INLET") = {4};
Physical Line("OUTLET") = {6};
Physical Line("FARFIELD") = {3,7};

