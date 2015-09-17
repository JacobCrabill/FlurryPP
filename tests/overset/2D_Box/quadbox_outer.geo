NX = 21;
NY = 21;

xmin = -5;
xmax = 5;
ymin = -5;
ymax = 5;

Point(1) = {xmin,ymin,0};
Point(2) = {xmax,ymin,0};
Point(3) = {xmax,ymax,0};
Point(4) = {xmin,ymax,0};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Transfinite Line{1,3} = NX;
Transfinite Line{2,4} = NY;

Line Loop(1) = {1,2,3,4};

Plane Surface(1) = {1};

Transfinite Surface{1};

Recombine Surface{1};

Physical Line("Left") = {4};
Physical Line("Right") = {2};
Physical Line("Top") = {3};
Physical Line("Bottom") = {1};
Physical Surface("FLUID") = {1};
