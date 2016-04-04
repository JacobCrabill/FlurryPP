xmin = -1;
xmax = 1;
ymin = 0;
ymax = 1;

NX = 4;
NY = 2;

Point (1) = {xmin, ymin, 0};
Point (2) = {xmax, ymin, 0};
Point (3) = {xmax, ymax, 0};
Point (4) = {xmin, ymax, 0};

Line (1) = {1,2};
Line (2) = {2,3};
Line (3) = {3,4};
Line (4) = {4,1};

Transfinite Line {1,-3} = NX+1;
Transfinite Line {2,-4} = NY+1;

Line Loop (1) = {1,2,3,4};

Plane Surface (1) = {1};

Transfinite Surface {1};
Recombine Surface {1};

Physical Surface ("FLUID") = {1};
Physical Line ("TOP") = {3};
Physical Line ("BOTTOM") = {1};
Physical Line ("LEFT") = {4};
Physical Line ("RIGHT") = {2};
