cl1 = .1;
cl2 = 2;

R1 = .5;
R2 = 3;
N_Tan = 20;
N_Out = 10;
Prog = 1.3;

Point(1) = {0,0,0,cl1};
Point(2) = {0,-R1,0,cl1};
Point(3) = {0,R1,0,cl1};
Point(4) = {0,-R2,0,cl2};
Point(5) = {0,R2,0,cl2};

Circle(1) = {2,1,3};
Circle(2) = {3,1,2};

Circle(3) = {4,1,5};
Circle(4) = {5,1,4};

Line(5) = {2,4};
Line(6) = {3,5};

Line Loop(1) = {5,3,-6,-1};
Line Loop(2) = {-2,6,4,-5};

Plane Surface(1) = {1};
Plane Surface(2) = {2};

Transfinite Line{1,2,3,4} = N_Tan;
Transfinite Line{5,6} = N_Out Using Progression Prog;

Transfinite Surface{1,2};
Recombine Surface {1,2};

Physical Surface("FLUID") = {1,2};
Physical Line("Wall") = {1,2};
Physical Line("Overset") = {3,4};
