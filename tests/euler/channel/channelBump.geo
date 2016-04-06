// Inviscid Channel-Flow with Smooth Gaussian Bump
H = 0.8;
L = 1.5;
NX = 6;
NY = 2;

nPoints = 250; // Geometry resolution for bottom bump surface

// Corner Points
Point(1) = {-L, 0, 0};
Point(2) = { L, 0, 0};
Point(3) = { L, H, 0};
Point(4) = {-L, H, 0};

// Bottom Surface: .0625*exp(-25*x^2)
pList[0] = 1;
x = -L;
dx = 2*L/(nPoints+1);
For i In {1:nPoints}
	x = x + dx;
	y = 0.0625*Exp(-25*x*x);
	pList[i] = newp;
	Point(pList[i]) = {x,y,0};
EndFor
pList[nPoints+1] = 2;
Spline(1) = pList[];

// Remaining Lines
Line(2) = {2,3};
Line(3) = {4,3};
Line(4) = {1,4};

Transfinite Line{1,3} = NX+1 Using Progression 1.0;
Transfinite Line{2,4} = NY+1 Using Progression 1.0;

Line Loop(1) = {1,2,-3,-4};
Plane Surface(1) = {1};
Transfinite Surface{1};
Recombine Surface{1};

Physical Line("Inlet") = {4};
Physical Line("Outlet") = {2};
Physical Line("Top") = {3};
Physical Line("Bottom") = {1};
Physical Surface("FLUID") = {1};

Mesh.ElementOrder = 10; // arbitrary orders available
