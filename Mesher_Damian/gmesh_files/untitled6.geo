cl=.05;
l=16.0;
h=10.0;

Point(1) = {0, 0, 0, cl};
Point(2) = {l, 0, 0, cl};
Point(3) = {l, h, 0, cl};
Point(4) = {0, h, 0, cl};
Point(5) = {6, 0, 0, cl};
Point(6) = {8, 2, 0, cl};
Point(7) = {10, 0, 0, cl};

Line(1) = {1, 5};
Line(2) = {5, 6};
Line(3) = {6, 7};
Line(4) = {7, 2};
Line(5) = {2, 3};
Line(6) = {3, 4};
Line(7) = {4, 1};

Line Loop(1) = {1, 2, 3, 4, 5, 6, 7};
Plane Surface(1) = {1};


//Transfinite Surface {1} Alternated;
Recombine Surface {1};

Physical Surface(1000) = {1};

Physical Line(1) = {6};
Physical Line(2) = {5, 7};

Mesh.SecondOrderIncomplete = 1;



