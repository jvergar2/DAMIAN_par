cli=.02;
cl=0.2;
l=20.0;
h=15.0;

Point(1) = {0, 0, 0, cli};
Point(2) = {2, 0, 0, cl};
Point(3) = {7.5, 0, 0, cl};
Point(4) = {10, 2.5, 0, cl};
Point(5) = {12.5, 0, 0, cl};
Point(6) = {18, 0, 0, cl};
Point(7) = {20, 0, 0, cli};
Point(8) = {20, 13, 0, cli};
Point(9) = {20, 15, 0, cli};
Point(10) = {18, 15, 0, cli};
Point(11) = {2, 15, 0, cli};
Point(12) = {0, 15, 0, cli};
Point(13) = {0, 13, 0, cli};
Point(14) = {2, 13, 0, cl};
Point(15) = {18, 13, 0, cl};


Line(1) = {2, 3};
Line(2) = {3, 4};
Line(3) = {4, 5};
Line(4) = {5, 6};
Line(5) = {6, 15};
Line(6) = {15, 14};
Line(7) = {14, 2};

Line(8) = {1, 2};
Line(9) = {14, 13};
Line(10) = {13, 1};

Line(11) = {6, 7};
Line(12) = {7, 8};
Line(13) = {8, 15};

Line(14) = {15, 10};
Line(15) = {10, 11};
Line(16) = {11, 14};

Line(17) = {11, 12};
Line(18) = {12, 13};

Line(19) = {8, 9};
Line(20) = {9, 10};

Line Loop(1) = {1, 2, 3, 4, 5, 6, 7};
Plane Surface(1) = {1};
Recombine Surface {1};

Line Loop(2) = {8, -7, 9, 10};
Plane Surface(2) = {2};
//Transfinite Surface {2} Alternated;
Recombine Surface {2};

Line Loop(3) = {11, 12, 13, -5};
Plane Surface(3) = {3};
//Transfinite Surface {3} Alternated;
Recombine Surface {3};

Line Loop(4) = {14, 15, 16, -6};
Plane Surface(4) = {4};
//Transfinite Surface {4} Alternated;
Recombine Surface {4};

Line Loop(5) = {17, 18, -9, -16};
Plane Surface(5) = {5};
//Transfinite Surface {5} Alternated;
Recombine Surface {5};

Line Loop(6) = {19, 20, -14, -13};
Plane Surface(6) = {6};
//Transfinite Surface {6} Alternated;
Recombine Surface {6};

Physical Surface(1000) = {1, 2, 3, 4, 5, 6};

Physical Line(1) = {20, 15, 17};
Physical Line(2) = {18, 10, 12, 19};

Mesh.SecondOrderIncomplete = 1;



