Point(14) = {0.0, 0.0, 0, 1e+22};
Point(15) = {1.0, 0.0, 0, 1e+22};
Point(16) = {1.0, 1.0, 0, 1e+22};
Point(17) = {0.0, 1.0, 0, 1e+22};
Line(19) = {14, 15};
Line(20) = {15, 16};
Line(21) = {16, 17};
Line(22) = {17, 14};
Line Loop(24) = {21, 22, 19, 20};
Plane Surface(24) = {24};
//stuff to make it a structured grid:
Transfinite Line "*" = 5;
Transfinite Surface "*";
Recombine Surface "*";
