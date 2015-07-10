Point(14) = {0.25, -0.005, 0, 1e+22};
Point(15) = {4.5, -0.005, 0, 1e+22};
Point(16) = {4.5, -0.145, 0, 1e+22};
Point(17) = {0.25, -0.145, 0, 1e+22};
Line(19) = {14, 15};
Line(20) = {15, 16};
Line(21) = {16, 17};
Line(22) = {17, 14};
Line Loop(24) = {21, 22, 19, 20};
Plane Surface(24) = {24};
//stuff for structured grid
Transfinite Line "*"=34;
Transfinite Surface "*";
Recombine Surface "*";
