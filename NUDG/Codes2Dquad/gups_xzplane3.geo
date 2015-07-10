Point(14) = {0.0024415, -1.5e-4, 0, 1e+22};
Point(15) = {4.99757, -1.5e-4, 0, 1e+22};
Point(16) = {4.99757, -0.1495, 0, 1e+22};
Point(17) = {0.0024415, -0.1495, 0, 1e+22};
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
