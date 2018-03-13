//h    =  1e-4;//doption(gmsh.hsize) ;//m
//hPCB =  13e-3;//doption(Geo.hPCB) ;//m
//ePCB =  2e-4;//doption(Geo.ePCB) ;//m
//hIC =  2e-3;//doption(Geo.hIC) ;//m
//eIC =  2e-4;//doption(Geo.eIC) ;//m
//h1   =  2e-3;//doption(Geo.h1) ;//m
//h2   =  1e-2;//doption(Geo.h2) ;//hPCB-hIC-h1;//m
//eAIR =  4e-4;//doption(Geo.eAIR) ;//m



Point(1)    = {0, 0, 0, h};
Point(2)    = {ePCB, 0, 0, h};
Point(3)    = {ePCB, h1, 0, h};
Point(4)    = {ePCB, h1+hIC, 0, h};
Point(5)    = {ePCB, h2, 0, h};
Point(6)    = {ePCB, h2+hIC, 0, h};
Point(7)    = {ePCB, hPCB, 0, h/2};
Point(8)    = {0, hPCB, 0, h};

// other corner of processors
Point(9)    = {ePCB+eIC, h1, 0, h};
Point(10)   = {ePCB+eIC, h1+hIC, 0, h};
Point(11)   = {ePCB+eIC, h2, 0, h};
Point(12)   = {ePCB+eIC, h2+hIC, 0, h};

// other corner of the areation
Point(13)   = {ePCB+eAIR, 0, 0, h};
Point(14)   = {ePCB+eAIR, hPCB, 0, h};
Point(15)   = {ePCB+eIC, 0, 0, h};

//______line __________________________________________
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 1};//PBC's cycle

Line(9) = {3, 9};
Line(10) = {9, 10};
Line(11) = {10, 4};//IC1's cycle

Line(12) = {5, 11};
Line(13) = {11, 12};
Line(14) = {12, 6};//IC2's cycle

Line(15) = {2, 15};
Line(16) = {15, 13};
Line(17) = {13, 14};
Line(18) = {14, 7};

//______surface________________________________________
Line Loop(20) = {1, 2, 3, 4, 5, 6, 7, 8};
Plane Surface(100) = {20};//motherBoard
Line Loop(21) = {9, 10, 11, -3};
Plane Surface(101) = {21};//processor IC1
Line Loop(22) = {12, 13, 14, -5};
Plane Surface(102) = {22};//processor IC2
Line Loop(23) = {15, 16, 17, 18, -6, -14, -13, -12, -4, -11, -10, -9, -2};
Plane Surface(103) = {23};//the conduct of areation

//______physical-line__________________________________
Physical Line("in1") = {15};//condition of dirichlet: T=T0, u_a=f_entre
Physical Line("in2") = {16};//condition of dirichlet: T=T0, u_a=f_entre
Physical Line("out") = {18};//condition of Robin
//Physical Line("wall") = {17};// paroi du conduit d'areation
Physical Line("borderfluid") = {2,9,10,11,4,12,13,14,6,17};// paroi du conduit d'areation interieur

//______physical-surface_______________________________
Physical Surface("PCB") = {100};//motherboard PCB
Physical Surface("IC1") = {101};//processor IC1
Physical Surface("IC2") = {102};//processor IC2
Physical Surface("AIR") = {103};//conduct of areation;
