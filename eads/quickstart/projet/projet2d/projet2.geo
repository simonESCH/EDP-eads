//parameter of the mesh

h    = .5;//mm

// dimension of the motherboard( PCB)
hPCB = 7;//mm
ePCB = 2;//mm

// dimension of the processor ( IC)
hIC = 2;//mm
eIC = 1;//mm

// placement of the processor
h1   = 1;//mm
h2   = hPCB-hIC-h1;//mm

//dimension of the areation( A)
eA   = 4;//mm





//______point__________________________________________
// corner of the motherboard
Point(1)    = {0, 0, 0, h};
Point(2)    = {ePCB, 0, 0, h};
Point(3)    = {ePCB, h1, 0, h};
Point(4)    = {ePCB, h1+hIC, 0, h};
Point(5)    = {ePCB, h2, 0, h};
Point(6)    = {ePCB, h2+hIC, 0, h};
Point(7)    = {ePCB, hPCB, 0, h};
Point(8)    = {0, hPCB, 0, h};

// other corner of processors
Point(9)    = {ePCB+eIC, h1, 0, h};
Point(10)   = {ePCB+eIC, h1+hIC, 0, h};
Point(11)   = {ePCB+eIC, h2, 0, h};
Point(12)   = {ePCB+eIC, h2+hIC, 0, h};

// other corner of the areation
Point(13)   = {ePCB+eA, 0, 0, h};
Point(14)   = {ePCB+eA, hPCB, 0, h};
Point(15)   = {ePCB+eIC, 0, 0, h};
Point(16)   = {ePCB+eIC, hPCB, 0, h};

//______line __________________________________________
Line(1) = {1, 8};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {15, 9};
Line(8) = {9, 10};
Line(9) = {10, 11};
Line(10) = {11, 12};
Line(11) = {12, 16};
Line(12) = {14, 13};
Line(13) = {13, 13};
Line(14) = {1, 2};
Line(15) = {2, 15};
Line(16) = {15, 13};
Line(17) = {3, 9};
Line(18) = {4, 10};
Line(19) = {5, 11};
Line(20) = {6, 12};
Line(21) = {8, 7};
Line(22) = {7, 16};
Line(23) = {16, 14};

//______surface________________________________________
Line Loop(24) = {1, 21, -6, -5, -4, -3, -2, -14};
Plane Surface(25) = {24};//motherBoard
Line Loop(26) = {7, 8, 9, 10, 11, 23, 12, -16};
Plane Surface(27) = {26};//air's passage
Line Loop(28) = {17, 8, -18, -3};
Plane Surface(29) = {28};//processor IC1
Line Loop(30) = {19, 10, -20, -5};
Plane Surface(31) = {30};//processor IC2
Line Loop(32) = {15, 16, -12, -23, -22, -6, 20, -10, -19, -4, 18, -8, -17, -2};
Plane Surface(33) = {32};//the conduct of areation

//______physical-line__________________________________
Physical Line("neumann") = {12, 14, 1, 21};//condition of neumann
Physical Line("in") = {15, 16};//condition of dirichlet: T=T0, u_a=f_entre
Physical Line("out") = {22, 23};//condition ... :

//______physical-surface______________________________
Physical Surface("PCB") = {25};//motherboard PCB
Physical Surface("IC1") = {29};//processor IC1
Physical Surface("IC2") = {31};//processor IC2
Physical Surface("POI") = {27};//passage of the air 
Physical Surface("AIR") = {33};//conduct of areation
