//__________dimension parameter_________________________________________________
h=.1;//mm

// dimension of the motherBoard PCB
ePCB    = 2;//mm
hPCB    = 100;//mm
lPCB    = 20;//mm

// dimension of the processor IC
eIC    = 2;//mm
hIC    = 10;//mm
lIC    = 5;//mm

// position of the processor 
h1      = 30;//mm
l1      = 7;//mm
h2      = 60;//mm
l2      = 7;//mm

// dimension of Air
eA      = 50;//mm

//__________creation of the Point_______________________________________________
Point(1)    = {0, 0, 0, h};
Point(2)    = {0, hPCB, 0, h};
Point(3)    = {0, hPCB, lPCB, h};
Point(4)    = {0, 0, lPCB, h};
Point(5)    = {ePCB, 0, 0, h};
Point(6)    = {ePCB, hPCB, 0, h};
Point(7)    = {ePCB, hPCB, lPCB, h};
Point(8)    = {ePCB, 0, lPCB, h};
Point(10)    = {ePCB, h1, l1, h};
Point(11)    = {ePCB, h1+hIC, l1, h};
Point(12)    = {ePCB, h1+hIC, l1+lIC, h};
Point(13)    = {ePCB, h1, l1+lIC, h};
Point(14)    = {ePCB+eIC, h1, l1, h};
Point(15)    = {ePCB+eIC, h1+hIC, l1, h};
Point(16)    = {ePCB+eIC, h1+hIC, l1+lIC, h};
Point(17)    = {ePCB+eIC, h1, l1+lIC, h};
Point(18)    = {ePCB, h2, l2, h};
Point(19)    = {ePCB, h2+hIC, l2, h};
Point(20)    = {ePCB, h2+hIC, l2+lIC, h};
Point(21)    = {ePCB, h2, l2+lIC, h};
Point(22)    = {ePCB+eIC, h2, l2, h};
Point(23)    = {ePCB+eIC, h2+hIC, l2, h};
Point(24)    = {ePCB+eIC, h2+hIC, l2+lIC, h};
Point(25)    = {ePCB+eIC, h2, l1+lIC, h};
Point(26)    = {ePCB+eA, 0, 0, h};
Point(27)    = {ePCB+eA, hPCB, 0, h};
Point(28)    = {ePCB+eA, hPCB, lPCB, h};
Point(29)    = {ePCB+eA, 0, lPCB, h};



//__________creation of the line________________________________________________
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};
Line(9) = {1, 5};
Line(10) = {2, 6};
Line(11) = {3, 7};
Line(12) = {4, 8};
Line(13) = {10, 11};
Line(14) = {11, 12};
Line(15) = {12, 13};
Line(16) = {13, 10};
Line(17) = {17, 14};
Line(18) = {14, 15};
Line(19) = {15, 16};
Line(20) = {16, 17};
Line(21) = {10, 14};
Line(22) = {11, 15};
Line(23) = {12, 16};
Line(24) = {13, 17};
Line(25) = {18, 19};
Line(26) = {19, 20};
Line(27) = {20, 21};
Line(28) = {21, 18};
Line(29) = {22, 23};
Line(30) = {23, 24};
Line(31) = {24, 25};
Line(32) = {25, 22};
Line(33) = {18, 22};
Line(34) = {19, 23};
Line(35) = {20, 24};
Line(36) = {21, 25};
Line(37) = {26, 27};
Line(38) = {27, 28};
Line(39) = {28, 29};
Line(40) = {29, 26};
Line(41) = {5, 26};
Line(42) = {6, 27};
Line(43) = {7, 28};
Line(44) = {8, 29};

//__________surface_____________________________________________________________
Line Loop(45) = {1, 2, 3, 4};
Plane Surface(46) = {45};
Line Loop(47) = {5, 6, 7, 8};
Line Loop(48) = {13, 14, 15, 16};
Line Loop(49) = {25, 26, 27, 28};
Plane Surface(50) = {47, 48, 49};
Line Loop(51) = {1, 10, -5, -9};
Plane Surface(52) = {51};
Line Loop(53) = {2, 11, -6, -10};
Plane Surface(54) = {53};
Line Loop(55) = {3, 12, -7, -11};
Plane Surface(56) = {55};
Line Loop(57) = {4, 9, -8, -12};
Plane Surface(58) = {57};
Plane Surface(59) = {48};
Line Loop(60) = {18, 19, 20, 17};
Plane Surface(61) = {60};
Line Loop(62) = {13, 22, -18, -21};
Plane Surface(63) = {62};
Line Loop(64) = {14, 23, -19, -22};
Plane Surface(65) = {64};
Line Loop(66) = {15, 24, -20, -23};
Plane Surface(67) = {66};
Line Loop(68) = {16, 21, -17, -24};
Plane Surface(69) = {68};
Plane Surface(70) = {49};
Line Loop(71) = {29, 30, 31, 32};
Plane Surface(72) = {71};
Line Loop(73) = {25, 34, -29, -33};
Plane Surface(74) = {73};
Line Loop(75) = {26, 35, -30, -34};
Plane Surface(76) = {75};
Line Loop(77) = {27, 36, -31, -35};
Plane Surface(78) = {77};
Line Loop(79) = {28, 33, -32, -36};
Plane Surface(80) = {79};
Line Loop(81) = {7, 44, -39, -43};
Plane Surface(82) = {81};
Line Loop(83) = {6, 43, -38, -42};
Plane Surface(84) = {83};
Line Loop(85) = {5, 42, -37, -41};
Plane Surface(86) = {85};
Line Loop(87) = {8, 41, -40, -44};
Plane Surface(88) = {87};
Line Loop(89) = {39, 40, 37, 38};
Plane Surface(90) = {89};

//__________volume______________________________________________________________
Surface Loop(91) = {50, 70, 59, 56, 46, 52, 54, 58};
Volume(92) = {91};
Surface Loop(93) = {78, 80, 74, 76, 72, 70};
Volume(94) = {93};
Surface Loop(95) = {67, 69, 63, 65, 61, 59};
Volume(96) = {95};
Surface Loop(97) = {86, 84, 82, 88, 90, 50, 61, 63, 65, 67, 69, 72, 74, 76, 78, 80};
Volume(98) = {97};

//__________physical surface____________________________________________________
Physical Surface("neumann") = {54, 52, 58, 56, 46, 90};
Physical Surface("entre") = {88};
Physical Surface("sortie") = {84, 86, 82};

//__________physical volume____________________________________________________
Physical Volume("PCB") = {92};
Physical Volume("IC1") = {96};
Physical Volume("IC2") = {94};
Physical Volume("Air") = {98};
