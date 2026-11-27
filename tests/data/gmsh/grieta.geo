// Placa de 2 x 1 con una grieta vertical de x = 1, y = 0 a y = 0.5 (curva
// incrustada en la superficie) para tests/test_interface_generation.py (ADR 0043).
DefineConstant[ h = 0.25 ];

Point(1) = {0, 0, 0, h};
Point(2) = {2, 0, 0, h};
Point(3) = {2, 1, 0, h};
Point(4) = {0, 1, 0, h};
Point(5) = {1, 0, 0, h};
Point(6) = {1, 0.5, 0, h};
Line(1) = {1, 5};
Line(2) = {5, 2};
Line(3) = {2, 3};
Line(4) = {3, 4};
Line(5) = {4, 1};
Line(6) = {5, 6};
Curve Loop(1) = {1, 2, 3, 4, 5};
Plane Surface(1) = {1};
Line{6} In Surface{1};

Physical Surface("Placa") = {1};
Physical Curve("Grieta") = {6};
Physical Curve("Abajo") = {1, 2};
