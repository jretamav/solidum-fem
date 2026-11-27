// Dos bloques de 2 x 1 apilados que se tocan en y = 1 sin compartir la curva
// (nodos coincidentes y distintos) para tests/test_interface_generation.py
// (ADR 0043, pair).
DefineConstant[ h = 0.5 ];

Point(1) = {0, 0, 0, h};
Point(2) = {2, 0, 0, h};
Point(3) = {2, 1, 0, h};
Point(4) = {0, 1, 0, h};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Point(5) = {0, 1, 0, h};
Point(6) = {2, 1, 0, h};
Point(7) = {2, 2, 0, h};
Point(8) = {0, 2, 0, h};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};
Curve Loop(2) = {5, 6, 7, 8};
Plane Surface(2) = {2};
Transfinite Curve{3, 5} = 5;

Physical Surface("Bloques") = {1, 2};
Physical Curve("TechoAbajo") = {3};
Physical Curve("PisoArriba") = {5};
Physical Curve("Base") = {1};
Physical Curve("Tapa") = {7};
Physical Point("Origen") = {1};
