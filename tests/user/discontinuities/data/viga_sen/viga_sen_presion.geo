// SIN placa inferior: carga 10P/11 como presión en la línea 3 (x = 0.23..0.25).
// Viga SEN de Schlangen (cortante en cuatro puntos), geometría del modelo FEAP
// de la tesis (Iviga4p2.dat), en metros. Entalla arriba (y = 0.08..0.10),
// cargas hacia arriba en la cara inferior (x = 0.02 y placa en x = 0.23..0.25),
// apoyos en la cara superior (placa en x = 0.19..0.21 y nodo en x = 0.42).
// h: tamaño de elemento en la zona central (x = 0.17..0.28), donde corre la grieta.
DefineConstant[ h = 0.010 ];
hout = Min(4*h, 0.020);

// Contorno de la viga (sentido antihorario)
Point(1)  = {0.0,    0.0,  0};
Point(2)  = {0.02,   0.0,  0};   // carga P/11
Point(3)  = {0.23,   0.0,  0};   // placa inferior
Point(4)  = {0.25,   0.0,  0};
Point(5)  = {0.44,   0.0,  0};
Point(6)  = {0.44,   0.10, 0};
Point(7)  = {0.42,   0.10, 0};   // apoyo derecho (ux, uy)
Point(8)  = {0.2225, 0.10, 0};   // boca de la entalla, lado derecho (nodo 170 de FEAP)
Point(9)  = {0.2225, 0.08, 0};
Point(10) = {0.2175, 0.08, 0};
Point(11) = {0.2175, 0.10, 0};   // boca de la entalla, lado izquierdo (nodo 162 de FEAP)
Point(12) = {0.21,   0.10, 0};   // placa superior
Point(13) = {0.19,   0.10, 0};
Point(14) = {0.0,    0.10, 0};
// Placa superior (apoyo uy en su centro)
Point(18) = {0.21, 0.11, 0};
Point(19) = {0.20, 0.11, 0};
Point(20) = {0.19, 0.11, 0};

Line(1)  = {1, 2};
Line(2)  = {2, 3};
Line(3)  = {3, 4};
Line(4)  = {4, 5};
Line(5)  = {5, 6};
Line(6)  = {6, 7};
Line(7)  = {7, 8};
Line(8)  = {8, 9};
Line(9)  = {9, 10};
Line(10) = {10, 11};
Line(11) = {11, 12};
Line(12) = {12, 13};
Line(13) = {13, 14};
Line(14) = {14, 1};
Curve Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14};
Plane Surface(1) = {1};


Line(19) = {12, 18};
Line(20) = {18, 19};
Line(21) = {19, 20};
Line(22) = {20, 13};
Curve Loop(3) = {-12, 19, 20, 21, 22};
Plane Surface(3) = {3};

// Tamaño: h en la zona central, hasta hout hacia los extremos
Field[1] = Box;
Field[1].VIn = h;
Field[1].VOut = hout;
Field[1].XMin = 0.17;
Field[1].XMax = 0.28;
Field[1].YMin = -0.02;
Field[1].YMax = 0.12;
Field[1].Thickness = 0.06;
Background Field = 1;
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;

Physical Surface("Concreto") = {1};
Physical Surface("Placas") = {3};
Physical Point("CargaIzq") = {2};
Physical Curve("CaraCarga") = {3};
Physical Point("ApoyoIzq") = {19};
Physical Point("ApoyoDer") = {7};
Physical Point("BocaIzq") = {11};
Physical Point("BocaDer") = {8};
