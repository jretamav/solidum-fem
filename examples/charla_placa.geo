// Placa cuadrada con agujero central para charla de IA.
// Variante de placa.geo con malla de cuadrilateros (Recombine + Frontal-Delaunay for Quads).

lc = 0.08;     // Tamano de los elementos de la malla
L = 2.0;       // Lado de la placa
r = 0.4;       // Radio del agujero central

// Contorno exterior
Point(1) = {-L/2, -L/2, 0, lc};
Point(2) = { L/2, -L/2, 0, lc};
Point(3) = { L/2,  L/2, 0, lc};
Point(4) = {-L/2,  L/2, 0, lc};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Contorno interior (agujero)
Point(5) = {0, 0, 0, lc};
Point(6) = {r, 0, 0, lc};
Point(7) = {0, r, 0, lc};
Point(8) = {-r, 0, 0, lc};
Point(9) = {0, -r, 0, lc};

Circle(5) = {6, 5, 7};
Circle(6) = {7, 5, 8};
Circle(7) = {8, 5, 9};
Circle(8) = {9, 5, 6};

Curve Loop(1) = {1, 2, 3, 4};
Curve Loop(2) = {5, 6, 7, 8};

Plane Surface(1) = {1, 2};

// Algoritmo Frontal-Delaunay for Quads + recombinacion -> Quad4.
Mesh.Algorithm = 8;          // Frontal-Delaunay for Quads
Mesh.RecombineAll = 1;       // Recombina triangulos en cuadrilateros
Mesh.RecombinationAlgorithm = 1; // Blossom (mejor calidad)
Recombine Surface{1};

// Grupos fisicos para BC en el YAML
Physical Curve("Apoyo_Izquierdo", 10) = {4};
Physical Curve("Carga_Derecha", 11) = {2};
Physical Surface("Dominio_Placa", 12) = {1};
