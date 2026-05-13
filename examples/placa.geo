// Parámetros de control
lc = 0.1;       // Tamaño de los elementos de la malla
L = 2.0;        // Lado de la placa
r = 0.4;        // Radio del agujero central

// --- Puntos del contorno exterior (Placa) ---
Point(1) = {-L/2, -L/2, 0, lc};
Point(2) = { L/2, -L/2, 0, lc};
Point(3) = { L/2,  L/2, 0, lc};
Point(4) = {-L/2,  L/2, 0, lc};

// --- Líneas del contorno exterior ---
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// --- Puntos del contorno interior (Agujero) ---
Point(5) = {0, 0, 0, lc};   // Centro
Point(6) = {r, 0, 0, lc};
Point(7) = {0, r, 0, lc};
Point(8) = {-r, 0, 0, lc};
Point(9) = {0, -r, 0, lc};

// --- Arcos para formar el círculo ---
Circle(5) = {6, 5, 7};
Circle(6) = {7, 5, 8};
Circle(7) = {8, 5, 9};
Circle(8) = {9, 5, 6};

// --- Definición de Bucles (Loops) ---
Curve Loop(1) = {1, 2, 3, 4};    // Perímetro exterior
Curve Loop(2) = {5, 6, 7, 8};    // Perímetro interior

// --- Definición de la Superficie (Resta de Loops) ---
Plane Surface(1) = {1, 2};

// --- Definición de Grupos Físicos (Para exportar) ---
Physical Curve("Apoyo_Izquierdo", 10) = {4};
Physical Curve("Carga_Derecha", 11) = {2};
Physical Surface("Dominio_Placa", 12) = {1};