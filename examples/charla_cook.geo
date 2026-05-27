// Cook's membrane — geometria canonica (Simo & Armero 1992).
// Trapecio: (0,0)-(48,44)-(48,60)-(0,44). Borde izquierdo empotrado,
// borde derecho cargado en cortante (desplazamiento vertical prescrito).
// Mallado transfinita estructurada -> Quad4 puro, sin triangulos residuales.

// Numero de divisiones en cada direccion (N+1 nodos por lado)
N_horiz = 24;   // divisiones a lo largo de los bordes superior e inferior
N_vert  = 20;   // divisiones a lo largo de los bordes izquierdo y derecho

// Vertices
Point(1) = { 0,  0, 0};   // esquina inferior izquierda (empotrada)
Point(2) = {48, 44, 0};   // esquina inferior derecha
Point(3) = {48, 60, 0};   // esquina superior derecha
Point(4) = { 0, 44, 0};   // esquina superior izquierda (empotrada)

// Lados
Line(1) = {1, 2};   // borde inferior
Line(2) = {2, 3};   // borde derecho (cargado)
Line(3) = {3, 4};   // borde superior
Line(4) = {4, 1};   // borde izquierdo (empotrado)

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// Mallado transfinita estructurada en quads
Transfinite Curve{1, 3} = N_horiz + 1;
Transfinite Curve{2, 4} = N_vert + 1;
Transfinite Surface{1};
Recombine Surface{1};

// Grupos fisicos para BC en el YAML
Physical Curve("Apoyo_Izquierdo", 10) = {4};
Physical Curve("Carga_Derecha",   11) = {2};
Physical Surface("Dominio_Cook",  12) = {1};
