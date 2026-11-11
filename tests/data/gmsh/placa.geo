// Placa de 2 x 1 para tests/test_mesh_from_geometry.py (ADR 0035).
// El tamaño de elemento h se sobrescribe desde el YAML (mesh.parameters).
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

Physical Surface("Placa") = {1};
Physical Curve("Izquierda") = {4};
Physical Curve("Derecha") = {2};
