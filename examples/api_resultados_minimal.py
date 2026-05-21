"""Ejemplo mínimo de consumo de la API pública de resultados (ADR 0002).

Muestra el patrón de uso esperado para consumidores externos como FenixBAR:
construir un Domain programáticamente, invocar ``solidum.run`` y leer el
``SolveResult`` inmutable resultante.

Se resuelve un marco 2D de una sola barra empotrada en un extremo con carga
transversal en el otro (cantilever clásico). Las fuerzas internas retornadas
siguen la convención §5 de Reglas.md: ``V`` horario positivo, ``M`` sagging
positivo.
"""

import os
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import solidum
from solidum import Domain, Frame2DEuler
from solidum.materials.elastic import Elastic1D
import numpy as np


# 1. Construir el dominio
domain = Domain()
n1 = domain.add_node(1, [0.0, 0.0])
n2 = domain.add_node(2, [2.0, 0.0])

material = Elastic1D(E=2.1e11)  # acero, Pa
elem = Frame2DEuler(1, [n1, n2], material, A=1e-3, I=8.33e-6)
domain.add_element(elem)

# 2. Condiciones de contorno: nodo 1 empotrado
n1.fix_dof('ux', 0.0)
n1.fix_dof('uy', 0.0)
n1.fix_dof('rz', 0.0)
domain.generate_equation_numbers()

# 3. Vector de cargas: 10 kN hacia abajo en el nodo 2
F = np.zeros(domain.total_dofs)
F[n2.dofs['uy']] = -10_000.0

# 4. Ejecutar — punto de entrada oficial para consumidores
result = solidum.run(domain, F_applied=F)

# 5. Leer resultados
print(f"Convergió: {result.converged}  (pasos: {result.num_steps})")
print(f"Desplazamiento vertical en nodo 2: {result.U[n2.dofs['uy']]:.6e} m")

print("\nReacciones en apoyos (DOFs restringidos):")
for node_id, dofs in result.reactions_by_node.items():
    for dof_name, value in dofs.items():
        print(f"  nodo {node_id}, {dof_name}: {value:+.4e}")

print("\nFuerzas internas por elemento (ejes locales, convención §5):")
for elem_id, ef in result.element_forces.items():
    print(f"  elem {elem_id}  kind={ef.kind}")
    for comp_name, arr in ef.components.items():
        print(f"    {comp_name}: i={arr[0]:+.4e}   j={arr[1]:+.4e}")

# 6. El resultado también queda cacheado en el dominio
assert domain.last_result is result
