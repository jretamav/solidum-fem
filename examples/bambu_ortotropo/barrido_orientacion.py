"""Barrido de la orientación de fibra sobre la tira de bambú.

Ejecuta `tira_45_grados.yaml` variando `theta` y reporta cómo cambia la
respuesta. Reproduce, desde un modelo FEM completo, la curva de constantes
aparentes que `tests/validation/test_off_axis_orthotropic.py` valida contra
las formas cerradas de Jones (1999) §2.8.

Uso
---
    python examples/bambu_ortotropo/barrido_orientacion.py

Sin argumentos. No modifica el YAML en disco: lo lee, sustituye `theta` en
memoria y escribe un temporal que borra al terminar.
"""
import os
import re
import sys
import tempfile

import numpy as np

sys.path.append(os.path.abspath(
    os.path.join(os.path.dirname(__file__), '..', '..')))

import solidum  # noqa: E402


AQUI = os.path.dirname(os.path.abspath(__file__))
MODELO = os.path.join(AQUI, 'tira_45_grados.yaml')

ANGULOS = (0.0, 15.0, 30.0, 45.0, 60.0, 75.0, 90.0)

# Del YAML. Se repiten aquí sólo para el cálculo del módulo aparente.
SIGMA_APLICADO = 10.0e6   # Pa
LONGITUD = 1.0            # m


def _correr(theta_deg: float):
    """Ejecuta el modelo con la orientación pedida. Devuelve (ux, uy) del
    nodo 25 — la esquina superior derecha, el punto más alejado del apoyo."""
    with open(MODELO, encoding='utf-8') as f:
        texto = f.read()
    texto = re.sub(r'theta: [0-9.]+', f'theta: {theta_deg}', texto)

    tmp = tempfile.NamedTemporaryFile(
        mode='w', suffix='.yaml', delete=False, encoding='utf-8')
    try:
        tmp.write(texto)
        tmp.close()
        resultado = solidum.run_yaml(tmp.name)
    finally:
        os.unlink(tmp.name)

    U = resultado.U
    return U[-2], U[-1]   # ux, uy del último nodo


def main() -> None:
    print()
    print("Tira de bambu traccionada — efecto de la orientacion de fibra")
    print("E1 = 15 GPa,  E2 = 0.8 GPa,  G12 = 0.7 GPa,  nu12 = 0.35")
    print("Traccion aplicada: 10 MPa")
    print()
    print("  theta    ux [mm]    uy [mm]    uy/ux     E_x [GPa]   E_x/E1")
    print("  " + "-" * 62)

    ux_ref = None
    for theta in ANGULOS:
        ux, uy = _correr(theta)
        E_x = SIGMA_APLICADO / (ux / LONGITUD)
        if ux_ref is None:
            ux_ref = ux
        E_ratio = ux_ref / ux
        print(f"  {theta:5.1f}  {ux * 1e3:9.4f}  {uy * 1e3:9.4f}  "
              f"{uy / ux:8.4f}  {E_x / 1e9:10.3f}  {E_ratio:8.4f}")

    print()
    print("Lectura de los resultados")
    print("-------------------------")
    print("* uy/ux en theta=0 y 90 grados es Poisson puro: -0.35 = -nu12 y")
    print("  -0.0187 = -nu21. La reciprocidad nu12/E1 = nu21/E2 obliga a que")
    print("  el segundo sea ~19 veces menor.")
    print()
    print("* Fuera de los ejes principales uy/ux se dispara: -2.25 a 15")
    print("  grados, muy por encima de cualquier Poisson admisible. Eso NO")
    print("  es contraccion transversal — es el acoplamiento traccion-")
    print("  cortante, que ningun material isotropo presenta: la tira se")
    print("  distorsiona ademas de alargarse.")
    print()
    print("* El maximo de distorsion relativa NO cae en 45 grados sino cerca")
    print("  de 15, donde la tira aun es rigida (E_x = 6.7 GPa) y el")
    print("  acoplamiento ya es fuerte. A 45 grados el material se ha")
    print("  ablandado tanto que ux crece y el cociente baja, aunque uy")
    print("  alcance ahi su maximo absoluto.")
    print()
    print("* E_x cae por debajo del 10% de E1 a 45 grados, pese a que")
    print("  E1/E2 = 18.75. El minimo NO se obtiene interpolando entre los")
    print("  dos modulos principales: el termino cruzado (1/G12 - 2*nu12/E1)")
    print("  domina la zona intermedia y con G12 pequeno hunde la curva.")
    print()


if __name__ == '__main__':
    main()
