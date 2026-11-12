# solidum_fem/examples/ejecutar_yaml.py
"""Ejecuta un modelo YAML y escribe sus resultados.

Uso: python examples/ejecutar_yaml.py ruta/al/modelo.yaml

Equivale a ``python -m solidum ruta/al/modelo.yaml``: valida el archivo,
resuelve el análisis declarado (el mismo despacho que ``solidum.run_yaml``,
con etapas, cargas de cuerpo y térmicas) y escribe junto al YAML lo que pide
su bloque ``output`` (ADR 0036): el modelo sin deformar, los VTK según
``results.frequency`` con su ``.pvd``, la historia ``history`` en CSV, la
salida de texto estilo FEAP y ``<nombre>_final.vtu``.
"""
import os
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from solidum.__main__ import main  # noqa: E402

if __name__ == "__main__":
    print("--- SOLIDUM FEM: MOTOR DE EJECUCIÓN YAML ---")
    sys.exit(main())
