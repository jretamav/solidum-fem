"""Ejecuta un modelo YAML y escribe sus resultados: ``solidum modelo.yaml``
(comando que crea la instalación) o, lo mismo, ``python -m solidum modelo.yaml``.

Lee y valida el archivo (todos los errores juntos), resuelve el análisis
declarado y escribe junto al YAML lo que pide su bloque ``output`` (ADR 0036):
VTK según ``results.frequency``, la historia ``history`` en CSV y la salida de
texto. Equivale a ``solidum.run_yaml(ruta, write_output=True)``.

Si el análisis se detiene (no converge, mecanismo, sistema mal planteado),
escribe los resultados del último paso convergido, explica la causa y
termina con el código 3.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path


def main(argv=None) -> int:
    from solidum.logging import utf8_output_when_redirected
    utf8_output_when_redirected()
    ap = argparse.ArgumentParser(prog="solidum",
                                 description="Resuelve un modelo YAML de Solidum FEM y escribe sus resultados.")
    ap.add_argument("modelo", help="archivo YAML del modelo")
    args = ap.parse_args(argv)
    ruta = Path(args.modelo)
    if not ruta.is_file():
        print(f"[!] No se encontró el archivo {ruta}", file=sys.stderr)
        return 1
    import solidum
    from solidum.math.solvers.diagnostics import (IllPosedSystemError, MechanismError,
                                                  SolverDivergedError)
    from solidum.utils.yaml_parser import YamlValidationError
    try:
        resultado = solidum.run_yaml(ruta, write_output=True)
    except YamlValidationError as exc:
        print(f"[!] El modelo tiene errores:\n{exc}", file=sys.stderr)
        return 2
    except (SolverDivergedError, MechanismError, IllPosedSystemError) as exc:
        # La nota la añade la salida sólo si escribió los resultados.
        notas = "".join(f"\n{n}" for n in getattr(exc, "__notes__", ()))
        print(f"[!] El análisis se detuvo:\n{exc}{notas}", file=sys.stderr)
        return 3
    print(f"Análisis completado ({type(resultado).__name__}); resultados en {ruta.resolve().parent}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
