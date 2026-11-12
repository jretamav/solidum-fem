"""Regenera las figuras de un modelo sin volver a calcular:
``python -m solidum.post modelo.yaml``.

Lee ``output.plots`` del YAML y dibuja con los resultados que escribió
``python -m solidum modelo.yaml`` (ADR 0039).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="python -m solidum.post",
                                 description="Dibuja las figuras de output.plots con los resultados "
                                             "ya escritos de un modelo YAML.")
    ap.add_argument("modelo", help="archivo YAML del modelo (ya resuelto)")
    args = ap.parse_args(argv)
    ruta = Path(args.modelo)
    if not ruta.is_file():
        print(f"[!] No se encontró el archivo {ruta}", file=sys.stderr)
        return 1
    import yaml
    from solidum.post import ResultFiles, make_plots, validate_plots_block
    with open(ruta, encoding="utf-8") as f:
        data = yaml.safe_load(f) or {}
    plots = (data.get("output") or {}).get("plots")
    if not plots:
        print(f"[!] {ruta.name} no declara output.plots.", file=sys.stderr)
        return 2
    errores = validate_plots_block(plots)
    if errores:
        print("[!] output.plots tiene errores:\n  " + "\n  ".join(errores), file=sys.stderr)
        return 2
    try:
        escritas = make_plots(ResultFiles.of_yaml(ruta), plots)
    except (FileNotFoundError, ValueError) as exc:
        print(f"[!] {exc}", file=sys.stderr)
        return 3
    for p in escritas:
        print(p)
    return 0


if __name__ == "__main__":
    sys.exit(main())
