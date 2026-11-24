"""Compilación común de los cuatro manuales con lualatex.

Los auxiliares (``.aux``, ``.log``, ``.out``, ``.toc``) no se escriben junto al
``.tex``: la carpeta del proyecto se sincroniza en la nube y esos archivos
cambian en cada compilación. Van a una carpeta de auxiliares, una subcarpeta
por manual, y sólo el PDF se copia junto al ``.tex``.

Carpeta de auxiliares: la variable de entorno ``SOLIDUM_LATEX_AUX`` si existe;
si no, ``<temporal del sistema>/solidum_manuals``. No se borra al terminar: el
registro (``<manual>.log``) se revisa después de compilar, y el builder imprime
su ruta.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from pathlib import Path


def aux_dir(tex_path: Path) -> Path:
    """Carpeta de auxiliares de un manual (se crea si no existe)."""
    base = os.environ.get("SOLIDUM_LATEX_AUX")
    raiz = Path(base) if base else Path(tempfile.gettempdir()) / "solidum_manuals"
    carpeta = raiz / tex_path.stem
    carpeta.mkdir(parents=True, exist_ok=True)
    return carpeta


def compile_pdf(tex_path: Path) -> bool:
    """Compila ``tex_path`` con dos pasadas de lualatex (para el índice).

    Se compila desde la carpeta del ``.tex`` (las rutas relativas de figuras y
    ``\\input`` no cambian) con ``-output-directory`` en la carpeta de
    auxiliares; al terminar bien, el PDF se copia junto al ``.tex``.
    """
    lualatex = shutil.which("lualatex")
    if lualatex is None:
        print("[!] lualatex no encontrado en PATH. Genera el .tex pero no compila.")
        return False
    aux = aux_dir(tex_path)
    for i in range(2):
        result = subprocess.run(
            [lualatex, "-interaction=nonstopmode", "-halt-on-error",
             f"-output-directory={aux}", tex_path.name],
            cwd=tex_path.parent,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
        )
        if result.returncode != 0:
            print(f"[!] lualatex (pasada {i+1}) fallo:")
            tail = result.stdout[-3000:]
            try:
                print(tail)
            except UnicodeEncodeError:
                # La consola (cp1252 en Windows) puede no aceptar algunos
                # glifos Unicode del registro; se degradan a ASCII.
                print(tail.encode("ascii", errors="replace").decode("ascii"))
            print(f"    Registro completo: {aux / (tex_path.stem + '.log')}")
            return False
    shutil.copy2(aux / (tex_path.stem + ".pdf"), tex_path.with_suffix(".pdf"))
    print(f"  -> registro de compilación: {aux / (tex_path.stem + '.log')}")
    return True
