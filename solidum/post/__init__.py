"""Post-proceso de los resultados de un modelo YAML (ADR 0039).

:func:`make_plots` dibuja las figuras de ``output.plots`` a partir de los
archivos que escribe Solidum; ``python -m solidum.post modelo.yaml`` las
regenera sin volver a calcular.
"""
from solidum.post.plots import PLOT_TYPES, ResultFiles, make_plots, validate_plots_block

__all__ = ["PLOT_TYPES", "ResultFiles", "make_plots", "validate_plots_block"]
