"""Infraestructura compartida del estudio de idealización isótropa en bambú.

Fase 0 del plan. Define los rangos de parámetros, las tres idealizaciones
isótropas candidatas y las métricas de error, de modo que las cuatro fases del
estudio usen exactamente las mismas definiciones.

Las decisiones metodológicas están fijadas en `../DISENO_DEL_ESTUDIO.md`; este
módulo es su traducción a código y no debe introducir criterios nuevos.
"""
from __future__ import annotations

import os
import subprocess
import sys
from dataclasses import dataclass

import numpy as np

sys.path.append(os.path.abspath(
    os.path.join(os.path.dirname(__file__), '..', '..', '..')))

from solidum.materials.elastic_2d import Elastic2D          # noqa: E402
from solidum.materials.orthotropic_2d import Orthotropic2D  # noqa: E402


# ----------------------------------------------------------------------
# Rangos de parámetros — §5 del diseño
# ----------------------------------------------------------------------

#: Módulo longitudinal de referencia. El estudio es invariante de escala en E_L
#: (todas las métricas son adimensionales), así que su valor sólo fija unidades.
E_L_REF = 15.0e9

#: Relación de anisotropía. Medidas en bambú: 5.12-8.02 según posición radial.
RANGO_ANISOTROPIA = (3.0, 12.0)

#: Cortante normalizado G_LT/E_L. Poco documentado — es el hueco que este
#: estudio pretende cuantificar en importancia (hipótesis H4).
RANGO_G_NORMALIZADO = (0.02, 0.10)

#: Coeficiente de Poisson mayor. Moso: 0.180-0.334.
RANGO_NU = (0.18, 0.34)

#: Valores centrales, usados cuando un barrido fija los demás parámetros.
CENTRO = {
    'E_L': E_L_REF,
    'anisotropia': 6.0,       # dentro del 5.12-8.02 medido
    'g_normalizado': 0.05,
    'nu_LT': 0.30,
}

#: Umbrales de admisibilidad — §7.3. Declarados ANTES de ver resultados.
UMBRAL_ADMISIBLE = 0.05
UMBRAL_RESERVAS = 0.15


def clasificar(error: float) -> str:
    """Traduce un error relativo a la escala de admisibilidad del §7.3."""
    if error < UMBRAL_ADMISIBLE:
        return "admisible"
    if error < UMBRAL_RESERVAS:
        return "con reservas"
    return "inadmisible"


# ----------------------------------------------------------------------
# Construcción de materiales
# ----------------------------------------------------------------------

@dataclass(frozen=True)
class ParametrosBambu:
    """Un punto del espacio de parámetros del estudio.

    Se parametriza por la RELACIÓN de anisotropía y el cortante NORMALIZADO en
    vez de por E_2 y G_12 absolutos, porque las métricas son adimensionales y
    así el barrido no depende de la escala de E_L.
    """
    E_L: float = E_L_REF
    anisotropia: float = 6.0      # E_L / E_T
    g_normalizado: float = 0.05   # G_LT / E_L
    nu_LT: float = 0.30

    @property
    def E_T(self) -> float:
        return self.E_L / self.anisotropia

    @property
    def G_LT(self) -> float:
        return self.g_normalizado * self.E_L

    @property
    def nu_TL(self) -> float:
        """Recíproco, por la relación de reciprocidad nu_LT/E_L = nu_TL/E_T."""
        return self.nu_LT * self.E_T / self.E_L

    def ortotropo(self, theta_deg: float = 0.0) -> Orthotropic2D:
        return Orthotropic2D(
            E1=self.E_L, E2=self.E_T, G12=self.G_LT,
            nu12=self.nu_LT, theta=theta_deg,
        )

    def es_admisible(self) -> bool:
        """El punto respeta la restricción |nu12| < sqrt(E1/E2)."""
        return abs(self.nu_LT) < np.sqrt(self.anisotropia)


# ----------------------------------------------------------------------
# Las tres idealizaciones isótropas — §6 del diseño
# ----------------------------------------------------------------------

def isotropo_I_L(p: ParametrosBambu) -> Elastic2D:
    """I-L: E = E_L, nu = nu_LT. La práctica habitual.

    Se mide el módulo longitudinal —el dato más disponible en la literatura de
    bambú— y se usa como si el material fuera isótropo.
    """
    return Elastic2D(E=p.E_L, nu=p.nu_LT, hypothesis='plane_stress')


def isotropo_I_T(p: ParametrosBambu) -> Elastic2D:
    """I-T: E = E_T, nu = nu_TL. La idealización conservadora."""
    return Elastic2D(E=p.E_T, nu=p.nu_TL, hypothesis='plane_stress')


# ----------------------------------------------------------------------
# Métricas de error — §7 del diseño
# ----------------------------------------------------------------------

def constantes_aparentes(material) -> tuple[float, float, float]:
    """(E_x, nu_xy, eta_xy_x) bajo tracción uniaxial en x.

    Se impone el estado de ESFUERZO y no el de deformación porque así se
    definen las constantes aparentes: un ensayo de tracción deja libres las
    demás componentes. Invirtiendo C se obtiene la flexibilidad S y de ahí:

        E_x      = 1 / S[0,0]
        nu_xy    = -S[1,0] / S[0,0]
        eta_xy,x =  S[2,0] / S[0,0]
    """
    S = np.linalg.inv(material.C)
    return 1.0 / S[0, 0], -S[1, 0] / S[0, 0], S[2, 0] / S[0, 0]


def m1_error_modulo(mat_iso, mat_orto) -> float:
    """M1 — error relativo en el módulo aparente E_x."""
    E_iso, _, _ = constantes_aparentes(mat_iso)
    E_orto, _, _ = constantes_aparentes(mat_orto)
    return abs(E_iso - E_orto) / E_orto


def m2_error_constitutiva(mat_iso, mat_orto) -> float:
    """M2 — error en norma de Frobenius de la matriz constitutiva.

    Captura el desajuste completo, no sólo el de una componente: un isótropo
    puede acertar E_x y seguir errando en el resto de la matriz.
    """
    dif = mat_iso.C - mat_orto.C
    return np.linalg.norm(dif, 'fro') / np.linalg.norm(mat_orto.C, 'fro')


def m3_acoplamiento_omitido(mat_orto) -> float:
    """M3 — acoplamiento tracción-cortante que la idealización descarta.

    NO se expresa como error relativo: el isótropo da eta = 0 para todo theta,
    de modo que el error relativo sería infinito. Se reporta el valor absoluto
    del acoplamiento omitido, que es lo físicamente interpretable.
    """
    _, _, eta = constantes_aparentes(mat_orto)
    return abs(eta)


def m4_error_desplazamiento(u_iso: np.ndarray, u_orto: np.ndarray) -> float:
    """M4 — error relativo en norma L2 del campo de desplazamientos."""
    return float(np.linalg.norm(u_iso - u_orto) / np.linalg.norm(u_orto))


def m5_error_esfuerzos(sig_iso: np.ndarray, sig_orto: np.ndarray,
                       pesos: np.ndarray | None = None) -> float:
    """M5 — error en esfuerzos en norma L2 ponderada por cuadratura.

    La ponderación por el peso de Gauss evita que la métrica dependa del número
    de elementos de la malla.
    """
    if pesos is None:
        pesos = np.ones(len(sig_orto))
    num = np.sum(pesos * np.sum((sig_iso - sig_orto) ** 2, axis=1))
    den = np.sum(pesos * np.sum(sig_orto ** 2, axis=1))
    return float(np.sqrt(num / den))


# ----------------------------------------------------------------------
# Trazabilidad
# ----------------------------------------------------------------------

def commit_actual() -> str:
    """Hash del commit con el que se generaron los resultados.

    Se serializa junto a los datos para que una figura sea siempre atribuible a
    una versión concreta del código.
    """
    try:
        r = subprocess.run(
            ['git', 'rev-parse', '--short', 'HEAD'],
            cwd=os.path.dirname(os.path.abspath(__file__)),
            capture_output=True, text=True, timeout=5,
        )
        if r.returncode == 0:
            return r.stdout.strip()
    except (FileNotFoundError, subprocess.TimeoutExpired):
        pass
    return "desconocido"


DIR_RESULTADOS = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                              'resultados')
DIR_FIGURAS = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           '..', 'figs')
