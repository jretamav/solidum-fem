# solidum_fem/solidum/materials/concrete_1d.py
"""``ManderConcrete1D`` — concreto uniaxial de daño-plasticidad con datos de
Mander. Ver ``docs/specs/ManderConcrete1D.md``.

Fenomenológico, con estructura de daño-plasticidad 1D: la envolvente de
compresión (Popovics), la deformación plástica (punto focal) y la
degradación de la tracción por compresión previa son de Mander, Priestley y
Park (1988); la descarga y la recarga son rectas con la rigidez dañada.
"""
from __future__ import annotations

import math

import numpy as np

from solidum.core.material import Material
from solidum.registry import MaterialRegistry

# Puntos de Gauss-Legendre para el área de la rama posterior al pico (la
# curva de Popovics es suave: 64 puntos la integran a precisión de máquina).
_N_GAUSS_POST = 64

_KEYS = ("eps_cmax", "eps_tmax", "h")


@MaterialRegistry.register
class ManderConcrete1D(Material):
    """Concreto uniaxial: envolvente de Mander, descarga de daño-plasticidad.

    Convención: ``σ > 0`` y ``ε > 0`` en tracción; los parámetros de
    compresión se dan en magnitud positiva.

    Parameters
    ----------
    fc : float
        Resistencia de pico en compresión (``f'cc`` en el núcleo confinado,
        ``f'c`` en el recubrimiento), > 0.
    eps_c : float
        Deformación de pico, > 0.
    E : float
        Módulo inicial, > ``fc/eps_c`` (Mander: ``5000·√f'c`` MPa).
    eps_cu : float
        Fin de la curva de Popovics, ≥ ``eps_c``.
    eps_sp : float, optional
        Deformación en que el esfuerzo llega a cero por una recta desde
        ``eps_cu`` (desprendimiento del recubrimiento), ≥ ``eps_cu``. Por
        omisión ``eps_cu``: el esfuerzo cae a cero en ``eps_cu``.
    ft : float, optional
        Resistencia a tracción, ≥ 0 (0: sin tracción).
    eps_tu : float, optional
        Deformación de tracción en que el esfuerzo llega a cero, ≥ ``ft/E``.
    G_f : float, optional
        Energía de fractura en tracción (N/m). Con ``ft > 0`` se da
        ``eps_tu`` o ``G_f``, no ambos. Con ``G_f``, ``ε_tu = 2·G_f/(ft·h)``
        con ``h`` la longitud característica del punto de integración.
    G_fc : float, optional
        Energía de fractura en compresión (N/m): la rama posterior al pico
        se estira para que su área valga ``G_fc/h``.
    density : float, optional
        Densidad (ADR 0008).
    confinement : dict, optional
        Argumentos de
        :func:`~solidum.materials.concrete_confinement.mander_confined_concrete`;
        calcula ``fc``, ``eps_c``, ``eps_cu`` y ``E`` (este último, si no se
        da) del núcleo confinado. Excluye ``fc``, ``eps_c`` y ``eps_cu``.

    Notes
    -----
    Variables internas: ``eps_cmax`` (máxima compresión alcanzada, ``e_m``),
    ``eps_tmax`` (máxima tracción alcanzada, medida desde el origen
    desplazado ``−e_p``) y ``h`` (longitud característica del punto, fija;
    0 si el elemento no la proporcionó). La deformación plástica ``e_p``,
    el daño ``d`` y el factor de tracción ``φ`` se deducen de ``e_m``.
    """
    STRAIN_DIM = 1
    PRIMARY_STATE_VAR = "eps_cmax"
    STATE_SCHEMA = {"eps_cmax": (), "eps_tmax": (), "h": ()}

    def __init__(self, fc: float | None = None, eps_c: float | None = None,
                 E: float | None = None, eps_cu: float | None = None,
                 eps_sp: float | None = None, ft: float = 0.0,
                 eps_tu: float | None = None, G_f: float | None = None,
                 G_fc: float | None = None, density: float | None = None,
                 confinement: dict | None = None):
        where = "ManderConcrete1D"
        if confinement is not None:
            given = [n for n, v in (("fc", fc), ("eps_c", eps_c), ("eps_cu", eps_cu))
                     if v is not None]
            if given:
                raise ValueError(
                    f"{where}: con 'confinement' no se declaran {given}: la función de "
                    f"confinamiento los calcula.")
            if not isinstance(confinement, dict):
                raise ValueError(f"{where}: 'confinement' debe ser un diccionario.")
            from solidum.materials.concrete_confinement import mander_confined_concrete
            core = mander_confined_concrete(**confinement)
            fc, eps_c, eps_cu = core["fc"], core["eps_c"], core["eps_cu"]
            if E is None:
                E = core["E"]
        missing = [n for n, v in (("fc", fc), ("eps_c", eps_c), ("E", E), ("eps_cu", eps_cu))
                   if v is None]
        if missing:
            raise ValueError(f"{where}: faltan {missing} (o 'confinement').")

        self.fc = _positive(fc, "fc")
        self.eps_c = _positive(eps_c, "eps_c")
        self.E = _positive(E, "E")
        if not self.E > self.fc / self.eps_c:
            raise ValueError(
                f"{where}: E = {self.E:.6g} debe ser mayor que la secante al pico "
                f"fc/eps_c = {self.fc / self.eps_c:.6g} (la curva de Popovics exige r > 1).")
        self.eps_cu = _positive(eps_cu, "eps_cu")
        if self.eps_cu < self.eps_c:
            raise ValueError(f"{where}: eps_cu = {self.eps_cu} < eps_c = {self.eps_c}.")
        self.eps_sp = self.eps_cu if eps_sp is None else _positive(eps_sp, "eps_sp")
        if self.eps_sp < self.eps_cu:
            raise ValueError(f"{where}: eps_sp = {self.eps_sp} < eps_cu = {self.eps_cu}.")

        self.ft = _nonnegative(ft, "ft")
        if self.ft == 0.0:
            if eps_tu is not None or G_f is not None:
                raise ValueError(
                    f"{where}: con ft = 0 (sin tracción) no se declaran eps_tu ni G_f.")
            self.eps_tu = None
            self.G_f = None
        else:
            if (eps_tu is None) == (G_f is None):
                raise ValueError(
                    f"{where}: con ft > 0 declare la rama descendente de tracción con "
                    f"eps_tu (deformación) o con G_f (energía de fractura), uno de los dos.")
            self.G_f = None if G_f is None else _positive(G_f, "G_f")
            self.eps_tu = None if eps_tu is None else _positive(eps_tu, "eps_tu")
            if self.eps_tu is not None and self.eps_tu < self.ft / self.E:
                raise ValueError(
                    f"{where}: eps_tu = {self.eps_tu} < ft/E = {self.ft / self.E:.6g}.")
        self.G_fc = None if G_fc is None else _positive(G_fc, "G_fc")
        if density is not None and not (math.isfinite(density) and density >= 0.0):
            raise ValueError(f"{where}: density = {density} debe ser ≥ 0.")
        self.density = density

        self._r = self.E / (self.E - self.fc / self.eps_c)
        self._eps_t = self.ft / self.E
        self._s_cu = self._popovics(self.eps_cu)[0]
        self._A_post = self._post_peak_area()
        if self.G_fc is not None and not self._A_post > 0.0:
            raise ValueError(
                f"{where}: G_fc no aplica: la envolvente no tiene rama posterior al "
                f"pico (eps_cu = eps_sp = eps_c).")

    # ------------------------------------------------------------------
    # Longitud característica (spec §5)
    # ------------------------------------------------------------------

    @property
    def uses_characteristic_length(self) -> bool:
        return self.G_f is not None or self.G_fc is not None

    def initial_state(self) -> dict:
        return {"eps_cmax": 0.0, "eps_tmax": 0.0, "h": 0.0}

    def initial_point_state(self, characteristic_length: float) -> dict:
        h = float(characteristic_length)
        if not (math.isfinite(h) and h > 0.0):
            raise ValueError(
                f"ManderConcrete1D: longitud característica h = {h} no válida (> 0).")
        if self.G_f is not None:
            h_max = 2.0 * self.G_f * self.E / self.ft ** 2
            if h >= h_max:
                raise ValueError(
                    f"ManderConcrete1D: la longitud característica h = {h:.4g} supera el "
                    f"límite 2·G_f·E/ft² = {h_max:.4g} (G_f = {self.G_f:.4g}, "
                    f"ft = {self.ft:.4g}, E = {self.E:.4g}): la rama de ablandamiento en "
                    f"tracción tendría retroceso en el propio material "
                    f"(ε_tu = 2·G_f/(ft·h) < ft/E). Use elementos más cortos.")
        return {"eps_cmax": 0.0, "eps_tmax": 0.0, "h": h}

    def _missing_h(self, energy: str) -> ValueError:
        return ValueError(
            f"ManderConcrete1D con {energy} necesita la longitud característica del "
            f"punto de integración, y el elemento que lo usa no la proporcionó. La "
            f"proporcionan las barras (Truss2D, Truss3D y sus versiones "
            f"corrotacionales) y los marcos con FiberSection; con otro elemento, dé la "
            f"rama descendente en deformaciones (eps_tu; eps_cu y eps_sp).")

    def _stretch(self, h: float) -> float:
        """Factor ``λ`` de la rama posterior al pico (1 sin ``G_fc``)."""
        if self.G_fc is None:
            return 1.0
        if not h > 0.0:
            raise self._missing_h("G_fc")
        return self.G_fc / (h * self._A_post)

    def tension_softening_end(self, h: float = 0.0) -> float:
        """``ε_tu`` efectiva: la dada, o ``2·G_f/(ft·h)`` con energía."""
        if self.G_f is None:
            return self.eps_tu if self.eps_tu is not None else self._eps_t
        if not h > 0.0:
            raise self._missing_h("G_f")
        return 2.0 * self.G_f / (self.ft * h)

    # ------------------------------------------------------------------
    # Envolventes (spec §2, §4, §5)
    # ------------------------------------------------------------------

    def _popovics(self, e: float):
        """``(s, ds/de)`` de la curva de Popovics (Mander 1988, ecs. 3–8)."""
        x = e / self.eps_c
        r = self._r
        if x <= 0.0:
            return 0.0, self.E
        xr = x ** r
        den = r - 1.0 + xr
        s = self.fc * x * r / den
        ds = self.fc * r * (r - 1.0) * (1.0 - xr) / (self.eps_c * den * den)
        return s, ds

    def _s_reference(self, e: float):
        """Envolvente de compresión sin regularizar: Popovics hasta
        ``eps_cu``, recta hasta cero en ``eps_sp`` y cero después."""
        if e <= self.eps_cu:
            return self._popovics(e)
        if e <= self.eps_sp and self.eps_sp > self.eps_cu:
            slope = -self._s_cu / (self.eps_sp - self.eps_cu)
            return self._s_cu + slope * (e - self.eps_cu), slope
        return 0.0, 0.0

    def compression_envelope(self, e: float, h: float = 0.0):
        """``(s, ds/de)`` de la envolvente de compresión, con la rama posterior
        al pico estirada por ``λ`` si hay ``G_fc`` (spec §5)."""
        if e <= self.eps_c:
            return self._popovics(e)
        lam = self._stretch(h)
        s, ds = self._s_reference(self.eps_c + (e - self.eps_c) / lam)
        return s, ds / lam

    def _post_peak_area(self) -> float:
        """``A_post = ∫ s_ref de`` desde el pico hasta ``eps_sp``."""
        area = 0.0
        if self.eps_cu > self.eps_c:
            xi, w = np.polynomial.legendre.leggauss(_N_GAUSS_POST)
            a, b = self.eps_c, self.eps_cu
            pts = 0.5 * (b - a) * xi + 0.5 * (b + a)
            area += 0.5 * (b - a) * sum(wi * self._popovics(p)[0] for p, wi in zip(pts, w))
        if self.eps_sp > self.eps_cu:
            area += 0.5 * self._s_cu * (self.eps_sp - self.eps_cu)
        return float(area)

    def plastic_strain(self, e_m: float, s_m: float) -> float:
        """Deformación plástica de compresión ``e_p`` del punto focal de Mander
        (ecs. 31–34) para la máxima compresión ``e_m`` y su esfuerzo ``s_m``."""
        if e_m <= 0.0:
            return 0.0
        a = max(self.eps_c / (self.eps_c + e_m), 0.09 * e_m / self.eps_c)
        e_a = a * math.sqrt(e_m * self.eps_c)
        return e_m - (e_m + e_a) * s_m / (s_m + self.E * e_a)

    def _tension_undamaged(self, t: float, h: float):
        """``(σ⁰_t, dσ⁰_t/dt)``: lineal hasta ft y recta hasta cero en ε_tu."""
        if t <= self._eps_t:
            return self.E * t, self.E
        eps_tu = self.tension_softening_end(h)
        if t < eps_tu:
            slope = -self.ft / (eps_tu - self._eps_t)
            return self.ft + slope * (t - self._eps_t), slope
        return 0.0, 0.0

    # ------------------------------------------------------------------
    # Contrato de Material
    # ------------------------------------------------------------------

    def admissibility_scale(self, state_vars=None) -> float:
        return self.fc

    def compute_state(self, strain, state_vars=None):
        try:
            eps = float(strain)
        except TypeError:                     # arreglo de un elemento
            eps = float(np.asarray(strain, dtype=float).reshape(-1)[0])
        if state_vars is None:
            em0, tm0, h = 0.0, 0.0, 0.0
        else:
            em0 = float(state_vars["eps_cmax"])
            tm0 = float(state_vars["eps_tmax"])
            h = float(state_vars.get("h", 0.0))

        e = -eps
        loading = e > em0
        em = e if loading else em0
        if em > 0.0:
            sm, dsm = self.compression_envelope(em, h)
            ep = self.plastic_strain(em, sm)
        else:
            sm, dsm, ep = 0.0, self.E, 0.0

        t = eps + ep
        tm = tm0
        if t <= 0.0:
            # Compresión, grieta cerrada: la rigidez (1 − d)·E no depende del
            # daño de tracción (recuperación total, spec §4).
            if loading:
                sigma, tangent = -sm, dsm
            elif em <= 0.0:
                sigma, tangent = self.E * eps, self.E
            elif sm <= 0.0:
                sigma, tangent = 0.0, 0.0              # aplastada
            else:
                E_r = sm / (em - ep)
                sigma, tangent = -E_r * (e - ep), E_r
        else:
            phi = max(0.0, 1.0 - ep / self.eps_c)
            if self.ft <= 0.0 or phi <= 0.0:
                sigma, tangent = 0.0, 0.0
                tm = max(tm0, t)
            elif t > tm0:
                tm = t
                s0, ds0 = self._tension_undamaged(t, h)
                sigma, tangent = phi * s0, phi * ds0
            elif tm0 <= self._eps_t:
                sigma, tangent = phi * self.E * t, phi * self.E
            else:
                secant = phi * self._tension_undamaged(tm0, h)[0] / tm0
                sigma, tangent = secant * t, secant
        return sigma, tangent, {"eps_cmax": em, "eps_tmax": tm, "h": h}

    def __repr__(self) -> str:
        return (f"ManderConcrete1D(fc={self.fc:.4g}, eps_c={self.eps_c:.4g}, E={self.E:.4g}, "
                f"eps_cu={self.eps_cu:.4g}, eps_sp={self.eps_sp:.4g}, ft={self.ft:.4g})")


def _positive(x, name: str) -> float:
    x = float(x)
    if not (math.isfinite(x) and x > 0.0):
        raise ValueError(f"ManderConcrete1D: {name} debe ser un número finito > 0 (recibido {x}).")
    return x


def _nonnegative(x, name: str) -> float:
    x = float(x)
    if not (math.isfinite(x) and x >= 0.0):
        raise ValueError(f"ManderConcrete1D: {name} debe ser un número finito ≥ 0 (recibido {x}).")
    return x
