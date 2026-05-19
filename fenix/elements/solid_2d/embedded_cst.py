"""``CST_Embedded2D`` — triángulo CST con discontinuidad interior embebida (KOS).

Fase 2 del ADR 0010. Cinemática enriquecida fiel al Cap. 2, 5, 7 de Retama
(2010); longitud efectiva ``l_d = (A_e/h)·cos(θ−α)`` del Cap. 6. Ver
``docs/specs/CST_Embedded2D.md`` para la formulación completa.
"""
from __future__ import annotations

from typing import List

import numpy as np

from fenix.core.cohesive_material import CohesiveMaterial
from fenix.core.discontinuity_state import DiscontinuityState
from fenix.core.element import Element
from fenix.core.material import Material
from fenix.core.node import Node
from fenix.elements.solid_2d._shared import (
    _compute_integrands,
    _compute_kinematics_tri3,
)
from fenix.registry import ElementRegistry


# Tolerancias del Newton local del salto (algoritmo de condensación).
# Centralizadas en ``fenix.constants`` (auditoría H-1.9) — re-exportadas
# como alias de módulo para preservar imports históricos.
from fenix.constants import (
    EMBEDDED_LOCAL_JUMP_MAX_ITER as _LOCAL_JUMP_MAX_ITER,
    EMBEDDED_LOCAL_JUMP_RTOL as _LOCAL_JUMP_RTOL,
)

# Bulks aceptados en fase 1 (ver punto abierto A de la spec, decisión cerrada
# 2026-05-18: la discrete approach presupone bulk elástico).
_ACCEPTED_BULKS = ('Elastic2D',)


@ElementRegistry.register
class CST_Embedded2D(Element):
    """
    Triángulo CST 2D con discontinuidad interior embebida (formulación KOS).

    En estado intacto se comporta exactamente como :class:`Tri3` (mismo `B`,
    misma rigidez, mismas fuerzas internas). Cuando un paso de carga cruza el
    umbral de tracción :math:`\\sigma_{t0}` del material cohesivo, el
    elemento activa una discontinuidad interna (irreversible) y a partir de
    ese paso enriquece su cinemática con un salto :math:`[[u]]` localizado
    en :math:`\\Gamma_d`. Los DOFs enriquecidos son elementales: se
    condensan estáticamente dentro de :meth:`compute_element_state` y nunca
    llegan al ensamblador global (ADR 0010 §3).

    Parameters
    ----------
    element_id : int
    nodes : list of Node
        Tres nodos del CST padre, en orden antihorario.
    material : Material
        Material elástico del bulk. **Fase 1: solo Elastic2D** (la discrete
        approach de Retama 2010 supone bulk elástico; mezclar plasticidad
        del bulk con cohesivo es out-of-scope).
    cohesive_material : CohesiveMaterial
        Material cohesivo *traction-jump* que gobierna :math:`\\Gamma_d` una
        vez activado. Debe declarar ``JUMP_DIM = 2``.
    thickness : float, default 1.0
        Espesor para estado plano.

    Notes
    -----
    La activación se evalúa en :meth:`prepare_step`, **al inicio de cada
    paso** del solver no lineal (ADR 0010 §5; evita chattering dentro del
    Newton). Una vez activada, la discontinuidad es irreversible.
    """

    DOF_NAMES = ['ux', 'uy']
    STRAIN_DIM = 3
    N_INTEGRATION_POINTS = 1

    def __init__(self, element_id: int, nodes: List[Node],
                 material: Material, cohesive_material: CohesiveMaterial,
                 thickness: float = 1.0):
        # Validación de bulk (fase 1)
        bulk_name = type(material).__name__
        if bulk_name not in _ACCEPTED_BULKS:
            raise ValueError(
                f"CST_Embedded2D(id={element_id}): bulk '{bulk_name}' no admitido "
                f"en fase 1. Bulks aceptados: {_ACCEPTED_BULKS}. La discrete "
                f"approach presupone bulk elástico (la disipación se concentra "
                f"en Γ_d). Otros bulks → out-of-scope en esta entrega."
            )
        # Validación de cohesivo
        if not isinstance(cohesive_material, CohesiveMaterial):
            raise TypeError(
                f"CST_Embedded2D(id={element_id}): cohesive_material debe ser "
                f"instancia de CohesiveMaterial; recibido "
                f"{type(cohesive_material).__name__}."
            )
        if getattr(cohesive_material, 'JUMP_DIM', None) != 2:
            raise ValueError(
                f"CST_Embedded2D(id={element_id}): cohesive_material "
                f"{type(cohesive_material).__name__} declara JUMP_DIM="
                f"{getattr(cohesive_material, 'JUMP_DIM', None)}; este elemento "
                f"2D requiere JUMP_DIM=2."
            )

        self.thickness = thickness
        self.cohesive_material = cohesive_material
        self.discontinuity_state: DiscontinuityState | None = None
        super().__init__(element_id, nodes, material)

    # ------------------------------------------------------------------
    # Hook de activación (ADR 0010 §5)
    # ------------------------------------------------------------------

    def prepare_step(self, U_committed: np.ndarray) -> None:
        """Chequea activación Rankine al inicio del paso con el estado convergido.

        Solo activa elementos intactos; una vez activado, nunca revierte.
        Usa el estado committed (no trial) para evitar chattering por el
        predictor lineal dentro del Newton (ADR 0010 §5).
        """
        if self.discontinuity_state is not None:
            return  # ya agrietado, irreversible

        u_e = self.get_local_displacements(U_committed)
        coords = self.get_coordinate_matrix(ndim=2)
        B, _ = _compute_kinematics_tri3(coords)
        strain = B @ u_e
        sigma, _, _ = self.material.compute_state(strain, self.state.vars[0])

        sigma_I, n_dir = _max_principal_stress_2d(sigma)
        sigma_t0 = getattr(self.cohesive_material, 'sigma_t0', None)
        if sigma_t0 is None:
            return  # cohesivo sin umbral declarado: no activable por Rankine

        if sigma_I > sigma_t0:
            self._activate(n_dir, coords)

    def _activate(self, n_dir: np.ndarray, coords: np.ndarray) -> None:
        """Materializa la DiscontinuityState al cruzar el umbral.

        Identifica el nodo solitario, orienta `n` hacia él, computa
        ``l_d = (A_e/h)·cos(θ−α)`` y crea el estado de la discontinuidad.
        """
        centroid = coords.mean(axis=0)

        # Productos escalares (x_i − centroid) · n para cada nodo
        rel = coords - centroid
        signs = rel @ n_dir
        n_positive = int(np.sum(signs > 0.0))

        # Convención: n apunta hacia el lado donde hay UN nodo (solitario).
        # Si la dirección elegida deja 2 nodos en lado positivo, invertirla.
        if n_positive == 2:
            n_dir = -n_dir
            signs = -signs
            n_positive = 1
        elif n_positive != 1:
            # Casos degenerados: línea perpendicular pasa exactamente por un nodo,
            # o todos los nodos están a un solo lado (geometría imposible si el
            # centroide es el promedio). Defensivo.
            raise RuntimeError(
                f"CST_Embedded2D(id={self.id}): activación degenerada — "
                f"{n_positive} nodos en lado positivo de Γ_d (esperado 1)."
            )

        solitary = int(np.argmax(signs))
        tangent = np.array([-n_dir[1], n_dir[0]])
        l_d = _compute_ld(coords, tangent, solitary)

        # Estado committed = jump cero, cohesivo virgen (se llenará en commit
        # tras la primera convergencia post-activación).
        self.discontinuity_state = DiscontinuityState(
            normal=np.asarray(n_dir, dtype=np.float64),
            tangent=np.asarray(tangent, dtype=np.float64),
            centroid=np.asarray(centroid, dtype=np.float64),
            solitary_node=solitary,
            l_d=l_d,
        )

    # ------------------------------------------------------------------
    # API principal: rigidez + fuerzas internas (estado intacto o agrietado)
    # ------------------------------------------------------------------

    def compute_element_state(self, u_e: np.ndarray):
        coords = self.get_coordinate_matrix(ndim=2)
        B, detJ = _compute_kinematics_tri3(coords)

        if self.discontinuity_state is None:
            return self._element_state_intact(u_e, B, detJ)
        return self._element_state_cracked(u_e, B, detJ)

    def _element_state_intact(self, u_e, B, detJ):
        """Bit-exact con Tri3 (mismo B, mismo material, mismo cuadratura)."""
        strain = B @ u_e
        sigma, C_tan, new_state = self.material.compute_state(strain, self.state.vars[0])
        self.state.vars_trial[0] = new_state
        self.state.stresses_trial[0] = sigma
        K_e, F_int_e = _compute_integrands(B, C_tan, sigma, detJ, 0.5, self.thickness)
        return K_e, F_int_e

    def _element_state_cracked(self, u_e, B_std, detJ):
        """Cinemática KOS + Newton local del salto + condensación estática."""
        ds = self.discontinuity_state
        G = np.column_stack([ds.normal, ds.tangent])              # 2×2
        i_star = ds.solitary_node
        B_phi = B_std[:, 2 * i_star:2 * i_star + 2]                 # 3×2 (B^φ)

        A_e = 0.5 * detJ                                            # área del CST
        vol = A_e * self.thickness                                  # factor de bulk
        C_e = self.material.C                                       # tensor elástico

        # Newton local sobre [[u]] hasta R^{[[u]]} = 0. CohesiveDamageIsotropic
        # tiene tangente piecewise constante (cargando vs descargando), así
        # que típicamente converge en 1-3 iteraciones. σ/t/R/K se evalúan
        # siempre con el ``jump`` actual y *antes* del update; salimos del
        # bucle sin actualizar (por convergencia o agotando max_iter) para
        # que las cantidades persistidas sean consistentes con ``jump_trial``.
        jump = np.copy(ds.jump_trial)
        coh_committed = ds.cohesive_state_committed or None

        for it in range(_LOCAL_JUMP_MAX_ITER):
            # Bulk: σ = C_e · (B_std·d − B^φ·G·[[u]])
            strain_bulk = B_std @ u_e - B_phi @ (G @ jump)
            sigma, C_tan_bulk, mat_state = self.material.compute_state(
                strain_bulk, self.state.vars[0]
            )

            # Cohesivo: t, T_tan en frame local
            t_local, T_coh, coh_state = self.cohesive_material.compute_traction(
                jump, coh_committed
            )

            # Residual y rigidez del salto:
            #   R^{[[u]]} = −G^T · (B^φ)^T · σ · vol + l_d · t · thickness
            #   K_jj    = G^T · (B^φ)^T · C_tan_bulk · B^φ · G · vol + l_d · T_coh · thickness
            # El factor ``thickness`` en los términos cohesivos cierra la
            # consistencia dimensional: el bulk lleva ``vol = A_e · thickness``
            # mientras que el cohesivo en 2D entrega ``t · l_d`` con unidades
            # N/m por unidad de espesor; multiplicarlo por ``thickness``
            # devuelve fuerza 3D N y rigidez N/m. Con ``thickness = 1.0`` el
            # bug es invisible (caso de todos los tests previos).
            R_jump = -G.T @ (B_phi.T @ sigma) * vol + ds.l_d * t_local * self.thickness
            K_jj = G.T @ (B_phi.T @ C_tan_bulk @ B_phi) @ G * vol + ds.l_d * T_coh * self.thickness

            scale = max(np.linalg.norm(t_local) * ds.l_d, 1.0)
            converged = np.linalg.norm(R_jump) < _LOCAL_JUMP_RTOL * scale
            if converged or it == _LOCAL_JUMP_MAX_ITER - 1:
                break

            jump = jump - np.linalg.solve(K_jj, R_jump)

        # Persistir trial
        ds.jump_trial = jump
        ds.cohesive_state_trial = coh_state
        self.state.vars_trial[0] = mat_state
        self.state.stresses_trial[0] = sigma

        # Matrices acopladas (vol ya aplicado):
        # K_dd = B_std^T · C_tan_bulk · B_std · vol         (6×6)
        # K_du = −B_std^T · C_tan_bulk · B^φ · G · vol      (6×2)
        K_dd = B_std.T @ C_tan_bulk @ B_std * vol
        K_du = -B_std.T @ (C_tan_bulk @ B_phi) @ G * vol

        # Condensación: K_cond = K_dd − K_du · K_jj^{-1} · K_du^T
        X = np.linalg.solve(K_jj, K_du.T)                   # (2×6)
        K_cond = K_dd - K_du @ X

        # F_int^cond = F_d − K_du · K_jj^{-1} · R_jump,
        # donde F_d = B_std^T · σ · vol.
        F_d = B_std.T @ sigma * vol
        F_int_cond = F_d - K_du @ np.linalg.solve(K_jj, R_jump)

        return K_cond, F_int_cond

    # ------------------------------------------------------------------
    # Commit del estado (incluye recuperación del salto convergido)
    # ------------------------------------------------------------------

    def commit_state(self) -> None:
        super().commit_state()
        if self.discontinuity_state is not None:
            self.discontinuity_state.commit()

    # ------------------------------------------------------------------
    # Post-procesamiento (compute_gauss_state extendido)
    # ------------------------------------------------------------------

    def compute_internal_forces(self, U_global: np.ndarray) -> dict:
        gs = self.compute_gauss_state(U_global)
        return {'stress': gs['stress'][0], 'strain': gs['strain'][0]}

    def compute_gauss_state(self, U_global: np.ndarray) -> dict:
        """σ, ε en el punto central + información de la grieta si está activada.

        En estado intacto coincide con :meth:`Tri3.compute_gauss_state`. En
        estado agrietado añade la clave ``'discontinuity'`` con la geometría
        de ``Γ_d`` (normal, tangente, centroide, l_d, nodo solitario) y el
        estado actual del salto (jump, tracción y damage en frame local).
        """
        u_e = self.get_local_displacements(U_global)
        coords = self.get_coordinate_matrix(ndim=2)
        B, _ = _compute_kinematics_tri3(coords)

        if self.discontinuity_state is None:
            strain = B @ u_e
            stress, _, _ = self.material.compute_state(strain, self.state.vars[0])
            centroid = coords.mean(axis=0)
            return {
                'points_natural': np.array([[1.0 / 3.0, 1.0 / 3.0]]),
                'points_global': centroid.reshape(1, 2),
                'strain': strain.reshape(1, 3),
                'stress': stress.reshape(1, 3),
            }

        # Estado agrietado: bulk evaluado con ε neta y discontinuity payload.
        ds = self.discontinuity_state
        G = np.column_stack([ds.normal, ds.tangent])
        i_star = ds.solitary_node
        B_phi = B[:, 2 * i_star:2 * i_star + 2]
        strain_bulk = B @ u_e - B_phi @ (G @ ds.jump_committed)
        stress, _, _ = self.material.compute_state(strain_bulk, self.state.vars[0])
        t_local, _, _ = self.cohesive_material.compute_traction(
            ds.jump_committed, ds.cohesive_state_committed or None
        )
        damage = (ds.cohesive_state_committed or {}).get('damage', 0.0)

        return {
            'points_natural': np.array([[1.0 / 3.0, 1.0 / 3.0]]),
            'points_global': ds.centroid.reshape(1, 2),
            'strain': strain_bulk.reshape(1, 3),
            'stress': stress.reshape(1, 3),
            'discontinuity': {
                'normal': np.copy(ds.normal),
                'tangent': np.copy(ds.tangent),
                'centroid': np.copy(ds.centroid),
                'solitary_node': ds.solitary_node,
                'l_d': ds.l_d,
                'jump': np.copy(ds.jump_committed),
                'traction': np.copy(t_local),
                'damage': float(damage),
            },
        }


# =========================================================================
# Helpers privados (módulo) — fuera de la clase para facilitar tests aislados
# =========================================================================


def _max_principal_stress_2d(sigma: np.ndarray):
    """Devuelve ``(σ_I, n)`` con ``σ_I`` la tensión principal mayor y ``n``
    su autovector unitario en el plano. ``sigma`` en Voigt ``[xx, yy, xy]``.
    """
    sxx, syy, sxy = float(sigma[0]), float(sigma[1]), float(sigma[2])
    mean = 0.5 * (sxx + syy)
    half_diff = 0.5 * (sxx - syy)
    radius = float(np.hypot(half_diff, sxy))
    sigma_I = mean + radius

    # Autovector: ángulo principal 2·θ_p = atan2(2·σ_xy, σ_xx − σ_yy)
    theta_p = 0.5 * np.arctan2(2.0 * sxy, sxx - syy)
    n = np.array([np.cos(theta_p), np.sin(theta_p)])
    return sigma_I, n


def _compute_ld(coords: np.ndarray, tangent: np.ndarray, solitary: int) -> float:
    """Longitud efectiva de Γ_d en el CST: ``l_d = (A_e/h)·|cos(θ − α)|``.

    Implementa la fórmula del Cap. 6 de Retama (2010). ``θ`` es el ángulo de
    la tangente a ``Γ_d`` y ``α`` el del lado opuesto al nodo solitario,
    ambos respecto a ``x`` global. ``A_e`` y ``h`` se calculan a partir de
    las coordenadas. El valor absoluto del coseno absorbe la convención de
    sentido (la longitud es positiva).
    """
    x0, x1, x2 = coords[0], coords[1], coords[2]
    A_e = 0.5 * abs((x1[0] - x0[0]) * (x2[1] - x0[1])
                    - (x2[0] - x0[0]) * (x1[1] - x0[1]))

    others = [i for i in range(3) if i != solitary]
    edge = coords[others[1]] - coords[others[0]]
    L_opp = float(np.hypot(edge[0], edge[1]))
    h = 2.0 * A_e / L_opp

    alpha = float(np.arctan2(edge[1], edge[0]))
    theta = float(np.arctan2(tangent[1], tangent[0]))
    return float((A_e / h) * abs(np.cos(theta - alpha)))
