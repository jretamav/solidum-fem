"""Kernels puntuales del camino por lotes (ADR 0014, fase 2).

Dos contratos, barridos sobre los registros completos para que un
componente nuevo entre solo:

1. **Materiales**: el kernel ``batch_kernel()`` (firma ``MAT_SIG``) debe
   reproducir ``compute_state`` — esfuerzo, tangente y estado — en régimen
   elástico, en régimen inelástico y en descarga desde un estado
   inelástico. Es la misma función compilada en ambos caminos, así que la
   comparación es exacta salvo en los lineales (``C·ε`` por BLAS frente a
   bucle explícito), donde se admite ``1e-14`` relativo.

2. **Elementos**: ``BATCH_KINEMATICS`` (firma ``KIN_SIG``) integra el
   volumen del elemento, produce jacobianos positivos y una ``B`` que
   anula las traslaciones de sólido rígido — las mismas propiedades que el
   barrido de contrato exige a ``K``.
"""
import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import solidum  # noqa: F401
from solidum.math.batch import StateSchema, element_is_batchable, material_is_batchable
from solidum.registry import ElementRegistry, MaterialRegistry, ThermalMaterialRegistry

from test_element_contract_sweep import (
    ESPESOR, FABRICAS, TERMICOS, _construir, nombres_programa_principal,
)
from test_material_contract_sweep import MUESTRAS

RNG = np.random.default_rng(20260922)


def _run_kernel(material, kernel, strain, S_old):
    schema = StateSchema(material.STATE_SCHEMA)
    S_new = np.zeros(schema.n_state)
    sigma = np.zeros(material.STRAIN_DIM)
    flag = np.zeros(1, dtype=np.int8)
    C = np.ascontiguousarray(material.batch_matrix(), dtype=np.float64)
    params = np.ascontiguousarray(material.batch_params(), dtype=np.float64).reshape(-1)
    C_out = np.zeros_like(C)
    Ct = kernel(np.ascontiguousarray(strain, dtype=np.float64), S_old, S_new, params, C, sigma, C_out, flag)
    return sigma, np.asarray(Ct), S_new, int(flag[0])


class TestKernelesDeMaterial(unittest.TestCase):

    def _strains(self, dim):
        small = np.zeros(dim)
        small[0] = 1.0e-8
        large = np.zeros(dim)
        large[0] = 4.0e-2
        large[-1] = 1.5e-2
        if dim == 6:
            large[1] = -1.0e-2
            large[4] = 0.5e-2
        unload = 0.5 * large
        return small, large, unload

    def test_todos_los_materiales_mecanicos_2d_3d_son_batchables(self):
        """Los materiales de sólidos (Voigt 3 y 6) declaran kernel; los 1D
        (fase 4 de la propuesta, diferida) sólo el esquema."""
        for nombre in MaterialRegistry.names():
            with self.subTest(material=nombre):
                material = MaterialRegistry.create(nombre, **MUESTRAS[nombre])
                if material.STRAIN_DIM in (3, 6):
                    self.assertTrue(material_is_batchable(material),
                                    f"{nombre} sin kernel por lotes")
                self.assertIsInstance(material.STATE_SCHEMA, dict)

    def test_kernel_reproduce_compute_state(self):
        for nombre in MaterialRegistry.names():
            material = MaterialRegistry.create(nombre, **MUESTRAS[nombre])
            if not material_is_batchable(material):
                continue
            with self.subTest(material=nombre):
                self._check_material(material, exact=material.PRIMARY_STATE_VAR is not None)

    def test_kernel_reproduce_compute_state_en_plane_stress(self):
        material = MaterialRegistry.create('VonMises2D', E=2.0e11, nu=0.3, sigma_y=2.5e8,
                                           H=1.0e9, hypothesis='plane_stress')
        self._check_material(material, exact=True)
        material = MaterialRegistry.create('IsotropicDamage2D', E=2.0e10, nu=0.2,
                                           kappa_0=1.0e-4, alpha=0.99, hypothesis='plane_strain')
        self._check_material(material, exact=True)
        material = MaterialRegistry.create('DruckerPrager2D', E=2.0e10, nu=0.2, cohesion=1.0e6,
                                           phi_deg=30.0, psi_deg=10.0, H=5.0e5)
        self._check_material(material, exact=True)

    def _check_material(self, material, exact):
        kernel = material.batch_kernel()
        schema = StateSchema(material.STATE_SCHEMA)
        dim = material.STRAIN_DIM
        S = np.zeros(schema.n_state)
        schema.pack(material.initial_state(), S)
        state = None
        for eps in self._strains(dim):
            sig_ref, C_ref, state = material.compute_state(eps.copy(), state)
            sig, Ct, S_new, _flag = _run_kernel(material, kernel, eps, S)
            if exact:
                np.testing.assert_array_equal(sig, np.asarray(sig_ref, float))
                np.testing.assert_array_equal(Ct, np.asarray(C_ref, float))
            else:
                escala = max(np.abs(np.asarray(sig_ref, float)).max(), 1e-300)
                np.testing.assert_allclose(sig, np.asarray(sig_ref, float),
                                           rtol=1.0e-14, atol=1.0e-14 * escala)
                np.testing.assert_array_equal(Ct, np.asarray(C_ref, float))
            fila_ref = np.zeros(schema.n_state)
            schema.pack(state, fila_ref)
            np.testing.assert_array_equal(S_new, fila_ref)
            S = S_new
        # Al menos un paso inelástico para los materiales con historia.
        if schema.n_state:
            self.assertTrue(np.any(S != 0.0) or material.PRIMARY_STATE_VAR == 'damage')

    def test_flag_del_newton_local_plane_stress(self):
        """Un predictor absurdamente lejano agota el Newton local: el kernel
        marca el punto y ``batch_report`` avisa una sola vez."""
        material = MaterialRegistry.create('VonMises2D', E=2.0e11, nu=0.3, sigma_y=1.0,
                                           H=0.0, hypothesis='plane_stress')
        kernel = material.batch_kernel()
        S = np.zeros(5)
        eps = np.array([1.0e3, 0.0, 0.0])
        _, _, _, flag = _run_kernel(material, kernel, eps, S)
        # Con o sin convergencia, el flag es coherente con compute_state:
        material._local_newton_warned = False
        material.compute_state(eps, None)
        self.assertEqual(flag, int(material._local_newton_warned))
        material._local_newton_warned = False
        material.batch_report(3)
        self.assertTrue(material._local_newton_warned)

    def test_material_termico_es_batchable(self):
        for nombre in ThermalMaterialRegistry.names():
            with self.subTest(material=nombre):
                mat = ThermalMaterialRegistry.create(nombre, k=50.0, dim=2)
                self.assertTrue(material_is_batchable(mat))
                grad = np.array([1.0, -2.0])
                q_ref, k_ref = mat.compute_flux(grad)
                sigma, Ct, _, _ = _run_kernel_thermal(mat, grad)
                np.testing.assert_allclose(sigma, -q_ref, rtol=1e-14)
                np.testing.assert_array_equal(Ct, k_ref)


def _run_kernel_thermal(mat, grad):
    kernel = mat.batch_kernel()
    C = np.ascontiguousarray(mat.batch_matrix(), dtype=np.float64)
    sigma = np.zeros(C.shape[0])
    flag = np.zeros(1, dtype=np.int8)
    S = np.zeros(0)
    Ct = kernel(np.ascontiguousarray(grad, float), S, S.copy(), np.zeros(0), C, sigma, np.zeros_like(C), flag)
    return sigma, np.asarray(Ct), None, int(flag[0])


class TestCinematicaDeElementos(unittest.TestCase):

    def _perturbed(self, elem):
        """Desplaza ligeramente los nodos para que el jacobiano no sea constante."""
        for node in elem.nodes:
            c = list(node.coordinates)
            for k in range(len(c)):
                c[k] += 0.03 * (RNG.random() - 0.5)
            node.coordinates = c
        return elem

    def test_elementos_batchables(self):
        """Todos los sólidos isoparamétricos y los térmicos declaran
        cinemática; los 1D no (quedan en el camino por elemento)."""
        esperados_no = {'Truss2D', 'Truss2DCorot', 'Truss3D',
                        'Truss3DCorot', 'Cable2DCorot', 'Cable3DCorot', 'Frame2DEuler',
                        'Frame2DEulerCorot', 'Frame2DTimoshenko', 'Frame3D'}
        for nombre in nombres_programa_principal(ElementRegistry):
            with self.subTest(elemento=nombre):
                elem = _construir(nombre)
                self.assertEqual(element_is_batchable(elem), nombre not in esperados_no)

    def test_cinematica_integra_el_volumen_y_anula_traslaciones(self):
        for nombre in nombres_programa_principal(ElementRegistry):
            elem = _construir(nombre)
            if not element_is_batchable(elem):
                continue
            with self.subTest(elemento=nombre):
                kin = type(elem).BATCH_KINEMATICS
                pts, w = elem.batch_quadrature()
                d = pts.shape[1]
                coords = elem.batch_reference_coordinates(d)
                n_dof = len(elem.DOF_NAMES) * len(elem.nodes)
                strain_dim = elem.STRAIN_DIM if elem.STRAIN_DIM is not None else elem.FLUX_DIM
                B = np.zeros((strain_dim, n_dof))
                volumen = 0.0
                for g in range(pts.shape[0]):
                    B.fill(np.nan)
                    detJ = kin(pts[g], coords, B)
                    self.assertGreater(detJ, 0.0)
                    self.assertTrue(np.all(np.isfinite(B)), f"{nombre}: B con NaN (fila sin rellenar)")
                    # Traslación uniforme en cada dirección: B·t = 0.
                    for eje in range(len(elem.DOF_NAMES)):
                        t = np.zeros(n_dof)
                        t[eje::len(elem.DOF_NAMES)] = 1.0
                        self.assertLess(np.abs(B @ t).max(), 1e-12 * max(np.abs(B).max(), 1.0))
                    volumen += detJ * w[g]
                esperado = FABRICAS[nombre][1]
                if nombre in TERMICOS:
                    esperado = 1.0 * (ESPESOR if d == 2 else 1.0)
                if esperado is not None:
                    escala = getattr(elem, type(elem).BATCH_SCALE) if type(elem).BATCH_SCALE else 1.0
                    self.assertAlmostEqual(volumen * escala / esperado, 1.0, places=12)

    def test_cinematica_por_lotes_coincide_con_compute_gauss_state(self):
        """Deformación por punto de Gauss: ``B_batch · u_e`` == la que reporta
        ``compute_gauss_state`` (camino por elemento), en geometría perturbada."""
        for nombre in nombres_programa_principal(ElementRegistry):
            if nombre in TERMICOS:
                continue
            elem = _construir(nombre)
            if not element_is_batchable(elem):
                continue
            with self.subTest(elemento=nombre):
                self._perturbed(elem)
                kin = type(elem).BATCH_KINEMATICS
                pts, _ = elem.batch_quadrature()
                d = pts.shape[1]
                coords = elem.batch_reference_coordinates(d)
                n_dof = len(elem.DOF_NAMES) * len(elem.nodes)
                u_e = 1.0e-3 * RNG.random(n_dof)
                # Campo global con sólo este elemento: DOFs numerados en orden.
                eq = 0
                for node in elem.nodes:
                    for name in node.dofs:
                        node.dofs[name] = eq
                        eq += 1
                U = np.zeros(eq)
                U[elem.get_global_dof_indices()] = u_e
                gs = elem.compute_gauss_state(U)
                B = np.zeros((elem.STRAIN_DIM, n_dof))
                for g in range(pts.shape[0]):
                    kin(pts[g], coords, B)
                    np.testing.assert_allclose(B @ u_e, gs['strain'][g], rtol=1e-13, atol=1e-18)

    def test_jacobiano_invertido_devuelve_detj_no_positivo_sin_lanzar(self):
        """Contrato ``KIN_SIG``: la cinemática por lotes **no lanza** (una
        excepción dentro de ``prange`` se perdería); señala el jacobiano
        degenerado devolviendo ``det J ≤ 0`` y sin escribir ``B``. Quien
        lanza es el ensamblador (``test_batch_assembly``)."""
        for nombre, dim, n_sig in (('Quad4', 2, 3), ('Hex8', 3, 6), ('Tet10', 3, 6), ('Tri6', 2, 3)):
            with self.subTest(elemento=nombre):
                elem = _construir(nombre)
                kin = type(elem).BATCH_KINEMATICS
                pts, _ = elem.batch_quadrature()
                coords = elem.batch_reference_coordinates(dim)
                # Reflejar el elemento invierte la orientación (det J < 0).
                coords = coords.copy()
                coords[:, 0] *= -1.0
                B = np.full((n_sig, dim * len(elem.nodes)), 7.0)
                detJ = kin(pts[0], coords, B)
                self.assertLessEqual(detJ, 0.0)
                self.assertTrue(np.all(B == 7.0))    # B intacta
                # La versión del camino por elemento sí lanza.
                with self.assertRaises(ValueError):
                    elem.get_coordinate_matrix = lambda ndim=dim, c=coords: c
                    elem.compute_gauss_state(np.zeros(elem.get_global_dof_indices().__len__()))


if __name__ == '__main__':
    unittest.main()
