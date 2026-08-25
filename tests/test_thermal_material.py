"""Tests del material térmico ``ThermalConduction`` (Etapa 8).

Cubre el bloque ``acceptance`` de ``docs/specs/ThermalConduction.md``.

Los coeficientes NO son unitarios a propósito: se usan valores físicos de
acero (k=45 W/(m·K), c=460 J/(kg·K), ρ=7850 kg/m³) y tensores con
componentes distintas entre sí, para que un error dimensional o un factor
perdido no quede enmascarado por multiplicaciones por 1.
"""
import numpy as np
import pytest

from solidum.materials.thermal_conduction import ThermalConduction

# Propiedades del acero al carbono — todas distintas de 1.
K_ACERO = 45.0
C_ACERO = 460.0
RHO_ACERO = 7850.0


class TestLeyDeFourier:
    """Verificación: q = -k·∇T exacto, isótropo y anisótropo."""

    def test_ley_de_fourier_isotropa_exacta(self):
        mat = ThermalConduction(k=K_ACERO, dim=2)
        grad = np.array([12.5, -3.75])
        q, k = mat.compute_flux(grad)

        np.testing.assert_allclose(q, -K_ACERO * grad, rtol=1e-14)
        np.testing.assert_allclose(k, K_ACERO * np.eye(2), rtol=1e-14)

    def test_ley_de_fourier_tensorial_exacta(self):
        k_tensor = np.array([[45.0, 7.0], [7.0, 12.0]])
        mat = ThermalConduction(k=k_tensor)
        grad = np.array([2.0, -5.0])
        q, k = mat.compute_flux(grad)

        np.testing.assert_allclose(q, -k_tensor @ grad, rtol=1e-14)
        np.testing.assert_allclose(k, k_tensor, rtol=1e-14)

    def test_flujo_opuesto_al_gradiente_isotropo(self):
        """q·∇T < 0 estricto: disipación termodinámica en toda dirección."""
        mat = ThermalConduction(k=K_ACERO, dim=3)
        rng = np.random.default_rng(20260825)
        for _ in range(20):
            grad = rng.normal(size=3)
            q, _ = mat.compute_flux(grad)
            assert float(q @ grad) < 0.0
            # Isótropo: además antiparalelo exacto.
            np.testing.assert_allclose(q, -K_ACERO * grad, rtol=1e-14)

    def test_anisotropia_desvia_el_flujo(self):
        """El test que distingue el contrato tensorial del escalar.

        Con k = diag(k1, k2) y k1 ≠ k2, el flujo deja de ser antiparalelo
        al gradiente: se desvía hacia el eje de mayor conductividad.
        """
        k1, k2 = 50.0, 5.0
        mat = ThermalConduction(k=np.diag([k1, k2]))
        grad = np.array([1.0, 1.0])  # a 45° de los ejes
        q, _ = mat.compute_flux(grad)

        np.testing.assert_allclose(q, np.array([-k1, -k2]), rtol=1e-12)

        # No antiparalelo: el ángulo con -∇T es estrictamente positivo.
        anti = -grad / np.linalg.norm(grad)
        q_dir = q / np.linalg.norm(q)
        assert not np.allclose(q_dir, anti, atol=1e-6)

        # Se desvía hacia el eje de mayor conductividad (x).
        assert abs(q[0]) > abs(q[1])

        # Sigue disipando pese a la desviación.
        assert float(q @ grad) < 0.0


class TestExpansionEscalarATensor:
    def test_expansion_2d(self):
        mat = ThermalConduction(k=K_ACERO, dim=2)
        np.testing.assert_allclose(
            mat.conductivity, K_ACERO * np.eye(2), atol=1e-14
        )
        assert mat.FLUX_DIM == 2
        assert mat.is_isotropic

    def test_expansion_3d(self):
        mat = ThermalConduction(k=K_ACERO, dim=3)
        np.testing.assert_allclose(
            mat.conductivity, K_ACERO * np.eye(3), atol=1e-14
        )
        assert mat.FLUX_DIM == 3
        assert mat.is_isotropic

    def test_tensor_3d_define_su_dimension(self):
        """Un k matricial manda sobre `dim`: el tamaño lo fija el tensor."""
        mat = ThermalConduction(k=np.diag([10.0, 20.0, 30.0]), dim=2)
        assert mat.FLUX_DIM == 3
        assert not mat.is_isotropic

    def test_tensor_multiplo_de_identidad_es_isotropo(self):
        mat = ThermalConduction(k=K_ACERO * np.eye(3))
        assert mat.is_isotropic


class TestSimetriaYPositividad:
    def test_simetria_exacta(self):
        for k in (K_ACERO, np.array([[45.0, 7.0], [7.0, 12.0]])):
            mat = ThermalConduction(k=k, dim=2)
            np.testing.assert_allclose(
                mat.conductivity, mat.conductivity.T, atol=1e-14
            )

    def test_autovalores_estrictamente_positivos(self):
        k_tensor = np.array([[45.0, 7.0], [7.0, 12.0]])
        mat = ThermalConduction(k=k_tensor)
        assert np.min(np.linalg.eigvalsh(mat.conductivity)) > 0.0


class TestDifusividad:
    def test_difusividad_derivada(self):
        mat = ThermalConduction(k=K_ACERO, c=C_ACERO, density=RHO_ACERO, dim=2)
        esperado = K_ACERO / (RHO_ACERO * C_ACERO)
        assert mat.thermal_diffusivity == pytest.approx(esperado, rel=1e-14)

    def test_capacidad_volumetrica(self):
        mat = ThermalConduction(k=K_ACERO, c=C_ACERO, density=RHO_ACERO)
        assert mat.volumetric_capacity() == pytest.approx(
            RHO_ACERO * C_ACERO, rel=1e-14
        )

    def test_difusividad_anisotropa_rechazada(self):
        """No es escalar cuando k es tensorial: sería un tensor."""
        mat = ThermalConduction(
            k=np.diag([50.0, 5.0]), c=C_ACERO, density=RHO_ACERO
        )
        with pytest.raises(ValueError, match="anisótropa"):
            _ = mat.thermal_diffusivity


class TestSinEstadoInterno:
    def test_material_es_sin_memoria(self):
        """Evaluaciones repetidas no alteran nada: modelo lineal sin historia."""
        mat = ThermalConduction(k=K_ACERO, dim=2)
        grad_a = np.array([10.0, 0.0])
        grad_b = np.array([0.0, -4.0])

        q_a1, _ = mat.compute_flux(grad_a)
        mat.compute_flux(grad_b)
        mat.compute_flux(grad_b * 100.0)
        q_a2, _ = mat.compute_flux(grad_a)

        np.testing.assert_allclose(q_a1, q_a2, rtol=1e-15)
        assert ThermalConduction.PRIMARY_STATE_VAR is None


class TestRechazoDeInputsInvalidos:
    def test_k_escalar_no_positivo(self):
        for k in (0.0, -45.0):
            with pytest.raises(ValueError, match="estrictamente positivo"):
                ThermalConduction(k=k, dim=2)

    def test_k_no_simetrico(self):
        with pytest.raises(ValueError, match="simétrico"):
            ThermalConduction(k=np.array([[1.0, 2.0], [3.0, 1.0]]))

    def test_k_no_definido_positivo(self):
        # Autovalores 3 y -1: viola la segunda ley.
        with pytest.raises(ValueError, match="definido positivo"):
            ThermalConduction(k=np.array([[1.0, 2.0], [2.0, 1.0]]))

    def test_k_no_cuadrado(self):
        with pytest.raises(ValueError, match="cuadrada"):
            ThermalConduction(k=np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]))

    def test_k_de_tamano_no_soportado(self):
        with pytest.raises(ValueError, match="no soportado"):
            ThermalConduction(k=np.eye(4))

    def test_dim_invalida(self):
        with pytest.raises(ValueError, match="dim=4"):
            ThermalConduction(k=K_ACERO, dim=4)

    def test_c_no_positivo(self):
        with pytest.raises(ValueError, match="c=.*estrictamente positivo"):
            ThermalConduction(k=K_ACERO, c=0.0, dim=2)

    def test_density_no_positiva(self):
        with pytest.raises(ValueError, match="density=.*estrictamente positiva"):
            ThermalConduction(k=K_ACERO, density=-1.0, dim=2)

    def test_gradiente_de_dimension_incompatible(self):
        mat = ThermalConduction(k=K_ACERO, dim=2)
        with pytest.raises(ValueError, match="FLUX_DIM=2"):
            mat.compute_flux(np.array([1.0, 2.0, 3.0]))


class TestCapacidadAusente:
    """El criterio del ADR 0008 trasladado al térmico, con mensaje accionable."""

    def test_estacionario_no_exige_capacidad(self):
        """Sólo con k, el material se construye y calcula flujo sin error."""
        mat = ThermalConduction(k=K_ACERO, dim=2)
        q, _ = mat.compute_flux(np.array([1.0, 0.0]))
        np.testing.assert_allclose(q, [-K_ACERO, 0.0], rtol=1e-14)
        assert mat.density is None
        assert mat.specific_heat is None

    def test_falla_al_pedir_capacidad_sin_declararla(self):
        mat = ThermalConduction(k=K_ACERO, dim=2)
        with pytest.raises(ValueError) as exc:
            mat.volumetric_capacity()
        assert "ThermalConduction" in str(exc.value)

    @pytest.mark.parametrize(
        "kwargs, esperados",
        [
            ({}, ["'density'", "'c'"]),
            ({"c": C_ACERO}, ["'density'"]),
            ({"density": RHO_ACERO}, ["'c'"]),
        ],
    )
    def test_mensaje_de_error_accionable(self, kwargs, esperados):
        """No basta con reportar el atributo ausente: debe decir qué hacer.

        Requisito fijado por el usuario al cerrar el §Diálogo de la spec.
        """
        mat = ThermalConduction(k=K_ACERO, dim=2, **kwargs)
        with pytest.raises(ValueError) as exc:
            mat.volumetric_capacity()
        msg = str(exc.value)

        # Nombra el material y cada parámetro que falta...
        assert "ThermalConduction" in msg
        for nombre in esperados:
            assert nombre in msg
        # ...y no nombra el que sí está.
        if "c" in kwargs:
            assert "'c'" not in msg
        if "density" in kwargs:
            assert "'density'" not in msg

        # ...e indica la salida, que es lo que lo hace accionable.
        assert "estacionario" in msg or "LinearSolver" in msg
        assert "Decláralo" in msg  # cubre "Decláralo" y "Decláralos"

    def test_consumidor_aparece_en_el_mensaje(self):
        """El mensaje sitúa el contexto en que se pidió la capacidad."""
        mat = ThermalConduction(k=K_ACERO, dim=2)
        with pytest.raises(ValueError, match="la matriz de capacidad"):
            mat.volumetric_capacity(consumer="la matriz de capacidad")


class TestRegistro:
    def test_esta_registrado(self):
        from solidum.registry import ThermalMaterialRegistry

        assert "ThermalConduction" in ThermalMaterialRegistry.names()
        assert ThermalMaterialRegistry.get("ThermalConduction") is ThermalConduction

    def test_no_contamina_el_registro_mecanico(self):
        """Familias paralelas: el térmico no aparece entre los mecánicos."""
        from solidum.registry import MaterialRegistry

        assert "ThermalConduction" not in MaterialRegistry.names()
