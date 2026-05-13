# ADR 0008 — Densidad como propiedad del material

- **Estado**: aceptado
- **Fecha**: 2026-05-13
- **Alcance**: clase base `Material` y todas las subclases; parser YAML; `Assembler`; catálogo de materiales; manual de usuario.

## Contexto

El cableado de fuerza de cuerpo introducido en los commits previos acepta un vector `body_force: [bx, by, bz]` global en YAML, aplicado uniformemente a todos los elementos. Esto sirve para estructuras monomaterial donde el usuario calcula `ρ·g` una vez y lo declara, pero **falla en estructuras multimaterial**: con un único vector global, no hay forma de distinguir el peso propio del acero (ρ=7850) del hormigón (ρ=2500) o de cables tensados (ρ=7850·η con factor de huecos).

Adicionalmente, la **densidad** es la magnitud física que abre la puerta a varios subsistemas pendientes:

- Matriz de masa elemental `Mₑ = ∫ρ·NᵀN dV`.
- Análisis modal (problema de autovalores generalizado `K·φ = λ·M·φ`).
- Integración temporal implícita (Newmark, HHT-α, Bossak-α).
- Integración temporal explícita (diferencias centradas, condicionalmente estable).

Sin densidad declarada por material, ninguno de estos subsistemas puede construirse. Introducirla ahora — con el peso propio multimaterial como primer consumidor — coloca la base que la dinámica reutilizará cuando llegue.

## Decisión

**La densidad es propiedad del material**, no del elemento. Se introduce como atributo `density` en la clase base `Material` con default `0.0`, y se acepta como parámetro opcional en el constructor de todas las subclases existentes. Físicamente, la densidad es una propiedad intrínseca del medio continuo (kg/m³), no de su discretización; vive en el mismo nivel jerárquico que `E`, `ν`, `σ_y`, `κ_0`, etc. — todas propiedades del material.

### YAML

Se introduce el bloque `gravity` como atajo declarativo para el caso particular `b = ρ·g` (peso propio), complementario al bloque general `body_force` (vector arbitrario uniforme):

```yaml
materials:
  - {id: 1, type: Elastic1D, E: 210.0e9, density: 7850.0}
  - {id: 2, type: Elastic1D, E:  30.0e9, density: 2500.0}

gravity: [0.0, -9.81]   # o [0.0, 0.0, -9.81] en 3D
```

El parser sigue aceptando `body_force` directo (compatibilidad con el cableado actual y útil para fuerzas de cuerpo que no son peso propio: campo electromagnético, fuerzas centrífugas precomputadas, etc.). Si el usuario declara ambos `gravity` y `body_force` en el mismo archivo, se lanza error de validación — el solapamiento es ambigüedad, no flexibilidad.

### Ensamblaje

`Assembler` incorpora un método paralelo a `assemble_body_load(b)`:

```python
def assemble_self_weight(self, g: np.ndarray) -> np.ndarray:
    """Acumula la integral consistente del peso propio, con b_element =
    element.material.density · g por elemento. Si algún material tiene
    density=0, log warning ruidoso (probable olvido del usuario)."""
```

Ambos métodos comparten el bucle interno sobre elementos (factorizado en un helper privado `_iterate_body_force(get_b_for_element)`). La diferencia es solo cómo se obtiene `b` para cada elemento: constante en `assemble_body_load`, derivado de la densidad en `assemble_self_weight`.

### Override a nivel de elemento

La densidad como propiedad del material es el contrato principal. No se introduce override a nivel de elemento por este ADR — cualquier especialización futura (perfiles HEB con masa efectiva tabulada, secciones no prismáticas, masas concentradas) será objeto de un ADR propio cuando llegue, no una promesa anticipada aquí.

## Contrato `Material.density`

```python
class Material(ABC):
    density: Optional[float] = None   # declarada solo si el análisis la requiere
```

- **Opcional al construir**: un material que se usa en análisis estático puro (sin peso propio ni matriz de masa) no necesita declarar densidad. La omisión es válida y el constructor no falla. Default `None` indica "no declarada", semánticamente distinto de "explícitamente cero".
- **Enforcement diferida al uso**: `Assembler.assemble_self_weight(g)` (y, en su día, `assemble_mass_matrix()`, `assemble_centrifugal()`, etc.) inspecciona `material.density` por elemento. Si encuentra `None`, **falla con `ValueError`** listando los materiales afectados por nombre. El fallo es ruidoso e inmediato en el momento en que se pide la magnitud que requiere masa, no antes y no después. Imposible propagar masa cero silenciosa a un análisis dinámico.
- **Valor `0.0` declarado explícitamente**: caso legítimo de material sin masa física (penalty, restricción, fixture). `assemble_self_weight` lo acepta, emite `WARNING` informativo, aporta cero al vector global. Distingue así "no he pensado en la masa" (`None` → error) de "este material no tiene masa por diseño" (`0.0` → warning).
- **Unidades del problema**: como con todas las propiedades materiales, el usuario es responsable de la consistencia de unidades. Si trabaja en N/m, las densidades deben ser kg/m³ con g en m/s²; si trabaja en N/mm, kg/mm³ con g en mm/s². El sistema de tolerancias (ADR 0007) ya garantiza invariancia bajo cambio de unidades.

## Consecuencias

**Inmediatas**:

- Peso propio físicamente correcto en estructuras multimaterial.
- API del YAML más intuitiva (`gravity: [0, -9.81]` versus `body_force: [0, -77008.5]` con el producto precalculado).
- Tests existentes con `body_force` siguen pasando sin modificación.

**Hacia adelante**:

- `Assembler.assemble_mass_matrix()` consume directamente `material.density` — patrón análogo al de `assemble_self_weight`, sin tocar el contrato del material.
- Análisis modal y dinámica transitoria heredan la propiedad sin rediseño.
- Cargas centrífugas (`b = ρ·ω²·r`) encajan en el mismo patrón.

**Deuda heredada**:

- Ninguna asimétrica al presente. La enforcement diferida (`None` al construir → `ValueError` al pedir peso propio o matriz de masa) cubre simultáneamente las dos exigencias en tensión: (a) no obligar a declarar densidad en análisis estáticos que no la usan, y (b) impedir fallos silenciosos de masa cuando entre dinámica.

## Alternativas consideradas

- **Densidad como propiedad del elemento** (`A`, `I`, `density` todos juntos). Rechazada: viola la separación elemento-material que es columna vertebral del proyecto (ver Reglas.md §3 sobre división de dominios). La densidad es física del medio, no de la malla.
- **`mass_per_length` directo en elementos como camino principal**. Más cercano a tablas de catálogo pero asume estructuralmente que todos los elementos son lineales/superficiales; no escala a sólidos 3D donde la masa intrínsecamente es `∫ρ dV`.
- **Mantener `body_force` como única vía**, exigiendo al usuario calcular `ρ·g` antes. Rechazada por el caso multimaterial y por bloquear el camino a dinámica.

## Referencias

- Bathe, K.-J. (2014), *Finite Element Procedures*, cap. 9 (matriz de masa).
- Cook, Malkus & Plesha, *Concepts and Applications of Finite Element Analysis*, §11.3.
- ANSYS Mechanical APDL — `MP, DENS` (densidad por material).
- Abaqus — `*DENSITY` keyword (densidad como atributo del material).
- ADR 0006 — Tolerancias del criterio de admisibilidad (patrón "fórmula centralizada, declarada por el material").
