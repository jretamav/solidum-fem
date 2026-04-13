from typing import Dict, Type, Any

class MaterialRegistry:
    """Registry para materiales. Mapea nombres a clases constructoras."""
    _materials: Dict[str, Type] = {}
    
    @classmethod
    def register(cls, name: str, material_class: Type) -> None:
        cls._materials[name] = material_class
        
    @classmethod
    def create(cls, name: str, **kwargs) -> Any:
        if name not in cls._materials:
            raise ValueError(f"Material '{name}' no registrado. Disponibles: {', '.join(cls._materials.keys())}")
        return cls._materials[name](**kwargs)

class ElementRegistry:
    """Registry para elementos. Mapea nombres a clases constructoras."""
    _elements: Dict[str, Type] = {}
    
    @classmethod
    def register(cls, name: str, element_class: Type) -> None:
        cls._elements[name] = element_class
        
    @classmethod
    def create(cls, name: str, **kwargs) -> Any:
        if name not in cls._elements:
            raise ValueError(f"Elemento '{name}' no registrado. Disponibles: {', '.join(cls._elements.keys())}")
        return cls._elements[name](**kwargs)

class SolverRegistry:
    """Registry para solucionadores."""
    _solvers: Dict[str, Type] = {}
    
    @classmethod
    def register(cls, name: str, solver_class: Type) -> None:
        cls._solvers[name] = solver_class
        
    @classmethod
    def create(cls, name: str, **kwargs) -> Any:
        if name not in cls._solvers:
            raise ValueError(f"Solucionador '{name}' no registrado. Disponibles: {', '.join(cls._solvers.keys())}")
        return cls._solvers[name](**kwargs)

class QuadratureRegistry:
    """Registry para reglas de integración (puntos, pesos)."""
    _rules: Dict[str, tuple] = {}
    
    @classmethod
    def register(cls, name: str, points: list, weights: list) -> None:
        cls._rules[name] = (points, weights)
        
    @classmethod
    def get(cls, name: str) -> tuple:
        if name not in cls._rules:
            raise ValueError(f"Cuadratura '{name}' no registrada. Disponibles: {', '.join(cls._rules.keys())}")
        return cls._rules[name]