"""Materiales de interfaz (ADR 0033): leyes tracción–desplazamiento relativo
de la familia ``InterfaceMaterialRegistry`` (sección YAML
``interface_materials``)."""
from solidum.materials.interface.bond_slip import BondSlipModelCode, BondSlipMultilinear
from solidum.materials.interface.cohesive_damage_isotropic import CohesiveDamageIsotropic

__all__ = ["BondSlipModelCode", "BondSlipMultilinear", "CohesiveDamageIsotropic"]
