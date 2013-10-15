"""A module for manipulating single material layers of structural parts.

Author: Perry Roth-Johnson
Last updated: October 11, 2013

"""


class Layer:
    """Define a layer made of a single material for a structural part.

    Parameters:
    -----------
    polygon : shapely.Polygon object, represents the boundary of this layer
    material : Material object, represents the material properties of this layer

    """
    def __init__(self, polygon, material, parent_part):
        self.polygon = polygon
        self.material = material
        self.parent_part = parent_part
        self.mass = self.polygon.area*self.material.rho  # mass per unit length

    def area_fraction(self):
        """Calculate the ratio of this part area to the cross-section area."""
        total_area = self.parent_part.parent_structure.area
        return self.polygon.area/total_area

    # def calculate_mass(self):
    #     """Calculate the mass (per unit length) of this part."""
    #     self.mass = self.polygon.area*self.material.rho
    #     return self.mass

    def mass_fraction(self):
        """Calculate the ratio of this part mass to the cross-section mass."""
        total_mass = self.parent_part.parent_structure.mass
        return self.mass/total_mass