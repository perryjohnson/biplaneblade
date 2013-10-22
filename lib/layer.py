"""A module for manipulating single material layers of structural parts.

Author: Perry Roth-Johnson
Last updated: October 11, 2013

"""


class Layer:
    """Define a layer made of a single material for a structural part.

    Parameters:
    -----------
    polygon : shapely.Polygon object, represents the boundary of this layer
    material : Material object, the material properties of this layer
    parent_part : structure.Part object, the part this layer belongs to
    mass : float, the mass per unit length of this layer

    """
    def __init__(self, polygon, material, parent_part):
        self.polygon = polygon
        self.material = material
        self.parent_part = parent_part
        self.mass = self.polygon.area*self.material.rho  # mass per unit length
        self.left = None  # saved later by <part>.get_and_save_edges()
        self.top = None  # saved later by <part>.get_and_save_edges()
        self.right = None  # saved later by <part>.get_and_save_edges()
        self.bottom = None  # saved later by <part>.get_and_save_edges()

    def area_fraction(self):
        """Calculate the ratio of this part area to the cross-section area."""
        total_area = self.parent_part.parent_structure.area
        return self.polygon.area/total_area

    def mass_fraction(self):
        """Calculate the ratio of this part mass to the cross-section mass."""
        total_mass = self.parent_part.parent_structure.mass
        return self.mass/total_mass

    def plot_edges(self, axes):
        """Plots the polygon edges of this layer."""
        if self.left is None:
            raise ValueError("Layer instance has attribute <layer>.left=None.\n  Try running <station>.structure.<part>.get_and_save_edges() first.")
        axes.plot(self.left[:,0], self.left[:,1], 'bo-')
        axes.plot(self.top[:,0], self.top[:,1], 'g^-')
        axes.plot(self.right[:,0], self.right[:,1], 'rs-')
        axes.plot(self.bottom[:,0], self.bottom[:,1], 'c*-')