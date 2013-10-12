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
    def __init__(self, polygon, material):
        self.polygon = polygon
        self.material = material
