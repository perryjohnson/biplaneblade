"""A module for organizing coordinate data for a blade station.

Author: Perry Roth-Johnson
Last updated: August 7, 2013

"""


class Coordinates:
    """Define the coordinates for a blade station."""
    def __init__(self, x1, x2, x3):
        self.x1 = x1   # units [m]
        self.x2 = x2   # units [m]
        self.x3 = x3   # units [m]
    def __str__(self):
        return """x1:  {0:7.4f} (meters)
x2:  {1:7.4f} (meters)
x3:  {2:7.4f} (meters)""".format(self.x1, self.x2, self.x3)
