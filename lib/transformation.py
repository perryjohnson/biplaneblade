"""A collection of functions for 2D transformations.

Author: Perry Roth-Johnson
Last updated: August 6, 2013

"""


import numpy as np


def rot_mat(theta):
    """Returns a 2D rotation matrix, given a rotation angle theta.

    Ref: http://en.wikipedia.org/wiki/Rotation_matrix

    """
    return np.array([[np.cos(theta), -np.sin(theta)],
                     [np.sin(theta),  np.cos(theta)]])

def rotate_coord_pair(x, y, t, degree_units=True):
    """Rotate a single (x,y) coordinate pair by theta degrees.

    Returns a new (x,y) coordinate pair that has been rotated.

    Parameters
    ----------
    x : float, x-coordinate (unrotated)
    y : float, y-coordinate (unrotated)
    t : float, rotation angle (default units: degrees)
    degree_units : boolean, True if t is in units of degrees, False if t is
        in units of radians

    """
    if degree_units:
        t = np.deg2rad(t)  # convert rotation angle to radians
    R = rot_mat(t)
    p_old = np.array([[x],
                      [y]])
    p_new = np.dot(R, p_old)
    return (float(p_new[0]), float(p_new[1]))
    