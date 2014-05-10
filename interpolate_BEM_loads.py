"""Interpolate BEM loads at station locations for Sandia and biplane blades.

Author: Perry Roth-Johnson
Last updated: May 9, 2014

"""


import scipy.interpolate as ipl


# BEM loads: thrust per length (y) at several spanwise locations (x)
f=ipl.interp1d(
    x=[0, 2.164, 6.611, 11.058, 16.617, 23.288, 29.958, 36.629, 43.300, 49.970,
    56.641, 63.312, 69.982, 76.653, 83.323, 88.882, 93.329, 97.776, 100],
    y=[0, 243.878, 299.154, 263.065, 2338.315, 3317.554, 4000.739, 4781.188,
    5851.984, 6812.364, 8020.803, 9184.607, 9815.452, 10768.313, 11609.353,
    12088.440, 11811.550, 8551.274, 0])
# Sandia blade, thrust per length (N/m) at each station's spanwise location
sandia = f([0, 0.5, 0.7, 0.9, 1.1, 1.3, 2.4, 2.6, 4.7, 6.8, 8.9, 11.4, 14.6,
    16.3, 17.9, 19.5, 22.2, 24.9, 27.6, 35.8, 43.9, 52, 60.2, 66.7, 68.3, 73.2,
    76.4, 84.6, 89.4, 94.3, 95.7, 97.2, 98.6, 100])
# biplane blade, thrust per length (N/m) at each station's spanwise location
biplane = f([0, 0.5, 0.7, 0.9, 1.1, 1.3, 2.4, 2.6, 4.7, 6.8, 8.9, 11.4, 14.6,
    16.3, 17.9, 19.5, 22.2, 24.9, 27.6, 30.333, 33.067, 35.8, 38.5, 41.2, 43.9,
    46.6, 49.3, 52, 60.2, 66.7, 68.3, 73.2, 76.4, 84.6, 89.4, 94.3, 95.7, 97.2,
    98.6, 100])
