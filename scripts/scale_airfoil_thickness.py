"""Scale existing airfoil coordinates to a new thickness.

The user can write new coordinates to a file, then plot the original and 
existing coordinates together to verify that scaling was successful. The user 
needs to manually input the original and new (desired) thickness-to-chord 
ratios for the original and new (scaled) airfoils.

Usage
-----
import scale_airfoil_thickness as sat
sat.write_and_plot(airfoil_dir='biplane_flap-sym_no-stagger/airfoils',
    original_path='SNL-100m-Transition51.txt',
    original_tc_ratio=0.51,
    new_path='SNL-100m-Transition51_53pt2.txt',
    new_tc_ratio=0.532)

Author: Perry Roth-Johnson
Last updated: August 16, 2013

"""


import numpy as np
import matplotlib.pyplot as plt
import os


def write_and_plot(airfoil_dir, original_path, original_tc_ratio, new_path,
    new_tc_ratio):
    """Write and plot the original and scaled airfoil coordinates.

    Parameters
    ----------
    airfoil_dir : str, path to directory containing 'original_path' and
        'new_path' airfoil coordinates
    original_path : str, path to file with (x,y) airfoil coordinates
    original_tc_ratio : float, original airfoil thickness-to-chord ratio
    new_path : str, path to new file with scaled (x,y) airfoil coordinates
    new_tc_ratio : float, new desired airfoil thickness-to-chord ratio

    Usage
    -----
    import scale_airfoil_thickness as sat
    sat.write_and_plot(airfoil_dir='biplane_flap-sym_no-stagger/airfoils',
        original_path='SNL-100m-Transition51.txt',
        original_tc_ratio=0.51,
        new_path='SNL-100m-Transition51_53pt2.txt',
        new_tc_ratio=0.532)

    """
    original_path = os.path.join(airfoil_dir, original_path)
    new_path = os.path.join(airfoil_dir, new_path)
    write_scaled_airfoil_coordinates(original_path, original_tc_ratio,
        new_path, new_tc_ratio)
    plot_scaled_airfoil_coordinates(original_path, original_tc_ratio, new_path,
        new_tc_ratio)

def write_scaled_airfoil_coordinates(original_path, original_tc_ratio,
    new_path, new_tc_ratio):
    """Write airfoil coordinates with a new thickness about the chord line.

    Parameters
    ----------
    original_path : str, path to file with (x,y) airfoil coordinates
    original_tc_ratio : float, original airfoil thickness-to-chord ratio
    new_path : str, path to new file with scaled (x,y) airfoil coordinates
    new_tc_ratio : float, new desired airfoil thickness-to-chord ratio

    Writes new airfoil coordinates to the 'new_path'.

    """
    a = np.loadtxt(original_path, dtype=[('x','f8'),('y','f8')])
    a['y'] = a['y']*(new_tc_ratio/original_tc_ratio)
    f = open(new_path, 'w')
    for point in a:
        f.write("{0:8.6f}\t{1:9.6f}\n".format(point['x'],point['y']))
    f.close()


def plot_scaled_airfoil_coordinates(original_path, original_tc_ratio, new_path, new_tc_ratio):
    """Plot original and new (scaled) airfoil coordinates.

    The plot shows a fine grid, which can be used to visually verify the 
    airfoil thicknesses.

    Parameters
    ----------
    path : str, path to file with (x,y) airfoil coordinates
    original_tc_ratio : float, original airfoil thickness-to-chord ratio
    new_path : str, path to new file with scaled (x,y) airfoil coordinates
    new_tc_ratio : float, new desired airfoil thickness-to-chord ratio

    """
    a = np.loadtxt(original_path, dtype=[('x','f8'),('y','f8')])
    b = np.loadtxt(new_path, dtype=[('x','f8'),('y','f8')])
    fig, axes = plt.subplots()
    axes.set_aspect('equal')
    axes.grid('on')
    axes.plot(a['x'],a['y'],'bo-',
        label='original t/c={0}'.format(original_tc_ratio))
    axes.plot(b['x'],b['y'],'rs-',label='new t/c={0}'.format(new_tc_ratio))
    axes.legend()
    plt.yticks(np.arange(-0.5,0.51,0.01))
    plt.xticks(np.arange(0.0,1.01,0.01))
    plt.show()

