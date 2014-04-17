"""Plot flapwise deflection vs. span for the Sandia blade.

Author: Perry Roth-Johnson
Last updated: April 10, 2014

"""


import numpy as np
import matplotlib.pyplot as plt
import lib.blade as bl
from scipy.interpolate import interp1d
from scipy import array


### parameters ###
glw=1.0            # global line width
gms=8.0           # global marker size
gmew=1.5           # global marker edge width
gmfc='None'        # global marker face color (none/empty)
mec_ms='red'       # marker edge color for monoplane spars


def extrap1d(interpolator):
    # copied from http://stackoverflow.com/questions/2745329/how-to-make-scipy-interpolate-give-a-an-extrapolated-result-beyond-the-input-ran
    xs = interpolator.x
    ys = interpolator.y
    def pointwise(x):
        if x < xs[0]:
            return ys[0]+(x-xs[0])*(ys[1]-ys[0])/(xs[1]-xs[0])
        elif x > xs[-1]:
            return ys[-1]+(x-xs[-1])*(ys[-1]-ys[-2])/(xs[-1]-xs[-2])
        else:
            return interpolator(x)
    def ufunclike(xs):
        return array(map(pointwise, array(xs)))
    return ufunclike

def plot_monoblade_displacement(component, skip_num=1, span=100.0, 
    filename='sandia_blade/beam_model/FIGURES/svy_disp_blade.mdt'):
    """Plot the displacement or rotation vs. span for a monoplane blade."""
    AB = np.loadtxt(filename)
    B = span
    AB[:,0] = AB[:,0]*B
    if component == 'flapwise':
        col=3
    else:
        raise NotImplementedError("`component` keyword must be 'flapwise'")
    plt.plot(AB[::skip_num,0], AB[::skip_num,col], 'rs--', markerfacecolor=gmfc,
        markersize=gms, linewidth=glw, markeredgewidth=gmew, 
        markeredgecolor=mec_ms, label='Sandia blade (beam model, DYMORE)',
        zorder=2)
    plt.xlabel('span [m]')
    plt.ylabel(component+' displacement [m]')

def plot_monoblade_force(component, x1, span=100.0,
    filename='sandia_blade/beam_model/FIGURES/svy_force_blade.mdt'):
    """Plot the force or moment resultant vs. span for a monoplane blade.

Input
-----
x1 <np.array>: x1-coordinates of all 24 spar stations

    """
    # read in all forces and moments calculated by DYMORE
    AB = np.loadtxt(filename)
    # the blade length, 100.0 meters
    B = span
    # multiply column 0 (eta-coordinates) by the span length, B
    AB[:,0] = AB[:,0]*B
    # plot the span (column 0) along the x-axis -------------------------------
    x = AB[:,0]
    if component == 'axial force':
        col = 1
    elif component == 'flapwise bending moment':
        col = 5
    else:
        raise NotImplementedError("`component` keyword must be 'axial force' or 'flapwise bending moment'")
    y = AB[:,col]
    # get force results at all the spar stations, using interpolation
    # (instead of using the default results at Gaussian integration points)
    f_i = interp1d(x,y)
    f_x = extrap1d(f_i)
    y1 = f_x(x1)  # use extrapolations function returned by 'extrap1d'
    # plot the results to the screen ------------------------------------------
    plt.plot(x1, y1/1000.0, 'rs--', markerfacecolor=gmfc, markersize=gms, 
        linewidth=glw, markeredgewidth=gmew, markeredgecolor=mec_ms, 
        label='monoplane blade', zorder=2)

# -----------------------------------------------------------------------------

plt.close('all')
# load the Sandia blade into memory
m = bl.MonoplaneBlade('Sandia blade SNL100-00', 'sandia_blade')
# x1 coordinates of all 24 spar stations
x1_stn = np.array([
     0.0,
     0.5,
     0.7,
     0.9,
     1.1,
     1.3,
     2.4,
     2.6,
     4.7,
     6.8,
     8.9,
    11.4,
    14.6,
    16.3,
    17.9,
    19.5,
    22.2,
    24.9,
    27.6,
    35.8,
    43.9,
    52.0,
    60.2,
    66.7,
    68.3,
    73.2,
    76.4,
    84.6,
    89.4,
    94.3,
    95.7,
    97.2,
    98.6,
   100.0
])

# plot flapwise deflection vs. span
plt.figure()
plot_monoblade_displacement(component='flapwise')
plt.legend(loc='upper left')
plt.show()
