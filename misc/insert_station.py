"""Calculate dimensions for a new station inserted between existing stations.

Inputs
------
A blade definition
Inboard endpoint (existing station 1 in blade definition)
Outboard endpoint (existing station 2 in blade definition)
x1-fraction distance to move between the endpoint stations (for example, a
    fraction distance of 0.5 selects a new station halfway between the selected
    endpoints)

Outputs (for new station)
-------
x1-location
twist
chord length
selected laminate thicknesses

Usage
-----
update the parameters in the 'SET THESE PARAMETERS' section below
run this script in an IPython terminal:
    |> %run insert_station
figure out the required thickness-to-chord ratio(s) at the new station
use scale_airfoils.py to create new airfoil profile(s) for the new station

Author: Perry Roth-Johnson
Last updated: August 16, 2013

"""


import scipy.interpolate as ipl
import lib.blade as bl


### SET THESE PARAMETERS ------------------------------------------------------
# import the blade definition
b1 = bl.BiplaneBlade(
    'biplane blade, flapwise symmetric, no stagger, rj/R=0.452, g/c=1.25',
    'biplane_flap-sym_no-stagger')
# select the inboard and outboard endpoints
i_endpt = 22
o_endpt = 24
# select a fraction distance to move between the endpoint stations
x_new_frac = 2.0/3.0
### ---------------------------------------------------------------------------

### CRUNCH THE NUMBERS
# grab the inboard and outboard endpoint stations
i_stn = b1._df.ix[i_endpt]
o_stn = b1._df.ix[o_endpt]
# collect x1, c, and t at endpts
x_i_stn = i_stn['x1']
c_i_stn = i_stn['chord']
t_i_stn = i_stn['thickness-to-chord ratio']*i_stn['chord']
x_o_stn = o_stn['x1']
c_o_stn = o_stn['chord']
t_o_stn = o_stn['thickness-to-chord ratio']*o_stn['chord']
# calculate the new x1-location
x_new_stn = x_new_frac*(x_o_stn-x_i_stn) + x_i_stn
# define eqn for airfoil thickness
def t_airfoil(t_to_c, c):
    """Airfoil thickness."""
    return t_to_c*c
# define interpolation equations for other blade quantities
chord = ipl.interp1d(x=[x_i_stn, x_o_stn], y=[c_i_stn, c_o_stn])
thickness = ipl.interp1d(x=[x_i_stn, x_o_stn], y=[t_i_stn, t_o_stn])
twist = ipl.interp1d(
    x=[i_stn['x1'], o_stn['x1']],
    y=[i_stn['twist'], o_stn['twist']])
spar_cap_height = ipl.interp1d(
    x=[i_stn['x1'], o_stn['x1']],
    y=[i_stn['spar cap height'], o_stn['spar cap height']])
SW_3_x2 = ipl.interp1d(
    x=[i_stn['x1'], o_stn['x1']],
    y=[i_stn['shear web 3 x2'], o_stn['shear web 3 x2']])
TE_reinf_height_uniax = ipl.interp1d(
    x=[i_stn['x1'], o_stn['x1']],
    y=[i_stn['TE reinf height uniax'], o_stn['TE reinf height uniax']])
TE_reinf_height_foam = ipl.interp1d(
    x=[i_stn['x1'], o_stn['x1']],
    y=[i_stn['TE reinf height foam'], o_stn['TE reinf height foam']])
LE_panel_height = ipl.interp1d(
    x=[i_stn['x1'], o_stn['x1']],
    y=[i_stn['LE panel height'], o_stn['LE panel height']])
aft_panel_height = ipl.interp1d(
    x=[i_stn['x1'], o_stn['x1']],
    y=[i_stn['aft panel height'], o_stn['aft panel height']])

# PRINT RESULTS
# print results to the screen, in the same order as in the CSV file
print ''
print 'x1 = {0:6.3f}'.format(x_new_stn)
print 'twist = {0:6.3f}'.format(float(twist(x_new_stn)))
print 'chord = {0:6.3f}'.format(float(chord(x_new_stn)))
print 'thickness-to-chord ratio = {0:6.3f}'.format(float(thickness(x_new_stn)/chord(x_new_stn)))
print 'spar cap height = {0:6.3f}'.format(float(spar_cap_height(x_new_stn)))
print 'shear web 3 x2 = {0:6.3f}'.format(float(SW_3_x2(x_new_stn)))
print 'TE reinf height uniax = {0:6.3f}'.format(float(TE_reinf_height_uniax(x_new_stn)))
print 'TE reinf height foam = {0:6.3f}'.format(float(TE_reinf_height_foam(x_new_stn)))
print 'LE panel height = {0:6.3f}'.format(float(LE_panel_height(x_new_stn)))
print 'aft panel height = {0:6.3f}'.format(float(aft_panel_height(x_new_stn)))
