"""Calculate dimensions for a new station inserted between existing stations.

Inputs
------
A blade definition
Inboard endpoint (existing station 1 in blade definition)
Outboard endpoint (existing station 2 in blade definition)
An airfoil profile (and its thickness-to-chord ratio) (new station)

Outputs (for new station)
-------
x1-location
twist
chord length
selected laminate thicknesses

Author: Perry Roth-Johnson
Last updated: August 15, 2013

"""


import scipy.interpolate as ipl
import scripts.blade as bl


### SET THESE PARAMETERS ------------------------------------------------------
# import the blade definition
b1 = bl.BiplaneBlade(
    'biplane blade, flapwise symmetric, no stagger, rj/R=0.452, g/c=1.25',
    'biplane_flap-sym_no-stagger')
# select the inboard and outboard endpoints
i_endpt = 21
o_endpt = 23
# select airfoil, with given t/c
af_name = 'DU91-W2-250'
af_t_to_c = 0.25
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
# define eqns for c_blade, t_blade, and t_af
c_blade = ipl.interp1d(x=[x_i_stn, x_o_stn],
                       y=[c_i_stn, c_o_stn])
t_blade = ipl.interp1d(x=[x_i_stn, x_o_stn],
                       y=[t_i_stn, t_o_stn])
def t_airfoil(t_to_c, c):
    return t_to_c*c
# set t_blade = t_airfoil, solve for x
def x_new(t_to_c, c_i, t_i, x_i, c_o, t_o, x_o):
    """Find the x1-location for a new station by setting t_blade = t_airfoil.

    Parameters (all floats)
    ----------
    t_to_c : thickness-to-chord ratio of the given airfoil
    c_i : chord of the inboard endpoint station
    t_i : thickness of the inboard endpoint station
    x_i : x1-location of the inboard endpoint station
    c_o : chord of the outboard endpoint station
    t_o : thickness of the outboard endpoint station
    x_o : x1-location of the outboard endpoint station

    """
    numer = t_to_c*c_i - t_i
    denom = ((t_o-t_i)/(x_o-x_i) - t_to_c*(c_o-c_i)/(x_o-x_i))
    return numer/denom + x_i
x_new_stn = x_new(
    t_to_c=af_t_to_c,
    c_i=c_i_stn,
    t_i=t_i_stn,
    x_i=x_i_stn,
    c_o=c_o_stn,
    t_o=t_o_stn,
    x_o=x_o_stn)
# now, given x_new_stn, calculate all the other parameters with lin interp
# (twist, laminate thicknesses, etc)
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
print 'chord = {0:6.3f}'.format(float(c_blade(x_new_stn)))
print 'spar cap height = {0:6.3f}'.format(float(spar_cap_height(x_new_stn)))
print 'shear web 3 x2 = {0:6.3f}'.format(float(SW_3_x2(x_new_stn)))
print 'TE reinf height uniax = {0:6.3f}'.format(float(TE_reinf_height_uniax(x_new_stn)))
print 'TE reinf height foam = {0:6.3f}'.format(float(TE_reinf_height_foam(x_new_stn)))
print 'LE panel height = {0:6.3f}'.format(float(LE_panel_height(x_new_stn)))
print 'aft panel height = {0:6.3f}'.format(float(aft_panel_height(x_new_stn)))
