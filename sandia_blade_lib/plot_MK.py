"""Plot the mass and stiffness data from Sandia and VABS.

First, data from mass and stiffness matrices for the Sandia blade are written
to the file 'sandia_blade/blade_props_from_VABS.csv'

Then, these data are plotted against published data from Griffith & Resor 2011.

Usage
-----
Open an IPython terminal and type:
|> %run plot_MK

Author: Perry Roth-Johnson
Last modified: April 9, 2014

"""

import lib.blade as bl
reload(bl)
import pandas as pd
import matplotlib.pyplot as plt


def rel_diff(vabs_data, sandia_data):
    """Calculate the percent relative difference."""
    return ((vabs_data-sandia_data)/sandia_data)*100.0

def prep_rel_diff_plot(axis,ymin=-100,ymax=100):
    """Prepare a relative difference plot."""
    axis2 = axis.twinx()
    axis2.set_ylabel('difference [%]', color='m')
    axis2.set_ylim([ymin,ymax])
    for tl in axis2.get_yticklabels():
        tl.set_color('m')
    axis2.grid('on')
    return axis2

# write all the mass and stiffness matrices from VABS to a csv file -----------
m = bl.MonoplaneBlade('Sandia blade SNL100-00', 'sandia_blade')
m.writecsv_mass_and_stiffness_props()

# plot VABS and Sandia datasets against one another ---------------------------
v=pd.DataFrame.from_csv('sandia_blade/blade_props_from_VABS.csv')
s=pd.DataFrame.from_csv('sandia_blade/blade_props_from_Sandia.csv')
plt.close('all')

# stiffness properties --------------------------------------------------------
f, axarr = plt.subplots(2,2, figsize=(12,8))
# ref for dual-axis plotting: http://matplotlib.org/examples/api/two_scales.html

# flapwise stiffness
twin_axis00 = prep_rel_diff_plot(axarr[0,0])
twin_axis00.plot(
    s['Blade Spanwise Coordinate'], rel_diff(v['K_55, EI_flap'],s['EI_flap']),
    'm^:', mec='m', mfc='None', mew=1, label='difference')
axarr[0,0].plot(s['Blade Spanwise Coordinate'],s['EI_flap'],'gx--',mfc='None',mew=1,label='Sandia (PreComp)')
axarr[0,0].plot(v['Blade Spanwise Coordinate'],v['K_55, EI_flap'],'ko-',mfc='None',mew=1,label='UCLA (VABS)')
axarr[0,0].set_xlabel('span [m]')
axarr[0,0].set_ylabel('flapwise stiffness [N*m^2]')
axarr[0,0].legend()
axarr[0,0].grid('on', axis='x')

# edgewise stiffness
twin_axis01 = prep_rel_diff_plot(axarr[0,1])
twin_axis01.plot(
    s['Blade Spanwise Coordinate'], rel_diff(v['K_66, EI_edge'],s['EI_edge']),
    'm^:', mec='m', mfc='None', mew=1, label='difference')
axarr[0,1].plot(s['Blade Spanwise Coordinate'],s['EI_edge'],'gx--',mfc='None',mew=1,label='Sandia (PreComp)')
axarr[0,1].plot(v['Blade Spanwise Coordinate'],v['K_66, EI_edge'],'ko-',mfc='None',mew=1,label='UCLA (VABS)')
axarr[0,1].set_xlabel('span [m]')
axarr[0,1].set_ylabel('edgewise stiffness [N*m^2]')
axarr[0,1].grid('on', axis='x')
axarr[0,1].legend()

# axial stiffness
twin_axis10 = prep_rel_diff_plot(axarr[1,0])
twin_axis10.plot(
    s['Blade Spanwise Coordinate'], rel_diff(v['K_11, EA_axial'],s['EA_axial']),
    'm^:', mec='m', mfc='None', mew=1, label='difference')
axarr[1,0].plot(s['Blade Spanwise Coordinate'],s['EA_axial'],'gx--',mfc='None',mew=1,label='Sandia (PreComp)')
axarr[1,0].plot(v['Blade Spanwise Coordinate'],v['K_11, EA_axial'],'ko-',mfc='None',mew=1,label='UCLA (VABS)')
axarr[1,0].set_xlabel('span [m]')
axarr[1,0].set_ylabel('axial stiffness [N]')
axarr[1,0].legend()
axarr[1,0].grid('on', axis='x')

# torsional stiffness
twin_axis11 = prep_rel_diff_plot(axarr[1,1])
twin_axis11.plot(
    s['Blade Spanwise Coordinate'], rel_diff(v['K_44, GJ_twist'],s['GJ_twist']),
    'm^:', mec='m', mfc='None', mew=1, label='difference')
axarr[1,1].plot(s['Blade Spanwise Coordinate'],s['GJ_twist'],'gx--',mfc='None',mew=1,label='Sandia (PreComp)')
axarr[1,1].plot(v['Blade Spanwise Coordinate'],v['K_44, GJ_twist'],'ko-',mfc='None',mew=1,label='UCLA (VABS)')
axarr[1,1].set_xlabel('span [m]')
axarr[1,1].set_ylabel('torsional stiffness [N*m^2]')
axarr[1,1].legend()
axarr[1,1].grid('on', axis='x')

plt.tight_layout()
plt.savefig('sandia_blade/Sandia_vs_VABS_stiffness_props.png')

# mass properties -------------------------------------------------------------
f2, axarr2 = plt.subplots(2,2, figsize=(12,8))

# mass density
twin_axis2_10 = prep_rel_diff_plot(axarr2[1,0])
twin_axis2_10.plot(
    s['Blade Spanwise Coordinate'], rel_diff(v['M_11, mu_mass'],s['mu_mass']),
    'm^:', mec='m', mfc='None', mew=1, label='difference')
axarr2[1,0].plot(s['Blade Spanwise Coordinate'],s['mu_mass'],'gx--',mfc='None',mew=1,label='Sandia (PreComp)')
axarr2[1,0].plot(v['Blade Spanwise Coordinate'],v['M_11, mu_mass'],'ko-',mfc='None',mew=1,label='UCLA (VABS)')
axarr2[1,0].set_xlabel('span [m]')
axarr2[1,0].set_ylabel('mass [kg/m]')
axarr2[1,0].legend()
axarr2[1,0].grid('on', axis='x')

# flapwise mass moment of inertia
twin_axis2_00 = prep_rel_diff_plot(axarr2[0,0])
twin_axis2_00.plot(
    s['Blade Spanwise Coordinate'], rel_diff(v['M_55, i22_flap'],s['i22_flap']),
    'm^:', mec='m', mfc='None', mew=1, label='difference')
axarr2[0,0].plot(s['Blade Spanwise Coordinate'],s['i22_flap'],'gx--',mfc='None',mew=1,label='Sandia (PreComp)')
axarr2[0,0].plot(v['Blade Spanwise Coordinate'],v['M_55, i22_flap'],'ko-',mfc='None',mew=1,label='UCLA (VABS)')
axarr2[0,0].set_xlabel('span [m]')
axarr2[0,0].set_ylabel('flapwise mass moment of inertia [kg*m]')
axarr2[0,0].legend()
axarr2[0,0].grid('on', axis='x')

# edgewise mass moment of inertia
twin_axis2_01 = prep_rel_diff_plot(axarr2[0,1])
twin_axis2_01.plot(
    s['Blade Spanwise Coordinate'], rel_diff(v['M_66, i33_edge'],s['i33_edge']),
    'm^:', mec='m', mfc='None', mew=1, label='difference')
axarr2[0,1].plot(s['Blade Spanwise Coordinate'],s['i33_edge'],'gx--',mfc='None',mew=1,label='Sandia (PreComp)')
axarr2[0,1].plot(v['Blade Spanwise Coordinate'],v['M_66, i33_edge'],'ko-',mfc='None',mew=1,label='UCLA (VABS)')
axarr2[0,1].set_xlabel('span [m]')
axarr2[0,1].set_ylabel('edgewise mass moment of inertia [kg*m]')
axarr2[0,1].legend()
axarr2[0,1].grid('on', axis='x')

plt.tight_layout()
plt.savefig('sandia_blade/Sandia_vs_VABS_mass_props.png')

plt.show()
