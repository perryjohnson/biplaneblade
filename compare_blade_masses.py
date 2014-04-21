"""A script to compare the Sandia and biplane blades' masses.

Author: Perry Roth-Johnson
Last updated: April 17, 2014

"""


import lib.blade as bl
reload(bl)
import lib.compare_blades as cb
reload(cb)
import matplotlib.pyplot as plt


biplane_flap_sym_no_stagger_flag = True
sandia_flag = True

# --- biplane blade, flapwise symmetric, no stagger----------------------------
if biplane_flap_sym_no_stagger_flag:
    b1 = bl.BiplaneBlade(
        'biplane blade, flapwise symmetric, no stagger, rj/R=0.452, g/c=1.25',
        'biplane_blade')
    b1.copy_all_airfoil_coords()

    # pre-process the airfoil coordinates
    for station in b1.list_of_stations:
        station.airfoil.create_polygon()
        station.structure.create_all_layers()
        station.structure.write_all_part_polygons()

# --- sandia blade ------------------------------------------------------------
if sandia_flag:
    m = bl.MonoplaneBlade('Sandia blade SNL100-00', 'sandia_blade')

    # pre-process the airfoil coordinates
    for station in m.list_of_stations:
        station.airfoil.create_polygon()
        station.structure.create_all_layers()
        station.structure.write_all_part_polygons()

# compare blade masses ----------------------------
plt.close('all')
cb.plot_mass_schedule(m, b1, show_stn_nums=True, blade1_stn_nums=[10,20],
    blade2_stn_nums=[10,24])
print 'stn   mass mono   mass bi   % diff'
print '---   ---------   -------   ------'
for stn in range(9,19):
    m_stn = m.list_of_stations[stn]
    b_stn = b1.list_of_stations[stn]
    pd = (b_stn.structure.mass - m_stn.structure.mass)/(m_stn.structure.mass)*100
    print '{0:3}   {1:9.0f}   {2:7.0f}   {3: 6.2f}'.format(
        m_stn.station_num, m_stn.structure.mass, b_stn.structure.mass, pd)
m.plot_percent_masses()
b1.plot_percent_masses()
# fig1, ax1 = plt.subplots()
# stn_to_plot = 12
# m.list_of_stations[stn_to_plot-1].plot_parts(ax1)
# fig2, ax2 = plt.subplots()
# b1.list_of_stations[stn_to_plot-1].plot_parts(ax2)
