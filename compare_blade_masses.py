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
    blade2_stn_nums=[10,24], print_flag=True)
m.plot_percent_masses()
b1.plot_percent_masses()
# m.plot_selected_cross_sections()
# b1.plot_selected_cross_sections(selected_stations=[1,7,11,13,16,22,25,29,32,34,36,39])
# b1.plot_selected_cross_sections(nrows=2, ncols=3, selected_stations=[10,13,17,20,22,24])
fig1, ax1 = plt.subplots()
m.list_of_stations[10-1].plot_parts(ax1)
fig2, ax2 = plt.subplots()
b1.list_of_stations[10-1].plot_parts(ax2)
