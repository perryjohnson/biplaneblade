"""A script to compare the Sandia and biplane blades' masses.

Author: Perry Roth-Johnson
Last updated: April 17, 2014

"""


import lib.blade as bl
reload(bl)
import lib.compare_blades as cb
reload(cb)


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
# cb.plot_mass_schedule(m,b1)
m.plot_percent_areas()
# m.plot_percent_masses()
b1.plot_percent_areas()
