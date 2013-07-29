"""A script to create the Sandia blade from its definition CSV file.

Author: Perry Roth-Johnson
Last updated: July 24, 2013

"""


import scripts.blade as bl

b = bl.Blade('Sandia blade SNL100-00', 'sandia_blade')
b.copy_all_airfoil_coords()

# make some plots of the chord and twist schedules
# b.plot_chord_schedule()
# b.plot_twist_schedule()

# create some airfoil plots
for station in b.list_of_stations:
    station.read_airfoil_coords()
    station.plot_airfoil_coords()