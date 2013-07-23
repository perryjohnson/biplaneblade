import scripts.blade as bld
reload(bld)

b = bld.Blade('Sandia blade SNL100-00', 'sandia_blade')
b.copy_all_airfoil_coords()

# make some plots of the chord and twist schedules
# b.plot_chord_schedule()
# b.plot_twist_schedule()
