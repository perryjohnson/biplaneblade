import os
import matplotlib.pyplot as plt
import numpy as np
import lib.blade as bl
reload(bl)

# load the sandia blade
m = bl.MonoplaneBlade('Sandia blade SNL100-00', 'sandia_blade',
    rotate_airfoil_coords=False)

# pre-process the airfoil coordinates
for station in m.list_of_stations:
    station.airfoil.create_polygon()
    station.structure.create_all_layers()
    station.structure.save_all_layer_edges()
    # station.plot_layer_edges()
    # station.structure.write_all_layer_edges()

# stn6 = m.list_of_stations[5]
# stn6.plot_layer_edges()
# stn6.plot_parts()

stn16 = m.list_of_stations[15]
# stn16.plot_layer_edges()
# stn16.plot_parts()

d = stn16.structure.write_all_layer_edges()

# f = open(os.path.join(stn16.station_path, 'curves.tg'), 'w')
# # write the curves of the LE panel into a file
# # left curve
# curve_num = 1
# f.write('curd {0} lp3\n'.format(curve_num))
# for coord_pair in stn16.structure.LE_panel.layer['foam'].left:
#     f.write('{0: .8f}  {1: .8f}  0.0\n'.format(coord_pair[0], coord_pair[1]))
# f.write(';;\n\n')
# # bottom curve
# curve_num += 1
# f.write('curd {0} lp3\n'.format(curve_num))
# for coord_pair in stn16.structure.LE_panel.layer['foam'].bottom:
#     f.write('{0: .8f}  {1: .8f}  0.0\n'.format(coord_pair[0], coord_pair[1]))
# f.write(';;\n\n')
# # right curve
# curve_num += 1
# f.write('curd {0} lp3\n'.format(curve_num))
# for coord_pair in stn16.structure.LE_panel.layer['foam'].right:
#     f.write('{0: .8f}  {1: .8f}  0.0\n'.format(coord_pair[0], coord_pair[1]))
# f.write(';;\n\n')
# # top curve
# curve_num += 1
# f.write('curd {0} lp3\n'.format(curve_num))
# for coord_pair in stn16.structure.LE_panel.layer['foam'].top:
#     f.write('{0: .8f}  {1: .8f}  0.0\n'.format(coord_pair[0], coord_pair[1]))
# f.write(';;\n\n')

# # write the curves of shear web 1 into the file
# # left biax layer -----
# # left curve
# curve_num += 1
# f.write('curd {0} lp3\n'.format(curve_num))
# for coord_pair in stn16.structure.shear_web_1.layer['biax, left'].left:
#     f.write('{0: .8f}  {1: .8f}  0.0\n'.format(coord_pair[0], coord_pair[1]))
# f.write(';;\n\n')
# # bottom curve
# curve_num += 1
# f.write('curd {0} lp3\n'.format(curve_num))
# for coord_pair in stn16.structure.shear_web_1.layer['biax, left'].bottom:
#     f.write('{0: .8f}  {1: .8f}  0.0\n'.format(coord_pair[0], coord_pair[1]))
# f.write(';;\n\n')
# # right curve
# curve_num += 1
# f.write('curd {0} lp3\n'.format(curve_num))
# for coord_pair in stn16.structure.shear_web_1.layer['biax, left'].right:
#     f.write('{0: .8f}  {1: .8f}  0.0\n'.format(coord_pair[0], coord_pair[1]))
# f.write(';;\n\n')
# # top curve
# curve_num += 1
# f.write('curd {0} lp3\n'.format(curve_num))
# for coord_pair in stn16.structure.shear_web_1.layer['biax, left'].top:
#     f.write('{0: .8f}  {1: .8f}  0.0\n'.format(coord_pair[0], coord_pair[1]))
# f.write(';;\n\n')
# # foam layer -----
# # left curve
# curve_num += 1
# f.write('curd {0} lp3\n'.format(curve_num))
# for coord_pair in stn16.structure.shear_web_1.layer['foam'].left:
#     f.write('{0: .8f}  {1: .8f}  0.0\n'.format(coord_pair[0], coord_pair[1]))
# f.write(';;\n\n')
# # bottom curve
# curve_num += 1
# f.write('curd {0} lp3\n'.format(curve_num))
# for coord_pair in stn16.structure.shear_web_1.layer['foam'].bottom:
#     f.write('{0: .8f}  {1: .8f}  0.0\n'.format(coord_pair[0], coord_pair[1]))
# f.write(';;\n\n')
# # right curve
# curve_num += 1
# f.write('curd {0} lp3\n'.format(curve_num))
# for coord_pair in stn16.structure.shear_web_1.layer['foam'].right:
#     f.write('{0: .8f}  {1: .8f}  0.0\n'.format(coord_pair[0], coord_pair[1]))
# f.write(';;\n\n')
# # top curve
# curve_num += 1
# f.write('curd {0} lp3\n'.format(curve_num))
# for coord_pair in stn16.structure.shear_web_1.layer['foam'].top:
#     f.write('{0: .8f}  {1: .8f}  0.0\n'.format(coord_pair[0], coord_pair[1]))
# f.write(';;\n\n')
# # right biax layer -----
# # left curve
# curve_num += 1
# f.write('curd {0} lp3\n'.format(curve_num))
# for coord_pair in stn16.structure.shear_web_1.layer['biax, right'].left:
#     f.write('{0: .8f}  {1: .8f}  0.0\n'.format(coord_pair[0], coord_pair[1]))
# f.write(';;\n\n')
# # bottom curve
# curve_num += 1
# f.write('curd {0} lp3\n'.format(curve_num))
# for coord_pair in stn16.structure.shear_web_1.layer['biax, right'].bottom:
#     f.write('{0: .8f}  {1: .8f}  0.0\n'.format(coord_pair[0], coord_pair[1]))
# f.write(';;\n\n')
# # right curve
# curve_num += 1
# f.write('curd {0} lp3\n'.format(curve_num))
# for coord_pair in stn16.structure.shear_web_1.layer['biax, right'].right:
#     f.write('{0: .8f}  {1: .8f}  0.0\n'.format(coord_pair[0], coord_pair[1]))
# f.write(';;\n\n')
# # top curve
# curve_num += 1
# f.write('curd {0} lp3\n'.format(curve_num))
# for coord_pair in stn16.structure.shear_web_1.layer['biax, right'].top:
#     f.write('{0: .8f}  {1: .8f}  0.0\n'.format(coord_pair[0], coord_pair[1]))
# f.write(';;\n\n')



# f.close()

# m.plot_selected_cross_sections()