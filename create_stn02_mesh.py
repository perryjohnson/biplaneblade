"""Write initial TrueGrid files for one Sandia blade station.

Usage
-----
start an IPython (qt)console with the pylab flag:
$ ipython qtconsole --pylab
or
$ ipython --pylab
Then, from the prompt, run this script:
|> %run create_stnXX_mesh
Once you are finished looking at the meshes, you can clean up extra files:
|> %run clean
(See the 'clean.py' script in this directory for details.)

Author: Perry Roth-Johnson
Last updated: March 13, 2014

"""


import matplotlib.pyplot as plt
import lib.blade as bl
import lib.poly_utils as pu
from shapely.geometry import Polygon


# SET THESE PARAMETERS -----------------
station_num = 2
cut_phase = True
# --------------------------------------
plt.close('all')

# load the Sandia blade
m = bl.MonoplaneBlade('Sandia blade SNL100-00', 'sandia_blade',
    rotate_airfoil_coords=False)

# pre-process the station dimensions
station = m.list_of_stations[station_num-1]
station.airfoil.create_polygon()
station.structure.create_all_layers()
station.structure.save_all_layer_edges()
station.structure.write_all_part_polygons()

# plot the parts
station.plot_parts(alternate_layers=False)

# access the structure for this station
st = station.structure

if cut_phase:
    # upper right -----------------------------------------------------------
    label = 'upper spar cap'

    # create the bounding polygon
    points = [
        (-0.75, 2.5),
        ( 0.75, 2.5),
        ( 0.75, 3.0),
        (-0.75, 3.0)
        ]
    bounding_polygon = Polygon(points)
    pu.plot_polygon(bounding_polygon, 'None', '#000000')

    # cut the new layer polygons
    pu.cut_plot_and_write_alt_layer(st.root_buildup, 'triax', label, 
        bounding_polygon)
    pu.cut_plot_and_write_alt_layer(st.external_surface, 'triax', label, 
        bounding_polygon)
    pu.cut_plot_and_write_alt_layer(st.external_surface, 'gelcoat', label, 
        bounding_polygon)
    pu.cut_plot_and_write_alt_layer(st.internal_surface_1, 'resin', label, 
        bounding_polygon)
    pu.cut_plot_and_write_alt_layer(st.internal_surface_1, 'triax', label, 
        bounding_polygon)
   

    # show the plot
    plt.show()

    # write the TrueGrid input file for mesh generation ---------------------
    st.write_truegrid_inputfile(interrupt_flag=True)
