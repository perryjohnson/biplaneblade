"""Write initial TrueGrid files for one Sandia blade station.

Usage
-----
start an IPython (qt)console with the pylab flag:
$ ipython qtconsole --pylab
or
$ ipython --pylab
Then, from the prompt, run this script:
|> %run create_stnXX_mesh

Author: Perry Roth-Johnson
Last updated: March 15, 2014

"""


import matplotlib.pyplot as plt
import lib.blade as bl
reload(bl)
import lib.poly_utils as pu
reload(pu)
from shapely.geometry import Polygon


# SET THESE PARAMETERS -----------------
station_num = 6
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

# upper spar cap -----------------------------------------------------------
label = 'upper spar cap'

# create the bounding polygon
points_usc = [
    (-0.75, 2.5),
    ( 0.75, 2.5),
    ( 0.75, 3.0),
    (-0.75, 3.0)
    ]
bounding_polygon = Polygon(points_usc)
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

# lower spar cap -----------------------------------------------------------
label = 'lower spar cap'

# create the bounding polygon
points_lsc = [
    (-0.75,-3.0),
    ( 0.75,-3.0),
    ( 0.75,-2.5),
    (-0.75,-2.5)
    ]
bounding_polygon = Polygon(points_lsc)
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

# TE reinforcement ---------------------------------------------------------
label = 'TE reinforcement'

# create the bounding polygon
points_te = [
    (1.84700000,   1.99925983),
    (1.95, 2.1),
    (3.0, 2.0),
    (3.0,-2.0),
    (1.95,-2.1),
    (1.84700000,  -1.99925983)
    ]
bounding_polygon = Polygon(points_te)
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

# LE panel -----------------------------------------------------------------
label = 'LE panel'

# create the bounding polygon
is1 = st.internal_surface_1.layer['resin']
points_le = [
    (-3.00,-3.0),
    (-0.75,-3.0),
    (-0.75000000, -2.58336990),
    # (-0.75994863,  -2.56986995),
    is1.polygon.interiors[0].coords[-1],
    is1.polygon.interiors[0].coords[-31],
    # (-0.75994863,   2.56986995),
    (-0.75000000,  2.58336990),
    (-0.75, 3.0),
    (-3.00, 3.0)
    ]
bounding_polygon = Polygon(points_le)
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

# upper aft panel 1 -------------------------------------------------------
label = 'upper aft panel 1'
# create the bounding polygon
points_ur = [
    (0.75, 2.8),
    (2.0,  2.8),
    (1.95, 2.1),
    points_te[0],
    (1.84700000,   1.99792189),
    (1.83700000,   1.99343115),
    (0.8, 2.0),
    (0.75949119,   2.56768163),
    (0.75000000,   2.58336990)
    ]
bounding_polygon = Polygon(points_ur)
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

# lower aft panel 1 -------------------------------------------------------
label = 'lower aft panel 1'
# create the bounding polygon
points_lr = [
    (0.75, -2.8),
    (2.0,  -2.8),
    (1.95, -2.1),
    points_te[-1],
    (1.84700000,   -points_ur[-5][1]),
    (1.83700000,   -points_ur[-4][1]),
    (0.8, -2.0),
    (0.75949119,  -2.56768163),
    (0.75000000,   -points_ur[-1][1])
    ]
bounding_polygon = Polygon(points_lr)
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
st.write_truegrid_inputfile(
    interrupt_flag=True, 
    additional_layers=[
        st.spar_cap.layer['upper'],
        st.spar_cap.layer['lower'],
        st.TE_reinforcement.layer['uniax'],
        st.aft_panel_1.layer['upper'],
        st.aft_panel_1.layer['lower'],
        st.LE_panel.layer['foam']
    ])
