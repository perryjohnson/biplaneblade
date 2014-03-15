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
import lib.poly_utils as pu
from shapely.geometry import Polygon


# SET THESE PARAMETERS -----------------
station_num = 2
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

# lower spar cap -----------------------------------------------------------
label = 'lower spar cap'

# create the bounding polygon
points = [
    (-0.75,-3.0),
    ( 0.75,-3.0),
    ( 0.75,-2.5),
    (-0.75,-2.5)
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

# TE reinforcement ---------------------------------------------------------
label = 'TE reinforcement'

# create the bounding polygon
points = [
    (1.84700000, 1.96787973),
    (1.98565, 2.10803),
    (3.0, 2.0),
    (3.0,-2.0),
    (1.98565,-2.10803),
    (1.84700000,-1.96787973)
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

# LE -----------------------------------------------------------------------
label = 'LE'

# create the bounding polygon
points = [
    (-3.00,-3.0),
    (-0.75,-3.0),
    (-0.75,-2.59318185),
    (-0.75999982, -2.57954345),
    (-0.75999982,  2.57954345),
    (-0.75, 2.59318185),
    (-0.75, 3.0),
    (-3.00, 3.0)
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

# upper right --------------------------------------------------------------
label = 'upper right'
# create the bounding polygon
points = [
    (0.75, 2.8),
    (2.0,  2.8),
    (1.98565, 2.10803),
    (1.84700000, 1.96787973),
    (1.83771501, 1.96282449),
    (0.8, 2.0),
    (0.75999982,   2.57954345),
    (0.75000000,   2.59318185)
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

# lower right --------------------------------------------------------------
label = 'lower right'
# create the bounding polygon
points = [
    (0.75, -2.8),
    (2.0,  -2.8),
    (1.98565, -2.10803),
    (1.84700000, -1.96787973),
    (1.83771501, -1.96282449),
    (0.8, -2.0),
    (0.75999982,   -2.57954345),
    (0.75000000,   -2.59318185)
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

# upper middle corner, internal surface -------------------------------------
label = 'upper middle corner'
# create the bounding polygon
points = [
    (0.75, 2.59318185),
    (0.75999982,   2.57954345),
    (0.77, 2.5),
    (0.75, 2.5)
    ]
bounding_polygon = Polygon(points)
pu.plot_polygon(bounding_polygon, 'None', '#000000')
# cut the new layer polygons
pu.cut_plot_and_write_alt_layer(st.internal_surface_1, 'resin', label, 
    bounding_polygon)
pu.cut_plot_and_write_alt_layer(st.internal_surface_1, 'triax', label, 
    bounding_polygon)

# lower middle corner, internal surface -------------------------------------
label = 'lower middle corner'
# create the bounding polygon
points = [
    (0.75, -2.59318185),
    (0.75999982,   -2.57954345),
    (0.77, -2.5),
    (0.75, -2.5)
    ]
bounding_polygon = Polygon(points)
pu.plot_polygon(bounding_polygon, 'None', '#000000')
# cut the new layer polygons
pu.cut_plot_and_write_alt_layer(st.internal_surface_1, 'resin', label, 
    bounding_polygon)
pu.cut_plot_and_write_alt_layer(st.internal_surface_1, 'triax', label, 
    bounding_polygon)

# upper left corner, internal surface -------------------------------------
label = 'upper left corner'
# create the bounding polygon
points = [
    (-0.75000000,   2.59318185),
    (-0.75999982,  2.57954345),
    (-0.77, 2.56),
    (-0.75, 2.56)
    ]
bounding_polygon = Polygon(points)
pu.plot_polygon(bounding_polygon, 'None', '#000000')
# cut the new layer polygons
pu.cut_plot_and_write_alt_layer(st.internal_surface_1, 'resin', label, 
    bounding_polygon)
pu.cut_plot_and_write_alt_layer(st.internal_surface_1, 'triax', label, 
    bounding_polygon)

# lower left corner, internal surface -------------------------------------
label = 'lower left corner'
# create the bounding polygon
points = [
    (-0.75000000,   -2.59318185),
    (-0.75999982,  -2.57954345),
    (-0.77, -2.56),
    (-0.75, -2.56)
    ]
bounding_polygon = Polygon(points)
pu.plot_polygon(bounding_polygon, 'None', '#000000')
# cut the new layer polygons
pu.cut_plot_and_write_alt_layer(st.internal_surface_1, 'resin', label, 
    bounding_polygon)
pu.cut_plot_and_write_alt_layer(st.internal_surface_1, 'triax', label, 
    bounding_polygon)

# upper right corner, internal surface -------------------------------------
label = 'upper right corner'
# create the bounding polygon
points = [
    (1.84700000,   1.96787973),
    (1.84700000,   1.95),
    (1.835, 1.95),
    (1.83771501, 1.96282449)
    ]
bounding_polygon = Polygon(points)
pu.plot_polygon(bounding_polygon, 'None', '#000000')
# cut the new layer polygons
pu.cut_plot_and_write_alt_layer(st.internal_surface_1, 'resin', label, 
    bounding_polygon)
pu.cut_plot_and_write_alt_layer(st.internal_surface_1, 'triax', label, 
    bounding_polygon)

# lower right corner, internal surface -------------------------------------
label = 'lower right corner'
# create the bounding polygon
points = [
    (1.84700000,   -1.96787973),
    (1.84700000,   -1.95),
    (1.835, -1.95),
    (1.83771501, -1.96282449)
    ]
bounding_polygon = Polygon(points)
pu.plot_polygon(bounding_polygon, 'None', '#000000')
# cut the new layer polygons
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
        st.TE_reinforcement.layer['uniax']
    ])
