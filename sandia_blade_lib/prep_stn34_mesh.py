"""Write initial TrueGrid files for one Sandia blade station.

Usage
-----
start an IPython (qt)console with the pylab flag:
$ ipython qtconsole --pylab
or
$ ipython --pylab
Then, from the prompt, run this script:
|> %run sandia_blade_lib/prep_stnXX_mesh.py
or
|> import sandia_blade_lib/prep_stnXX_mesh

Author: Perry Roth-Johnson
Last updated: April 10, 2014

"""


import matplotlib.pyplot as plt
import lib.blade as bl
import lib.poly_utils as pu
from shapely.geometry import Polygon


# SET THESE PARAMETERS -----------------
station_num = 34
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
station.plot_parts()

# access the structure for this station
st = station.structure

# upper --------------------------------------------------------------------
label = 'upper'

# create the bounding polygon
is1 = st.internal_surface_1.layer['resin']
points_u = [
    (is1.polygon.interiors[0].coords[0][0],0.014),   # InternalSurface1_resin.txt
    is1.polygon.interiors[0].coords[0],   # InternalSurface1_resin.txt
    is1.polygon.interiors[0].coords[50-40],   # InternalSurface1_resin.txt
    (is1.polygon.interiors[0].coords[50-40][0],0.014)   # InternalSurface1_resin.txt
    ]
bounding_polygon = Polygon(points_u)
pu.plot_polygon(bounding_polygon, 'None', '#000000')

# cut the new layer polygons
pu.cut_plot_and_write_alt_layer(st.external_surface, 'triax', label, 
    bounding_polygon)
pu.cut_plot_and_write_alt_layer(st.external_surface, 'gelcoat', label, 
    bounding_polygon)
pu.cut_plot_and_write_alt_layer(st.internal_surface_1, 'resin', label, 
    bounding_polygon)
pu.cut_plot_and_write_alt_layer(st.internal_surface_1, 'triax', label, 
    bounding_polygon)

# lower --------------------------------------------------------------------
label = 'lower'

# create the bounding polygon
points_l = [
    (is1.polygon.interiors[0].coords[0][0],-0.0063),   # InternalSurface1_resin.txt
    is1.polygon.interiors[0].coords[0],   # InternalSurface1_resin.txt
    is1.polygon.interiors[0].coords[50-40],   # InternalSurface1_resin.txt
    (is1.polygon.interiors[0].coords[50-40][0],-0.0063)   # InternalSurface1_resin.txt
    ]
bounding_polygon = Polygon(points_l)
pu.plot_polygon(bounding_polygon, 'None', '#000000')

# cut the new layer polygons
pu.cut_plot_and_write_alt_layer(st.external_surface, 'triax', label, 
    bounding_polygon)
pu.cut_plot_and_write_alt_layer(st.external_surface, 'gelcoat', label, 
    bounding_polygon)
pu.cut_plot_and_write_alt_layer(st.internal_surface_1, 'resin', label, 
    bounding_polygon)
pu.cut_plot_and_write_alt_layer(st.internal_surface_1, 'triax', label, 
    bounding_polygon)

# aft, upper 1 -------------------------------------------------------------
label = 'aft, upper 1'

# create the bounding polygon
is1t = st.internal_surface_1.layer['triax']
points_au1 = [
    points_u[-1],
    points_u[-2],
    is1t.polygon.interiors[0].coords[72-55],  # InternalSurface1_triax.txt
    is1t.polygon.exterior.coords[24-3],  # InternalSurface1_triax.txt
    (is1t.polygon.exterior.coords[24-3][0], 0.014) # InternalSurface1_triax.txt
    ]
bounding_polygon = Polygon(points_au1)
pu.plot_polygon(bounding_polygon, 'None', '#000000')

# cut the new layer polygons
pu.cut_plot_and_write_alt_layer(st.external_surface, 'triax', label, 
    bounding_polygon)
pu.cut_plot_and_write_alt_layer(st.external_surface, 'gelcoat', label, 
    bounding_polygon)
pu.cut_plot_and_write_alt_layer(st.internal_surface_1, 'resin', label, 
    bounding_polygon)
pu.cut_plot_and_write_alt_layer(st.internal_surface_1, 'triax', label, 
    bounding_polygon)

# aft, lower 1 -------------------------------------------------------------
label = 'aft, lower 1'

# create the bounding polygon
points_al1 = [
    (points_au1[0][0], -0.006),
    points_au1[1],
    points_au1[2],
    points_au1[3],
    (points_au1[3][0], -0.006)
    ]
bounding_polygon = Polygon(points_al1)
pu.plot_polygon(bounding_polygon, 'None', '#000000')

# cut the new layer polygons
pu.cut_plot_and_write_alt_layer(st.external_surface, 'triax', label, 
    bounding_polygon)
pu.cut_plot_and_write_alt_layer(st.external_surface, 'gelcoat', label, 
    bounding_polygon)
pu.cut_plot_and_write_alt_layer(st.internal_surface_1, 'resin', label, 
    bounding_polygon)
pu.cut_plot_and_write_alt_layer(st.internal_surface_1, 'triax', label, 
    bounding_polygon)

# aft, upper 2 ------------------------------------------------
label = 'aft, upper 2'

# create the bounding polygon
est = st.external_surface.layer['triax']
esg = st.external_surface.layer['gelcoat']
points_au2 = [
    points_au1[-1],
    points_au1[-2],
    est.polygon.exterior.coords[0],
    (esg.polygon.exterior.coords[0][0], 0.0001),
    (esg.polygon.exterior.coords[-2][0], 0.014)
    ]
bounding_polygon = Polygon(points_au2)
pu.plot_polygon(bounding_polygon, 'None', '#000000')

# cut the new layer polygons
pu.cut_plot_and_write_alt_layer(st.external_surface, 'triax', label, 
    bounding_polygon)
pu.cut_plot_and_write_alt_layer(st.external_surface, 'gelcoat', label, 
    bounding_polygon)

# aft, lower 2 ------------------------------------------------
label = 'aft, lower 2'

# create the bounding polygon
points_al2 = [
    (points_au2[1][0],-0.006),
    points_au2[1],
    points_au2[2],
    points_au2[3],
    (points_au2[4][0],-0.006)
    ]
bounding_polygon = Polygon(points_al2)
pu.plot_polygon(bounding_polygon, 'None', '#000000')

# cut the new layer polygons
pu.cut_plot_and_write_alt_layer(st.external_surface, 'triax', label, 
    bounding_polygon)
pu.cut_plot_and_write_alt_layer(st.external_surface, 'gelcoat', label, 
    bounding_polygon)

# fore, upper 1 -------------------------------------------------------------
label = 'fore, upper 1'

# create the bounding polygon
points_fu1 = [
    points_u[0],
    points_u[1],
    is1t.polygon.interiors[0].coords[0],  # InternalSurface1_triax.txt
    is1t.polygon.exterior.coords[0],  # InternalSurface1_triax.txt
    (is1t.polygon.exterior.coords[0][0], 0.014) # InternalSurface1_triax.txt
    ]
bounding_polygon = Polygon(points_fu1)
pu.plot_polygon(bounding_polygon, 'None', '#000000')

# cut the new layer polygons
pu.cut_plot_and_write_alt_layer(st.external_surface, 'triax', label, 
    bounding_polygon)
pu.cut_plot_and_write_alt_layer(st.external_surface, 'gelcoat', label, 
    bounding_polygon)
pu.cut_plot_and_write_alt_layer(st.internal_surface_1, 'resin', label, 
    bounding_polygon)
pu.cut_plot_and_write_alt_layer(st.internal_surface_1, 'triax', label, 
    bounding_polygon)

# fore, lower 1 -------------------------------------------------------------
label = 'fore, lower 1'

# create the bounding polygon
points_fl1 = [
    (points_fu1[1][0],-0.006),
    points_fu1[1],
    points_fu1[2],
    points_fu1[3],
    (points_fu1[3][0],-0.006)
    ]
bounding_polygon = Polygon(points_fl1)
pu.plot_polygon(bounding_polygon, 'None', '#000000')

# cut the new layer polygons
pu.cut_plot_and_write_alt_layer(st.external_surface, 'triax', label, 
    bounding_polygon)
pu.cut_plot_and_write_alt_layer(st.external_surface, 'gelcoat', label, 
    bounding_polygon)
pu.cut_plot_and_write_alt_layer(st.internal_surface_1, 'resin', label, 
    bounding_polygon)
pu.cut_plot_and_write_alt_layer(st.internal_surface_1, 'triax', label, 
    bounding_polygon)

# LE, upper -------------------------------------------------------------
label = 'LE, upper'

# create the bounding polygon
points_leu = [
    is1t.polygon.exterior.coords[0],  # InternalSurface1_triax.txt
    (is1t.polygon.exterior.coords[0][0], 0.014), # InternalSurface1_triax.txt
    (-0.04,0.014),
    (-0.04,0.0)
    ]
bounding_polygon = Polygon(points_leu)
pu.plot_polygon(bounding_polygon, 'None', '#000000')

# cut the new layer polygons
pu.cut_plot_and_write_alt_layer(st.external_surface, 'triax', label, 
    bounding_polygon)
pu.cut_plot_and_write_alt_layer(st.external_surface, 'gelcoat', label, 
    bounding_polygon)

# LE, lower -------------------------------------------------------------
label = 'LE, lower'

# create the bounding polygon
points_lel = [
    is1t.polygon.exterior.coords[0],  # InternalSurface1_triax.txt
    (is1t.polygon.exterior.coords[0][0], -0.006), # InternalSurface1_triax.txt
    (-0.04,-0.006),
    (-0.04,0.0)
    ]
bounding_polygon = Polygon(points_lel)
pu.plot_polygon(bounding_polygon, 'None', '#000000')

# cut the new layer polygons
pu.cut_plot_and_write_alt_layer(st.external_surface, 'triax', label, 
    bounding_polygon)
pu.cut_plot_and_write_alt_layer(st.external_surface, 'gelcoat', label, 
    bounding_polygon)



# show the plot
plt.show()

# write the TrueGrid input file for mesh generation ---------------------
st.write_truegrid_inputfile(
    interrupt_flag=True, 
    soft_warning=False)
