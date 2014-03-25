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
Last updated: March 24, 2014

"""


import matplotlib.pyplot as plt
import lib.blade as bl
reload(bl)
import lib.poly_utils as pu
reload(pu)
from shapely.geometry import Polygon


# SET THESE PARAMETERS -----------------
station_num = 31
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
usc = st.spar_cap.layer['upper']
is1 = st.internal_surface_1.layer['resin']
points_usc = [
    tuple(usc.left[0]),              # SparCap_upper.txt
    (usc.left[0][0], 0.1),
    is1.polygon.interiors[0].coords[519-263],   # InternalSurface1_resin.txt
    tuple(usc.right[1]),             # SparCap_upper.txt
    (usc.right[1][0], 0.3),
    (usc.left[0][0], 0.3)
    ]
bounding_polygon = Polygon(points_usc)
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

# lower spar cap -----------------------------------------------------------
label = 'lower spar cap'

# create the bounding polygon
lsc = st.spar_cap.layer['lower']
points_lsc = [
    tuple(lsc.left[1]),
    (lsc.left[1][0], 0.0),
    is1.polygon.interiors[0].coords[313-263],   # InternalSurface1_resin.txt
    tuple(lsc.right[0]),        # SparCap_lower.txt
    (lsc.right[0][0], -0.15),
    (lsc.left[1][0], -0.15)
    ]
bounding_polygon = Polygon(points_lsc)
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

# TE reinforcement, upper 1 ------------------------------------------------
label = 'TE reinforcement, upper 1'

# create the bounding polygon
ter = st.TE_reinforcement.layer['foam']
points_teu1 = [
    (ter.top[0][0], 0.3),              # TE_Reinforcement_foam.txt
    tuple(ter.top[0]),                  # TE_Reinforcement_foam.txt
    (0.6, 0.15),
    is1.polygon.interiors[0].coords[478-263],  # InternalSurface1_resin.txt
    (is1.polygon.interiors[0].coords[478-263][0],   0.3)  # InternalSurface1_resin.txt
    ]
bounding_polygon = Polygon(points_teu1)
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
pu.cut_plot_and_write_alt_layer(st.TE_reinforcement, 'foam', label,
    bounding_polygon)
pu.cut_plot_and_write_alt_layer(st.TE_reinforcement, 'uniax', label, 
    bounding_polygon)

# TE reinforcement, lower 1 ------------------------------------------------
label = 'TE reinforcement, lower 1'

# create the bounding polygon
points_tel1 = [
    (ter.bottom[0][0], -0.15),              # TE_Reinforcement_foam.txt
    tuple(ter.bottom[1]),                  # TE_Reinforcement_foam.txt
    (0.6, -0.02),
    (1.0, 0.05),
    points_teu1[-2],         # InternalSurface1_resin.txt
    (points_teu1[-1][0],   -0.15)                # InternalSurface1_resin.txt
    ]
bounding_polygon = Polygon(points_tel1)
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
pu.cut_plot_and_write_alt_layer(st.TE_reinforcement, 'foam', label,
    bounding_polygon)
pu.cut_plot_and_write_alt_layer(st.TE_reinforcement, 'uniax', label, 
    bounding_polygon)

# TE reinforcement, upper 2 ------------------------------------------------
label = 'TE reinforcement, upper 2'

# create the bounding polygon
is1t = st.internal_surface_1.layer['triax']
points_teu2 = [
    points_teu1[-1],
    points_teu1[-2],
    is1t.polygon.interiors[0].coords[376-176],  # InternalSurface1_triax.txt
    is1t.polygon.exterior.coords[44-3],  # InternalSurface1_triax.txt
    (is1t.polygon.exterior.coords[44-3][0],   0.3) # InternalSurface1_triax.txt
    ]
bounding_polygon = Polygon(points_teu2)
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
pu.cut_plot_and_write_alt_layer(st.TE_reinforcement, 'foam', label,
    bounding_polygon)
pu.cut_plot_and_write_alt_layer(st.TE_reinforcement, 'uniax', label, 
    bounding_polygon)

# TE reinforcement, lower 2 ------------------------------------------------
label = 'TE reinforcement, lower 2'

# create the bounding polygon
points_tel2 = [
    (points_teu2[0][0], -0.1),
    points_teu2[1],
    points_teu2[2],
    points_teu2[3],
    (points_teu2[3][0], -0.1)
    ]
bounding_polygon = Polygon(points_tel2)
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
pu.cut_plot_and_write_alt_layer(st.TE_reinforcement, 'foam', label,
    bounding_polygon)
pu.cut_plot_and_write_alt_layer(st.TE_reinforcement, 'uniax', label, 
    bounding_polygon)

# TE reinforcement, upper 3 ------------------------------------------------
label = 'TE reinforcement, upper 3'

# create the bounding polygon
teru = st.TE_reinforcement.layer['uniax']
est = st.external_surface.layer['triax']
esg = st.external_surface.layer['gelcoat']
points_teu3 = [
    points_teu2[-1],
    points_teu2[-2],
    ter.polygon.exterior.coords[0],
    teru.polygon.exterior.coords[0],
    (est.polygon.exterior.coords[-1][0], 0.0025),
    est.polygon.exterior.coords[-2],
    esg.polygon.exterior.coords[-2],
    (esg.polygon.exterior.coords[-2][0], 0.3)
    ]
bounding_polygon = Polygon(points_teu3)
pu.plot_polygon(bounding_polygon, 'None', '#000000')

# cut the new layer polygons
pu.cut_plot_and_write_alt_layer(st.external_surface, 'triax', label, 
    bounding_polygon)
pu.cut_plot_and_write_alt_layer(st.external_surface, 'gelcoat', label, 
    bounding_polygon)
pu.cut_plot_and_write_alt_layer(st.TE_reinforcement, 'foam', label,
    bounding_polygon)
pu.cut_plot_and_write_alt_layer(st.TE_reinforcement, 'uniax', label, 
    bounding_polygon)

# TE reinforcement, lower 3 ------------------------------------------------
label = 'TE reinforcement, lower 3'

# create the bounding polygon
points_tel3 = [
    (points_teu3[0][0], -0.1),
    points_teu3[1],
    points_teu3[2],
    points_teu3[3],
    points_teu3[4],
    est.polygon.exterior.coords[-1],
    esg.polygon.exterior.coords[-1],
    (points_teu3[4][0], -0.1)
    ]
bounding_polygon = Polygon(points_tel3)
pu.plot_polygon(bounding_polygon, 'None', '#000000')

# cut the new layer polygons
pu.cut_plot_and_write_alt_layer(st.external_surface, 'triax', label, 
    bounding_polygon)
pu.cut_plot_and_write_alt_layer(st.external_surface, 'gelcoat', label, 
    bounding_polygon)
pu.cut_plot_and_write_alt_layer(st.TE_reinforcement, 'foam', label,
    bounding_polygon)
pu.cut_plot_and_write_alt_layer(st.TE_reinforcement, 'uniax', label, 
    bounding_polygon)

# LE panel -----------------------------------------------------------------
label = 'LE panel'

# create the bounding polygon
lep = st.LE_panel.layer['foam']
is1 = st.internal_surface_1.layer['resin']
points_le = [
    (-1.0,-0.15),
    (lep.bottom[0][0],-0.15),
    (lep.bottom[0][0],0.3),
    (-1.0, 0.3)
    ]
bounding_polygon = Polygon(points_le)
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




# show the plot
plt.show()

# write the TrueGrid input file for mesh generation ---------------------
st.write_truegrid_inputfile(
    interrupt_flag=True, 
    additional_layers=[
        st.spar_cap.layer['upper'],
        st.spar_cap.layer['lower'],
        st.LE_panel.layer['foam']
    ],
    alt_TE_reinforcement=True,
    soft_warning=False)
