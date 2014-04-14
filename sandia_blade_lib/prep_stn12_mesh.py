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
station_num = 12
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

# upper spar cap -----------------------------------------------------------
label = 'upper spar cap'

# create the bounding polygon
usc = st.spar_cap.layer['upper']
is2 = st.internal_surface_2.layer['resin']
points_usc = [
    (-0.75, usc.left[0][1]),  # SparCap_upper.txt
    is2.polygon.interiors[0].coords[0],  # InternalSurface2_resin.txt
    ( 0.74, 1.91931038),  # InternalSurface2_resin.txt
    ( 0.75, usc.right[1][1]),  # SparCap_upper.txt
    ( 0.75, 2.1),
    (-0.75, 2.1)
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
pu.cut_plot_and_write_alt_layer(st.internal_surface_2, 'resin', label, 
    bounding_polygon)
pu.cut_plot_and_write_alt_layer(st.internal_surface_2, 'triax', label, 
    bounding_polygon)

# lower spar cap -----------------------------------------------------------
label = 'lower spar cap'

# create the bounding polygon
lsc = st.spar_cap.layer['lower']
points_lsc = [
    (-0.75,-2.1),
    ( 0.75,-2.1),
    ( 0.75000000,  lsc.right[0][1]),  # SparCap_lower.txt
    ( 0.74000000,  -1.92424072),  # InternalSurface2_resin.txt
    is2.polygon.interiors[0].coords[1],  # InternalSurface2_resin.txt
    (-0.75000000,  lsc.left[1][1])   # SparCap_lower.txt
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
pu.cut_plot_and_write_alt_layer(st.internal_surface_2, 'resin', label, 
    bounding_polygon)
pu.cut_plot_and_write_alt_layer(st.internal_surface_2, 'triax', label, 
    bounding_polygon)

# TE reinforcement ---------------------------------------------------------
label = 'TE reinforcement'

# create the bounding polygon
ter = st.TE_reinforcement.layer['foam']
points_te = [
    (4.0,-2.0),
    (ter.bottom[1][0], -2.0),          # TE_Reinforcement_foam.txt
    tuple(ter.bottom[1]),  # TE_Reinforcement_foam.txt
    (2.86832393,  -1.16833178),  # InternalSurface3_resin.txt
    (3.0, 0.0),
    (2.86824864,   1.20797337),  # InternalSurface3_resin.txt
    tuple(ter.top[0]),  # TE_Reinforcement_foam.txt
    (ter.top[0][0], 2.0),           # TE_Reinforcement_foam.txt
    (4.0, 2.0)
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
pu.cut_plot_and_write_alt_layer(st.internal_surface_3, 'resin', label, 
    bounding_polygon)
pu.cut_plot_and_write_alt_layer(st.internal_surface_3, 'triax', label, 
    bounding_polygon)

# LE panel -----------------------------------------------------------------
label = 'LE panel'

# create the bounding polygon
lep = st.LE_panel.layer['foam']
is1 = st.internal_surface_1.layer['resin']
points_le = [
    (-3.00,-2.0),
    (-0.836,-2.0),
    tuple(lep.bottom[0]),  # LE_Panel_foam.txt
    is1.polygon.interiors[0].coords[0],  # InternalSurface1_resin.txt
    (-1.5, 0.0),
    is1.polygon.interiors[0].coords[1],  # InternalSurface1_resin.txt
    tuple(lep.top[1]),  # LE_Panel_foam.txt
    (-0.836, 2.0),
    (-3.00, 2.0)
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
ap1u = st.aft_panel_1.layer['upper']
is3 = st.internal_surface_3.layer['resin']
points_ur = [
    (0.836, 2.1),
    (ap1u.right[1][0], 2.1),           # AftPanel1_upper.txt
    tuple(ap1u.right[1]),  # AftPanel1_upper.txt
    (2.865, 1.18),                # by eye
    (1.2, 1.0),
    is3.polygon.interiors[0].coords[-2],       # InternalSurface3_resin.txt
    tuple(ap1u.left[0])        # AftPanel1_upper.txt
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
pu.cut_plot_and_write_alt_layer(st.internal_surface_3, 'resin', label, 
    bounding_polygon)
pu.cut_plot_and_write_alt_layer(st.internal_surface_3, 'triax', label, 
    bounding_polygon)

# lower aft panel 1 -------------------------------------------------------
label = 'lower aft panel 1'

# create the bounding polygon
ap1l = st.aft_panel_1.layer['lower']
points_lr = [
    (0.836, -2.1),
    (ap1l.right[0][0], -2.1),          # AftPanel1_lower.txt
    tuple(ap1l.right[0]),  # AftPanel1_lower.txt
    (2.865, -1.15),               # by eye
    (1.2, -1.0),
    is3.polygon.interiors[0].coords[-1],       # InternalSurface3_resin.txt
    tuple(ap1l.left[1])        # AftPanel1_lower.txt
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
pu.cut_plot_and_write_alt_layer(st.internal_surface_3, 'resin', label, 
    bounding_polygon)
pu.cut_plot_and_write_alt_layer(st.internal_surface_3, 'triax', label, 
    bounding_polygon)

# above shear web 1 ----------------------------------------------------------
label = 'above shear web 1'

# create the bounding polygon
points_asw1 = [
    (-0.75, 2.1),
    (-0.75, 1.0),
    (-0.836, 1.0),
    (-0.836, 2.1)
    ]
bounding_polygon = Polygon(points_asw1)
pu.plot_polygon(bounding_polygon, 'None', '#000000')

# cut the new layer polygons
pu.cut_plot_and_write_alt_layer(st.root_buildup, 'triax', label, 
    bounding_polygon)
pu.cut_plot_and_write_alt_layer(st.external_surface, 'triax', label, 
    bounding_polygon)
pu.cut_plot_and_write_alt_layer(st.external_surface, 'gelcoat', label, 
    bounding_polygon)

# below shear web 1 ----------------------------------------------------------
label = 'below shear web 1'

# create the bounding polygon
points_bsw1 = [
    (-0.75, -2.1),
    (-0.75, -1.0),
    (-0.836, -1.0),
    (-0.836, -2.1)
    ]
bounding_polygon = Polygon(points_bsw1)
pu.plot_polygon(bounding_polygon, 'None', '#000000')

# cut the new layer polygons
pu.cut_plot_and_write_alt_layer(st.root_buildup, 'triax', label, 
    bounding_polygon)
pu.cut_plot_and_write_alt_layer(st.external_surface, 'triax', label, 
    bounding_polygon)
pu.cut_plot_and_write_alt_layer(st.external_surface, 'gelcoat', label, 
    bounding_polygon)

# above shear web 2 ----------------------------------------------------------
label = 'above shear web 2'

# create the bounding polygon
points_asw2 = [
    (0.75, 2.1),
    (0.75, 1.0),
    (0.836, 1.0),
    (0.836, 2.1)
    ]
bounding_polygon = Polygon(points_asw2)
pu.plot_polygon(bounding_polygon, 'None', '#000000')

# cut the new layer polygons
pu.cut_plot_and_write_alt_layer(st.root_buildup, 'triax', label, 
    bounding_polygon)
pu.cut_plot_and_write_alt_layer(st.external_surface, 'triax', label, 
    bounding_polygon)
pu.cut_plot_and_write_alt_layer(st.external_surface, 'gelcoat', label, 
    bounding_polygon)

# below shear web 2 ----------------------------------------------------------
label = 'below shear web 2'

# create the bounding polygon
points_bsw2 = [
    (0.75, -2.1),
    (0.75, -1.0),
    (0.836, -1.0),
    (0.836, -2.1)
    ]
bounding_polygon = Polygon(points_bsw2)
pu.plot_polygon(bounding_polygon, 'None', '#000000')

# cut the new layer polygons
pu.cut_plot_and_write_alt_layer(st.root_buildup, 'triax', label, 
    bounding_polygon)
pu.cut_plot_and_write_alt_layer(st.external_surface, 'triax', label, 
    bounding_polygon)
pu.cut_plot_and_write_alt_layer(st.external_surface, 'gelcoat', label, 
    bounding_polygon)

# left of shear web 1 -------------------------------------------------------
label = 'left of shear web 1'

# create the bounding polygon
points_lsw1 = points_le[2:-2]
bounding_polygon = Polygon(points_lsw1)
pu.plot_polygon(bounding_polygon, 'None', '#000000')
# cut the new layer polygons
pu.cut_plot_and_write_alt_layer(st.internal_surface_1, 'resin', label, 
    bounding_polygon)
pu.cut_plot_and_write_alt_layer(st.internal_surface_1, 'triax', label, 
    bounding_polygon)

# right of shear web 1 -------------------------------------------------------
label = 'right of shear web 1'

# create the bounding polygon
points_rsw1 = [
    points_usc[0],
    points_usc[1],
    (0.0, 0.0),
    points_lsc[-2],
    points_lsc[-1]
    ]
bounding_polygon = Polygon(points_rsw1)
pu.plot_polygon(bounding_polygon, 'None', '#000000')
# cut the new layer polygons
pu.cut_plot_and_write_alt_layer(st.internal_surface_2, 'resin', label, 
    bounding_polygon)
pu.cut_plot_and_write_alt_layer(st.internal_surface_2, 'triax', label, 
    bounding_polygon)

# left of shear web 2 -------------------------------------------------------
label = 'left of shear web 2'

# create the bounding polygon
points_lsw2 = [
    points_usc[3],
    points_usc[2],
    (0.0, 0.0),
    points_lsc[3],
    points_lsc[2]
    ]
bounding_polygon = Polygon(points_lsw2)
pu.plot_polygon(bounding_polygon, 'None', '#000000')
# cut the new layer polygons
pu.cut_plot_and_write_alt_layer(st.internal_surface_2, 'resin', label, 
    bounding_polygon)
pu.cut_plot_and_write_alt_layer(st.internal_surface_2, 'triax', label, 
    bounding_polygon)

# right of shear web 2 -------------------------------------------------------
label = 'right of shear web 2'

# create the bounding polygon
points_rsw2 = [
    points_ur[-1],
    points_ur[-2],
    (1.5, 0.0),
    points_lr[-2],
    points_lr[-1]
    ]
bounding_polygon = Polygon(points_rsw2)
pu.plot_polygon(bounding_polygon, 'None', '#000000')
# cut the new layer polygons
pu.cut_plot_and_write_alt_layer(st.internal_surface_3, 'resin', label, 
    bounding_polygon)
pu.cut_plot_and_write_alt_layer(st.internal_surface_3, 'triax', label, 
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
        st.TE_reinforcement.layer['foam'],
        st.aft_panel_1.layer['upper'],
        st.aft_panel_1.layer['lower'],
        st.LE_panel.layer['foam'],
        st.shear_web_1.layer['biax, left'],
        st.shear_web_1.layer['foam'],
        st.shear_web_1.layer['biax, right'],
        st.shear_web_2.layer['biax, left'],
        st.shear_web_2.layer['foam'],
        st.shear_web_2.layer['biax, right']
    ],
    soft_warning=False)
