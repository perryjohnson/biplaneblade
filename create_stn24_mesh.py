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
station_num = 24
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
is2 = st.internal_surface_2.layer['resin']
points_usc = [
    (-0.75, usc.left[0][1]),              # SparCap_upper.txt
    is2.polygon.interiors[0].coords[0],   # InternalSurface2_resin.txt
    is2.polygon.interiors[0].coords[48-34],   # InternalSurface2_resin.txt
    ( 0.75, usc.right[1][1]),             # SparCap_upper.txt
    ( 0.75, 0.7),
    (-0.75, 0.7)
    ]
bounding_polygon = Polygon(points_usc)
pu.plot_polygon(bounding_polygon, 'None', '#000000')

# cut the new layer polygons
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
    (-0.75,-0.35),
    ( 0.75,-0.35),
    (0.75000000,  lsc.right[0][1]),      # SparCap_lower.txt
    is2.polygon.interiors[0].coords[47-34],   # InternalSurface2_resin.txt
    is2.polygon.interiors[0].coords[1],   # InternalSurface2_resin.txt
    (-0.75000000,  lsc.left[1][1])        # SparCap_lower.txt
    ]
bounding_polygon = Polygon(points_lsc)
pu.plot_polygon(bounding_polygon, 'None', '#000000')

# cut the new layer polygons
pu.cut_plot_and_write_alt_layer(st.external_surface, 'triax', label, 
    bounding_polygon)
pu.cut_plot_and_write_alt_layer(st.external_surface, 'gelcoat', label, 
    bounding_polygon)
pu.cut_plot_and_write_alt_layer(st.internal_surface_2, 'resin', label, 
    bounding_polygon)
pu.cut_plot_and_write_alt_layer(st.internal_surface_2, 'triax', label, 
    bounding_polygon)

# TE reinforcement, upper 1 ------------------------------------------------
label = 'TE reinforcement, upper 1'

# create the bounding polygon
ter = st.TE_reinforcement.layer['foam']
is3 = st.internal_surface_3.layer['resin']
points_teu1 = [
    (ter.top[0][0], 0.35),              # TE_Reinforcement_foam.txt
    tuple(ter.top[0]),                  # TE_Reinforcement_foam.txt
    is3.polygon.interiors[0].coords[436-186],  # InternalSurface3_resin.txt
    (2.5, 0.1),
    is3.polygon.interiors[0].coords[424-186],  # InternalSurface3_resin.txt
    (is3.polygon.interiors[0].coords[424-186][0],   0.35)                # InternalSurface3_resin.txt
    ]
bounding_polygon = Polygon(points_teu1)
pu.plot_polygon(bounding_polygon, 'None', '#000000')

# cut the new layer polygons
pu.cut_plot_and_write_alt_layer(st.external_surface, 'triax', label, 
    bounding_polygon)
pu.cut_plot_and_write_alt_layer(st.external_surface, 'gelcoat', label, 
    bounding_polygon)
pu.cut_plot_and_write_alt_layer(st.internal_surface_3, 'resin', label, 
    bounding_polygon)
pu.cut_plot_and_write_alt_layer(st.internal_surface_3, 'triax', label, 
    bounding_polygon)
pu.cut_plot_and_write_alt_layer(st.TE_reinforcement, 'foam', label,
    bounding_polygon)
pu.cut_plot_and_write_alt_layer(st.TE_reinforcement, 'uniax', label, 
    bounding_polygon)

# TE reinforcement, lower 1 ------------------------------------------------
label = 'TE reinforcement, lower 1'

# create the bounding polygon
points_tel1 = [
    (ter.bottom[0][0], -0.1),              # TE_Reinforcement_foam.txt
    tuple(ter.bottom[1]),                  # TE_Reinforcement_foam.txt
    is3.polygon.interiors[0].coords[295-186],  # InternalSurface3_resin.txt
    points_teu1[-3],
    points_teu1[-2],         # InternalSurface3_resin.txt
    (points_teu1[-1][0],   -0.1)                # InternalSurface3_resin.txt
    ]
bounding_polygon = Polygon(points_tel1)
pu.plot_polygon(bounding_polygon, 'None', '#000000')

# cut the new layer polygons
pu.cut_plot_and_write_alt_layer(st.external_surface, 'triax', label, 
    bounding_polygon)
pu.cut_plot_and_write_alt_layer(st.external_surface, 'gelcoat', label, 
    bounding_polygon)
pu.cut_plot_and_write_alt_layer(st.internal_surface_3, 'resin', label, 
    bounding_polygon)
pu.cut_plot_and_write_alt_layer(st.internal_surface_3, 'triax', label, 
    bounding_polygon)
pu.cut_plot_and_write_alt_layer(st.TE_reinforcement, 'foam', label,
    bounding_polygon)
pu.cut_plot_and_write_alt_layer(st.TE_reinforcement, 'uniax', label, 
    bounding_polygon)

# TE reinforcement, upper 2 ------------------------------------------------
label = 'TE reinforcement, upper 2'

# create the bounding polygon
is3t = st.internal_surface_3.layer['triax']
points_teu2 = [
    points_teu1[-1],
    points_teu1[-2],
    is3t.polygon.interiors[0].coords[259-127],  # InternalSurface3_triax.txt
    is3t.polygon.exterior.coords[40-3],  # InternalSurface3_triax.txt
    (is3t.polygon.exterior.coords[40-3][0],   0.35) # InternalSurface3_triax.txt
    ]
bounding_polygon = Polygon(points_teu2)
pu.plot_polygon(bounding_polygon, 'None', '#000000')

# cut the new layer polygons
pu.cut_plot_and_write_alt_layer(st.external_surface, 'triax', label, 
    bounding_polygon)
pu.cut_plot_and_write_alt_layer(st.external_surface, 'gelcoat', label, 
    bounding_polygon)
pu.cut_plot_and_write_alt_layer(st.internal_surface_3, 'resin', label, 
    bounding_polygon)
pu.cut_plot_and_write_alt_layer(st.internal_surface_3, 'triax', label, 
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
pu.cut_plot_and_write_alt_layer(st.internal_surface_3, 'resin', label, 
    bounding_polygon)
pu.cut_plot_and_write_alt_layer(st.internal_surface_3, 'triax', label, 
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
    (est.polygon.exterior.coords[-1][0], 0.005),
    est.polygon.exterior.coords[-2],
    esg.polygon.exterior.coords[-2],
    (esg.polygon.exterior.coords[-2][0], 0.35)
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
    (-3.00,-1.6),
    (-0.836,-1.6),
    tuple(lep.bottom[0]),  # LE_Panel_foam.txt
    is1.polygon.interiors[0].coords[-2],  # InternalSurface1_resin.txt
    (-1.5, 0.0),
    is1.polygon.interiors[0].coords[-1],  # InternalSurface1_resin.txt
    tuple(lep.top[1]),  # LE_Panel_foam.txt
    (-0.836, 1.3),
    (-3.00, 1.3)
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

# upper aft panel 1 -------------------------------------------------------
label = 'upper aft panel 1'

# create the bounding polygon
ap1u = st.aft_panel_1.layer['upper']
is3 = st.internal_surface_3.layer['resin']
points_ap1u = [
    (0.836, 1.3),
    (ap1u.right[1][0], 1.3),              # AftPanel1_upper.txt
    tuple(ap1u.right[1]),                 # AftPanel1_upper.txt
    (2.14, 0.20),
    (1.2, 0.2),
    is3.polygon.interiors[0].coords[-2],  # InternalSurface3_resin.txt
    tuple(ap1u.left[0])                   # AftPanel1_upper.txt
    ]
bounding_polygon = Polygon(points_ap1u)
pu.plot_polygon(bounding_polygon, 'None', '#000000')
# cut the new layer polygons
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
points_ap1l = [
    (0.836, -1.6),
    (ap1l.right[0][0], -1.6),             # AftPanel1_lower.txt
    tuple(ap1l.right[0]),                 # AftPanel1_lower.txt
    (2.14, 0.1),
    (1.2, 0.1),
    is3.polygon.interiors[0].coords[-1],  # InternalSurface3_resin.txt
    tuple(ap1l.left[1])                   # AftPanel1_lower.txt
    ]
bounding_polygon = Polygon(points_ap1l)
pu.plot_polygon(bounding_polygon, 'None', '#000000')
# cut the new layer polygons
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
    (-0.75, 0.1),
    (-0.836, 0.1),
    (-0.836, 2.1)
    ]
bounding_polygon = Polygon(points_asw1)
pu.plot_polygon(bounding_polygon, 'None', '#000000')

# cut the new layer polygons
pu.cut_plot_and_write_alt_layer(st.external_surface, 'triax', label, 
    bounding_polygon)
pu.cut_plot_and_write_alt_layer(st.external_surface, 'gelcoat', label, 
    bounding_polygon)

# below shear web 1 ----------------------------------------------------------
label = 'below shear web 1'

# create the bounding polygon
points_bsw1 = [
    (-0.75, -2.1),
    (-0.75, -0.1),
    (-0.836, -0.1),
    (-0.836, -2.1)
    ]
bounding_polygon = Polygon(points_bsw1)
pu.plot_polygon(bounding_polygon, 'None', '#000000')

# cut the new layer polygons
pu.cut_plot_and_write_alt_layer(st.external_surface, 'triax', label, 
    bounding_polygon)
pu.cut_plot_and_write_alt_layer(st.external_surface, 'gelcoat', label, 
    bounding_polygon)

# above shear web 2 ----------------------------------------------------------
label = 'above shear web 2'

# create the bounding polygon
points_asw2 = [
    (0.75, 2.1),
    (0.75, 0.1),
    (0.836, 0.1),
    (0.836, 2.1)
    ]
bounding_polygon = Polygon(points_asw2)
pu.plot_polygon(bounding_polygon, 'None', '#000000')

# cut the new layer polygons
pu.cut_plot_and_write_alt_layer(st.external_surface, 'triax', label, 
    bounding_polygon)
pu.cut_plot_and_write_alt_layer(st.external_surface, 'gelcoat', label, 
    bounding_polygon)

# below shear web 2 ----------------------------------------------------------
label = 'below shear web 2'

# create the bounding polygon
points_bsw2 = [
    (0.75, -2.1),
    (0.75, -0.1),
    (0.836, -0.1),
    (0.836, -2.1)
    ]
bounding_polygon = Polygon(points_bsw2)
pu.plot_polygon(bounding_polygon, 'None', '#000000')

# cut the new layer polygons
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
    points_ap1u[-1],
    points_ap1u[-2],
    (1.5, 0.0),
    points_ap1l[-2],
    points_ap1l[-1]
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
        # st.TE_reinforcement.layer['uniax'],
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
    alt_TE_reinforcement=True,
    soft_warning=False)
