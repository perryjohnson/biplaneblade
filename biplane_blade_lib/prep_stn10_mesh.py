"""Write initial TrueGrid files for one biplane blade station.

Usage
-----
start an IPython (qt)console with the pylab flag:
$ ipython qtconsole --pylab
or
$ ipython --pylab
Then, from the prompt, run this script:
|> %run biplane_blade_lib/prep_stnXX_mesh.py
or
|> import biplane_blade_lib/prep_stnXX_mesh

Author: Perry Roth-Johnson
Last updated: April 23, 2014

"""


import matplotlib.pyplot as plt
import lib.blade as bl
reload(bl)
import lib.poly_utils as pu
reload(pu)
from shapely.geometry import Polygon
from shapely.affinity import translate


# SET THESE PARAMETERS -----------------
station_num = 10
# --------------------------------------
plt.close('all')

# load the biplane blade
b1 = bl.BiplaneBlade(
    'biplane blade, flapwise symmetric, no stagger, rj/R=0.452, g/c=1.25',
    'biplane_blade')

# pre-process the station dimensions
station = b1.list_of_stations[station_num-1]
station.airfoil.create_polygon()
station.structure.create_all_layers()
station.structure.save_all_layer_edges()
station.structure.write_all_part_polygons()

# plot the parts
station.plot_parts()

# access the structure and airfoil for this station
st = station.structure
af = station.airfoil
x3_off = af.lower_chord * af.gap_to_chord_ratio * af.gap_fraction


# upper spar cap -----------------------------------------------------------
label = 'upper spar cap'

# create the bounding polygon
usc = st.lower_spar_cap.layer['upper']
is2 = st.lower_internal_surface_2.layer['resin']
points_usc = [
    (-0.75, usc.left[0][1]),             # lower_SparCap_upper.txt
    is2.polygon.interiors[0].coords[0],  # lower_InternalSurface2_resin.txt
    ( 0.742, -0.69618860),               # lower_InternalSurface2_resin.txt
    ( 0.75, usc.right[1][1]),            # lower_SparCap_upper.txt
    ( 0.75, 0.0),
    (-0.75, 0.0)
    ]
bounding_polygon = Polygon(points_usc)
pu.plot_polygon(bounding_polygon, 'None', '#000000')

# cut the new layer polygons
pu.cut_plot_and_write_alt_layer(st.lower_root_buildup, 'triax', label, 
    bounding_polygon, airfoil='lower')
pu.cut_plot_and_write_alt_layer(st.lower_external_surface, 'triax', label, 
    bounding_polygon, airfoil='lower')
pu.cut_plot_and_write_alt_layer(st.lower_external_surface, 'gelcoat', label, 
    bounding_polygon, airfoil='lower')
pu.cut_plot_and_write_alt_layer(st.lower_internal_surface_2, 'resin', label, 
    bounding_polygon, airfoil='lower')
pu.cut_plot_and_write_alt_layer(st.lower_internal_surface_2, 'triax', label, 
    bounding_polygon, airfoil='lower')

# lower spar cap -----------------------------------------------------------
label = 'lower spar cap'

# create the bounding polygon
lsc = st.lower_spar_cap.layer['lower']
points_lsc = [
    (-0.75,-4.0),
    ( 0.75,-4.0),
    ( 0.75000000,  lsc.right[0][1]),     # lower_SparCap_lower.txt
    ( 0.74200000,  -3.70894871),         # lower_InternalSurface2_resin.txt
    is2.polygon.interiors[0].coords[1],  # lower_InternalSurface2_resin.txt
    (-0.75000000,  lsc.left[1][1])       # lower_SparCap_lower.txt
    ]
bounding_polygon = Polygon(points_lsc)
pu.plot_polygon(bounding_polygon, 'None', '#000000')

# cut the new layer polygons
pu.cut_plot_and_write_alt_layer(st.lower_root_buildup, 'triax', label, 
    bounding_polygon, airfoil='lower')
pu.cut_plot_and_write_alt_layer(st.lower_external_surface, 'triax', label, 
    bounding_polygon, airfoil='lower')
pu.cut_plot_and_write_alt_layer(st.lower_external_surface, 'gelcoat', label, 
    bounding_polygon, airfoil='lower')
pu.cut_plot_and_write_alt_layer(st.lower_internal_surface_2, 'resin', label, 
    bounding_polygon, airfoil='lower')
pu.cut_plot_and_write_alt_layer(st.lower_internal_surface_2, 'triax', label, 
    bounding_polygon, airfoil='lower')

# TE reinforcement ---------------------------------------------------------
label = 'TE reinforcement'

# create the bounding polygon
ter = st.lower_TE_reinforcement.layer['uniax']
points_te = [
    (3.5,-4.0),
    (ter.bottom[1][0], -4.0),          # lower_TE_Reinforcement_foam.txt
    tuple(ter.bottom[1]),              # lower_TE_Reinforcement_foam.txt
    (2.36172800,  -3.14392266),        # lower_InternalSurface3_resin.txt
    (3.0, -x3_off),
    (2.36172800,  -1.20243664),        # lower_InternalSurface3_resin.txt
    tuple(ter.top[0]),                 # lower_TE_Reinforcement_foam.txt
    (ter.top[0][0], 0.0),              # lower_TE_Reinforcement_foam.txt
    (3.5, 0.0)
    ]
bounding_polygon = Polygon(points_te)
pu.plot_polygon(bounding_polygon, 'None', '#000000')

# cut the new layer polygons
pu.cut_plot_and_write_alt_layer(st.lower_root_buildup, 'triax', label, 
    bounding_polygon, airfoil='lower')
pu.cut_plot_and_write_alt_layer(st.lower_external_surface, 'triax', label, 
    bounding_polygon, airfoil='lower')
pu.cut_plot_and_write_alt_layer(st.lower_external_surface, 'gelcoat', label, 
    bounding_polygon, airfoil='lower')
pu.cut_plot_and_write_alt_layer(st.lower_internal_surface_3, 'resin', label, 
    bounding_polygon, airfoil='lower')
pu.cut_plot_and_write_alt_layer(st.lower_internal_surface_3, 'triax', label, 
    bounding_polygon, airfoil='lower')

# LE panel -----------------------------------------------------------------
label = 'LE panel'

# create the bounding polygon
lep = st.lower_LE_panel.layer['foam']
is1 = st.lower_internal_surface_1.layer['resin']
points_le = [
    (-3.00,-4.0),
    (-0.836,-4.0),
    tuple(lep.bottom[0]),                 # lower_LE_Panel_foam.txt
    is1.polygon.interiors[0].coords[-2],  # lower_InternalSurface1_resin.txt
    (-1.5, -x3_off),
    is1.polygon.interiors[0].coords[-1],  # lower_InternalSurface1_resin.txt
    tuple(lep.top[1]),                    # lower_LE_Panel_foam.txt
    (-0.836, 0.0),
    (-3.00, 0.0)
    ]
bounding_polygon = Polygon(points_le)
pu.plot_polygon(bounding_polygon, 'None', '#000000')
# cut the new layer polygons
pu.cut_plot_and_write_alt_layer(st.lower_root_buildup, 'triax', label, 
    bounding_polygon, airfoil='lower')
pu.cut_plot_and_write_alt_layer(st.lower_external_surface, 'triax', label, 
    bounding_polygon, airfoil='lower')
pu.cut_plot_and_write_alt_layer(st.lower_external_surface, 'gelcoat', label, 
    bounding_polygon, airfoil='lower')
pu.cut_plot_and_write_alt_layer(st.lower_internal_surface_1, 'resin', label, 
    bounding_polygon, airfoil='lower')
pu.cut_plot_and_write_alt_layer(st.lower_internal_surface_1, 'triax', label, 
    bounding_polygon, airfoil='lower')

# upper aft panel 1 -------------------------------------------------------
label = 'upper aft panel 1'
# create the bounding polygon
ap1u = st.lower_aft_panel_1.layer['upper']
is3 = st.lower_internal_surface_3.layer['resin']
points_ap1u = [
    (0.836, 0.0),
    (ap1u.right[1][0], 0.0),              # lower_AftPanel1_upper.txt
    tuple(ap1u.right[1]),                 # lower_AftPanel1_upper.txt
    (2.35714654,  -1.24251551),           # lower_InternalSurface3_resin.txt
    (2.0, -1.5),
    is3.polygon.interiors[0].coords[-2],  # lower_InternalSurface3_resin.txt
    tuple(ap1u.left[0])                   # lower_AftPanel1_upper.txt
    ]
bounding_polygon = Polygon(points_ap1u)
pu.plot_polygon(bounding_polygon, 'None', '#000000')
# cut the new layer polygons
pu.cut_plot_and_write_alt_layer(st.lower_root_buildup, 'triax', label, 
    bounding_polygon, airfoil='lower')
pu.cut_plot_and_write_alt_layer(st.lower_external_surface, 'triax', label, 
    bounding_polygon, airfoil='lower')
pu.cut_plot_and_write_alt_layer(st.lower_external_surface, 'gelcoat', label, 
    bounding_polygon, airfoil='lower')
pu.cut_plot_and_write_alt_layer(st.lower_internal_surface_3, 'resin', label, 
    bounding_polygon, airfoil='lower')
pu.cut_plot_and_write_alt_layer(st.lower_internal_surface_3, 'triax', label, 
    bounding_polygon, airfoil='lower')

# lower aft panel 1 -------------------------------------------------------
label = 'lower aft panel 1'
# create the bounding polygon
ap1l = st.lower_aft_panel_1.layer['lower']
points_ap1l = [
    (0.836, -4.0),
    (ap1l.right[0][0], -4.0),             # lower_AftPanel1_lower.txt
    tuple(ap1l.right[0]),                 # lower_AftPanel1_lower.txt
    (2.35717756,  -3.10392041),           # lower_InternalSurface3_resin.txt
    (2.0, -2.5),
    is3.polygon.interiors[0].coords[-1],  # lower_InternalSurface3_resin.txt
    tuple(ap1l.left[1])                   # lower_AftPanel1_lower.txt
    ]
bounding_polygon = Polygon(points_ap1l)
pu.plot_polygon(bounding_polygon, 'None', '#000000')
# cut the new layer polygons
pu.cut_plot_and_write_alt_layer(st.lower_root_buildup, 'triax', label, 
    bounding_polygon, airfoil='lower')
pu.cut_plot_and_write_alt_layer(st.lower_external_surface, 'triax', label, 
    bounding_polygon, airfoil='lower')
pu.cut_plot_and_write_alt_layer(st.lower_external_surface, 'gelcoat', label, 
    bounding_polygon, airfoil='lower')
pu.cut_plot_and_write_alt_layer(st.lower_internal_surface_3, 'resin', label, 
    bounding_polygon, airfoil='lower')
pu.cut_plot_and_write_alt_layer(st.lower_internal_surface_3, 'triax', label, 
    bounding_polygon, airfoil='lower')

# above shear web 1 ----------------------------------------------------------
label = 'above shear web 1'

# create the bounding polygon
points_asw1 = [
    (-0.75, 0.0),
    (-0.75, -1.5),
    (-0.836, -1.5),
    (-0.836, 0.0)
    ]
bounding_polygon = Polygon(points_asw1)
pu.plot_polygon(bounding_polygon, 'None', '#000000')

# cut the new layer polygons
pu.cut_plot_and_write_alt_layer(st.lower_root_buildup, 'triax', label, 
    bounding_polygon, airfoil='lower')
pu.cut_plot_and_write_alt_layer(st.lower_external_surface, 'triax', label, 
    bounding_polygon, airfoil='lower')
pu.cut_plot_and_write_alt_layer(st.lower_external_surface, 'gelcoat', label, 
    bounding_polygon, airfoil='lower')

# below shear web 1 ----------------------------------------------------------
label = 'below shear web 1'

# create the bounding polygon
points_bsw1 = [
    (-0.75, -4.0),
    (-0.75, -2.5),
    (-0.836, -2.5),
    (-0.836, -4.0)
    ]
bounding_polygon = Polygon(points_bsw1)
pu.plot_polygon(bounding_polygon, 'None', '#000000')

# cut the new layer polygons
pu.cut_plot_and_write_alt_layer(st.lower_root_buildup, 'triax', label, 
    bounding_polygon, airfoil='lower')
pu.cut_plot_and_write_alt_layer(st.lower_external_surface, 'triax', label, 
    bounding_polygon, airfoil='lower')
pu.cut_plot_and_write_alt_layer(st.lower_external_surface, 'gelcoat', label, 
    bounding_polygon, airfoil='lower')

# above shear web 2 ----------------------------------------------------------
label = 'above shear web 2'

# create the bounding polygon
points_asw2 = [
    (0.75, 0.0),
    (0.75, -1.5),
    (0.836, -1.5),
    (0.836, 0.0)
    ]
bounding_polygon = Polygon(points_asw2)
pu.plot_polygon(bounding_polygon, 'None', '#000000')

# cut the new layer polygons
pu.cut_plot_and_write_alt_layer(st.lower_root_buildup, 'triax', label, 
    bounding_polygon, airfoil='lower')
pu.cut_plot_and_write_alt_layer(st.lower_external_surface, 'triax', label, 
    bounding_polygon, airfoil='lower')
pu.cut_plot_and_write_alt_layer(st.lower_external_surface, 'gelcoat', label, 
    bounding_polygon, airfoil='lower')

# below shear web 2 ----------------------------------------------------------
label = 'below shear web 2'

# create the bounding polygon
points_bsw2 = [
    (0.75, -4.0),
    (0.75, -2.5),
    (0.836, -2.5),
    (0.836, -4.0)
    ]
bounding_polygon = Polygon(points_bsw2)
pu.plot_polygon(bounding_polygon, 'None', '#000000')

# cut the new layer polygons
pu.cut_plot_and_write_alt_layer(st.lower_root_buildup, 'triax', label, 
    bounding_polygon, airfoil='lower')
pu.cut_plot_and_write_alt_layer(st.lower_external_surface, 'triax', label, 
    bounding_polygon, airfoil='lower')
pu.cut_plot_and_write_alt_layer(st.lower_external_surface, 'gelcoat', label, 
    bounding_polygon, airfoil='lower')

# left of shear web 1 -------------------------------------------------------
label = 'left of shear web 1'

# create the bounding polygon
points_lsw1 = points_le[2:-2]
bounding_polygon = Polygon(points_lsw1)
pu.plot_polygon(bounding_polygon, 'None', '#000000')
# cut the new layer polygons
pu.cut_plot_and_write_alt_layer(st.lower_internal_surface_1, 'resin', label, 
    bounding_polygon, airfoil='lower')
pu.cut_plot_and_write_alt_layer(st.lower_internal_surface_1, 'triax', label, 
    bounding_polygon, airfoil='lower')

# right of shear web 1 -------------------------------------------------------
label = 'right of shear web 1'

# create the bounding polygon
points_rsw1 = [
    points_usc[0],
    points_usc[1],
    (0.0, -x3_off),
    points_lsc[-2],
    points_lsc[-1]
    ]
bounding_polygon = Polygon(points_rsw1)
pu.plot_polygon(bounding_polygon, 'None', '#000000')
# cut the new layer polygons
pu.cut_plot_and_write_alt_layer(st.lower_internal_surface_2, 'resin', label, 
    bounding_polygon, airfoil='lower')
pu.cut_plot_and_write_alt_layer(st.lower_internal_surface_2, 'triax', label, 
    bounding_polygon, airfoil='lower')

# left of shear web 2 -------------------------------------------------------
label = 'left of shear web 2'

# create the bounding polygon
points_lsw2 = [
    points_usc[3],
    points_usc[2],
    (0.0, -x3_off),
    points_lsc[3],
    points_lsc[2]
    ]
bounding_polygon = Polygon(points_lsw2)
pu.plot_polygon(bounding_polygon, 'None', '#000000')
# cut the new layer polygons
pu.cut_plot_and_write_alt_layer(st.lower_internal_surface_2, 'resin', label, 
    bounding_polygon, airfoil='lower')
pu.cut_plot_and_write_alt_layer(st.lower_internal_surface_2, 'triax', label, 
    bounding_polygon, airfoil='lower')

# right of shear web 2 -------------------------------------------------------
label = 'right of shear web 2'

# create the bounding polygon
points_rsw2 = [
    points_ap1u[-1],
    points_ap1u[-2],
    (1.5, -x3_off),
    points_ap1l[-2],
    points_ap1l[-1]
    ]
bounding_polygon = Polygon(points_rsw2)
pu.plot_polygon(bounding_polygon, 'None', '#000000')
# cut the new layer polygons
pu.cut_plot_and_write_alt_layer(st.lower_internal_surface_3, 'resin', label, 
    bounding_polygon, airfoil='lower')
pu.cut_plot_and_write_alt_layer(st.lower_internal_surface_3, 'triax', label, 
    bounding_polygon, airfoil='lower')

# -----------------------------------------------------------------------------

list_of_mesh_layers = []
# translate all the alt layers
for (name, layer) in st.lower_root_buildup.alt_layer.items():
    layer.move(x3_off, alt_layer=True)
    list_of_mesh_layers.append(layer)
for (name, layer) in st.lower_external_surface.alt_layer.items():
    layer.move(x3_off, alt_layer=True)
    list_of_mesh_layers.append(layer)
for (name, layer) in st.lower_internal_surface_1.alt_layer.items():
    layer.move(x3_off, alt_layer=True)
    list_of_mesh_layers.append(layer)
for (name, layer) in st.lower_internal_surface_2.alt_layer.items():
    layer.move(x3_off, alt_layer=True)
    list_of_mesh_layers.append(layer)
for (name, layer) in st.lower_internal_surface_3.alt_layer.items():
    layer.move(x3_off, alt_layer=True)
    list_of_mesh_layers.append(layer)
# translate all the remaining regular layers
st.lower_spar_cap.layer['upper'].move(x3_off)
st.lower_spar_cap.layer['lower'].move(x3_off)
st.lower_TE_reinforcement.layer['uniax'].move(x3_off)
st.lower_aft_panel_1.layer['upper'].move(x3_off)
st.lower_aft_panel_1.layer['lower'].move(x3_off)
st.lower_LE_panel.layer['foam'].move(x3_off)
st.lower_shear_web_1.layer['biax, left'].move(x3_off)
st.lower_shear_web_1.layer['foam'].move(x3_off)
st.lower_shear_web_1.layer['biax, right'].move(x3_off)
st.lower_shear_web_2.layer['biax, left'].move(x3_off)
st.lower_shear_web_2.layer['foam'].move(x3_off)
st.lower_shear_web_2.layer['biax, right'].move(x3_off)

list_of_mesh_layers.append(st.lower_spar_cap.layer['upper'])
list_of_mesh_layers.append(st.lower_spar_cap.layer['lower'])
list_of_mesh_layers.append(st.lower_TE_reinforcement.layer['uniax'])
list_of_mesh_layers.append(st.lower_aft_panel_1.layer['upper'])
list_of_mesh_layers.append(st.lower_aft_panel_1.layer['lower'])
list_of_mesh_layers.append(st.lower_LE_panel.layer['foam'])
list_of_mesh_layers.append(st.lower_shear_web_1.layer['biax, left'])
list_of_mesh_layers.append(st.lower_shear_web_1.layer['foam'])
list_of_mesh_layers.append(st.lower_shear_web_1.layer['biax, right'])
list_of_mesh_layers.append(st.lower_shear_web_2.layer['biax, left'])
list_of_mesh_layers.append(st.lower_shear_web_2.layer['foam'])
list_of_mesh_layers.append(st.lower_shear_web_2.layer['biax, right'])


# plot the lower airfoil in the local beam coordinate system
#  (translate it up by the appropriate gap distance: x3_off)
fig,ax = plt.subplots()
fmt1 = "Station #{0}, {1}, {2}% span\n"
fmt2 = "lower airfoil in local beam coordinate system (x3-offset = {3:+.4f})"
fmt = fmt1 + fmt2
ax.set_title(fmt.format(station.station_num, station.airfoil.name, 
    station.coords.x1, x3_off))
lp2 = translate(af.lower_polygon, yoff=x3_off)
(minx, miny, maxx, maxy) = lp2.bounds
ax.set_xlim([minx*1.2,maxx*1.2])
ax.set_ylim([miny*1.2,maxy*1.2])
plt.grid('on')
ax.set_xlabel('x2 [meters]')
ax.set_ylabel('x3 [meters]')
ax.set_aspect('equal')
for layer in list_of_mesh_layers:
    station.plot_polygon(layer.polygon, ax, layer.face_color, layer.edge_color,
        alpha=0.8)


# show the plots
plt.show()


# write the TrueGrid input file for mesh generation ---------------------
st.write_truegrid_inputfile(
    interrupt_flag=True, 
    additional_layers=[
        st.lower_spar_cap.layer['upper'],
        st.lower_spar_cap.layer['lower'],
        st.lower_TE_reinforcement.layer['uniax'],
        st.lower_aft_panel_1.layer['upper'],
        st.lower_aft_panel_1.layer['lower'],
        st.lower_LE_panel.layer['foam'],
        st.lower_shear_web_1.layer['biax, left'],
        st.lower_shear_web_1.layer['foam'],
        st.lower_shear_web_1.layer['biax, right'],
        st.lower_shear_web_2.layer['biax, left'],
        st.lower_shear_web_2.layer['foam'],
        st.lower_shear_web_2.layer['biax, right']
    ],
    soft_warning=False)
