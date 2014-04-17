"""Write initial TrueGrid files for one biplane blade station.

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
Last updated: April 17, 2014

"""


import matplotlib.pyplot as plt
import lib.blade as bl
import lib.poly_utils as pu
from shapely.geometry import Polygon


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
# station.structure.save_all_layer_edges()
# station.structure.write_all_part_polygons()

# plot the parts
station.plot_parts()

# access the structure for this station
st = station.structure

# # upper spar cap -----------------------------------------------------------
# label = 'upper spar cap'

# # create the bounding polygon
# is2 = st.internal_surface_2.layer['resin']
# points_usc = [
#     (-0.75, 2.26509618),
#     (-0.74, 2.25715046),
#     is2.polygon.interiors[0].coords[11],
#     ( 0.75, 2.30577725),
#     ( 0.75, 2.6),
#     (-0.75, 2.6)
#     ]
# bounding_polygon = Polygon(points_usc)
# pu.plot_polygon(bounding_polygon, 'None', '#000000')

# # cut the new layer polygons
# pu.cut_plot_and_write_alt_layer(st.root_buildup, 'triax', label, 
#     bounding_polygon)
# pu.cut_plot_and_write_alt_layer(st.external_surface, 'triax', label, 
#     bounding_polygon)
# pu.cut_plot_and_write_alt_layer(st.external_surface, 'gelcoat', label, 
#     bounding_polygon)
# pu.cut_plot_and_write_alt_layer(st.internal_surface_2, 'resin', label, 
#     bounding_polygon)
# pu.cut_plot_and_write_alt_layer(st.internal_surface_2, 'triax', label, 
#     bounding_polygon)

# # lower spar cap -----------------------------------------------------------
# label = 'lower spar cap'

# # create the bounding polygon
# points_lsc = [
#     (-0.75,-2.6),
#     ( 0.75,-2.6),
#     (0.75000000,  -2.30355608),
#     (0.74000000,  -2.29511264),
#     (-0.74000000,  -2.27314235),
#     (-0.75000000,  -2.28110958)
#     ]
# bounding_polygon = Polygon(points_lsc)
# pu.plot_polygon(bounding_polygon, 'None', '#000000')

# # cut the new layer polygons
# pu.cut_plot_and_write_alt_layer(st.root_buildup, 'triax', label, 
#     bounding_polygon)
# pu.cut_plot_and_write_alt_layer(st.external_surface, 'triax', label, 
#     bounding_polygon)
# pu.cut_plot_and_write_alt_layer(st.external_surface, 'gelcoat', label, 
#     bounding_polygon)
# pu.cut_plot_and_write_alt_layer(st.internal_surface_2, 'resin', label, 
#     bounding_polygon)
# pu.cut_plot_and_write_alt_layer(st.internal_surface_2, 'triax', label, 
#     bounding_polygon)

# # TE reinforcement ---------------------------------------------------------
# label = 'TE reinforcement'

# # create the bounding polygon
# points_te = [
#     (3.8,-2.0),
#     (2.35372800, -2.6),
#     (2.35372800,  -1.64819531),
#     (2.36372800,  -1.62867451),
#     (3.0, 0.0),
#     (2.36372800,   1.64660095),
#     (2.35372800,   1.66616240),
#     (2.35372800, 2.6),
#     (3.8, 2.0)
#     ]
# bounding_polygon = Polygon(points_te)
# pu.plot_polygon(bounding_polygon, 'None', '#000000')

# # cut the new layer polygons
# pu.cut_plot_and_write_alt_layer(st.root_buildup, 'triax', label, 
#     bounding_polygon)
# pu.cut_plot_and_write_alt_layer(st.external_surface, 'triax', label, 
#     bounding_polygon)
# pu.cut_plot_and_write_alt_layer(st.external_surface, 'gelcoat', label, 
#     bounding_polygon)
# pu.cut_plot_and_write_alt_layer(st.internal_surface_3, 'resin', label, 
#     bounding_polygon)
# pu.cut_plot_and_write_alt_layer(st.internal_surface_3, 'triax', label, 
#     bounding_polygon)

# # LE panel -----------------------------------------------------------------
# label = 'LE panel'

# # create the bounding polygon
# points_le = [
#     (-3.00,-2.6),
#     (-0.836,-2.6),
#     (-0.836, -2.24087898),
#     (-0.846, -2.22832757),
#     (-1.5, 0.0),
#     (-0.846,  2.21202909),
#     (-0.836,  2.22461448),
#     (-0.836, 2.6),
#     (-3.00, 2.6)
#     ]
# bounding_polygon = Polygon(points_le)
# pu.plot_polygon(bounding_polygon, 'None', '#000000')
# # cut the new layer polygons
# pu.cut_plot_and_write_alt_layer(st.root_buildup, 'triax', label, 
#     bounding_polygon)
# pu.cut_plot_and_write_alt_layer(st.external_surface, 'triax', label, 
#     bounding_polygon)
# pu.cut_plot_and_write_alt_layer(st.external_surface, 'gelcoat', label, 
#     bounding_polygon)
# pu.cut_plot_and_write_alt_layer(st.internal_surface_1, 'resin', label, 
#     bounding_polygon)
# pu.cut_plot_and_write_alt_layer(st.internal_surface_1, 'triax', label, 
#     bounding_polygon)

# # upper aft panel 1 -------------------------------------------------------
# label = 'upper aft panel 1'
# # create the bounding polygon
# points_ur = [
#     (0.836, 2.8),
#     (2.35372800, 2.6),
#     (2.35372800,   1.56461914),
#     (2.38, 1.52),
#     (1.2, 2.0),
#     (0.846,   2.20753447),
#     (0.836,   2.21969129)
#     ]
# bounding_polygon = Polygon(points_ur)
# pu.plot_polygon(bounding_polygon, 'None', '#000000')
# # cut the new layer polygons
# pu.cut_plot_and_write_alt_layer(st.root_buildup, 'triax', label, 
#     bounding_polygon)
# pu.cut_plot_and_write_alt_layer(st.external_surface, 'triax', label, 
#     bounding_polygon)
# pu.cut_plot_and_write_alt_layer(st.external_surface, 'gelcoat', label, 
#     bounding_polygon)
# pu.cut_plot_and_write_alt_layer(st.internal_surface_3, 'resin', label, 
#     bounding_polygon)
# pu.cut_plot_and_write_alt_layer(st.internal_surface_3, 'triax', label, 
#     bounding_polygon)

# # lower aft panel 1 -------------------------------------------------------
# label = 'lower aft panel 1'
# # create the bounding polygon
# points_lr = [
#     (0.836, -2.8),
#     (2.35372800, -2.6),
#     (2.35372800,  -1.54681749),
#     (2.4, -1.45),
#     (1.2, -2.0),
#     (0.846,  -2.20350638),
#     (0.836,  -2.21585048)
#     ]
# bounding_polygon = Polygon(points_lr)
# pu.plot_polygon(bounding_polygon, 'None', '#000000')
# # cut the new layer polygons
# pu.cut_plot_and_write_alt_layer(st.root_buildup, 'triax', label, 
#     bounding_polygon)
# pu.cut_plot_and_write_alt_layer(st.external_surface, 'triax', label, 
#     bounding_polygon)
# pu.cut_plot_and_write_alt_layer(st.external_surface, 'gelcoat', label, 
#     bounding_polygon)
# pu.cut_plot_and_write_alt_layer(st.internal_surface_3, 'resin', label, 
#     bounding_polygon)
# pu.cut_plot_and_write_alt_layer(st.internal_surface_3, 'triax', label, 
#     bounding_polygon)

# # above shear web 1 ----------------------------------------------------------
# label = 'above shear web 1'

# # create the bounding polygon
# points_asw1 = [
#     (-0.75, 2.7),
#     (-0.75, 2.0),
#     (-0.836, 2.0),
#     (-0.836, 2.7)
#     ]
# bounding_polygon = Polygon(points_asw1)
# pu.plot_polygon(bounding_polygon, 'None', '#000000')

# # cut the new layer polygons
# pu.cut_plot_and_write_alt_layer(st.root_buildup, 'triax', label, 
#     bounding_polygon)
# pu.cut_plot_and_write_alt_layer(st.external_surface, 'triax', label, 
#     bounding_polygon)
# pu.cut_plot_and_write_alt_layer(st.external_surface, 'gelcoat', label, 
#     bounding_polygon)

# # below shear web 1 ----------------------------------------------------------
# label = 'below shear web 1'

# # create the bounding polygon
# points_bsw1 = [
#     (-0.75, -2.7),
#     (-0.75, -2.0),
#     (-0.836, -2.0),
#     (-0.836, -2.7)
#     ]
# bounding_polygon = Polygon(points_bsw1)
# pu.plot_polygon(bounding_polygon, 'None', '#000000')

# # cut the new layer polygons
# pu.cut_plot_and_write_alt_layer(st.root_buildup, 'triax', label, 
#     bounding_polygon)
# pu.cut_plot_and_write_alt_layer(st.external_surface, 'triax', label, 
#     bounding_polygon)
# pu.cut_plot_and_write_alt_layer(st.external_surface, 'gelcoat', label, 
#     bounding_polygon)

# # above shear web 2 ----------------------------------------------------------
# label = 'above shear web 2'

# # create the bounding polygon
# points_asw2 = [
#     (0.75, 2.7),
#     (0.75, 2.0),
#     (0.836, 2.0),
#     (0.836, 2.7)
#     ]
# bounding_polygon = Polygon(points_asw2)
# pu.plot_polygon(bounding_polygon, 'None', '#000000')

# # cut the new layer polygons
# pu.cut_plot_and_write_alt_layer(st.root_buildup, 'triax', label, 
#     bounding_polygon)
# pu.cut_plot_and_write_alt_layer(st.external_surface, 'triax', label, 
#     bounding_polygon)
# pu.cut_plot_and_write_alt_layer(st.external_surface, 'gelcoat', label, 
#     bounding_polygon)

# # below shear web 2 ----------------------------------------------------------
# label = 'below shear web 2'

# # create the bounding polygon
# points_bsw2 = [
#     (0.75, -2.7),
#     (0.75, -2.0),
#     (0.836, -2.0),
#     (0.836, -2.7)
#     ]
# bounding_polygon = Polygon(points_bsw2)
# pu.plot_polygon(bounding_polygon, 'None', '#000000')

# # cut the new layer polygons
# pu.cut_plot_and_write_alt_layer(st.root_buildup, 'triax', label, 
#     bounding_polygon)
# pu.cut_plot_and_write_alt_layer(st.external_surface, 'triax', label, 
#     bounding_polygon)
# pu.cut_plot_and_write_alt_layer(st.external_surface, 'gelcoat', label, 
#     bounding_polygon)

# # left of shear web 1 -------------------------------------------------------
# label = 'left of shear web 1'

# # create the bounding polygon
# points_lsw1 = points_le[2:-2]
# bounding_polygon = Polygon(points_lsw1)
# pu.plot_polygon(bounding_polygon, 'None', '#000000')
# # cut the new layer polygons
# pu.cut_plot_and_write_alt_layer(st.internal_surface_1, 'resin', label, 
#     bounding_polygon)
# pu.cut_plot_and_write_alt_layer(st.internal_surface_1, 'triax', label, 
#     bounding_polygon)

# # right of shear web 1 -------------------------------------------------------
# label = 'right of shear web 1'

# # create the bounding polygon
# points_rsw1 = [
#     points_usc[0],
#     points_usc[1],
#     (0.0, 0.0),
#     points_lsc[-2],
#     points_lsc[-1]
#     ]
# bounding_polygon = Polygon(points_rsw1)
# pu.plot_polygon(bounding_polygon, 'None', '#000000')
# # cut the new layer polygons
# pu.cut_plot_and_write_alt_layer(st.internal_surface_2, 'resin', label, 
#     bounding_polygon)
# pu.cut_plot_and_write_alt_layer(st.internal_surface_2, 'triax', label, 
#     bounding_polygon)

# # left of shear web 2 -------------------------------------------------------
# label = 'left of shear web 2'

# # create the bounding polygon
# points_lsw2 = [
#     points_usc[3],
#     points_usc[2],
#     (0.0, 0.0),
#     points_lsc[3],
#     points_lsc[2]
#     ]
# bounding_polygon = Polygon(points_lsw2)
# pu.plot_polygon(bounding_polygon, 'None', '#000000')
# # cut the new layer polygons
# pu.cut_plot_and_write_alt_layer(st.internal_surface_2, 'resin', label, 
#     bounding_polygon)
# pu.cut_plot_and_write_alt_layer(st.internal_surface_2, 'triax', label, 
#     bounding_polygon)

# # right of shear web 2 -------------------------------------------------------
# label = 'right of shear web 2'

# # create the bounding polygon
# points_rsw2 = [
#     points_ur[-1],
#     points_ur[-2],
#     (1.5, 0.0),
#     points_lr[-2],
#     points_lr[-1]
#     ]
# bounding_polygon = Polygon(points_rsw2)
# pu.plot_polygon(bounding_polygon, 'None', '#000000')
# # cut the new layer polygons
# pu.cut_plot_and_write_alt_layer(st.internal_surface_3, 'resin', label, 
#     bounding_polygon)
# pu.cut_plot_and_write_alt_layer(st.internal_surface_3, 'triax', label, 
#     bounding_polygon)




# show the plot
plt.show()

# # write the TrueGrid input file for mesh generation ---------------------
# st.write_truegrid_inputfile(
#     interrupt_flag=True, 
#     additional_layers=[
#         st.spar_cap.layer['upper'],
#         st.spar_cap.layer['lower'],
#         st.TE_reinforcement.layer['uniax'],
#         st.aft_panel_1.layer['upper'],
#         st.aft_panel_1.layer['lower'],
#         st.LE_panel.layer['foam'],
#         st.shear_web_1.layer['biax, left'],
#         st.shear_web_1.layer['foam'],
#         st.shear_web_1.layer['biax, right'],
#         st.shear_web_2.layer['biax, left'],
#         st.shear_web_2.layer['foam'],
#         st.shear_web_2.layer['biax, right']
#     ],
#     soft_warning=False)
