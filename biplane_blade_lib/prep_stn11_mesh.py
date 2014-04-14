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
station_num = 11
# --------------------------------------
plt.close('all')

# load the biplane blade
b1 = bl.BiplaneBlade(
    'biplane blade, flapwise symmetric, no stagger, rj/R=0.452, g/c=1.25',
    'biplane_blade',
    rotate_airfoil_coords=False)

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
# points_usc = [
#     (-0.75, 2.09763311),  # SparCap_upper.txt
#     (-0.74, 2.08939698),  # InternalSurface2_resin.txt
#     ( 0.74, 2.13635553),  # InternalSurface2_resin.txt
#     ( 0.75, 2.14526482),  # SparCap_upper.txt
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
#     ( 0.75000000,  -2.14489854),  # SparCap_lower.txt
#     ( 0.74000000,  -2.13619881),  # InternalSurface2_resin.txt
#     (-0.74000000,  -2.11461209),  # InternalSurface2_resin.txt
#     (-0.75000000,  -2.12288194)   # SparCap_lower.txt
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
# is3 = st.internal_surface_3.layer['resin']
# points_te = [
#     (3.8,-2.0),
#     (2.58339700, -2.6),          # TE_Reinforcement_foam.txt
#     (2.58339700,  -1.40409979),  # TE_Reinforcement_foam.txt
#     (2.59334086,  -1.38461963),  # InternalSurface3_resin.txt
#     (3.0, 0.0),
#     is3.polygon.interiors[0].coords[73],
#     (2.58339700,   1.43169225),  # TE_Reinforcement_foam.txt
#     (2.58339700, 2.6),           # TE_Reinforcement_foam.txt
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
#     (-0.836, -2.09486467),  # LE_Panel_foam.txt
#     (-0.846, -2.08215284),  # InternalSurface1_resin.txt
#     (-1.5, 0.0),
#     (-0.846,  2.05679990),  # InternalSurface1_resin.txt
#     (-0.836,  2.06949738),  # LE_Panel_foam.txt
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
#     (0.836, 2.6),
#     (2.58339700, 2.6),           # AftPanel1_upper.txt
#     (2.58339700,   1.41312361),  # AftPanel1_upper.txt
#     (2.59, 1.39),                # by eye
#     (1.2, 1.9),
#     (0.846,   2.07117275),       # InternalSurface3_resin.txt
#     (0.836,   2.08281631)        # AftPanel1_upper.txt
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
#     (0.836, -2.6),
#     (2.58339700, -2.6),          # AftPanel1_lower.txt
#     (2.58339700,  -1.38561726),  # AftPanel1_lower.txt
#     (2.59, -1.36),                # by eye
#     (1.2, -1.9),
#     (0.846,  -2.06820336),       # InternalSurface3_resin.txt
#     (0.836,  -2.08012677)        # AftPanel1_lower.txt
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
#         st.TE_reinforcement.layer['foam'],
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
