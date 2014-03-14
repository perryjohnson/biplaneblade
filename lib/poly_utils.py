"""Preprocessing tools to cut polygons and write their new coordinates.

Author: Perry Roth-Johnson
Last modified: March 14, 2014

"""


import matplotlib.pyplot as plt
from shapely.geometry import Polygon
from descartes import PolygonPatch


# cut up the layer polygons to prepare for grid generation
def cut_polygon(original, bounding):
    """Cut the original layer polygon with the bounding polygon."""
    return original.intersection(bounding)

def plot_polygon(p, face_color, edge_color='r'):
    """Plot a polygon on the current axes."""
    # get the current axes, so we can add polygons to the plot
    ax = plt.gcf().gca()
    patch = PolygonPatch(p, fc=face_color, ec=edge_color, alpha=0.8)
    ax.add_patch(patch)

def cut_plot_and_write_alt_layer(part, material, ext_label, b_polygon):
    """Cut, plot, and write a polygon for an alternate layer."""
    l = part.layer[material]
    # cut polygon
    p_new = cut_polygon(l.polygon, b_polygon)
    # check if extra negligibly small polygons were created
    if p_new.geom_type != 'Polygon':
        print '  Found a non-Polygon made of {0} polygons.'.format(len(p_new.geoms))
        good_poly_index = None
        for i,p in enumerate(p_new.geoms):
            print '  polygon[{0}]: area={1}, centroid={2}'.format(i,p.area,p.centroid.xy)
            if p.area > 1.0e-08:
                # only keep polygons with significant area
                good_poly_index = i
                print '  ...keep polygon[{0}]!'.format(i)
            else:
                # throw out any polygon with insignificant area
                print '  ...throw out polygon[{0}]'.format(i)
        # overwrite p_new with the good polygon
        p_new = p_new.geoms[good_poly_index]
    # plot polygon
    plot_polygon(p_new, l.face_color)
    new_layer_name = material + ', ' + ext_label
    l.parent_part.add_new_layer(new_layer_name, p_new, material)
    # write polygon
    part.alt_layer[new_layer_name].write_polygon_edges()