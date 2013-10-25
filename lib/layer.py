"""A module for manipulating single material layers of structural parts.

Author: Perry Roth-Johnson
Last updated: October 11, 2013

"""


import numpy as np
from shapely.geometry import asLineString


class Layer:
    """Define a layer made of a single material for a structural part.

    Parameters:
    -----------
    polygon : shapely.Polygon object, represents the boundary of this layer
    material : Material object, the material properties of this layer
    parent_part : structure.Part object, the part this layer belongs to
    name : str, the name of this layer
    mass : float, the mass per unit length of this layer
    left : np.array, coords for left edge, saved later by get_and_save_edges()
    top : np.array, coords for top edge, saved later by get_and_save_edges()
    right : np.array, coords for right edge, saved later by
        get_and_save_edges()
    bottom : np.array, coords for bottom edge, saved later by
        get_and_save_edges()

    """
    def __init__(self, polygon, material, parent_part, name=''):
        self.polygon = polygon
        self.material = material
        self.parent_part = parent_part
        self.name = name
        self.mass = self.polygon.area*self.material.rho  # mass per unit length
        self.left = None  # saved later by <part>.get_and_save_edges()
        self.top = None  # saved later by <part>.get_and_save_edges()
        self.right = None  # saved later by <part>.get_and_save_edges()
        self.bottom = None  # saved later by <part>.get_and_save_edges()

    def area_fraction(self):
        """Calculate the ratio of this part area to the cross-section area."""
        total_area = self.parent_part.parent_structure.area
        return self.polygon.area/total_area

    def mass_fraction(self):
        """Calculate the ratio of this part mass to the cross-section mass."""
        total_mass = self.parent_part.parent_structure.mass
        return self.mass/total_mass

    def plot_edges(self, axes):
        """Plots the polygon edges of this layer."""
        if self.left is None:
            raise ValueError("Layer instance has attribute <layer>.left=None.\n  Try running <station>.structure.<part>.get_and_save_edges() first.")
        axes.plot(self.left[:,0], self.left[:,1], 'bo-')
        axes.plot(self.top[:,0], self.top[:,1], 'g^-')
        axes.plot(self.right[:,0], self.right[:,1], 'rs-')
        axes.plot(self.bottom[:,0], self.bottom[:,1], 'c*-')

    def get_edges(self):
        """Returns 4 arrays of coords for each edge of this layer."""
        p = self.polygon  # get the polygon for this layer
        # store the polygon exterior coords as a numpy array
        a = np.array(p.exterior.coords)
        # get the x- and y-coordinates of the polygon exterior
        x = a[:,0]
        y = a[:,1]
        if self.parent_part.__class__.__name__ == 'RootBuildup':
            # find the indices where the x-coord is equal to zero
            match_x = np.nonzero(x==0.0)[0]
            # find the indices where the y-coord is equal to the right edge
            match_y = np.nonzero(y==0.0)[0]
            # group all the indices together in a sorted array
            match = np.append(match_x, match_y)
        elif self.parent_part.__class__.__name__ == 'LE_Panel':
            # get the coordinates for the right edge of the LE panel
            r = self.parent_part.right
            # find the indices where the x-coordinate is equal to the right edge
            #   of the LE panel (the "corners" of the LE panel)
            match = np.nonzero(x==r)[0]
        elif self.parent_part.__class__.__name__ == 'SparCap':
            # get the coordinates for the left and right edges
            l = self.parent_part.left
            r = self.parent_part.right
            # find the indices where the x-coord is equal to the left edge
            match_l = np.nonzero(x==l)[0]
            # find the indices where the x-coord is equal to the right edge
            match_r = np.nonzero(x==r)[0]
            # group all the indices together in a sorted array
            match = np.append(match_l, match_r)
        elif self.parent_part.__class__.__name__ == 'AftPanel':
            # get the coordinates for the left and right edges
            l = self.parent_part.left
            r = self.parent_part.right
            # find the indices where the x-coord is equal to the left edge
            match_l = np.nonzero(x==l)[0]
            # find the indices where the x-coord is equal to the right edge
            match_r = np.nonzero(x==r)[0]
            # group all the indices together in a sorted array
            match = np.append(match_l, match_r)
        elif self.parent_part.__class__.__name__ == 'TE_Reinforcement':
            # get the coordinates for the left edge
            l = self.parent_part.left
            # find the indices where the x-coord is equal to the left edge
            match = np.nonzero(x==l)[0]
        elif self.parent_part.__class__.__name__ == 'ShearWeb':
            # get the coordinates for the left and right edges
            if self.name == 'biax, left':
                l = self.parent_part.left
                r = self.parent_part.left + self.parent_part.base_biax
            elif self.name == 'foam':
                l = self.parent_part.left + self.parent_part.base_biax
                r = self.parent_part.right - self.parent_part.base_biax
            elif self.name == 'biax, right':
                l = self.parent_part.right - self.parent_part.base_biax
                r = self.parent_part.right
            # find the indices where the x-coord is equal to the left edge
            match_l = np.nonzero(x==l)[0]
            # find the indices where the x-coord is equal to the right edge
            match_r = np.nonzero(x==r)[0]
            # group all the indices together in a sorted array
            match = np.append(match_l, match_r)
        match.sort()
        # split the polygon up at each of the corners into 4 "edges"
        edge1 = a[match[0]:match[1]+1,:]
        edge2 = a[match[1]:match[2]+1,:]
        edge3 = a[match[2]:match[3]+1,:]
        try:
            edge4 = a[match[3]:match[4]+1,:]
        except IndexError:
            edge4 = np.append(a[match[3]:,:],a[1:match[0]+1,:],axis=0)
        return (edge1, edge2, edge3, edge4)

    def get_and_save_edges(self):
        """Identifies and saves the left, top, right, and bottom (LTRB) edges.

        This method saves the LTRB edges as attributes within the layer object.

        self.left : np.array, coords for left edge
        self.top : np.array, coords for top edge
        self.right : np.array, coords for right edge
        self.bottom : np.array, coords for bottom edge

        """
        edges = self.get_edges()
        # get centroids
        centroids = []
        for edge in edges:
            centroids.append(asLineString(edge).centroid)
        # determine which edges are top, bottom, left, and right
        l = range(4)  # list of indices, one for each edge
        c = np.array([[centroids[0].x, centroids[0].y],
                      [centroids[1].x, centroids[1].y],
                      [centroids[2].x, centroids[2].y],
                      [centroids[3].x, centroids[3].y]])
        cx = c[:,0]
        cy = c[:,1]
        x_max_ind = np.nonzero(cx==cx.max())[0][0]
        x_min_ind = np.nonzero(cx==cx.min())[0][0]
        y_max_ind = np.nonzero(cy==cy.max())[0][0]
        y_min_ind = np.nonzero(cy==cy.min())[0][0]
        if self.parent_part.__class__.__name__ == 'RootBuildup':
            # find centroid at x=0
            ind_x = np.nonzero(cx==0.0)[0][0]
            l.remove(ind_x)  # remove the index for the right edge
            # find centroid at y=0
            ind_y = np.nonzero(cy==0.0)[0][0]
            l.remove(ind_y)  # remove the index for the left edge
            if self.name == 'triax, lower left':
                self.right = edges[ind_x]  # right edge saved!
                self.left = edges[ind_y]  # left edge saved!
            elif self.name == 'triax, lower right':
                self.left = edges[ind_x]  # left edge saved!
                self.right = edges[ind_y]  # right edge saved!
            elif self.name == 'triax, upper right':
                self.left = edges[ind_x]  # left edge saved!
                self.right = edges[ind_y]  # right edge saved!
            elif self.name == 'triax, upper left':
                self.right = edges[ind_x]  # right edge saved!
                self.left = edges[ind_y]  # left edge saved!
            # find top and bottom edges
            if centroids[l[0]].y > centroids[l[1]].y:
                self.top = edges[l[0]]     # top edge saved!
                self.bottom = edges[l[1]]  # bottom edge saved!
            else:
                self.top = edges[l[1]]     # top edge saved!
                self.bottom = edges[l[0]]  # bottom edge saved!
        elif self.parent_part.__class__.__name__ == 'LE_Panel':
            l.remove(x_min_ind)  # remove the index for the left edge
            self.left = edges[x_min_ind]  # left edge saved!
            l.remove(y_max_ind)  # remove the index for the top edge
            self.top = edges[y_max_ind]  # top edge saved!
            l.remove(y_min_ind)  # remove the index for the bottom edge
            self.bottom = edges[y_min_ind]  # bottom edge saved!
            self.right = edges[l[0]]  # right edge saved!
        elif self.parent_part.__class__.__name__ == 'SparCap':
            l.remove(x_min_ind)  # remove the index for the left edge
            self.left = edges[x_min_ind]  # left edge saved!
            l.remove(x_max_ind)  # remove the index for the right edge
            self.right = edges[x_max_ind]  # right edge saved!
            if centroids[l[0]].y > centroids[l[1]].y:
                self.top = edges[l[0]]     # top edge saved!
                self.bottom = edges[l[1]]  # bottom edge saved!
            else:
                self.top = edges[l[1]]     # top edge saved!
                self.bottom = edges[l[0]]  # bottom edge saved!
        elif self.parent_part.__class__.__name__ == 'AftPanel':
            l.remove(x_min_ind)  # remove the index for the left edge
            self.left = edges[x_min_ind]  # left edge saved!
            l.remove(x_max_ind)  # remove the index for the right edge
            self.right = edges[x_max_ind]  # right edge saved!
            if centroids[l[0]].y > centroids[l[1]].y:
                self.top = edges[l[0]]     # top edge saved!
                self.bottom = edges[l[1]]  # bottom edge saved!
            else:
                self.top = edges[l[1]]     # top edge saved!
                self.bottom = edges[l[0]]  # bottom edge saved!
        elif self.parent_part.__class__.__name__ == 'TE_Reinforcement':
            l.remove(x_max_ind)  # remove the index for the right edge
            self.right = edges[x_max_ind]  # right edge saved!
            l.remove(y_max_ind)  # remove the index for the top edge
            self.top = edges[y_max_ind]  # top edge saved!
            l.remove(y_min_ind)  # remove the index for the bottom edge
            self.bottom = edges[y_min_ind]  # bottom edge saved!
            self.left = edges[l[0]]  # left edge saved!
        elif self.parent_part.__class__.__name__ == 'ShearWeb':
            l.remove(x_min_ind)  # remove the index for the left edge
            self.left = edges[x_min_ind]  # left edge saved!
            l.remove(x_max_ind)  # remove the index for the right edge
            self.right = edges[x_max_ind]  # right edge saved!
            if centroids[l[0]].y > centroids[l[1]].y:
                self.top = edges[l[0]]     # top edge saved!
                self.bottom = edges[l[1]]  # bottom edge saved!
            else:
                self.top = edges[l[1]]     # top edge saved!
                self.bottom = edges[l[0]]  # bottom edge saved!
