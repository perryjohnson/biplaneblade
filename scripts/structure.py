"""A module for organizing structural part data for a blade station.

Author: Perry Roth-Johnson
Last updated: August 7, 2013

"""


import numpy as np
from math import isnan


class Part:
    """Define the dimensions of a structural part."""
    def __init__(self, base, height):
        self.base = base
        self.height = height
    def __str__(self):
        return """base:    {0} (meters)
height:  {1} (meters)""".format(self.base, self.height)
    def exists(self):
        """Checks if a structural part exists at this station."""
        if isnan(self.base) and isnan(self.height):
            return False
        else:
            return True


class ShearWeb(Part):
    """Define the biax (skin) and foam (core) dimensions of a shear web.

    Parameters
    ----------
    base_biax : float, base of biax material in shear web
    base_foam : float, base of foam material in shear web
    x2 : float, chordwise distance from pitch axis to edge of shear web
    height : NaN (DO NOT SPECIFY), height of shear web

    """
    def __init__(self, base_biax, base_foam, x2, height=np.nan):
        Part.__init__(self, base=(2.0*base_biax+base_foam), height=height)
        self.base_biax = base_biax
        self.base_foam = base_foam
        self.x2 = x2
    def __str__(self):
        return """base:    {0:6.4f} (meters)
|-> base_biax:  {1:6.4f} (meters)
|-> base_foam:  {2:6.4f} (meters)
height:  {3} (meters)
x2:      {4:6.4f} (meters)""".format(self.base, self.base_biax,
    self.base_foam, self.height, self.x2)


class TE_Reinforcement(Part):
    """Define uniax and foam dimensions of a trailing edge reinforcement."""
    def __init__(self, base, height_uniax, height_foam):
        Part.__init__(self, base, height=(height_uniax+height_foam))
        self.height_uniax = height_uniax
        self.height_foam = height_foam
    def __str__(self):
        return """base:    {0:6.4f} (meters)
height:  {1:6.4f} (meters)
|-> height_uniax:  {2:6.4f} (meters)
|-> height_foam:   {3:6.4f} (meters)""".format(self.base, self.height,
    self.height_uniax, self.height_foam)


class InternalSurface(Part):
    """Define triax and resin dimensions of the internal surface."""
    def __init__(self, base, height_triax, height_resin):
        Part.__init__(self, base, height=(height_triax+height_resin))
        self.height_triax = height_triax
        self.height_resin = height_resin
    def __str__(self):
        return """base:    {0:6.4f} (meters)
height:  {1:6.4f} (meters)
|-> height_triax:  {2:6.4f} (meters)
|-> height_resin:  {3:6.4f} (meters)""".format(self.base, self.height,
    self.height_triax, self.height_resin)


class ExternalSurface(Part):
    """Define triax and gelcoat dimensions of the external surface."""
    def __init__(self, base, height_triax, height_gelcoat):
        Part.__init__(self, base, height=(height_triax+height_gelcoat))
        self.height_triax = height_triax
        self.height_gelcoat = height_gelcoat
    def __str__(self):
        return """base:    {0:6.4f} (meters)
height:  {1:6.4f} (meters)
|-> height_triax:  {2:6.4f} (meters)
|-> height_gelcoat:  {3:6.4f} (meters)""".format(self.base, self.height,
    self.height_triax, self.height_gelcoat)


class MonoplaneStructure:
    """Define the monoplane laminate schedule (internal dimensions)."""
    def __init__(self, h_RB, b_SC, h_SC, b_SW1_biax, b_SW1_foam, x2_SW1,
                 b_SW2_biax, b_SW2_foam, x2_SW2, b_SW3_biax, b_SW3_foam,
                 x2_SW3, b_TE_reinf, h_TE_reinf_uniax, h_TE_reinf_foam,
                 h_LE_panel, h_aft_panel, h_int_surf_triax, h_int_surf_resin,
                 h_ext_surf_triax, h_ext_surf_gelcoat):
        self.root_buildup = Part(np.nan, h_RB)
        self.spar_cap = Part(b_SC, h_SC)
        self.shear_web_1 = ShearWeb(b_SW1_biax, b_SW1_foam, x2_SW1)
        self.shear_web_2 = ShearWeb(b_SW2_biax, b_SW2_foam, x2_SW2)
        self.shear_web_3 = ShearWeb(b_SW3_biax, b_SW3_foam, x2_SW3)
        self.TE_reinforcement = TE_Reinforcement(b_TE_reinf, h_TE_reinf_uniax, 
                                                 h_TE_reinf_foam)
        self.LE_panel = Part(np.nan, h_LE_panel)
        self.aft_panel = Part(np.nan, h_aft_panel)
        self.internal_surface = InternalSurface(np.nan, h_int_surf_triax, h_int_surf_resin)
        self.external_surface = ExternalSurface(np.nan, h_ext_surf_triax, h_ext_surf_gelcoat)

    def __str__(self):
        """Returns a string of all the internal dimensions for this structure."""
        s = ''
        s += "--- ROOT BUILDUP ---\n"
        s += str(self.root_buildup) + '\n'
        s += "--- SPAR CAP ---\n"
        s += str(self.spar_cap) + '\n'
        s += "--- SHEAR WEB 1 ---\n"
        s += str(self.shear_web_1) + '\n'
        s += "--- SHEAR WEB 2 ---\n"
        s += str(self.shear_web_2) + '\n'
        s += "--- SHEAR WEB 3 ---\n"
        s += str(self.shear_web_3) + '\n'
        s += "---TE REINFORCEMENT ---\n"
        s += str(self.TE_reinforcement) + '\n'
        s += "--- LE PANEL ---\n"
        s += str(self.LE_panel) + '\n'
        s += "--- AFT PANEL ---\n"
        s += str(self.aft_panel) + '\n'
        s += "--- INTERNAL SURFACE ---\n"
        s += str(self.internal_surface) + '\n'
        s += "--- EXTERNAL SURFACE ---\n"
        s += str(self.external_surface) + '\n'
        return s

    def which_parts_exist(self):
        """Check which structural parts exist at this monoplane station.

        Returns a dictionary of booleans.

        """
        d = {'root buildup': self.root_buildup.exists(),
             'spar cap': self.spar_cap.exists(),
             'shear web 1': self.shear_web_1.exists(),
             'shear web 2': self.shear_web_2.exists(),
             'shear web 3': self.shear_web_3.exists(),
             'TE reinforcement': self.TE_reinforcement.exists(),
             'LE panel': self.LE_panel.exists(),
             'aft panel': self.aft_panel.exists(),
             'internal surface': self.internal_surface.exists(),
             'external surface': self.external_surface.exists()}
        return d


class BiplaneStructure:
    """Define the biplane laminate schedule (internal dimensions)."""
    def __init__(self, h_RB, b_SC, h_SC, b_SW1_biax, b_SW1_foam, x2_SW1,
                 b_SW2_biax, b_SW2_foam, x2_SW2, b_SW3_biax, b_SW3_foam,
                 x2_SW3, b_TE_reinf, h_TE_reinf_uniax, h_TE_reinf_foam,
                 h_LE_panel, h_aft_panel, h_int_surf_triax, h_int_surf_resin,
                 h_ext_surf_triax, h_ext_surf_gelcoat, h_RB_u, b_SC_u, h_SC_u,
                 b_SW1_biax_u, b_SW1_foam_u, x2_SW1_u, b_SW2_biax_u,
                 b_SW2_foam_u, x2_SW2_u, b_SW3_biax_u, b_SW3_foam_u, x2_SW3_u,
                 b_TE_reinf_u, h_TE_reinf_uniax_u, h_TE_reinf_foam_u,
                 h_LE_panel_u, h_aft_panel_u, h_int_surf_triax_u,
                 h_int_surf_resin_u, h_ext_surf_triax_u, h_ext_surf_gelcoat_u):
        self.lower_root_buildup = Part(np.nan, h_RB)
        self.lower_spar_cap = Part(b_SC, h_SC)
        self.lower_shear_web_1 = ShearWeb(b_SW1_biax, b_SW1_foam, x2_SW1)
        self.lower_shear_web_2 = ShearWeb(b_SW2_biax, b_SW2_foam, x2_SW2)
        self.lower_shear_web_3 = ShearWeb(b_SW3_biax, b_SW3_foam, x2_SW3)
        self.lower_TE_reinforcement = TE_Reinforcement(b_TE_reinf, h_TE_reinf_uniax, 
                                                 h_TE_reinf_foam)
        self.lower_LE_panel = Part(np.nan, h_LE_panel)
        self.lower_aft_panel = Part(np.nan, h_aft_panel)
        self.lower_internal_surface = InternalSurface(np.nan, h_int_surf_triax, h_int_surf_resin)
        self.lower_external_surface = ExternalSurface(np.nan, h_ext_surf_triax, h_ext_surf_gelcoat)
        self.upper_root_buildup = Part(np.nan, h_RB_u)
        self.upper_spar_cap = Part(b_SC_u, h_SC_u)
        self.upper_shear_web_1 = ShearWeb(b_SW1_biax_u, b_SW1_foam_u, x2_SW1_u)
        self.upper_shear_web_2 = ShearWeb(b_SW2_biax_u, b_SW2_foam_u, x2_SW2_u)
        self.upper_shear_web_3 = ShearWeb(b_SW3_biax_u, b_SW3_foam_u, x2_SW3_u)
        self.upper_TE_reinforcement = TE_Reinforcement(b_TE_reinf_u,
            h_TE_reinf_uniax_u, h_TE_reinf_foam_u)
        self.upper_LE_panel = Part(np.nan, h_LE_panel_u)
        self.upper_aft_panel = Part(np.nan, h_aft_panel_u)
        self.upper_internal_surface = InternalSurface(np.nan,
            h_int_surf_triax_u, h_int_surf_resin_u)
        self.upper_external_surface = ExternalSurface(np.nan,
            h_ext_surf_triax_u, h_ext_surf_gelcoat_u)

    def __str__(self):
        """Returns a string of all the internal dimensions for this structure."""
        s = ''
        s += "****** Lower Airfoil ******\n"
        s += "  --- ROOT BUILDUP ---\n"
        s += "  " + str(self.lower_root_buildup) + '\n'
        s += "  --- SPAR CAP ---\n"
        s += "  " + str(self.lower_spar_cap) + '\n'
        s += "  --- SHEAR WEB 1 ---\n"
        s += "  " + str(self.lower_shear_web_1) + '\n'
        s += "  --- SHEAR WEB 2 ---\n"
        s += "  " + str(self.lower_shear_web_2) + '\n'
        s += "  --- SHEAR WEB 3 ---\n"
        s += "  " + str(self.lower_shear_web_3) + '\n'
        s += "  ---TE REINFORCEMENT ---\n"
        s += "  " + str(self.lower_TE_reinforcement) + '\n'
        s += "  --- LE PANEL ---\n"
        s += "  " + str(self.lower_LE_panel) + '\n'
        s += "  --- AFT PANEL ---\n"
        s += "  " + str(self.lower_aft_panel) + '\n'
        s += "  --- INTERNAL SURFACE ---\n"
        s += "  " + str(self.lower_internal_surface) + '\n'
        s += "  --- EXTERNAL SURFACE ---\n"
        s += "  " + str(self.lower_external_surface) + '\n'
        s += "****** Upper Airfoil ******\n"
        s += "  --- ROOT BUILDUP ---\n"
        s += "  " + str(self.upper_root_buildup) + '\n'
        s += "  --- SPAR CAP ---\n"
        s += "  " + str(self.upper_spar_cap) + '\n'
        s += "  --- SHEAR WEB 1 ---\n"
        s += "  " + str(self.upper_shear_web_1) + '\n'
        s += "  --- SHEAR WEB 2 ---\n"
        s += "  " + str(self.upper_shear_web_2) + '\n'
        s += "  --- SHEAR WEB 3 ---\n"
        s += "  " + str(self.upper_shear_web_3) + '\n'
        s += "  ---TE REINFORCEMENT ---\n"
        s += "  " + str(self.upper_TE_reinforcement) + '\n'
        s += "  --- LE PANEL ---\n"
        s += "  " + str(self.upper_LE_panel) + '\n'
        s += "  --- AFT PANEL ---\n"
        s += "  " + str(self.upper_aft_panel) + '\n'
        s += "  --- INTERNAL SURFACE ---\n"
        s += "  " + str(self.upper_internal_surface) + '\n'
        s += "  --- EXTERNAL SURFACE ---\n"
        s += "  " + str(self.upper_external_surface) + '\n'
        return s

    def which_parts_exist(self):
        """Check which structural parts exist at this biplane station.

        Returns a dictionary of booleans.

        """
        d = {'lower root buildup': self.lower_root_buildup.exists(),
             'lower spar cap': self.lower_spar_cap.exists(),
             'lower shear web 1': self.lower_shear_web_1.exists(),
             'lower shear web 2': self.lower_shear_web_2.exists(),
             'lower shear web 3': self.lower_shear_web_3.exists(),
             'lower TE reinforcement': self.lower_TE_reinforcement.exists(),
             'lower LE panel': self.lower_LE_panel.exists(),
             'lower aft panel': self.lower_aft_panel.exists(),
             'lower internal surface': self.lower_internal_surface.exists(),
             'lower external surface': self.lower_external_surface.exists(),
             'upper root buildup': self.upper_root_buildup.exists(),
             'upper spar cap': self.upper_spar_cap.exists(),
             'upper shear web 1': self.upper_shear_web_1.exists(),
             'upper shear web 2': self.upper_shear_web_2.exists(),
             'upper shear web 3': self.upper_shear_web_3.exists(),
             'upper TE reinforcement': self.upper_TE_reinforcement.exists(),
             'upper LE panel': self.upper_LE_panel.exists(),
             'upper aft panel': self.upper_aft_panel.exists(),
             'upper internal surface': self.upper_internal_surface.exists(),
             'upper external surface': self.upper_external_surface.exists()}
        return d