"""A module for organizing material property data for a blade.

Author: Perry Roth-Johnson
Last updated: October 11, 2013

"""


class _Material:
    """Define a material."""
    logfile_name = 'material.log'
    number_of_materials = 0
    def __init__(self, name):
        _Material.number_of_materials += 1
        self.material_num = _Material.number_of_materials
        self.name = name
    def __del__(self):
        _Material.number_of_materials -= 1
        print " Material deleted, and now _Material.number_of_materials = {0}".format(_Material.number_of_materials)


class IsotropicMaterial(_Material):
    """Define the properties of an isotropic material."""
    def __init__(self, name, E, nu, rho):
        _Material.__init__(self, name)
        self.type = 'isotropic'
        self.E = float(E)       # Young's modulus
        self.nu = float(nu)     # Poisson's ratio
        self.rho = float(rho)   # density
        # calculate the shear modulus from E and nu, ref:
        # wikipedia.org/wiki/Young%27s_modulus#Relation_among_elastic_constants
        self.G = float(E)/(2.0*(1.0+float(nu)))
        self.logf = open(_Material.logfile_name, "a")
        self.logf.write("............(Created material #{0})............\n".format(self.material_num))
        print " Created material #{0}, {1}".format(self.material_num, self.name)
        self.logf.write(str(self) + '\n')
        self.logf.flush()
        self.logf.close()
    def __str__(self):
        return """Isotropic Material ---
Name: {0}
E:    {1:8.2e} (Pa)
G:    {2:8.2e} (Pa)
nu:   {3:4.2f} (-)
rho:  {4:6.1f} (kg/m^3)""".format(self.name, self.E, self.G, self.nu, self.rho)


class OrthotropicMaterial(_Material):
    """Define the properties of an orthotropic material."""
    def __init__(self, name, E1, E2, E3, G12, G13, G23, nu12, nu13, nu23, rho):
        _Material.__init__(self, name)
        self.type = 'orthotropic'
        # Young's modulus
        self.E1 = float(E1)
        self.E2 = float(E2)
        self.E3 = float(E3)
        # shear modulus
        self.G12 = float(G12)
        self.G13 = float(G13)
        self.G23 = float(G23)
        # Poisson's ratio
        self.nu12 = float(nu12)
        self.nu13 = float(nu13)
        self.nu23 = float(nu23)
        # density
        self.rho = float(rho)
        # apply the reciprocity relations for the other Poisson's ratios
        self.nu21 = float(nu12) * (float(E2)/float(E1))
        self.nu31 = float(nu13) * (float(E3)/float(E1))
        self.nu32 = float(nu23) * (float(E3)/float(E2))
        self.logf = open(_Material.logfile_name, "a")
        self.logf.write("............(Created material #{0})............\n".format(self.material_num))
        print " Created material #{0}, {1}".format(self.material_num, self.name)
        self.logf.write(str(self) + '\n')
        self.logf.flush()
        self.logf.close()
    def __str__(self):
        return """Orthotropic Material ---
Name:  {0}
E1:    {1:8.2e} (Pa)
E2:    {2:8.2e} (Pa)
E3:    {3:8.2e} (Pa)
G12:   {4:8.2e} (Pa)
G13:   {5:8.2e} (Pa)
G23:   {6:8.2e} (Pa)
nu12:  {7:4.2f} (-)
nu13:  {8:4.2f} (-)
nu23:  {9:4.2f} (-)
nu21:  {10:4.2f} (-)
nu31:  {11:4.2f} (-)
nu32:  {12:4.2f} (-)
rho:   {13:6.1f} (kg/m^3)""".format(self.name, self.E1, self.E2, self.E3,
    self.G12, self.G13, self.G23, self.nu12, self.nu13, self.nu23, self.nu21,
    self.nu31, self.nu32, self.rho)
