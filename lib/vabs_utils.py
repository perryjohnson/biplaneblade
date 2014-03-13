"""
A module to translate data from a 2D cross-section grid file (from TrueGrid)
into a VABS input file.

Functions to update
-------------------
self._write_element_layers()
self._write_layers()
self._write_materials()

Author: Perry Roth-Johnson
Last updated: March 12, 2014

"""


import numpy as np
import abaqus_utils2 as au
reload(au)


class VabsInputFile:
    """The VabsInputFile class contains methods for translating data from an
    AbaqusGrid object (from TrueGrid) into a VABS input file.

    Usage:
    import lib.vabs_utils as vu
    f = vu.VabsInputFile(
    vabs_filename='sandia_blade/mesh_stn01.vabs',
    grid=g,
    debug_flag=True)

    """
    def __init__(self, vabs_filename, grid, debug_flag=False):
        self.vabs_filename = vabs_filename
        self.grid = grid
        self.set_flags()
        self._write_input_file(debug_flag=debug_flag)

    # Rewrite this! User should set flags from the __init__() method.
    def set_flags(self):
        self.format_flag = 1
        self.number_of_layers = 1
        self.Timoshenko_flag = 1
        self.recover_flag = 0
        self.thermal_flag = 0
        self.curve_flag = 0
        self.oblique_flag = 0
        self.trapeze_flag = 0
        self.Vlasov_flag = 0
        if self.curve_flag:
            self.k1 = 0
            self.k2 = 0
            self.k3 = 0

    def _write_header(self):
        flag1_fmt = '{0:d} {1:d}\n'
        flag2_comments = '# Timoshenko_flag  recover_flag  thermal_flag'
        flag2_fmt = '{0:d} {1:d} {2:d}    ' + flag2_comments + '\n'
        flag3_comments = '# curve_flag  oblique_flag  trapeze_flag  Vlasov_flag'
        flag3_fmt = '{0:d} {1:d} {2:d} {3:d}  ' + flag3_comments + '\n\n'
        num_comments = '# nnode  nelem  nmate'
        num_fmt = '{0:d} {1:d} {2:d}   ' + num_comments + '\n\n'
        self.vabs_file.write(flag1_fmt.format(self.format_flag,
                                              self.number_of_layers))
        self.vabs_file.write(flag2_fmt.format(self.Timoshenko_flag,
                                              self.recover_flag,
                                              self.thermal_flag))
        self.vabs_file.write(flag3_fmt.format(self.curve_flag,
                                              self.oblique_flag,
                                              self.trapeze_flag,
                                              self.Vlasov_flag))
        self.vabs_file.write(num_fmt.format(self.grid.number_of_nodes,
                                            self.grid.number_of_elements,
                                            1))
                                            # self.grid.number_of_materials))

    def _write_nodes(self):
        n = str(len(str(self.grid.number_of_nodes)))
        fmt = '{0:>'+n+'d}' + 5*' ' + '{1:> 10.8f}' + 2*' ' + '{2:> 10.8f}\n'
        for node in self.grid.list_of_nodes:
            self.vabs_file.write(fmt.format(node.node_num,
                                            node.x2,
                                            node.x3))
        self.vabs_file.write('\n')

    def _write_element_connectivity(self):
        nn = str(len(str(self.grid.number_of_nodes)))
        nfmt = '>' + nn + 'd'
        ne = str(len(str(self.grid.number_of_elements)))
        efmt = '>' + ne + 'd'
        fmt = '{0:'+efmt+'}     {1:'+nfmt+'} {2:'+nfmt+'} {3:'+nfmt+'} {4:'+nfmt+'} {5:'+nfmt+'} {6:'+nfmt+'} {7:'+nfmt+'} {8:'+nfmt+'} {9:'+nfmt+'}\n'
        for element in self.grid.list_of_elements:
            self.vabs_file.write(fmt.format(element.elem_num,
                                            element.node1.node_num,
                                            element.node2.node_num,
                                            element.node3.node_num,
                                            element.node4.node_num,
                                            element.node5.node_num,
                                            element.node6.node_num,
                                            element.node7.node_num,
                                            element.node8.node_num,
                                            0))
                                            # element['node9']))
        self.vabs_file.write('\n')

    def _write_element_layers(self):
        n = str(len(str(self.grid.number_of_elements)))
        fmt = '{0:>'+n+'d}' + 5*' ' + '{1:d} {2:>7.2f}\n'
        for element in self.grid.list_of_elements:
            self.vabs_file.write(fmt.format(element.elem_num,1,element.theta1))
        self.vabs_file.write('\n')

    def _write_layers(self):
        self.vabs_file.write('1    1     0.00\n\n')

    def _write_materials(self):
        # write material properties  # HARDCODED FOR FOAM! CHANGE LATER!
        mat_id = 1
        orth = 0
        name = 'Foam'
        E = 2.56000e+08  # Young's modulus = 0.256 GPa
        nu = 3.00000e-01  # Poisson's ratio = 0.3
        rho = 2.00000e+02  # density = 200 kg/m^3
        isotropic_ln1 = '{0:<4d}{1:<4d}# {2}\n'
        isotropic_ln2 = '{3:>11.5e}   {4:>11.5e}\n'
        isotropic_ln3 = '{5:>11.5e}\n'
        isotropic_fmt = isotropic_ln1 + isotropic_ln2 + isotropic_ln3
        self.vabs_file.write(isotropic_fmt.format(mat_id, orth, name,
                                                  E, nu,
                                                  rho))
        # self.vabs_file.write('1   0   # Foam\n')
        # self.vabs_file.write('2.56000e+08   3.00000e-01\n')
        # self.vabs_file.write('2.00000e+02\n')

    def _write_input_file(self, debug_flag=False):
        """Writes the VABS input file.

        This non-public method is automatically run when a new VabsInputFile
        instance is created.

        """
        if debug_flag:
            print 'VABS input file: ' + self.vabs_filename
        # open the input file
        self.vabs_file = open(self.vabs_filename, 'w+')
        # write to the input file
        self._write_header()
        self._write_nodes()
        self._write_element_connectivity()
        self._write_element_layers()
        self._write_layers()
        self._write_materials()
        # close the input file
        self.vabs_file.close()
