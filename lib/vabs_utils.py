"""A module to translate data from an AbaqusGrid object (from TrueGrid) into a 
VABS input file.

Functions to update
-------------------
self._write_element_layers()
self._write_layers()
self._write_materials()

Author: Perry Roth-Johnson
Last updated: March 12, 2014

"""


import numpy as np
import pandas as pd
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
    material_filename='sandia_blade/materials.csv',
    layer_filename='sandia_blade/layers.csv',
    debug_flag=True)

    """
    def __init__(self, vabs_filename, grid, material_filename, layer_filename,
        debug_flag=False,
        flags={
            'format'           : 1,
            'Timoshenko'       : 1,
            'recover'          : 0,
            'thermal'          : 0,
            'curve'            : 0,
            'k1'               : 0,
            'k2'               : 0,
            'k3'               : 0,
            'oblique'          : 0,
            'trapeze'          : 0,
            'Vlasov'           : 0
        }):
        self.vabs_filename = vabs_filename
        self.grid = grid
        # read material file to determine the number of materials
        self.material_filename = material_filename
        self._mf = pd.read_csv(self.material_filename)
        self.number_of_materials = len(self._mf)
        # read layer file to determine the number of layers
        self.layer_filename = layer_filename
        self._lf = pd.read_csv(self.layer_filename)
        self.number_of_layers = len(self._lf)
        self.flags = flags
        self._write_input_file(debug_flag=debug_flag)

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

    def _write_header(self):
        flag1_fmt = '{0:d} {1:d}\n'
        flag2_comments = '# Timoshenko_flag  recover_flag  thermal_flag'
        flag2_fmt = '{0:d} {1:d} {2:d}    ' + flag2_comments + '\n'
        flag3_comments = '# curve_flag  oblique_flag  trapeze_flag  Vlasov_flag'
        flag3_fmt = '{0:d} {1:d} {2:d} {3:d}  ' + flag3_comments + '\n\n'
        flag3_alt_fmt1 = '{0:d} {1:d} {2:d} {3:d}  ' + flag3_comments + '\n'
        flag3_alt_fmt2 = '{4:6.8f} {5:6.8f} {6:6.8f}\n\n'
        flag3_alt_fmt = flag3_alt_fmt1 + flag3_alt_fmt2
        num_comments = '# nnode  nelem  nmate'
        num_fmt = '{0:d} {1:d} {2:d}   ' + num_comments + '\n\n'
        self.vabs_file.write(flag1_fmt.format(self.flags['format'],
                                              self.number_of_layers))
        self.vabs_file.write(flag2_fmt.format(self.flags['Timoshenko'],
                                              self.flags['recover'],
                                              self.flags['thermal']))
        if self.flags['curve'] == 1:
            self.vabs_file.write(flag3_alt_fmt.format(self.flags['curve'],
                                                      self.flags['oblique'],
                                                      self.flags['trapeze'],
                                                      self.flags['Vlasov'],
                                                      self.flags['k1'],
                                                      self.flags['k2'],
                                                      self.flags['k3']))
        else:
            self.vabs_file.write(flag3_fmt.format(self.flags['curve'],
                                                  self.flags['oblique'],
                                                  self.flags['trapeze'],
                                                  self.flags['Vlasov']))
        self.vabs_file.write(num_fmt.format(self.grid.number_of_nodes,
                                            self.grid.number_of_elements,
                                            self.number_of_materials))

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
                                            element.node9.node_num))
        self.vabs_file.write('\n')

    def _write_element_layers(self):
        n = str(len(str(self.grid.number_of_elements)))
        fmt = '{0:>'+n+'d}' + 5*' ' + '{1:d} {2:>7.2f}\n'
        for element in self.grid.list_of_elements:
            self.vabs_file.write(fmt.format(
                element.elem_num,
                element.layer_num,
                element.theta1)
            )
        self.vabs_file.write('\n')

    def _write_layers(self):
        fmt = '{0}        {1}    {2:3.2f}   # {3}\n'
        for l in range(self.number_of_layers):
            self.vabs_file.write(fmt.format(
                self._lf['layer number'][l],
                self._lf['material number'][l],
                self._lf['layup orientation angle'][l],
                self._lf['layer name'][l]))
        self.vabs_file.write('\n')

    def _write_materials(self):
        # text formatting for isotropic materials
        isotropic_ln1 = '{0:<4d}{1:<4d}# {2}\n'
        isotropic_ln2 = '{3:>11.5e}   {4:>11.5e}\n'
        isotropic_ln3 = '{5:>11.5e}\n\n'
        isotropic_fmt = isotropic_ln1 + isotropic_ln2 + isotropic_ln3
        # text formatting for orthotropic materials
        orthotropic_ln1 = '{0:<4d}{1:<4d}# {2}\n'
        orthotropic_ln2 = '{3:>11.5e}   {4:>11.5e}   {5:>11.5e}\n'
        orthotropic_ln3 = '{6:>11.5e}   {7:>11.5e}   {8:>11.5e}\n'
        orthotropic_ln4 = '{9:>11.5e}   {10:>11.5e}   {11:>11.5e}\n'
        orthotropic_ln5 = '{12:>11.5e}\n\n'
        orthotropic_fmt = (orthotropic_ln1 + orthotropic_ln2 +
            orthotropic_ln3 + orthotropic_ln4 + orthotropic_ln5)
        for m in range(self.number_of_materials):
            if self._mf['type'][m] == 'isotropic':
                orth = 0
                line = isotropic_fmt.format(
                    self._mf['number'][m], orth, self._mf['name'][m], 
                    self._mf['E1'][m], self._mf['nu12'][m],
                    self._mf['rho'][m])
                self.vabs_file.write(line)
            elif self._mf['type'][m] == 'orthotropic':
                orth = 1
                line = orthotropic_fmt.format(
                    self._mf['number'][m], orth, self._mf['name'][m], 
                    self._mf['E1'][m], self._mf['E2'][m], self._mf['E3'][m],
                    self._mf['G12'][m], self._mf['G13'][m], self._mf['G23'][m],
                    self._mf['nu12'][m], self._mf['nu13'][m], self._mf['nu23'][m],
                    self._mf['rho'][m])
                self.vabs_file.write(line)
            else:
                raise Warning("The material type {0} is undefined!".format(material['type']))
