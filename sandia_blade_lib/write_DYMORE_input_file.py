"""Write code blocks for a DYMORE input file.

Usage
-----
start an IPython (qt)console with the pylab flag:
$ ipython qtconsole --pylab
or
$ ipython --pylab
Then, from the prompt, run this script:
|> %run write_DYMORE_input_file
Finally, open the newly created '*.dat' files, and manually copy the contents
into your DYMORE input file.

Author: Perry Roth-Johnson
Last updated: April 10, 2014

"""


import lib.blade as bl
import lib.dymore_utils as du


# load the Sandia blade
m = bl.MonoplaneBlade('Sandia blade SNL100-00', 'sandia_blade')

# orientation distribution ----------------------------------------------------
# open a new file
ODfile = 'sandia_blade_OD.dat'
f = open(ODfile, 'w')
# write the header
f.write("""@ORIENTATION_DISTRIBUTION_DEFINITION {
  @ORIENTATION_DISTRIBUTION_NAME { orientationBlade } {
    @ORIENTATION_DEFINITION_TYPE { TWIST_ANGLE }
    @COORDINATE_TYPE { ETA_COORDINATE }

""")

print " Writing @ORIENTATION_DISTRIBUTION_DEFINITION {...}"
for station in m.list_of_stations:
    print "   Station #{0}...".format(station.station_num)
    # write the ETA_COORDINATE and TWIST_ANGLE in the DYMORE input file format
    eta = station.coords.x1/m.list_of_stations[-1].coords.x1
    fmt1 = '    @ETA_COORDINATE {{{0:8.5f}}}\n'.format(eta)
    fmt2 = '    @TWIST_ANGLE    {{{0:8.5f}}}\n\n'.format(station.airfoil.twist)
    f.write(fmt1+fmt2)
# write the footer
f.write("""  }
}
""")
# close the text file
f.close()
print " See '{0}' for the @ORIENTATION_DISTRIBUTION_DEFINITION block.".format(ODfile)


# mass and stiffness matrices -------------------------------------------------
# open a new file
MKfile = 'sandia_blade_MK.dat'
f = open(MKfile, 'w')
# write the header
f.write("""@BEAM_PROPERTY_DEFINITION {
  @BEAM_PROPERTY_NAME {propBlade} {
    @PROPERTY_DEFINITION_TYPE { 6X6_MATRICES }
    @COORDINATE_TYPE { ETA_COORDINATE }

""")

print " Writing @BEAM_PROPERTY_DEFINITION {...}"
for station in m.list_of_stations:
    print "   Station #{0}...".format(station.station_num)
    # get the VABS output filename
    vabsMK = 'sandia_blade/stn{0:02d}/mesh_stn{0:02d}.vabs.K'.format(
        station.station_num)
    # write the mass and stiffness matrices in the DYMORE input file format
    du.writeMKmatrices(f, vabsMK,
        {'eta': station.coords.x1/m.list_of_stations[-1].coords.x1},
        CoordType='ETA_COORDINATE', debug_flag=False)

# write the footer
f.write("""  }
}
""")
# close the text file
f.close()
print " See '{0}' for the @BEAM_PROPERTY_DEFINITION block.".format(MKfile)
