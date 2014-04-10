"""A module to translate the outputs from VABS to inputs for DYMORE.

For example, if VABS was used to calculate the mass and stiffness matrices for
a series of cross-sections in a beam, this module can construct the DYMORE
input file for a beam with those cross-sectional properties.

Author: Perry Roth-Johnson
Last modified: August 1, 2013

"""

import numpy as np


def readFile(filestr):
    """
    Read a file into memory.

    Parameters
    ----------
    filestr : <string>
        The filename of the file to be read into memory.

    Returns
    -------
    fileLines : <list of strings>
        A list of strings, each of which represents one line in the file.
    """

    f = open(filestr, 'r')  # read the file
    fileLines = f.readlines()  # parse all characters in file into one long string
    return fileLines


def defineRegularExpressions():
    """
    Define regular expressions to search for mass and stiffness matrices in a 
    VABS-formatted output file.

    Parameters
    ----------
    <none>

    Returns
    -------
    cm_x2_pat : <pattern object>
        The regex pattern for the x2-coordinate of the center of mass.
    cm_x3_pat : <pattern object>
        The regex pattern for the x3-coordinate of the center of mass.
    mpus_pat : <pattern object>
        The regex pattern for the mass per unit span.
    i1_pat : <pattern object>
        The regex pattern for the moment of inertia about the x1-axis.
    i2_pat : <pattern object>
        The regex pattern for the moment of inertia about the x2-axis.
    i3_pat : <pattern object>
        The regex pattern for the moment of inertia about the x3-axis.
    K_head_pat : <pattern object>
        The regex pattern for the Timoshenko Stiffness Matrix header.
    """

    import re
    cm_x2_pat = re.compile(r'\s\sXm2 =.+')
    cm_x3_pat = re.compile(r'\s\sXm3 =.+')
    mpus_pat = re.compile(r'\sMass Per Unit Span\s+=.+')
    i1_pat = re.compile(r'\sMass Moments of Intertia about x1 axis\s+=.+')
    i2_pat = re.compile(r'\sMass Moments of Intertia about x2 axis\s+=.+')
    i3_pat = re.compile(r'\sMass Moments of Intertia about x3 axis\s+=.+')
    K_head_pat = re.compile(r'\sTimoshenko Stiffness Matrix.+')
    
    return (cm_x2_pat, cm_x3_pat, mpus_pat, i1_pat, i2_pat, i3_pat, K_head_pat)


def pullMKmatrices(MKlines, print_flag=False):
    """
    Pull the numerical values of the mass and stiffness matrices from the VABS 
    output file (*.dat.K).

    Parameters
    ----------
    MKlines : <list of strings>
        The contents of the VABS output file (*.dat.K).
        Each string is a line from the VABS output file.
    print_flag : <logical>
        Set to True to print out extra debugging information to the screen.

    Returns
    -------
    cm_x2 : <float>
        The x2-coordinate of the center of mass.
    cm_x3 : <float>
        The x3-coordinate of the center of mass.
    mpus : <float>
        The mass per unit span.
    i1 : <float>
        The moment of inertia about the x1-axis.
    i2 : <float>
        The moment of inertia about the x2-axis.
    i3 : <float>
        The moment of inertia about the x3-axis.
    K : <np.array>
        The Timoshenko stiffness matrix.
    """

    (cm_x2_pat, cm_x3_pat, mpus_pat, i1_pat, i2_pat, i3_pat, K_head_pat) = defineRegularExpressions()

    for i in range(len(MKlines)):
        cm_x2_match = cm_x2_pat.match(MKlines[i])
        cm_x3_match = cm_x3_pat.match(MKlines[i])
        mpus_match = mpus_pat.match(MKlines[i])
        i1_match = i1_pat.match(MKlines[i])
        i2_match = i2_pat.match(MKlines[i])
        i3_match = i3_pat.match(MKlines[i])
        K_head_match = K_head_pat.match(MKlines[i])

        if cm_x2_match:
            cm_x2 = float(MKlines[i].replace('Xm2 =', ''))
        elif cm_x3_match:
            cm_x3 = float(MKlines[i].replace('Xm3 =', ''))
        elif mpus_match:
            mpus = float(MKlines[i].replace('Mass Per Unit Span                     =', ''))
        elif i1_match:
            i1 = float(MKlines[i].replace('Mass Moments of Intertia about x1 axis =', ''))
        elif i2_match:
            i2 = float(MKlines[i].replace('Mass Moments of Intertia about x2 axis =', ''))
        elif i3_match:
            i3 = float(MKlines[i].replace('Mass Moments of Intertia about x3 axis =', ''))
        elif K_head_match:
            K_head_line = i

    K_row_1 = np.fromstring(MKlines[K_head_line+3], sep=' ')
    K_row_2 = np.fromstring(MKlines[K_head_line+4], sep=' ')
    K_row_3 = np.fromstring(MKlines[K_head_line+5], sep=' ')
    K_row_4 = np.fromstring(MKlines[K_head_line+6], sep=' ')
    K_row_5 = np.fromstring(MKlines[K_head_line+7], sep=' ')
    K_row_6 = np.fromstring(MKlines[K_head_line+8], sep=' ')
    K = np.vstack((K_row_1,K_row_2,K_row_3,K_row_4,K_row_5,K_row_6))

    if print_flag:
        print "center of mass, x2 = " + str(cm_x2)
        print "center of mass, x3 = " + str(cm_x3)
        print "mass per unit span = " + str(mpus)
        print "moment of inertia about x1 axis = " + str(i1)
        print "moment of inertia about x2 axis = " + str(i2)
        print "moment of inertia about x3 axis = " + str(i3)
        print ""
        print "Stiffness matrix:"
        print K

    return (cm_x2, cm_x3, mpus, i1, i2, i3, K)


def makeFile(dymoreFileName):
    """
    Make a temporary file on the hard disk to store DYMORE-formatted data.

    Parameters
    ----------
    dymoreFileName : <string>
        The filename of the file to write to.
       
    Returns
    -------
    tempFile : <file object>
        A file handle to the temporary file.
    """

    tempFile = open(dymoreFileName, 'w+')
    return tempFile


def writeDymoreMK(f, CoordType, coord, cm_x2, cm_x3, mpus, i1, i2, i3, K):
    """
    Description.

    Parameters
    ----------
    f : <file object>
        The file handle that data will be written to.
    CoordType : <string>
        Acceptable values are: 'ETA_COORDINATE',
                               'CURVILINEAR_COORDINATE', or
                               'AXIAL_COORDINATE'
    coord : <float>
        The spanwise coordinate of this cross-section.
        This coordinate should match the CoordType specified above.
    cm_x2 : <float>
        The x2-coordinate of the center of mass.
    cm_x3 : <float>
        The x3-coordinate of the center of mass.
    mpus : <float>
        The mass per unit span.
    i1 : <float>
        The moment of inertia about the x1-axis.
    i2 : <float>
        The moment of inertia about the x2-axis.
    i3 : <float>
        The moment of inertia about the x3-axis.
    K : <np.array>
        The Timoshenko stiffness matrix.

    Returns
    -------
    <none>

    Example output (written to a file)
    ----------------------------------
    @ETA_COORDINATE {0.00000e+00} {
      @STIFFNESS_MATRIX { 7.6443255182E+09,   -3.5444961981E-04,   -1.5092432335E-03,    3.3599749794E+06,    2.7710007447E-01,    4.1602501550E-02,
                                               2.8284702841E+08,   -2.8863166160E+01,   -4.5930836014E-01,   -3.3517886643E+05,    3.4162114776E-03,
                                                                    3.5606703330E+08,   -4.0749872012E-01,    3.6079611429E-02,   -4.2577508629E+04,
                                                                                         8.8773955810E+08,    1.8897378940E-03,    8.3869473951E-04,
                                                                                                              4.5282893600E+10,   -5.7686739280E-02,
                                                                                                                                   2.2625281359E+09}
      @MASS_PER_UNIT_SPAN {7.1224712000E+02}
      @MOMENTS_OF_INERTIA {3.9569408290E+03,
                           3.6961203640E+03,
                           2.6082046495E+02}
      @CENTRE_OF_MASS_LOCATION {-1.6834618673E-17,
                                -1.1472480873E-16}
    }
    """

    tab = '  '
    if CoordType == 'ETA_COORDINATE':
        f.write(tab*2 + '@ETA_COORDINATE {' + ('%11.5e' % coord) + '} {\n')
    elif CoordType == 'CURVILINEAR_COORDINATE':
        # f.write(tab*2 + '@CURVILINEAR_COORDINATE {' + ('%11.5e' % coord) + '} {\n')
        raise NotImplementedError("CURVILINEAR_COORDINATE feature is not yet supported.")
    elif CoordType == 'AXIAL_COORDINATE':
        f.write(tab*2 + '@AXIAL_COORDINATE {' + ('%11.5e' % coord) + '} {\n')
    f.write(tab*3 +   '@STIFFNESS_MATRIX {' + ('%17.10e' % K[0,0]) + ',' + ('%20.10e' % K[0,1]) + ',' + ('%20.10e' % K[0,2]) + ',' + ('%20.10e' % K[0,3]) + ',' + ('%20.10e' % K[0,4]) + ',' + ('%20.10e' % K[0,5]) + ',' + '\n')
    f.write(tab*3 + ' '*37 +                                               ('%20.10e' % K[1,1]) + ',' + ('%20.10e' % K[1,2]) + ',' + ('%20.10e' % K[1,3]) + ',' + ('%20.10e' % K[1,4]) + ',' + ('%20.10e' % K[1,5]) + ',' + '\n')
    f.write(tab*3 + ' '*(37+21*1) +                                                                     ('%20.10e' % K[2,2]) + ',' + ('%20.10e' % K[2,3]) + ',' + ('%20.10e' % K[2,4]) + ',' + ('%20.10e' % K[2,5]) + ',' + '\n')
    f.write(tab*3 + ' '*(37+21*2) +                                                                                                  ('%20.10e' % K[3,3]) + ',' + ('%20.10e' % K[3,4]) + ',' + ('%20.10e' % K[3,5]) + ',' + '\n')
    f.write(tab*3 + ' '*(37+21*3) +                                                                                                                               ('%20.10e' % K[4,4]) + ',' + ('%20.10e' % K[4,5]) + ',' + '\n')
    f.write(tab*3 + ' '*(37+21*4) +                                                                                                                                                            ('%20.10e' % K[5,5]) + '}' + '\n')
    f.write(tab*3 +   '@MASS_PER_UNIT_SPAN {' + ('%17.10e' % mpus) + '}\n')
    f.write(tab*3 +   '@MOMENTS_OF_INERTIA {' + ('%17.10e' % i1) + ',\n')
    f.write(tab*3 + ' '*21 +                    ('%17.10e' % i2) + ',\n')
    f.write(tab*3 + ' '*21 +                    ('%17.10e' % i3) + '}\n')
    f.write(tab*3 +   '@CENTRE_OF_MASS_LOCATION {' + ('%17.10e' % cm_x2) + ',\n')
    f.write(tab*3 + ' '*26 +                         ('%17.10e' % cm_x3) + '}\n')
    f.write(tab*2 + '}\n')
    f.write(tab*2 + '\n')

    return


def writeMKmatrices(DYMOREfileHandle, vabsMKfilepath, station_data, CoordType='ETA_COORDINATE', debug_flag=False):
    """
    Write the mass and stiffness matrices for one cross-section to a file.

    Parameters
    ----------
    DYMOREfileHandle : <file object>
        The file handle that data will be written to.
    vabsMKfilepath : <string>
        The path to the VABS file that contains the mass and stiffness matrices\
        for this cross-section.
        Example: 'VABS/M_and_K_matrices/spar_station_04.dat.K'
    station_data : <dictionary>
        The dictionary of layup parameters for this spar station.
    CoordType : <string>
        Acceptable values are: 'ETA_COORDINATE',
                               'CURVILINEAR_COORDINATE', or
                               'AXIAL_COORDINATE'
    debug_flag : <logical>
        Set to True to print out extra debugging information to the screen.

    Returns
    -------
    <none>
    """

    if debug_flag:
        print "CoordType = " + CoordType
        if CoordType == 'ETA_COORDINATE':
            print "eta =", station_data['eta']
        elif CoordType == 'CURVILINEAR_COORDINATE':
            # print "s =", station_data['s']
            raise NotImplementedError("CURVILINEAR_COORDINATE feature is not yet supported.")
        elif CoordType == 'AXIAL_COORDINATE':
            # print "x1 =", station_data['x1']
            raise NotImplementedError("AXIAL_COORDINATE feature is not yet supported.")
    MKlines = readFile(vabsMKfilepath)
    (cm_x2, cm_x3, mpus, i1, i2, i3, K) = pullMKmatrices(MKlines, print_flag=debug_flag)
    if CoordType == 'ETA_COORDINATE':
        writeDymoreMK(DYMOREfileHandle, CoordType, station_data['eta'], cm_x2, cm_x3, mpus, i1, i2, i3, K)
    elif CoordType == 'CURVILINEAR_COORDINATE':
        # writeDymoreMK(DYMOREfileHandle, CoordType, station_data['s'], cm_x2, cm_x3, mpus, i1, i2, i3, K)
        raise NotImplementedError("CURVILINEAR_COORDINATE feature is not yet supported.")
    elif CoordType == 'AXIAL_COORDINATE':
        # writeDymoreMK(DYMOREfileHandle, CoordType, station_data['x1'], cm_x2, cm_x3, mpus, i1, i2, i3, K)
        raise NotImplementedError("AXIAL_COORDINATE feature is not yet supported.")

    return


def formatComments(comments):
    """
    Format the comments for a Dymore code block. Maximum comment length is 5 lines of 120 characters each

    Parameters
    ----------
    comments : <string>
        The original, unformatted comment.

    Returns
    -------
    comments : <string>
        The comment, reformatted for Dymore.
    """

    if len(comments) > 120:
        if len(comments) < 120*2:    # break the comment up into 2 lines
            line1 = comments[0:120]
            line2 = comments[120:]
            comments = line1 + '/' + line2
        elif len(comments) < 120*3:  # break the comment up into 3 lines
            line1 = comments[0:120]
            line2 = comments[120:240]
            line3 = comments[240:]
            comments = line1 + '/' + line2 + '/' + line3
        elif len(comments) < 120*4:  # break the comment up into 4 lines
            line1 = comments[0:120]
            line2 = comments[120:240]
            line3 = comments[240:360]
            line4 = comments[360:]
            comments = line1 + '/' + line2 + '/' + line3 + '/' + line4
        elif len(comments) < 120*5:  # break the comment up into 5 lines
            line1 = comments[0:120]
            line2 = comments[120:240]
            line3 = comments[240:360]
            line4 = comments[360:480]
            line5 = comments[480:]
            comments = line1 + '/' + line2 + '/' + line3 + '/' + line4 + '/' + line5
        else:
            print "***WARNING*** the comment was too long and was truncated to 600 characters!"
            line1 = comments[0:120]
            line2 = comments[120:240]
            line3 = comments[240:360]
            line4 = comments[360:480]
            line5 = comments[480:600]
            comments = line1 + '/' + line2 + '/' + line3 + '/' + line4 + '/' + line5

    return comments


if __name__ == '__main__':  #run this code if DYMOREutilities is called directly from the command line (good for debugging)
    print "hello world"
    
