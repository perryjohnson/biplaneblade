"""Wrapper file to run VABS on input files for a station cross-section.

Author: Perry Roth-Johnson
Last modified: April 25, 2014

"""


import os


# -----------------------------------------------
# update these parameters!
# list_of_station_nums = range(1,40+1)
list_of_station_nums = [1,2,3,4,5,6,7,8,9,10,25,28,29,30,31,32,33,34,35,36,37,38,39,40]
# -----------------------------------------------

for station_num in list_of_station_nums:
    stn_str = 'stn{0:02d}'.format(station_num)
    path_to_VABS_exe = 'D:\\Programs\\VABS\\vabs_3-7'
    relative_path_to_VABS_input_file = 'biplane_blade\\{0}'.format(stn_str)
    VABS_input_filename = 'mesh_{0}.vabs'.format(stn_str)

    absolute_path_to_VABS_input_file = os.path.join(os.getcwd(), 
        relative_path_to_VABS_input_file, VABS_input_filename)

    if not os.path.exists(path_to_VABS_exe):
        raise ValueError("The path '{0}' to the VABS executable does not exist!".format(path_to_VABS_exe))
    if not os.path.exists(absolute_path_to_VABS_input_file):
        raise ValueError("The path '{0}' to the VABS input file does not exist!".format(absolute_path_to_VABS_input_file))

    cwd = os.getcwd()
    os.chdir(path_to_VABS_exe)
    print "RUNNING VABS....."
    vabs_command = 'VABSIII.exe ' + absolute_path_to_VABS_input_file
    os.system(vabs_command)

    if os.path.exists(absolute_path_to_VABS_input_file + '.K'):
        print "generated {0}.K with mass and stiffness matrices!".format(VABS_input_filename)
    os.chdir(cwd)
