"""Delete leftover station paths in ./sandia_blade and in ./biplane_blade"""


import os
import shutil


# --- sandia blade ------------------------------------------------------------
# generate list of station paths
list_of_station_paths = []
root_path = r".\sandia_blade"
for station in range(1,34+1):
    d = 'stn{0:02d}'.format(station)
    p = os.path.join(root_path, d)
    list_of_station_paths.append(p)

# delete each station path
for station_path in list_of_station_paths:
    if os.path.exists(station_path):
        # delete the station path, even if it has contents
        shutil.rmtree(station_path)
        print " [Deleted station path] {0}".format(station_path)

# delete 'selected_cross-sections.png'
scs_file = os.path.join(root_path, 'selected_cross-sections.png')
if os.path.isfile(scs_file):
    os.remove(scs_file)
    print " [Deleted selected cross-sections] {0}".format(scs_file)

# --- biplane blade 1 ---------------------------------------------------------
# generate list of station paths
list_of_station_paths = []
root_path = r".\biplane_flap-sym_no-stagger"
for station in range(1,34+1):
    d = 'stn{0:02d}'.format(station)
    p = os.path.join(root_path, d)
    list_of_station_paths.append(p)

# delete each station path
for station_path in list_of_station_paths:
    if os.path.exists(station_path):
        # delete the station path, even if it has contents
        shutil.rmtree(station_path)
        print " [Deleted station path] {0}".format(station_path)

# delete 'selected_cross-sections.png'
scs_file = os.path.join(root_path, 'selected_cross-sections.png')
if os.path.isfile(scs_file):
    os.remove(scs_file)
    print " [Deleted selected cross-sections] {0}".format(scs_file)

# --- biplane blade 2 ---------------------------------------------------------
# generate list of station paths
list_of_station_paths = []
root_path = r".\biplane_flap-sym_stagger"
for station in range(1,34+1):
    d = 'stn{0:02d}'.format(station)
    p = os.path.join(root_path, d)
    list_of_station_paths.append(p)

# delete each station path
for station_path in list_of_station_paths:
    if os.path.exists(station_path):
        # delete the station path, even if it has contents
        shutil.rmtree(station_path)
        print " [Deleted station path] {0}".format(station_path)

# delete 'selected_cross-sections.png'
scs_file = os.path.join(root_path, 'selected_cross-sections.png')
if os.path.isfile(scs_file):
    os.remove(scs_file)
    print " [Deleted selected cross-sections] {0}".format(scs_file)

# --- biplane blade 3 ---------------------------------------------------------
# generate list of station paths
list_of_station_paths = []
root_path = r".\biplane_flap-asym_no-stagger"
for station in range(1,34+1):
    d = 'stn{0:02d}'.format(station)
    p = os.path.join(root_path, d)
    list_of_station_paths.append(p)

# delete each station path
for station_path in list_of_station_paths:
    if os.path.exists(station_path):
        # delete the station path, even if it has contents
        shutil.rmtree(station_path)
        print " [Deleted station path] {0}".format(station_path)

# delete 'selected_cross-sections.png'
scs_file = os.path.join(root_path, 'selected_cross-sections.png')
if os.path.isfile(scs_file):
    os.remove(scs_file)
    print " [Deleted selected cross-sections] {0}".format(scs_file)

# --- biplane blade 4 ---------------------------------------------------------
# generate list of station paths
list_of_station_paths = []
root_path = r".\biplane_flap-asym_stagger"
for station in range(1,34+1):
    d = 'stn{0:02d}'.format(station)
    p = os.path.join(root_path, d)
    list_of_station_paths.append(p)

# delete each station path
for station_path in list_of_station_paths:
    if os.path.exists(station_path):
        # delete the station path, even if it has contents
        shutil.rmtree(station_path)
        print " [Deleted station path] {0}".format(station_path)

# delete 'selected_cross-sections.png'
scs_file = os.path.join(root_path, 'selected_cross-sections.png')
if os.path.isfile(scs_file):
    os.remove(scs_file)
    print " [Deleted selected cross-sections] {0}".format(scs_file)

# clean up unused variables
del list_of_station_paths
del root_path
del d
del p
del station
del station_path
del os
del shutil