"""Delete leftover files and paths in each blade_path"""


import os
import shutil


def clean_blade_path(blade_path):
    """Removes all files except blade defintion and airfoil path."""
    clean_list = os.listdir(blade_path)
    clean_list.remove('blade_definition.csv')
    try:
        clean_list.remove('materials.csv')
    except ValueError:
        pass
    clean_list.remove('airfoils')
    for item in clean_list:
        item = os.path.join(blade_path,item)
        if os.path.exists(item):
            if os.path.isfile(item):
                os.remove(item)
                print " [Deleted file] {0}".format(item)
            else:
                shutil.rmtree(item)
                print " [Deleted path] {0}".format(item)

clean_blade_path('sandia_blade')
clean_blade_path('biplane_flap-sym_no-stagger')
clean_blade_path('biplane_flap-sym_stagger')
clean_blade_path('biplane_flap-asym_no-stagger')
clean_blade_path('biplane_flap-asym_stagger')
