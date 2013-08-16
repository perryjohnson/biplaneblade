import scripts.scale_airfoil_thickness as sat


airfoil_dir = 'biplane_flap-sym_no-stagger/airfoils'
airfoil_list = []

airfoil_list.append({
    'original_path' : 'SNL-100m-Transition51.txt',
    'original_tc_ratio' : 0.51,
    'new_path' : 'SNL-100m-Transition51_53pt2.txt',
    'new_tc_ratio' : 0.532
    })
airfoil_list.append({
    'original_path' : 'DU99-W-405.txt',
    'original_tc_ratio' : 0.405,
    'new_path' : 'DU99-W-405_40pt8.txt',
    'new_tc_ratio' : 0.408
    })
airfoil_list.append({
    'original_path' : 'DU91-W2-250_26.txt',
    'original_tc_ratio' : 0.26,
    'new_path' : 'DU91-W2-250_26_25pt5.txt',
    'new_tc_ratio' : 0.255
    })
airfoil_list.append({
    'original_path' : 'DU93-W-210_23.txt',
    'original_tc_ratio' : 0.23,
    'new_path' : 'DU93-W-210_23_23pt5.txt',
    'new_tc_ratio' : 0.235
    })
airfoil_list.append({
    'original_path' : 'DU93-W-210.txt',
    'original_tc_ratio' : 0.21,
    'new_path' : 'DU93-W-210_21pt75.txt',
    'new_tc_ratio' : 0.2175
    })
airfoil_list.append({
    'original_path' : 'DU93-W-210.txt',
    'original_tc_ratio' : 0.21,
    'new_path' : 'DU93-W-210_20pt25.txt',
    'new_tc_ratio' : 0.2025
    })
airfoil_list.append({
    'original_path' : 'DU93-W-210.txt',
    'original_tc_ratio' : 0.21,
    'new_path' : 'DU93-W-210_20pt4.txt',
    'new_tc_ratio' : 0.204
    })

for airfoil in airfoil_list:
    sat.write_and_plot(
        airfoil_dir,
        airfoil['original_path'],
        airfoil['original_tc_ratio'],
        airfoil['new_path'],
        airfoil['new_tc_ratio']
        )
