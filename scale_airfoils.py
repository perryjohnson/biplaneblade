import scripts.scale_airfoil_thickness as sat


airfoil_dir = 'biplane_flap-sym_no-stagger/airfoils'
airfoil_list = []

airfoil_list.append({
    'original_path' : 'DU91-W2-250_26.txt',
    'original_tc_ratio' : 0.26,
    'new_path' : 'DU91-W2-250_25pt1.txt',
    'new_tc_ratio' : 0.251
    })
airfoil_list.append({
    'original_path' : 'DU91-W2-250_26.txt',
    'original_tc_ratio' : 0.26,
    'new_path' : 'DU91-W2-250_24pt1.txt',
    'new_tc_ratio' : 0.241
    })

for airfoil in airfoil_list:
    sat.write_and_plot(
        airfoil_dir,
        airfoil['original_path'],
        airfoil['original_tc_ratio'],
        airfoil['new_path'],
        airfoil['new_tc_ratio']
        )
