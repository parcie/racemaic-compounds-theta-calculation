import gemmi
import json
from group_decomposition import decompose
import threading
from pathlib import Path

def run_decompose_with_timeout(group_name,time_limit):
    result = [None]

    def wrapper():
        try:
            result[0] = decompose(group_name)
        except Exception as e:
            print(f"An error occurred in thread for group {group_name}: {e}")
        finally:
            print(f"Thread for group {group_name} has finished.")
    calculation_thread = threading.Thread(target=wrapper)
    calculation_thread.start()
    calculation_thread.join(timeout=time_limit)
    if calculation_thread.is_alive():
        print(f"Time limit exceeded for group {group_name}. Skipping this group.")
        calculation_thread.join(timeout=time_limit)  # Ensure thread is cleaned up properly
        return None
    else:
        return result[0]

current_path = Path(__file__).resolve()
project_root = current_path.parent.parent

file_folder = project_root / 'file_in_process'
blank_table = file_folder / 'blank_table.json'
with open(blank_table,'r+') as file1:
    sg_info = json.load(file1)

time_limit = 600

sg_table = {}
for i in range(1,231):
    key = f'{i}'
    sg_this_number = sg_info[key]
    sg_this_number_copy = sg_this_number.copy()
    sg_standard = sg_this_number['standard_name']

    isomorphic_groups_dict = sg_this_number['isomorphic']
    isomorphic_groups = list(isomorphic_groups_dict.keys())
    print(isomorphic_groups)
    isomorphic_groups_dict = {}

    for j in range(len(isomorphic_groups)):
        a_group_name = isomorphic_groups[j]

        a_group_object = run_decompose_with_timeout(a_group_name, time_limit)
        if a_group_object is None:
            continue

        same_chirality_ops = a_group_object.chirality_retain_positions
        same_chirality_ops_string = []
        if same_chirality_ops is None:
            same_chirality_ops_string = []
        else:

            for y in range(len(same_chirality_ops)):
                ops_j = same_chirality_ops[y]
                ops_j_string = ops_j.triplet()
                same_chirality_ops_string.append(ops_j_string)

            different_chirality_ops = a_group_object.chirality_change_positions
            different_chirality_ops_string = []
            if different_chirality_ops is None:
                different_chirality_ops_string=[]
            else:
                for x in range(len(different_chirality_ops)):
                    ops_x = different_chirality_ops[x]
                    ops_x_string = ops_x.triplet()
                    different_chirality_ops_string.append(ops_x_string)
        points_info = {'chirality_retain_positions': same_chirality_ops_string,
                       'chirality_change_positions': different_chirality_ops_string}

        isomorphic_groups_dict[a_group_name] = points_info
        if a_group_name == sg_standard:
            sg_this_number['max_sohncke_subgroup'] = a_group_object.max_sohncke_subgroup
        else:
            continue

    sg_this_number_copy['isomorphic'] =isomorphic_groups_dict
    sg_table[i] = sg_this_number_copy

    file2_path = project_root / 'file_in_process' /'decomposition_result.json'
    with open('/Users/chenfangyi/Documents/はたらく/final_cal/sg_test/430_table.json', 'w') as file2:
        file2.truncate(0)
        json.dump(sg_table,file2,indent=4)

print('over')
file3_path = project_root / 'file_in_process' /'decomposition_result2.json'
with open(file3_path, 'w') as file3:
        json.dump(sg_table,file3,indent=4)












