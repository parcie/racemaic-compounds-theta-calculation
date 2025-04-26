import json
import gemmi
from group_decomposition import decompose
from pathlib import Path

all_sgs_table = gemmi.spacegroup_table()
sgs_by_number = {}
sgs_by_number_dict = {}
for i in range(1,231):
    sgs_this_number = {}
    sgs_by_number[i] = sgs_this_number
    sgs_standard = gemmi.find_spacegroup_by_number(i)
    sgs_standard_name = sgs_standard.hm
    sgs_by_number_dict[i]= {'standard_name':sgs_standard_name}

for sg in all_sgs_table:
    sg_number = sg.number
    sg_name = sg.hm
    sgs_this_number = sgs_by_number[sg_number]
    sgs_this_number[sg_name] = None
    sgs_this_number_dict = sgs_by_number_dict[sg_number]
    sgs_this_number_dict["max_sohncke_subgroup"] = None
    sgs_this_number_dict['isomorphic'] = sgs_this_number
    sgs_by_number_dict[sg_number] = sgs_this_number_dict

current_path = Path(__file__).resolve()
project_root = current_path.parent.parent
file_folder = project_root / 'file_in_process'
file_name = file_folder/'blank_table.json'
with open(file_name,'w') as file1:
    json.dump(sgs_by_number_dict,file1,indent = 4)