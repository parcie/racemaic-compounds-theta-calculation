import gemmi
from itertools import combinations
import numpy as np
from fractions import Fraction

class decompose:

    def __init__(self,sg_name):
        sg = gemmi.SpaceGroup(sg_name)
        self.sg = sg
        self.sg_name = sg_name
        ops_list = sg.operations()

        cen_ops_list = ops_list.cen_ops
        sym_ops_list = ops_list.sym_ops
        cen_generators = self.centering_adding(cen_ops_list)
        all_generators = self.all_generators(cen_generators,sym_ops_list)
        self.all_ops = all_generators
        print('all_ops',all_generators)
        self.sym_ops = sym_ops_list
        self.cen_generators = cen_generators
        crystal_system_str = sg.crystal_system_str()
        self.crystal_system_str = crystal_system_str

        if sg.is_sohncke() is True:
            result = 'already sohncke'
            chirality_retain_positions = all_generators
            chirality_change_positions = None
            max_sohncke_subgroup = sg_name
        else:
            all_info = self.non_sohncke_sgs()
            result = all_info['result']
            chirality_retain_positions = all_info['chirality_retain_positions']
            chirality_change_positions = all_info['chirality_change_positions']
            max_sohncke_subgroup = all_info['max_sohncke_subgroup']

        self.result = result
        self.chirality_retain_positions = chirality_retain_positions
        self.chirality_change_positions = chirality_change_positions
        self.max_sohncke_subgroup = max_sohncke_subgroup
        return


    def non_sohncke_sgs(self):
        sg = self.sg
        sg_name = self.sg_name
        sg_number = sg.number
        sym_ops_list = self.sym_ops
        cen_generators = self.cen_generators
        crystal_system_str = self.crystal_system_str
        all_generators = self.all_ops

        point_group_minus42m = [115,116,117,118,119,120]

        # for space group with only 2-order axis, a shortcut
        if crystal_system_str in ['monoclinic','orthorhombic']:
            points_classified_result = self.orthorhomobic_group_tools()
            chirality_retain_positions = points_classified_result[0]
            chirality_change_positions = points_classified_result[1]
            subgroup_name = self.find_group_by_ops_list(chirality_retain_positions)
            if subgroup_name is None:
                max_sohncke_subgroup_name = None
            else:
                max_sohncke_subgroup_name = subgroup_name
            result = 'success'

        # for the -4m2 point group ###################
        elif sg_number in point_group_minus42m:
            if sg_number in [119, 120]:
                points_classified_result = self.minus42m_pointgroup_method(sym_ops_list,cen_generators)
                chirality_retain = points_classified_result[0]
                chirality_retain_positions = self.all_generators(cen_generators,chirality_retain)
                chirality_change = points_classified_result[1]
                chirality_change_positions = self.all_generators(cen_generators,chirality_change)
                max_sohncke_subgroup_name = 'F222'

            else:
                points_classified_result = self.minus42m_pointgroup_method(all_generators, cen_generators)
                chirality_retain = points_classified_result[0]
                chirality_retain_positions = self.all_generators(cen_generators, chirality_retain)
                chirality_change = points_classified_result[1]
                chirality_change_positions = self.all_generators(cen_generators, chirality_change)
                max_sohncke_subgroup_name = 'C222'

            result = '-4m2 success'
        elif sg_number in [81,82]:
            points_classified_result = self.minus4_pointgroup_method(cen_generators)
            chirality_retain_positions = points_classified_result[0]
            chirality_change_positions = points_classified_result[1]
            if sg_number == 81:
                max_sohncke_subgroup_name = 'P2'
            else:
                max_sohncke_subgroup_name = 'I2'
            result = '-4 success'

        # directly decompose with transformation matrix = E
        else:
            print('do it directly')
            direct_decompose_subgroup_result = self.find_subgroup(all_generators)
            # print(direct_decompose_subgroup_result)

            direct_decompose_subgroup_name = direct_decompose_subgroup_result[0]
            direct_subgroup_ops_list = direct_decompose_subgroup_result[1]

            # using the symmorphoc spacegroup
            if direct_decompose_subgroup_name == 'decompose error':
                print('move into symmorphic tools')
                points_classified = self.symmorphic_group_tool()
                chirality_retain_positions = points_classified[0]
                chirality_change_positions = points_classified[1]
                # when the symmorphic spacegroup cannot be decomposed by T=E
                #######################
                if chirality_retain_positions is None:
                    #try to use point group method
                    print('point group method')
                    points_classified = self.primitive_lattice_method(sym_ops_list,cen_generators)
                    chirality_retain_positions = points_classified[0]
                    chirality_change_positions = points_classified[1]
                    if chirality_retain_positions is None:
                        max_sohncke_subgroup_name == None
                        result = 'failed'
                    else:
                       result = 'through primitive lattice,success'

                else:
                    result = 'success'
                    max_sohncke_subgroup_name = self.find_group_by_ops_list(chirality_retain_positions)

            elif direct_decompose_subgroup_name == 'this spaecegroup is already a sohncke spacegroup':
                chirality_retain_positions = all_generators
                chirality_change_positions = None
                max_sohncke_subgroup_name = sg_name
                result = 'already sohncke'

            elif direct_decompose_subgroup_name == 'cannot decompose by index 2':
                chirality_change_positions = None
                chirality_retain_positions = None
                max_sohncke_subgroup_name = None
                result = 'failed'

            # successfully decomposed by direct decompose with T=E
            else:
                chirality_retain_positions = direct_subgroup_ops_list
                chirality_change_positions = [item for item in all_generators if item not in chirality_retain_positions]
                max_sohncke_subgroup_name = direct_decompose_subgroup_name
                result = 'success'

        chirality_retain_positions = chirality_retain_positions
        chirality_change_positions = chirality_change_positions
        max_sohncke_subgroup = max_sohncke_subgroup_name
        return_dict = {'result':result,'chirality_retain_positions':chirality_retain_positions,
                       'chirality_change_positions':chirality_change_positions,'max_sohncke_subgroup':max_sohncke_subgroup}
        return return_dict

    @staticmethod
    def centering_adding(cen_ops_list):
        cen_generators = []
        for i in range(len(cen_ops_list)):
            cen_ops = cen_ops_list[i]
            if cen_ops == [0,0,0]:
                cen_generators.append(gemmi.Op('x,y,z'))
            elif cen_ops == [0,12,12]:
                cen_generators.append(gemmi.Op('x,y+1/2,z+1/2'))
            elif cen_ops == [12,0,12]:
                cen_generators.append(gemmi.Op('x+1/2,y,z+1/2'))
            elif cen_ops == [12,12,0]:
                cen_generators.append(gemmi.Op('x+1/2,y+1/2,z'))
            elif cen_ops == [12,12,12]:
                cen_generators.append(gemmi.Op('x+1/2,y+1/2,z+1/2'))
            elif cen_ops == [16,8,8]:
                cen_generators.append(gemmi.Op('x+2/3,y+1/3,z+1/3'))
            elif cen_ops == [8,16,16]:
                cen_generators.append(gemmi.Op('x+1/3,y+2/3,z+2/3'))
            else:continue
        return cen_generators

    @staticmethod
    def all_generators(cen_generators,sym_generators):
        all_generators = []
        for i in range(len(cen_generators)):
            center = cen_generators[i]
            for j in range(len(sym_generators)):
                op = sym_generators[j]
                all_op = center*op
                all_generators.append(all_op)
        return all_generators

    @staticmethod
    def find_subgroup(ops_list):
        group_size = len(ops_list)
        ops_list_copy = list(ops_list)

        group_ops = gemmi.GroupOps(ops_list)
        group = gemmi.find_spacegroup_by_ops(group_ops)

        if group_size % 2 != 0:
            return ["cannot decompose by index 2",None]

        elif group.is_sohncke() == True:
            return ["this spaecegroup is already a sohncke spacegroup",None]

        else:
            subgroup_size_minus1 = group_size / 2 - 1
            first_point = gemmi.Op('x,y,z')
            del ops_list_copy[0]
            sub_sg = None

            for subgroup_ops_tuple in combinations(ops_list_copy, int(subgroup_size_minus1)):
                subgroup_ops_list = list(subgroup_ops_tuple)
                subgroup_ops_list.append(first_point)
                subgroup_ops = gemmi.GroupOps(subgroup_ops_list)  # only sym_ops in subgroup
                sub_sg = gemmi.find_spacegroup_by_ops(subgroup_ops)
                if sub_sg is None:
                    continue
                elif sub_sg.is_sohncke() == True:
                    return [sub_sg.hm, subgroup_ops_list]
                else:
                    continue

            if sub_sg is None:
                return ['decompose error',None]
            else:
                return [sub_sg,'not sohncke']


    def symmorphic_group_tool(self):
        ops = self.all_ops
        groupops = gemmi.GroupOps(ops)
        symmorphic_groupops = groupops.derive_symmorphic()
        symmorphic_sym_ops = symmorphic_groupops.sym_ops
        symmorphic_cen = symmorphic_groupops.cen_ops
        symmorphic_cen_ops = self.centering_adding(symmorphic_cen)

        symmorphic_all_ops = self.all_generators(symmorphic_cen_ops,symmorphic_sym_ops)
        translation_ops = self.record_translation(ops,symmorphic_all_ops)

        symmorphic_subgroup_name, symmorphic_subgroup_ops_list = self.find_subgroup(symmorphic_all_ops)

        error_list = ['decompose_error','this spaecegroup is already a sohncke spacegroup','cannot decompose by index 2']
        if symmorphic_subgroup_name in error_list:
            return [None,None]
        else:
            subgroup_ops_list = self.re_add_translation(symmorphic_all_ops,translation_ops,symmorphic_subgroup_ops_list)
            chirality_retain_positions = subgroup_ops_list

            chirality_change_positions = [item for item in ops if item not in chirality_retain_positions]
            return [chirality_retain_positions, chirality_change_positions]

    def sym_ops_classification_method(self,sym_ops_list):
        subgroup_name = self.find_subgroup(sym_ops_list)

        subgroup = gemmi.find_spacegroup_by_name(subgroup_name)
        subgroup_ops = subgroup.operations()

        if subgroup is not None:
            sym_ops_retain_chirality = subgroup_ops.sym_ops
            sym_ops_change_chirality = [item for item in sym_ops_list if item not in sym_ops_retain_chirality]

        else:
            sym_ops_retain_chirality = None
            sym_ops_change_chirality = None
        return [chirality_retain_positions, chirality_change_positions]


    @staticmethod
    def find_group_by_ops_list(ops_list):
        groupops = gemmi.GroupOps(ops_list)
        spacegroup = gemmi.find_spacegroup_by_ops(groupops)
        if spacegroup is None:
            return None
        else:
            return spacegroup.hm

    @staticmethod
    def operation_division(dividend_op,divisor_op):
        '''
        op1*op2=op3, op3 is dividend, op1 is divisor ,get the op2
        '''
        if dividend_op.triplet() == divisor_op.triplet():
            return gemmi.Op('x,y,z')
        else:
            dividend_seitz = dividend_op.seitz()
            divisor_seitz = divisor_op.seitz()

            dividend_seitz_array = np.array([[float(element) for element in row] for row in dividend_seitz])

            divisor_seitz_array = np.array([[float(element) for element in row] for row in divisor_seitz])
            divisor_inv = np.linalg.inv(divisor_seitz_array)
            remainder_seitz_array = np.dot(dividend_seitz_array,divisor_inv)
            remainder_seitz_list = list(remainder_seitz_array)
            remainder_op = gemmi.seitz_to_op(remainder_seitz_list)
            return remainder_op

    def record_translation(self,ops,symmorphic_ops):
        translation_ops = []

        for i in range(len(ops)):
            op_i = ops[i]
            symmorphic_op_i = symmorphic_ops[i]
            removed_translation_i = self.operation_division(op_i,symmorphic_op_i)
            translation_ops.append(removed_translation_i)
        return translation_ops

    @staticmethod
    def re_add_translation(symmorphic_all_ops,translation_ops,symmorphic_subgroup_ops_list):
        symmorphic_all_ops_dict = {}
        for i in range(len(symmorphic_all_ops)):
            symmorphic_all_ops_dict[i] = symmorphic_all_ops[i]

        translation_ops_dict = {}
        for i in range(len(translation_ops)):
            translation_ops_dict[i] = translation_ops[i]

        re_add_translation_ops = []
        for i in range(len(symmorphic_subgroup_ops_list)):
            subgroup_op = symmorphic_subgroup_ops_list[i]
            key_i = next((k for k ,v in symmorphic_all_ops_dict.items() if v == subgroup_op),None)
            translation_i = translation_ops_dict[key_i]
            group_op = translation_i * subgroup_op
            re_add_translation_ops.append(group_op)

        symmetry_retain_ops = re_add_translation_ops
        return symmetry_retain_ops

    @staticmethod
    def symbol_cal(string1):
        symbol=list(string1)
        length=len(symbol)
        minus=0
        for i in range(length):
            if symbol[i]=="-":
                if 119<ord(symbol[i+1])<123:
                    minus=minus+1
        return minus

    def orthorhomobic_group_tools(self):
        operation_list = self.all_ops
        print('ops_list',operation_list)
        chirality_retain_op_list = []
        chirality_change_op_list = []

        for i in range(len(operation_list)):
            temp_op = operation_list[i]
            temp_str = temp_op.triplet()
            minus = self.symbol_cal(temp_str)
            if minus in [0,2]:
                chirality_retain_op_list.append(temp_op)
            else:
                chirality_change_op_list.append(temp_op)
        print('chirality_retain',chirality_retain_op_list,'chriality_change',chirality_change_op_list)
        return [chirality_retain_op_list,chirality_change_op_list]


    def minus42m_pointgroup_method(self,ops,cen_generators):
        groupops = gemmi.GroupOps(self.all_ops)
        symmorphic_groupops = groupops.derive_symmorphic()
        symmorphic_sym_ops = symmorphic_groupops.sym_ops
        symmorphic_cen = symmorphic_groupops.cen_ops
        symmorphic_cen_ops = self.centering_adding(symmorphic_cen)

        symmorphic_all_ops = self.all_generators(symmorphic_cen_ops,symmorphic_sym_ops)
        translation_ops = self.record_translation(ops, symmorphic_all_ops)

        symmorphic_subgroup_sym_ops = [gemmi.Op('x,y,z'),gemmi.Op('-x,-y,z'),gemmi.Op('y,x,-z'),gemmi.Op('-y,-x,-z')]
        subgroup_ops_list = self.re_add_translation(symmorphic_all_ops, translation_ops, symmorphic_subgroup_sym_ops)

        chirality_retain_positions = subgroup_ops_list
        chirality_change_positions = [item for item in ops if item not in chirality_retain_positions]
        return [chirality_retain_positions, chirality_change_positions]

    def minus4_pointgroup_method(self,cen_generators):
        ops = self.all_ops
        subgroup_sym_ops = [gemmi.Op('x,y,z'),gemmi.Op('-x,-y,z')]
        subgroup_ops_list = self.all_generators(cen_generators,subgroup_sym_ops)
        chirality_retain_positions = subgroup_ops_list
        chirality_change_positions = [item for item in ops if item not in chirality_retain_positions]
        return [chirality_retain_positions, chirality_change_positions]

    # furtherly turn the symmorphic sg to its corresponding primitive lattice
    def primitive_lattice_method(self,sym_ops_list,cen_generators):
        groupops = gemmi.GroupOps(sym_ops_list)
        symmorphic_groupops = groupops.derive_symmorphic()
        symmorphic_sym_ops = symmorphic_groupops.sym_ops
        translation_ops = self.record_translation(sym_ops_list, symmorphic_all_ops)

        subgroup_result = self.find_subgroup(symmorphic_sym_ops)
        subgroup_name = subgroup_result[0]
        error_type = ['decompose_error','this spaecegroup is already a sohncke spacegroup','cannot decompose by index 2']
        if subgroup_name in error_type:
            return [None,None]
        else:
            subgroup_ops_list = subgroup_result[1]

            chirality_retain_positions_pre = self. self.re_add_translation(symmorphic_all_ops, translation_ops, )
            chirality_retain_positions = self.all_generators(cen_generators,chirality_retain_positions_pre)
            chirality_change_positions_primitive = [item for item in sym_ops_list if item not in chirality_retain_positions]
            chirality_change_positions = self.all_generators(cen_generators,chirality_change_positions_primitive)
            return [chirality_retain_positions, chirality_change_positions]




