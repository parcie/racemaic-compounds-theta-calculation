import gemmi
import ccdc
from ccdc import io
from ccdc import crystal
from communal_tools import math_method

class chiral_line:

    def __init__(self,ccdc_code,same_ops_list,point_center):
        '''
        fix the type of same_ops_list element
        :param ccdc_code:
        :param same_ops_list: list of ops,ops are in gemmi.Op , these general positions are of same chirality
        :param point_center:
        '''
        self.ccdc_code = ccdc_code
        csd_reader = io.MoleculeReader('CSD')
        crys = csd_reader.crystal(ccdc_code)
        sg_name = crys.spacegroup_symbol
        sg = gemmi.SpaceGroup(sg_name)
        sg_number = sg.number

        self.crystal = crys
        self.first_point = point_center
        sg_number_and_setting = crys.spacegroup_number_and_setting
        spacegroup_setting = sg_number_and_setting[1]  # judge whether is in standard stype
        if spacegroup_setting == 1:
            self.standard_sg = True
        else:
            self.standard_sg = False

        list_minus4_axis = [81,82,111,112,113,114,115,116,117,118,119,120,121,122]  #with -4 axis
        points_number_incell = len(same_ops_list)
        self.points_number_incell = points_number_incell
        group_ops = gemmi.GroupOps(same_ops_list)
        same_ops_list_sym = group_ops.sym_ops

        min_distance = None
        axis_type = None
        if sg_number == 2:
            line_direction = self.GetLinedirection3(ccdc_code)  # line direction in crystallographic coordinate
            basic_op = None
            packing_number = None
        elif sg_number in [6, 7, 8, 9]:
            list_info = self.GetLinedirection2(ccdc_code)
            line_direction = list_info[0]
            basic_op = list_info[1]
            packing_number = 2
        elif sg_number in list_minus4_axis: # with -4 axis
            list_info = self.method_for_minus_4(same_ops_list_sym)
            line_direction = list_info[0]
            basic_op = list_info[1]
            packing_number = 2
            axis_type = 'rotoinversion'

        else:
            list_info = self.GeneralMethod(same_ops_list_sym)
            line_direction = list_info['line_direction']
            basic_op = list_info['basis_op']
            min_distance = list_info['min_distance']
            packing_number = list_info['max_axis_order']
            axis_type = list_info['axis_type']

        self.line_direction = line_direction
        self.basic_op = basic_op
        self.packing_number = packing_number
        self.axis_type = axis_type
        self.min_distance = min_distance

    @property
    def PackingNumber(self):
        return self.packing_number

    @property
    def Direction(self):
        return self.line_direction

    @property
    def DirectionOrtho(self):
        frac_direction = self.line_direction
        ccdc_code = self.ccdc_code
        csd_reader = io.MoleculeReader('CSD')
        crys = csd_reader.crystal(ccdc_code)
        a = crys.cell_lengths.a
        b = crys.cell_lengths.b
        c = crys.cell_lengths.c
        alpha = crys.cell_angles.alpha
        beta = crys.cell_angles.beta
        gamma = crys.cell_angles.gamma
        ortho_direction = math_method.Orthogonalization(frac_direction,a,b,c,alpha,beta,gamma)
        return ortho_direction

    @property
    def PackingOp(self):
        return self.basic_op

    @property
    def AxisOrder(self):
        return self.packing_number

    @property
    def PackingPairDistance(self):
        return self.min_distance

    @property
    def AxidType(self):
        return self.axis_type

    def GeneralMethod(self,same_ops_list):
        highest_axis_info = self.FindHighestOrderAxisDirection(same_ops_list, False)
        is_unique = highest_axis_info['unique']
        point_center = self.first_point
        if is_unique is True:
            line_direction = highest_axis_info['line_direction']
            basic_op = highest_axis_info['basis_op']
            type = highest_axis_info['axis_type']
            packing_number = highest_axis_info['max_axis_order']
            min_distance = None

        elif is_unique is False:
            line_directions = highest_axis_info['line_direction']
            basic_ops = highest_axis_info['basis_op']
            type = highest_axis_info['axis_type']
            packing_number = highest_axis_info['max_axis_order']
            list_info = self.DistanceMethod(basic_ops,line_directions,point_center)
            line_direction = list_info['line_direction']
            basic_op = list_info['basis_op']
            min_distance = list_info['min_distance']
        else:
            line_direction = None
            basic_op = None
            min_distance = None
            packing_number = None
            type = None

        list_info = {'line_direction':line_direction,'basis_op':basic_op,
                       'min_distance':min_distance,"max_axis_order":packing_number,
                       'axis_type':type}
        return list_info


    @staticmethod
    def DistanceMethod(basic_ops,line_directions,molecule_center):
        op_E = gemmi.Op('x,y,z')
        start_point = op_E.apply_to_xyz(molecule_center)

        distances = []
        for i in range(len(basic_ops)):
            opi_string = basic_ops[i]
            opi = gemmi.Op(opi_string)
            pointi = opi.apply_to_xyz(molecule_center)
            distance_i = math_method.PPdistance(start_point,pointi)
            distances.append(distance_i)

        min_index = distances.index(min(distances))
        line_direction = line_directions[min_index]
        basic_op = basic_ops[min_index]
        min_distance = min(distances)
        return_info = {'line_direction':line_direction,'basis_op':basic_op,'min_distance':min_distance}
        return return_info

    def method_for_minus_4(self,same_ops_list):
        '''
        direction correspond to -4axis
        2-fold rotation

        :param same_ops_list:
        :return:
        '''
        crystal = self.crystal
        is_standard_type = self.standard_sg
        if is_standard_type is True:
            line_direction = (0,0,1)
            axis_type = "2-fold rotation"
            basic_op = self.FindBasisop(axis_type,line_direction,same_ops_list)
            return [line_direction,basic_op]
            #packed in the 2-fold rotation
        else:
            operators = crystal.symmetry_operators
            highest_axis_info = self.FindHighestOrderAxisDirection(operators,True)
            line_direction = highest_axis_info['line_direction']
            axis_type = "2-fold rotation"
            basic_op = self.FindBasisop(axis_type, line_direction, same_ops_list)
            return [line_direction, basic_op]


    def FindHighestOrderAxisDirection(self,operators,choose_rotoinverse: bool=False):
        '''

        :param operators:operators can be whether all or just same ops list, element in operators must be strings of operation's representation
        :param choose_rotoinverse:
        :return: info_to-return['unique']:whether the highest order ond priority axis
        '''
        crystal = self.crystal

        max_order =1
        for i in range(len(operators)):
            opi_gemmi = operators[i]
            opi = opi_gemmi.triplet()
            description = crystal.symmetry_operator_description(opi)
            order_string = description[0]
            if order_string.isdigit() is True:
                order = int(order_string)
                if order>max_order:
                    max_order = order
                else:
                    continue

        line_direction_list = []
        highest_order_ops = [] # axis list for all of the same order,may include different type
        for i in range(len(operators)):
            op1_gemmi = operators[i]
            op1 = op1_gemmi.triplet()
            description = crystal.symmetry_operator_description(op1)
            axis_order_str = description[0]
            if axis_order_str.isdigit() == True:
                if int(axis_order_str) == max_order:
                    direction_i = self.GetOpDirection(op1)
                    line_direction_list.append(direction_i)
                    highest_order_ops.append(op1)

        if choose_rotoinverse is True:  # choose roto-inversion as the line direction
            line_direction = line_direction_list[0]
            basisop = None
            info_to_return = {'unique':True,'line_direction':line_direction,'axis_type':"rotoinversion","basis_op":basisop,"max_axis_order":max_order}
            return info_to_return

        screw_axis = {} # op in string :direction in tuple
        rotation_axis = {}

        for i in range(len(highest_order_ops)):
            op1 = highest_order_ops[i]
            description = crystal.symmetry_operator_description(op1)
            if 'screw' in description:
                this_line_direction = line_direction_list[i]
                screw_axis[op1] = this_line_direction

            elif 'rotation' in description:
                this_line_direction = line_direction_list[i]
                rotation_axis[op1] = this_line_direction

        if len(screw_axis) == 1:
            basis_op_in_list = list(screw_axis.keys())
            line_direction_in_list = list(screw_axis.values())
            info_to_return = {'unique':True,'line_direction':line_direction_in_list[0],'axis_type':"screw",
                                "basis_op":basis_op_in_list[0],"max_axis_order":max_order}
            return info_to_return
        elif len(screw_axis) > 1:
            basis_op_in_list = list(screw_axis.keys())
            line_direction_in_list = list(screw_axis.values())
            info_to_return = {'unique':False,'line_direction':line_direction_in_list,'axis_type':"screw",
                              "basis_op":basis_op_in_list,"max_axis_order":max_order}
            return info_to_return
        elif len(rotation_axis) == 1:
            basis_op_in_list = list(rotation_axis.keys())
            line_direction_in_list = list(rotation_axis.values())
            info_to_return = {'unique': True, 'line_direction': line_direction_in_list[0], 'axis_type': "rotation",
                              "basis_op": basis_op_in_list[0], "max_axis_order": max_order}
            return info_to_return
        elif len(rotation_axis) > 1:
            basis_op_in_list = list(rotation_axis.keys())
            line_direction_in_list = list(rotation_axis.values())
            info_to_return = {'unique':False,'line_direction':line_direction_in_list,'axis_type':"screw",
                              "basis_op":basis_op_in_list,"max_axis_order":max_order}
            return info_to_return

        else:
            return "no-rotational symmetry"


    def FindBasisop(self,axis_type:type==str,line_direction,ops):
        '''
        find an operation that fits the axis type and the axis's direction
        return th operation in string
        '''

        def get_screw_component(string):
            screw_component = string[-9:]
            screw_square = 0
            for i in range(len(screw_component)):
                chari = screw_component[i]
                if chari.isdigit() is True:
                    numberi = float(chari)
                    screw_square = screw_square+numberi*numberi
                else:
                    continue
            screw_length = screw_square**0.5
            return screw_length

        crys = self.crystal
        line_direction_list = list(line_direction)
        line_direction_string = str(line_direction_list)

        target_ops = []
        for i in range(len(ops)):
            op = ops[i]
            op_string = op.triplet()
            description = crys.symmetry_operator_description(op_string)
            if axis_type in description:
                if line_direction_string in description:
                    target_ops.append(op)
                else:
                    continue
            else:
                continue
        if len(target_ops) == 1:
            return target_ops[0]
        elif screw in axis_type:
            screw_component_max = 1
            target_op = target_ops[0]
            for i in range(len(target_ops)):
                opi = target_ops[i]
                descriptioni = crys.symmetry_operator_description(opi)
                screw_component = get_screw_component(descriptioni)
                if screw_component<screw_component_max:
                    screw_component_max = screw_component
                    target_op = opi
                else:
                    continue
            return target_op
        else:
            return target_ops[0]


    def GetOpDirection(self, op):
        '''
        :param op: the string of the operation
        :return: if it relates with any rotational axis(include rotoinversion and screw), return the axis' direction
        '''
        # get axis direction from op description
        crys = self.crystal
        description = crys.symmetry_operator_description(op)

        index_direction = description.find('axis with direction [')  # str_len 21
        start_index = index_direction + 21
        str_direction = description[start_index:start_index + 7]  # get the str
        dir_num_list = {}
        a = 0
        for i in range(len(str_direction)):
            t_str = str_direction[i]
            if t_str.isnumeric() == True:
                dir_num_list[a] = int(str_direction[i])
                a = a + 1
            else:
                continue
        direction = tuple(dir_num_list.values())
        return direction


    @staticmethod
    def GetLinedirection3(ccdc_code):
        '''
        this is for the chiral line direction in P-1 space group
        :return: line direction in tuple, which must be parallel with crystallographic coordinate
        '''
        csd_reader = io.MoleculeReader('CSD')
        crys = csd_reader.crystal(ccdc_code)

        a = crys.cell_lengths.a
        b = crys.cell_lengths.b
        c = crys.cell_lengths.c

        y = min(a, b, c)
        if y == a:
            line_direction = (1, 0, 0)
        elif y == b:
            line_direction = (0, 1, 0)
        elif y == c:
            line_direction = (0, 0, 1)
        else:
            line_direction = None

        return line_direction

    @staticmethod
    def GetLinedirection2(ccdc_code):
        '''
        this function is for the space groups without axis order higher than 1, but hve glide plane
        the space group number include N0(6,7,8,9)
        :param line_direction: the direction vector(3 dimensional) in tuple,like (0,0,1), which represent the direction perpendicular to the glide plane
        :param base_op: the corresponding operation that represent glide operation
        '''

        csd_reader = io.MoleculeReader('CSD')
        crys = csd_reader.crystal(ccdc_code)
        operators = crys.symmetry_operators
        dir_num_list = []
        base_op = ''
        str_direction = ''
        for i in range(len(operators)):
            op1 = operators[i]
            description = crys.symmetry_operator_description(op1)
            index1 = description.find('Glide plane perpendicular to [')
            if index1 != -1:
                start_index = index1 + 30
                str_direction = description[start_index:start_index + 7]
                base_op = op1
                break
            else:
                continue

        for i in range(len(str_direction)):
            t_str = str_direction[i]
            if t_str.isnumeric() is True:
                dir_num_list.append(t_str)
        direction = tuple(dir_num_list)
        direction_int_tuple = tuple(map(int, direction))

        return [direction_int_tuple, base_op]  # return direction in tuple and base_op in string


