from communal_tools import math_method
from line_direction import chiral_line
from molecule import molecule
from point_packing import packing_points
from model_forming2_mnselecting import folded_plane
from model_forming2_mnselecting import overlapped_folded_plane
from fractions import Fraction
import json
from ccdc import io
import gemmi
import math
import numpy as np

class a_racemic_sample:

    def __init__(self,ccdc_code,sg_info_table):
        self.sg_info_table = sg_info_table
        csd_reader = io.MoleculeReader('CSD')
        crys = csd_reader.crystal(ccdc_code)

        sg_name = crys.spacegroup_symbol

        sg = gemmi.SpaceGroup(sg_name)
        self.sg_name = sg.hm
        base_op = sg.basisop
        self.base_op = base_op

        sg_setting = crys.spacegroup_number_and_setting
        self.sg_number_and_setting = sg_setting
        sg_number = sg_setting[0]

        points_set = self.general_position_classification()
        chirality_retain_positions = points_set[0]
        chirality_change_positions = points_set[1]
        self.chirality_retain_positions = chirality_retain_positions
        self.chirality_change_positions = chirality_change_positions

        a = crys.cell_lengths.a
        b = crys.cell_lengths.b
        c = crys.cell_lengths.c
        alpha = crys.cell_angles.alpha
        beta = crys.cell_angles.beta
        gamma = crys.cell_angles.gamma

        self.cell_parameters = {'a': a, 'b': b, 'c': c, 'alpha': alpha, 'beta': beta, 'gamma': gamma}

        first_molecule = molecule(ccdc_code)
        first_molecule_center_frac = first_molecule.first_fractional_mol_center
        first_molecule_center_ortho = first_molecule.first_orthogonal_mol_center
        z_prime = first_molecule.z_prime


        try:
            first_chiral_line = chiral_line(ccdc_code,chirality_retain_positions,first_molecule_center_ortho)
        except Exception as e:
            print(f'there is an error in conforming chiral lines,error type {e}')

        line_direction_frac = first_chiral_line.Direction
        self.line_direction = line_direction_frac
        packing_basic_op = first_chiral_line.PackingOp
        self.chiral_line = first_chiral_line  # a chiral line object
        packing_number = first_chiral_line.PackingNumber


        try:
            packing1 = packing_points(sg_number,
                                      packing_number,
                                      packing_basic_op,
                                      chirality_retain_positions)
        except Exception as e:
            print(f'there is an error in points packing,error type {e}')

        same_chirality_packed_ops = packing1.after_packing
        different_chirality_packed_ops = self.inverse_centers(same_chirality_packed_ops)

        try:
            folded_plane1 = folded_plane(sg_number,first_molecule_center_frac,same_chirality_packed_ops,different_chirality_packed_ops,
                                     line_direction_frac,a,b,c,alpha,beta,gamma)
        except Exception as e:
            print(f'there is an error in conforming folded plane,error type {e}')

        m = folded_plane1.m
        m_is_zero = math.isclose(0, m, rel_tol=1e-5)
        if m_is_zero:  # why overlapped
            self.is_overlap = True

            folded_plane2 = overlapped_folded_plane(sg_number,first_molecule_center_frac,same_chirality_packed_ops,
                                                    line_direction_frac,a,b,c,alpha,beta,gamma)
            self.theta = folded_plane2.theta
            self.sigma = folded_plane2.sigma
            self.m = folded_plane2.m
            self.n = folded_plane2.n
            self.x = folded_plane2.x
            self.y = folded_plane2.y
            self.k = folded_plane2.k
            self.x_vector = folded_plane2.x_vector
            self.folded_plane_object = folded_plane2
            print(ccdc_code)
        else:
            self.is_overlap = False
            self.theta = folded_plane1.theta
            self.sigma = folded_plane1.sigma
            self.m = folded_plane1.m
            self.n = folded_plane1.n
            self.x = folded_plane1.x
            self.y = folded_plane1.y
            self.k = folded_plane1.k
            self.x_vector = folded_plane1.x_vector
            self.folded_plane_object = folded_plane1
            print(ccdc_code)
            return


    @property
    def its_position_sets(self):
        return self.general_position_classification()

    @property
    def its_folded_plane(self):
        return self.folded_plane_object

    def general_position_classification(self):
        '''
        :return:
        '''
        sg_info_table = self.sg_info_table
        sg_number_and_setting = self.sg_number_and_setting
        sg_number = sg_number_and_setting[0]
        sg_setting = sg_number_and_setting[1]
        sg_name = self.sg_name

        key1 = str(sg_number)

        sgs_by_numebr = sg_info_table[key1]
        isomorphic_sgs = sgs_by_numebr['isomorphic']

        isomorphic_sgs_names = list(isomorphic_sgs.keys())

        if sg_name not in isomorphic_sgs_names:
            return "this space group symbol is too irregular to be recognized "
        else:
            this_sg = isomorphic_sgs[sg_name]
            chirality_retian_positions_str = this_sg['chirality_retain_positions']
            chirality_change_positions_str = this_sg['chirality_change_positions']
            chirality_retian_positions = []
            chirality_change_positions = []
            for i in range(len(chirality_change_positions_str)):
                op_retain_str = chirality_retian_positions_str[i]
                op_change_str = chirality_change_positions_str[i]
                op_retain = gemmi.Op(op_retain_str)
                op_change = gemmi.Op(op_change_str)
                chirality_retian_positions.append(op_retain)
                chirality_change_positions.append(op_change)

            return [chirality_retian_positions,chirality_change_positions]

    def inverse_centers(self,ops_centers):
        def RotTransSeparation(seitz):
            rot_part = seitz[0:3, 0:3]
            x_trans = seitz[0][3]
            y_trans = seitz[1][3]
            z_trans = seitz[2][3]
            trans_matrix = np.array([[x_trans,0,0],[0,y_trans,0],[0,0,z_trans]])
            trans_list = [x_trans,y_trans,z_trans]
            return rot_part,trans_matrix,trans_list

        def TransSeitzToList(matrix):
            x_trans = matrix[0,0]
            y_trans = matrix[1,1]
            z_trans = matrix[2,2]
            trans_list = np.array([x_trans,y_trans,z_trans])
            return trans_list

        def CombineTwoPart(rot,trans):
            whole_seitz = np.column_stack((rot,trans))
            last_row = np.array([0,0,0,1])
            whole_seitz = np.vstack((whole_seitz,last_row))
            return whole_seitz

        chirality_change_positions = self.chirality_change_positions
        # define the relationship between two different chirlaity sets
        inverse_op = gemmi.Op('-x,-y,-z')
        reflection_opb = gemmi.Op('-x,y,-z')
        reflection_opa = gemmi.Op('x,-y,-z')
        reflection_opc = gemmi.Op('-x,-y,z')

        if inverse_op in chirality_change_positions:
            rot_relation_op_seitz = np.array([[-1,0,0],[0,-1,0],[0,0,-1]])
            trans_list = np.array([0, 0, 0])
        elif reflection_opb in chirality_change_positions:
            rot_relation_op_seitz = np.array([[-1,0,0],[0,1,0],[0,0,-1]])
            trans_list = np.array([0, 0, 0])
        elif reflection_opa in chirality_change_positions:
            rot_relation_op_seitz = np.array([[1,0,0],[0,-1,0],[0,0,-1]])
            trans_list = np.array([0, 0, 0])
        elif reflection_opc in chirality_change_positions:
            rot_relation_op_seitz = np.array([[-1,0,0],[0,-1,0],[0,0,1]])
            trans_list = np.array([0, 0, 0])
        else:
            relation_op = chirality_change_positions[0]
            relation_op_seitz = relation_op.seitz()
            relation_op_seitz = np.array(relation_op_seitz, dtype=object)
            rot_relation_op_seitz,trans_relation_op_seitz,trans_list = RotTransSeparation(relation_op_seitz)

        changed_chirality_centers = []
        for i in range(len(ops_centers)):
            centeri = ops_centers[i]
            centeri_seitz = centeri.seitz()
            centeri_seitz = np.array(centeri_seitz)
            centeri_rot_seitz,centeri_trans_seitz,centeri_trans_list = RotTransSeparation(centeri_seitz)
            rot_changed_seitz = np.dot(centeri_rot_seitz,rot_relation_op_seitz)
            rot_changed_trans = np.dot(centeri_trans_seitz,rot_relation_op_seitz)
            rot_changed_trans_list = TransSeitzToList(rot_changed_trans)
            all_trans = np.array(rot_changed_trans_list)+np.array(trans_list)
            changed_whole_seitz = CombineTwoPart(rot_changed_seitz,all_trans)
            changed_whole_seitz = list(changed_whole_seitz)
            op_center = gemmi.seitz_to_op(changed_whole_seitz)
            changed_chirality_centers.append(op_center)
        return changed_chirality_centers


def sg_table_read(file_name):
    '''
    :param file_name:
    :return: the table of spacegroup info of points classification
    '''
    try:
        with open(file_name) as file:
            sg_info_table = json.load(file)
            return sg_info_table
    except FileNotFoundError:
        print("The file was not found.")
    except json.JSONDecodeError:
        print("The file contains invalid JSON.")

