'''
what info that is available
line direction and its corresponding operation
packing number
'''
import gemmi
from communal_tools import math_method
import numpy as np

class packing_points:

    def __init__(self,sg_number,packing_number,basic_op,same_chiral_points):
        '''

        :param sg_number: the index number of the space group
        :param packing_number: number of points in a pack, highest axis order or 2(for monoclinic)
        :param basic_op: the intrinsic operation in the packing
        :param same_chiral_points: in gemmi.Op type
        '''
        self.crys_points = same_chiral_points
        self.packing_number = packing_number
        self.basisop = basic_op
        same_chiral_points_groupops = gemmi.GroupOps(same_chiral_points)
        same_chiral_points_sym = same_chiral_points_groupops.sym_ops

        same_chiral_points_center = same_chiral_points_groupops.cen_ops

        self.crys_points_sym = same_chiral_points_sym
        self.crys_points_cen = same_chiral_points_center

        if sg_number == 2:
            packed_points = same_chiral_points
        elif len(same_chiral_points) == 1:
            packed_points = same_chiral_points
        elif sg_number in [7,9]:
            packed_points = self.NonAxisPacking()
        else:
            packed_points = self.AxisPacking()
        self.after_packing = packed_points
        return

    @staticmethod
    def AddingLattice(sym_ops,cen_ops):
        all_ops = []
        for i in range(len(sym_ops)):
            opi = sym_ops[i]
            for j in range(len(cen_ops)):
                cenj = cen_ops[j]
                if cenj == [0,0,0]:
                    cenj_op = gemmi.Op('x,y,z')
                elif cenj == [12,12,0]:
                    cenj_op = gemmi.Op('x+1/2,y+1/2,z')
                elif cenj == [12,0,12]:
                    cenj_op = gemmi.Op('x+1/2,y,z+1/2')
                elif cenj == [0,12,12]:
                    cenj_op = gemmi.Op('x,y+1/2,z+1/2')
                elif cenj == [12,12,12]:
                    cenj_op = gemmi.Op('x+1/2,y+1/2,z+1/2')
                elif cenj == [16,8,8]:
                    cenj_op = gemmi.Op('x+2/3,y+1/3,z+1/3')
                elif cenj == [8,16,16]:
                    cenj_op = gemmi.Op('x+1/3,y+2/3,z+2/3')
                else:
                    print('lattice_type error')
                cen_add_op = cenj_op * opi
                all_ops.append(cen_add_op)
        return all_ops

    @property
    def GeneralPositions(self):
        return self.after_packing

    def NonAxisPacking(self):
        '''
        this is for space group Cm,Cc,
        return: the center ops
        '''
        points = self.crys_points_sym
        points_center_ops = self.GetOpsCenter(points)
        points_center_op_list = [points_center_ops]
        cen_ops = self.crys_points_cen
        all_packed_ops = self.AddingLattice(points_center_op_list,cen_ops)
        return all_packed_ops

    def AxisPacking(self):
        points = self.crys_points_sym
        lattice_centers = self.crys_points_cen
        packing_number = self.packing_number
        points_number = len(self.crys_points)
        packed_ops_centers = []
        if packing_number == points_number:
            packed_ops_center = self.GetOpsCenter(points)
            packed_ops_centers.append(packed_ops_center)
        else:
            points_set = self.SubsetPoints()
            packed_ops_centers = self.SubsetPacking(points_set)
        all_packed_ops_centers = self.AddingLattice(packed_ops_centers,lattice_centers)

        return all_packed_ops_centers

    def RelationsBetweenSubsets(self,subsets):
        recorded_relations = []
        for i in range(len(subsets)):
            subset = subsets[i]
            if gemmi.Op('x,y,z') not in subset:
                recorded_relations.append(subset[0])
            else:
                center = self.GetOpsCenter(subset)
        recorded_relations.insert(0,gemmi.Op('x,y,z'))
        return recorded_relations,center

    def SubsetPacking(self,points_in_set):
        recorded_relations,first_ops_center = self.RelationsBetweenSubsets(points_in_set)
        all_ops_centers = []
        for i in range(len(recorded_relations)):
            relation_i = recorded_relations[i]
            temp_center = first_ops_center*relation_i
            all_ops_centers.append(temp_center)
        return all_ops_centers

    def SubsetPoints(self):
        points = self.crys_points_sym
        packing_number = self.packing_number
        basic_op = self.basisop
        set_number = len(points)/packing_number
        remain_points = points.copy()
        remain_points.remove(gemmi.Op('x,y,z'))
        remain_points.insert(0,gemmi.Op('x,y,z'))

        points_in_sets = []
        if int(set_number)%1 != 0:
            return "missing general position"
        else:
            while len(remain_points) != 0:
                p1 =remain_points[0]
                set1 = self.GrouopSet(p1,basic_op)
                points_in_sets.append(set1)
                remain_points = [item for item in remain_points if item not in set1]

        return points_in_sets

    @staticmethod
    def GrouopSet(p1, op):
        start_point = p1
        process_point = p1 * op
        points_set = [process_point]
        while process_point != start_point:
            process_point = process_point * op
            points_set.append(process_point)
        return points_set

    @staticmethod
    def GetOpsCenter(ops):
        number = len(ops)
        ops_seitz_array = []
        for i in range(number):
            op = ops[i]
            op_seitz = op.seitz()
            op_seitz_array = np.array(op_seitz)
            ops_seitz_array.append(op_seitz_array)

        sum_array = ops_seitz_array[0]
        for i in range(1,number):
            array_i = ops_seitz_array[i]
            sum_array = sum_array+array_i
        center_array = sum_array/number
        center_op = gemmi.seitz_to_op(center_array)
        return center_op

'''
sg_number = 12
packing_number = 2
basic_op = gemmi.Op('-x,y,-z')
same_chiral_points = [gemmi.Op('x,y,z'),gemmi.Op('-x,y,-z'),gemmi.Op('x+1/2,y+1/2,z'),gemmi.Op('-x+1/2,y+1/2,-z')]
packed_points = packing_points(sg_number,packing_number,basic_op,same_chiral_points)
print(packed_points.after_packing)
'''











