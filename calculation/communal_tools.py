import numpy as np
import os
import math
import gemmi
from itertools import combinations
import sys

class math_method:

    def __init__(self):
        return

    @staticmethod
    def Orthogonalization(coor, a, b, c, alpha, beta, gamma):
        '''
        to transform point's coordinates under crystallographic coordinate system
        to its coordinates in orthogonal 3D coordinate system
        :param: coor: list of its crystallographic coordinate
        :param a/b/c: float,cell parameters under the unit A
        :param alpha/beta/gamma: float,cell parameters in degree
        :return: return an array
        '''
        alpha_rad = np.radians(alpha)
        beta_rad = np.radians(beta)
        gamma_rad = np.radians(gamma)
        cos_alpha = np.cos(alpha_rad)
        cos_beta = np.cos(beta_rad)
        cos_gamma = np.cos(gamma_rad)
        sin_gamma = np.sin(gamma_rad)
        if math.isclose(gamma,120,rel_tol = 1e-9) is True:
            T = np.array([[a, b * np.cos(gamma_rad), 0],
                  [0, b * np.sin(gamma_rad), 0],
                  [0, 0, c]])
        else:
            T = np.array([[a, b * cos_gamma, c * cos_beta],
                          [0, b * sin_gamma, c * (cos_alpha - cos_beta * cos_gamma) / sin_gamma],
                          [0, 0, c * np.sqrt(1 - cos_alpha ** 2 - cos_beta ** 2 - cos_gamma ** 2 + 2 * cos_alpha * cos_beta * cos_gamma)]])
        coor_array = np.array(coor)
        after_ortho = np.dot(T,coor_array)
        result = np.around(after_ortho,3)
        return list(result)

    @staticmethod
    def Deorthogonalization(cartesian_coords, a, b, c, alpha, beta, gamma):
        '''
        to transform point's coordinates from orthogoanl 3d system to crystallographic coordinate system
        '''
        alpha_rad = np.radians(alpha)
        beta_rad = np.radians(beta)
        gamma_rad = np.radians(gamma)
        cos_alpha = np.cos(alpha_rad)
        cos_beta = np.cos(beta_rad)
        cos_gamma = np.cos(gamma_rad)
        sin_gamma = np.sin(gamma_rad)
        if math.isclose(gamma, 120, rel_tol=1e-9) is True:
            T = np.array([[a, b * np.cos(gamma_rad), 0],
                          [0, b * np.sin(gamma_rad), 0],
                          [0, 0, c]])
        else:
            T = np.array([[a, b * cos_gamma, c * cos_beta],
                          [0, b * sin_gamma, c * (cos_alpha - cos_beta * cos_gamma) / sin_gamma],
                          [0, 0, c * np.sqrt(
                              1 - cos_alpha ** 2 - cos_beta ** 2 - cos_gamma ** 2 + 2 * cos_alpha * cos_beta * cos_gamma)]])
        T_inv = np.linalg.inv(T)
        cartesian_coords = np.array(cartesian_coords)
        fractional_coords = np.dot(T_inv, cartesian_coords)
        result = np.around(fractional_coords,3)
        return list(result)

    @staticmethod
    def RemoveDuplicate(dict1):  #to remove duplicates in a dictionary
        unique_dict = {}
        for key, value in dict1.items():
            if value not in unique_dict.values():
                unique_dict[key] = value
        return unique_dict

    @staticmethod
    def PPdistance(point1,point2):
        '''
        :param point1: type = list
        :param point2: type = list
        :return: the distance between two points, float type
        '''
        x1 = point1[0]
        y1 = point1[1]
        z1 = point1[2]
        x2 = point2[0]
        y2 = point2[1]
        z2 = point2[2]
        PPdistance = (((x1 - x2) ** 2) + ((y1 - y2) ** 2) + ((z1 - z2) ** 2)) + +0.5
        final_distance  = np.around(PPdistance,3)
        return final_distance

    @staticmethod
    def PointLineDistance(point,line_point1,line_point2):
        '''
        :param point:
        :param line_point1: point in line, these two points define a line in space
        :param line_point2: point in line
        :return: the distance between one point and a line in space, float type
        '''

        point = np.array(point)
        line_point1 = np.array(line_point1)
        line_point2 = np.array(line_point2)
        line_vector = line_point2 - line_point1
        point_vector = point - line_point1
        projection = np.dot(point_vector, line_vector) / np.dot(line_vector, line_vector)
        projected_point = line_point1 + projection * line_vector
        distance = np.linalg.norm(point - projected_point)
        final_distance = np.around(distance, 3)
        return final_distance


    def LineLineDistance(self,line1_point1, line2_point1, direction):
        '''
        :param line1_point1:type=list
        :param line2_point1:type=list
        :param direction: a list that describes a direction's vector
        :return: two parallel lines‘ distance, type = float
        '''

        line1_point1 = np.array(line1_point1)
        line2_point1 = np.array(line2_point1)
        line2_point2 = np.array(line2_point1 + direction)
        direction = np.array(direction)
        cross_product = np.cross(direction, direction)

        # If the cross product is a zero vector, the lines are parallel
        if np.linalg.norm(cross_product) == 0:
            return self.PointLineDistance(line1_point1, line2_point1, line2_point2)

        connecting_vector = line2_point1 - line1_point1
        shortest_distance = abs(np.dot(connecting_vector, cross_product)) / np.linalg.norm(cross_product)
        final_distance  = np.around(shortest_distance,3)
        return final_distance

    @staticmethod
    def GetCentroid(atom,n_atom):
        '''
        :param atom: atoms' coordinates(list), type=list
        :param n_atom: atom number in a molecule
        :return: type=list
        '''
        centroid = []
        x_sum = 0
        y_sum = 0
        z_sum = 0
        for i in range(n_atom):
            a1 = atom[i]
            xi = a1[0]
            yi = a1[1]
            zi = a1[2]
            x_sum = x_sum + xi
            y_sum = y_sum + yi
            z_sum = z_sum + zi
        x_centroid = np.around(x_sum / n_atom,3)
        y_centroid = np.around(y_sum / n_atom,3)
        z_centroid = np.around(z_sum / n_atom,3)
        centroid = [x_centroid, y_centroid, z_centroid]
        return centroid

    @staticmethod
    def WithinCell(point):  # translate point（crystallographic coordinate) into incell
        x = point[0]
        y = point[1]
        z = point[2]
        if x < 0:
            num = 0 - x
            x_in_cell = x + math.ceil(num)
        elif x > 1:
            num = x - 1
            x_in_cell = x - math.ceil(num)
        else:
            x_in_cell = x
        if y < 0:
            num = 0 - y
            y_in_cell = y + math.ceil(num)
        elif y > 1:
            num = y - 1
            y_in_cell = y - math.ceil(num)
        else:
            y_in_cell = y
        if z < 0:
            num = 0 - z
            z_in_cell = z + math.ceil(num)
        elif z > 1:
            num = z - 1
            z_in_cell = z - math.ceil(num)
        else:
            z_in_cell = z
        point_in_cell = [x_in_cell, y_in_cell, z_in_cell]
        return point_in_cell


class cell_method:

    def __int__(self):
        return

    @staticmethod
    def FindNearestPair(same_chirality_points, different_chirality_points):
        '''
        find the nearest two different chirality points to get the enlarging direction
        :param same_chirality_points: orthogonlaized
        :param different_chirality_points:
        :return: a list of point pair
        '''
        min_distance = float('inf')
        nearest_pair = None
        for point1 in same_chirality_points:
            for point2 in different_chirality_points:
                distance = np.linalg.norm(np.array(point1) - np.array(point2))
                if distance < min_distance:
                    min_distance = distance
                    nearest_pair = [point1, point2]
        return nearest_pair

    @staticmethod
    def PointEnlarging(points, EnlargeMatrix):
        '''
        to enalrge points in 2*2*2 cells for calculation preparation
        :param points: packed points' coordinates in crystallographic fractional system
        :param EnlargeMatrix: the matrix that define which direction to enlarge,
        tp make sure all possible surrounding different chirality lines are found
        :return: enlarged relatively different chirality points list
        '''
        len_points = len(points)

        enlarged_points = {}

        vector_x = EnlargeMatrix[0]
        vector_y = EnlargeMatrix[1]
        vector_z = EnlargeMatrix[2]

        enlarge_x = np.array(vector_x)
        enlarge_y = np.array(vector_y)
        enlarge_z = np.array(vector_z)

        for i in range(len_points):
            temp_point = np.array(points[i])
            enlarged_points[1 * len_points + i] = list(temp_point + enlarge_x)
            enlarged_points[2 * len_points + i] = list(temp_point + enlarge_y)
            enlarged_points[3 * len_points + i] = list(temp_point + enlarge_z)
            enlarged_points[4 * len_points + i] = list(temp_point + enlarge_x + enlarge_y)
            enlarged_points[5 * len_points + i] = list(temp_point + enlarge_x + enlarge_z)
            enlarged_points[6 * len_points + i] = list(temp_point + enlarge_y + enlarge_z)
            enlarged_points[7 * len_points + i] = list(temp_point + enlarge_x + enlarge_y + enlarge_z)

        enlarged_points_final = math_method.RemoveDuplicate(enlarged_points)
        enlarged_points_list = list(enlarged_points_final.values())
        return enlarged_points_list


    @staticmethod
    def EnlargeMatrix(point1, point2):
        '''
        point1 and point 2 are two different chirality points, this is to define the point 1's enlarge direction
        enlarge is to expand its cell from 1*1*1 to 2*2*1, the line direction axes not expand
        :return: a 3*3 matrix
        '''
        delta_x = point1[0] - point2[0]
        delta_y = point1[1] - point2[1]
        delta_z = point1[2] - point2[2]

        if delta_x > 0:
            dx = 1
        else:
            dx = -1
        if delta_y > 0:
            dy = 1
        else:
            dy = -1
        if delta_z > 0:
            dz = 1
        else:
            dz = -1

        EnlargeMatrix_dict = {}
        EnlargeMatrix_dict[0] = [-dx, 0, 0]
        EnlargeMatrix_dict[1] = [0, -dy, 0]
        EnlargeMatrix_dict[2] = [0, 0, -dz]

        EnlargeMatrix = list(EnlargeMatrix_dict.values())

        return EnlargeMatrix

    @staticmethod
    def EnlargeMatrix_2D(point1, point2, line_direction_frac):
        delta_x = point1[0] - point2[0]
        delta_y = point1[1] - point2[1]
        delta_z = point1[2] - point2[2]

        if delta_x > 0:
            dx = 1
        else:
            dx = -1
        if delta_y > 0:
            dy = 1
        else:
            dy = -1
        if delta_z > 0:
            dz = 1
        else:
            dz = -1

        line_direction_x = line_direction_frac[0]
        line_direction_y = line_direction_frac[1]
        line_direction_z = line_direction_frac[2]
        if line_direction_x != 0:
            dx = 0
        if line_direction_y != 0:
            dy = 0
        if line_direction_z != 0:
            dz = 0

        EnlargeMatrix_dict = {}
        EnlargeMatrix_dict[0] = [-dx,0,0]
        EnlargeMatrix_dict[1] = [0,-dy,0]
        EnlargeMatrix_dict[2] = [0,0,-dz]
        EnlargeMatrix = list(EnlargeMatrix_dict.values())

        return EnlargeMatrix
