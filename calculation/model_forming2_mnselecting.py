import numpy as np
from itertools import combinations
from communal_tools import math_method
from communal_tools import cell_method
import math


class folded_plane:

    def __init__(self,sg_number,molecule_frac,packed_same_chirality_points,packed_different_chirality_points,
                 line_direction_frac,a,b,c,alpha,beta,gamma):
        '''
        :param sg_number: the index of the space group
        :param molecule_frac: list of molecule coordinate
        :param packed_same_chirality_points: gemmi.op
        :param packed_different_chirality_points: gemmi.op
        :param line_direction_frac: in the type of [0,0,1]
        :param a,b,c,alpha,beta,gamma:cell parameters, angle in radius,length in atrome
        '''
        self.sg_number = sg_number
        self.cell_parameters = {'a':a,"b":b,'c':c,'alpha':alpha,'beta':beta,'gamma':gamma}

        #points in gemmi.ops type
        self.same_chirlaity_general_points = packed_same_chirality_points
        self.different_chirality_general_points = packed_different_chirality_points

        chirality_retain_points_frac = self.GetFracPointsCoordinates(packed_same_chirality_points,molecule_frac)
        chirality_change_points_frac = self.GetFracPointsCoordinates(packed_different_chirality_points,molecule_frac)
        self.same_chirlaity_points_frac = chirality_retain_points_frac
        self.different_chirality_points_frac = chirality_change_points_frac

        self.line_direction_frac = list(line_direction_frac)
        line_direction_ortho = math_method.Orthogonalization(line_direction_frac,a,b,c,alpha,beta,gamma)
        self.line_direction_ortho = list(line_direction_ortho)

        incell_same_points_number = len(packed_same_chirality_points)
        if incell_same_points_number == 1:
            point1 = chirality_retain_points_frac[0]
            point2 = chirality_change_points_frac[0]
            print(point1,point2)
            info = self.OnePairModel(point1,point2)
        elif incell_same_points_number == 2:
            whether_collinear_inset, y_inset_frac_vector = self.WhetherAllCollinear(chirality_retain_points_frac,line_direction_frac)
            if whether_collinear_inset is True:
                info = self.CollinearTwoLinesMethod(y_inset_frac_vector,chirality_retain_points_frac,chirality_change_points_frac)
            else:  # two same chirality points represent different lines
                info = self.MultiPairMethod(line_direction_frac,chirality_retain_points_frac,chirality_change_points_frac)
        # if there are more than 2 points in a set
        else:  # if all points are collinear/ some of them are collinear/ non-collinear
            whether_all_collinear, y_inset_vector = self.WhetherAllCollinear(chirality_retain_points_frac,line_direction_frac)
            if whether_all_collinear is True:
                info = self.CollinearTwoLinesMethod(y_inset_frac_vector, chirality_retain_points_frac,
                                                    chirality_change_points_frac)
            else:
                whether_collinear_exist,y_inset_vector,one_line_one_point1 = self.WhetherCollinearExists(chirality_retain_points_frac,line_direction_frac)
                if whether_collinear_exist is False:
                    info = self.MultiPairMethod(line_direction_frac, chirality_retain_points_frac, chirality_change_points_frac)
                else:
                    points1 = one_line_one_point1
                    points2 = self.RemoveDuplicateCollinearPoints(chirality_change_points_frac,line_direction_frac)
                    info = self.MultiPairMethod(y_inset_vector,points1,points2)

        theta_angle = info['theta']
        sigma_angle = info['sigma']
        m_value = info['m']
        n_value = info['n']
        x_value = info['x']
        y_value = info['y']
        x_vector_frac = info['x_vector']
        k_value = info['k']

        self.x_vector = x_vector_frac
        self.theta = theta_angle
        self.sigma = sigma_angle
        self.m = m_value
        self.n = n_value
        self.x = x_value
        self.y = y_value
        self.k = k_value
        return


    def MultiPairMethod(self,y_frac,points1,points2):
        '''
        each point in point set represent a line and  there are no collinear situation
        points should be as near as possible to make sure find the true x direction
        :param y_frac:
        :param points1: in frac coordinate, more than 1
        :param points2:in frac coordinate,more than 1
        :return:
        '''
        line_direction_frac = self.line_direction_frac
        line_direction_ortho = self.line_direction_ortho
        cell_parameters = self.cell_parameters
        print('1')
        y_ortho = math_method.Orthogonalization(y_frac,cell_parameters['a'],cell_parameters['b'],cell_parameters['c'],
                                                      cell_parameters['alpha'],cell_parameters['beta'],cell_parameters['gamma'])
        y_value = np.linalg.norm(y_ortho)

        center_points1 = math_method.GetCentroid(points1,len(points1))
        center_points2 = math_method.GetCentroid(points2, len(points2))
        enlarge_matrix = cell_method.EnlargeMatrix_2D(center_points2,center_points1,line_direction_frac)
        enlarged_points2 = cell_method.PointEnlarging(points2,enlarge_matrix)

        points1_ortho = self.GetOrthoPointsCoordinate(points1,cell_parameters['a'],cell_parameters['b'],cell_parameters['c'],
                                                      cell_parameters['alpha'],cell_parameters['beta'],cell_parameters['gamma'])

        enlarged_points2_ortho = self.GetOrthoPointsCoordinate(enlarged_points2,cell_parameters['a'],cell_parameters['b'],cell_parameters['c'],
                                                      cell_parameters['alpha'],cell_parameters['beta'],cell_parameters['gamma'])
        # select one point1 and two point23 which in the nearset situation
        all_mns = []
        all_mns_lines = []
        for i in range(len(points1_ortho)):
            point1 = points1_ortho[i]
            two_nearest_lines_info = self.TwoNearestLines(point1,enlarged_points2_ortho,line_direction_ortho)
            two_distances = two_nearest_lines_info['distances']
            two_points2 = two_nearest_lines_info['lines']
            all_mns.append(two_distances)
            all_mns_lines.append(two_points2)

        minimum_mn_pair = min(all_mns,key=lambda x:(x[0],x[1]))
        m = minimum_mn_pair[0]
        n = minimum_mn_pair[1]
        index_point1 = all_mns.index(minimum_mn_pair)
        selected_point1 = points1_ortho[index_point1]
        two_points2 = all_mns_lines[index_point1]

        #
        x_pair_points_ortho,x_cell_direction = self.XPointsPairFix(cell_parameters,two_points2,line_direction_frac)

        x_direction_ortho = list(np.array(x_pair_points_ortho[0]) - np.array(x_pair_points_ortho[1]))
        x = np.linalg.norm(x_direction_ortho)

        sigma_radius, sigma_label = self.GetSigma(x_direction_ortho, line_direction_ortho, cell_parameters['alpha'],
                                                  cell_parameters['beta'], cell_parameters['gamma'])

        if sigma_label != "None":
            sigma_angle = cell_parameters[sigma_label]
        else:
            sigma_angle = (sigma_radius * 180) / np.pi

        sg_number = self.sg_number
        theta_radius = self.GetTheta(sg_number, x, y_value, m, n, sigma_radius)
        theta_angle = (theta_radius * 180) / np.pi
        theta_angle = round(theta_angle, 3)
        k = (m+n)/(x*np.sin(sigma_radius))
        output_dict = {'theta': theta_angle, 'sigma': sigma_angle, 'm': m, 'n': n, 'x': x, 'y': y_value,'x_vector':x_cell_direction,'k':k}
        return output_dict

    @staticmethod
    def XPointsPairFix(cell_parameters,two_points_ortho,line_direction_frac):
        '''
        to make sure that the x vector get the proper two points
        :param two_points_ortho:
        :return:
        '''
        def find_intersection(P1,P2,v1,v2): # find the cross point, P1 with line direction v1
            A = np.array([
                [v1[0], -v2[0]],  # Coefficients for x equation
                [v1[1], -v2[1]],  # Coefficients for y equation
                [v1[2], -v2[2]]  # Coefficients for z equation
            ])
            B = np.array([
                P2[0] - P1[0],  # x difference between points P1 and P2
                P2[1] - P1[1],  # y difference between points P1 and P2
                P2[2] - P1[2]  # z difference between points P1 and P2
            ])
            result, residuals, rank, s_values = np.linalg.lstsq(A, B, rcond=None)
            t, s = result
            intersection_point = P1 + t * np.array(v1)
            return intersection_point

        def is_parallel(v1,v2s):
            v1_array = np.array(v1)
            for v2 in v2s:
                v2_array = np.array(v2)
                cross_vector = np.cross(v1_array,v2_array)
                cross_value = np.linalg.norm(cross_vector)
                cross_value = np.around(cross_value,1)
                if math.isclose(cross_value,0,rel_tol=1e-6) is True:
                    return True
                else:
                    continue
            return False

        point21_ortho = two_points_ortho[0]
        point22_ortho = two_points_ortho[1]
        point21_frac = math_method.Deorthogonalization(point21_ortho,cell_parameters['a'],cell_parameters['b'],cell_parameters['c'],
                                                      cell_parameters['alpha'],cell_parameters['beta'],cell_parameters['gamma'])
        point22_frac = math_method.Deorthogonalization(point22_ortho,cell_parameters['a'],cell_parameters['b'],cell_parameters['c'],
                                                      cell_parameters['alpha'],cell_parameters['beta'],cell_parameters['gamma'])
        initial_x_vector_array = np.array(point21_frac)-np.array(point22_frac)
        initial_x_vector_array2 = np.around(initial_x_vector_array,2)
        initial_x_vector = list(initial_x_vector_array2)

        potential_x_direction = [[1,0,0],[0,1,0],[0,0,1],
                                 [1,1,0],[0,1,1],[1,0,1],
                                 [-1,1,0],[0,-1,1],[-1,0,1],
                                 [1,1,1],[1,1,-1],[-1,1,1],[1,-1,1],
                                 [-1,-1,1],[-1,1,-1],[1,-1,-1],
                                 [2,1,1],[1,2,1],[1,1,2],
                                 [2,-1,1],[1,-2,1],[1,-1,2],
                                 [-2,1,1],[-1,2,1],[-1,1,2],
                                 [2,1,-1],[1,2,-1],[1,1,-2],
                                 [1,2,2],[2,1,2],[2,2,1],
                                 [-1,2,2],[-2,1,2],[-2,2,1],
                                 [1,2,-2],[2,1,-2],[2,2,-1]]

        potential_x_direction.remove(line_direction_frac)

        whether_parallel = is_parallel(initial_x_vector,potential_x_direction)

        if whether_parallel is True:
            return [point21_ortho,point22_ortho],initial_x_vector

        refined_point22_frac = point21_frac
        point22_frac_array = np.array(point22_frac)
        output_x_direction = None
        for i in range(len(potential_x_direction)):
            x_direction = potential_x_direction[i]
            cross_point_array = find_intersection(point21_frac,point22_frac,x_direction,line_direction_frac)
            distance = np.linalg.norm(cross_point_array-point22_frac_array)
            are_equal_y1 = math.isclose(distance,1,rel_tol = 1e-6)
            are_equal_y2 = math.isclose(distance,2,rel_tol = 1e-6)
            are_itself = math.isclose(distance,0,rel_tol = 1e-6)

            if are_equal_y1 or are_equal_y2 is True:
                refined_point22_frac = list(cross_point_array)
                output_x_direction = x_direction
            elif are_itself is True:
                return two_points_ortho,x_direction
            else:
                continue
        refined_point22_ortho = math_method.Orthogonalization(refined_point22_frac,cell_parameters['a'],cell_parameters['b'],cell_parameters['c'],
                                                      cell_parameters['alpha'],cell_parameters['beta'],cell_parameters['gamma'])

        return [point21_ortho,refined_point22_ortho],output_x_direction

    @staticmethod
    def TwoNearestLines(point1_ortho,points2_ortho,line_direction_ortho):
        '''
        point1 represent a line. points2 represent a series of lines, find the nearest two lines in points2 to the point1
        :param point1_ortho:
        :param points2_ortho:
        :param line_direction_ortho:
        :return: {'lines':selected 2 points,'distances':two minimum values}
        '''
        line_line_distance_set = []
        math_object = math_method()
        for i in range(len(points2_ortho)):
            point2_i = points2_ortho[i]
            distance_i = math_object.LineLineDistance(point1_ortho,point2_i,line_direction_ortho)
            line_line_distance_set.append(distance_i)

        min1_index = line_line_distance_set.index(min(line_line_distance_set))
        temp = line_line_distance_set[min1_index]
        line_line_distance_set[min1_index] = float('inf')
        min2_index = line_line_distance_set.index(min(line_line_distance_set))
        line_line_distance_set[min1_index] = temp
        point21 = points2_ortho[min1_index]
        point22 = points2_ortho[min2_index]
        min_distance1 = line_line_distance_set[min1_index]
        min_distance2 = line_line_distance_set[min2_index]
        return {'lines':[point21,point22],'distances':[min_distance1,min_distance2]}

    @staticmethod
    def GetFracPointsCoordinates(points_ops,molecule_frac):
        points_frac_coor = []
        for i in range(len(points_ops)):
            ops = points_ops[i]
            coor = ops.apply_to_xyz(molecule_frac)
            points_frac_coor.append(coor)
        return points_frac_coor

    @staticmethod
    def GetOrthoPointsCoordinate(points,a,b,c,alpha,beta,gamma):
        points_ortho_coor = []
        for i in range(len(points)):
            temp_point = points[i]
            temp_ortho = math_method.Orthogonalization(temp_point,a,b,c,alpha,beta,gamma)
            points_ortho_coor.append(temp_ortho)
        return points_ortho_coor

    def OnePairModel(self,point1,point2):  # if there is only one point in a homo-chiral domain
        '''
        in this situation,there is only a pair of enantiomer points in a cell
        :param point1: in fractional coordinate, list
        :param point2: in fractional coordinate,list
        :return:
        '''

        line_direction_frac1 = self.line_direction_frac
        line_direction_ortho = self.line_direction_ortho
        y_value = np.linalg.norm(line_direction_ortho)

        enlarge_matrix = cell_method.EnlargeMatrix_2D(point1,point2,line_direction_frac1)
        point1_list = [point1]
        enlarged_points1 = cell_method.PointEnlarging(point1_list,enlarge_matrix)

        cell_parameters = self.cell_parameters
        a = cell_parameters['a']
        b = cell_parameters['b']
        c = cell_parameters['c']
        alpha = cell_parameters['alpha']
        beta = cell_parameters['beta']
        gamma = cell_parameters['gamma']
        enlarged_points1_ortho = self.GetOrthoPointsCoordinate(enlarged_points1,a,b,c,alpha,beta,gamma)
        point2_ortho = math_method.Orthogonalization(point2,a,b,c,alpha,beta,gamma)

        mn_output = self.GetMN(same_point=point2_ortho,different_points=enlarged_points1_ortho,line_direction=line_direction_ortho)
        m = mn_output['m']
        n = mn_output['n']
        selected_two_lines = mn_output['selected_lines']

        point11 = selected_two_lines[0]
        point12 = selected_two_lines[1]
        x_vector = np.array(point11) - np.array(point12)
        x = np.linalg.norm(x_vector)

        cell_parameters = self.cell_parameters
        sigma_radius,sigma_label = self.GetSigma(x_vector,line_direction_ortho,cell_parameters['alpha'],
                                     cell_parameters['beta'],cell_parameters['gamma'])

        x_vector_frac = math_method.Deorthogonalization(x_vector,cell_parameters['a'],cell_parameters['b'],cell_parameters['c'],cell_parameters['alpha'],
                                     cell_parameters['beta'],cell_parameters['gamma'])

        if sigma_label != 'None':
            sigma_angle = cell_parameters[sigma_label]
        else:
            sigma_angle = (sigma_radius*180)/np.pi

        sg_number = self.sg_number
        theta_radius = self.GetTheta(sg_number,x,y_value,m,n,sigma_radius)
        theta_angle = (theta_radius*180)/np.pi
        theta_angle = round(theta_angle,3)
        k = (m+n)/(x*np.sin(sigma_radius))
        output_dict = {'theta':theta_angle,'sigma':sigma_angle,'m':m,'n':n,'x':x,'y':y_value,'x_vector':x_vector_frac,'k':k}
        return output_dict

    @staticmethod
    def GetTheta(space_group_number,x,y,m,n,sigma):
        '''
        :param space_group_number: the index of the crystal's space group among 230
        :return: theta in radius system
        '''
        print('-----this part is theta calculation----')
        print('m',m,'n',n,'x',x,'y',y,'sigma',sigma)
        k = (m + n) / (x * np.sin(sigma))

        tan_theta = (k * x * np.sin(sigma)) / (y + np.abs(x * np.cos(sigma)))
        theta_radius = np.arctan(tan_theta)
        theta_radius = round(theta_radius, 4)

        return theta_radius

    @staticmethod
    def GetSigma(x_direction,line_direction_ortho,alpha,beta,gamma):
        '''
        :param x_direction: the vector's direction in list. the x vector is from a chiral(nearest different chiral line)
         line's point to the nearest point in its same chiral line(next nearest to different chiral line)
        :param line_direction_ortho,line direction after orthogonalization
        :param alpha, beta, gamma:crystal's cell parameters
        :return: sigma in radius system
        '''

        alpha_radius = (alpha / 180) * np.pi
        beta_radius = (beta / 180) * np.pi
        gamma_radius = (gamma / 180) * np.pi

        length_x = np.linalg.norm(x_direction)

        length_line = np.linalg.norm(line_direction_ortho)
        vector_inner = np.dot(x_direction, line_direction_ortho)

        cos_sigma = vector_inner / (length_x * length_line)
        cos_sigma = round(cos_sigma,3)

        cos_alpha = round(np.cos(alpha_radius),3)
        cos_beta = round(np.cos(beta_radius), 3)
        cos_gamma = round(np.cos(gamma_radius), 3)

        angle_cos_dict = {'alpha':cos_alpha,'beta':cos_beta,'gamma':cos_gamma}
        angle_list = {'alpha':alpha,'beta':beta,'gamma':gamma}
        cos_list = [cos_alpha,cos_beta,cos_gamma]

        if cos_sigma in cos_list:
            for key,value in angle_cos_dict.items():
                if value == cos_sigma:
                    sigma_label = key
                    sigma_theta = angle_list[sigma_label]

            sigma_radius = (sigma_theta*np.pi)/180
        else:
            sigma_radius = np.arccos(cos_sigma)
            sigma_label = 'None'
        print('sigma info',sigma_radius,sigma_label)
        return sigma_radius,sigma_label

    @staticmethod
    def GetMN(same_point,different_points,line_direction):
        '''
        :param same_point: in ortho coordinate
        :param different_points: in ortho coordinate
        :param line_direction: in ortho coordinates
        :return:lines are represented by points in list type
        '''
        def multi_same_distance_selection(different_points,line_direction):
            maths = math_method()
            nearest_pair = None
            min_distance = float('inf')
            for point1,point2 in combinations(different_points,2):
                ll_distance = maths.LineLineDistance(point1,point2,line_direction)
                if ll_distance < min_distance:
                    min_distance = ll_distance
                    nearest_pair = [point1,point2]
            return nearest_pair

        line_line_distance = []
        math_object = math_method()
        for i in range(len(different_points)):
            distance_i = math_object.LineLineDistance(different_points[i], same_point, line_direction)
            line_line_distance.append(distance_i)
        print('lldistance in GetMN',line_line_distance)
        sorted_line_line_distance = line_line_distance.copy()
        sorted_line_line_distance.sort()
        m = sorted_line_line_distance[0]
        n = sorted_line_line_distance[1]  #get m n by distance
        if m != n:
            index_m = line_line_distance.index(m)
            index_n = line_line_distance.index(n)
            selected_different_chirality_lines = [different_points[index_m], different_points[index_n]]
        else:  # in this case
            indexes = [index for index,element in enumerate(line_line_distance) if element == m]
            if len(indexes) >2:
                all_surrounding_different_points = [different_points[i] for i in indexes]
                selected_different_chirality_lines = multi_same_distance_selection(all_surrounding_different_points,line_direction)
            elif len(indexes) == 2:
                index_m = indexes[0]
                index_n = indexes[1]
                selected_different_chirality_lines = [different_points[index_m],different_points[index_n]]
            else:
                print('there is an error in finding nearest two different chirality points')

        output = {'m':m,'n':n,'selected_lines':selected_different_chirality_lines}

        return output

    @staticmethod
    def WhetherAllCollinear(points,line_direction):
        '''
        points can be both ortho or fractional, it decided the output vector
        to judge whether all these points are in the line with given line direction
        if True, return the interval vector
        else False,return False with None
        '''

        line_direction_array = np.array(line_direction)
        is_collinear = True
        vector_norm = []
        vectors = []

        for point_pair in combinations(points,2):
            point1 = np.array(point_pair[0])
            point2 = np.array(point_pair[1])
            vector = point1 - point2
            cross_value = np.cross(vector,line_direction_array)
            cross_length = np.linalg.norm(cross_value)
            if cross_length != 0: # if non-collinear exists, in the case that all points are incell
                is_collinear = False
                break
            else:
                distance = np.linalg.norm(vector)
                vector_norm.append(distance)
                vectors.append(vector)

        if is_collinear is True:
            y_value = min(vector_norm)
            y_index = vector_norm.index(y_value)
            y_vector = vectors[y_index]
        else:
            y_vector = [None]

        return is_collinear,list(y_vector)

    @staticmethod
    def WhetherCollinearExists(points,line_direction):
        '''
        to judge whether there are collinear points, if there are collinear, return the minimum distance vector and the non-collinear points
        :param points: usually fractional, but also can be ortho
        :return:collinear judgement(bool); y vector if Ture(list), list of points without collinear
        '''
        def PointsInLineSet(input_points,line_direction_array1):
            '''
            put the collinear points into a list to represent a line, return the lines' list
            '''
            lines = []
            remaining_points = input_points[:]  # Copy the original points list
            while remaining_points:
                point1 = remaining_points[0]
                one_line = [point1]
                points_on_line = []  # Track points to remove later
                for point2 in remaining_points[1:]:  # Iterate over remaining points (excluding point1)
                    vector1 = np.array(point1) - np.array(point2)
                    cross_value1 = np.cross(vector1, line_direction_array1)
                    if np.linalg.norm(cross_value1) < 1e-9:  # Use a small tolerance instead of == 0
                        one_line.append(point2)
                        points_on_line.append(point2)  # Mark point2 for removal
                remaining_points.remove(point1)
                for p in points_on_line:
                    remaining_points.remove(p)
                lines.append(one_line)
            return lines

        def FindYVector(points_in_set):
            def find_nearest_two_points(ps):
                min_distance = float('inf')
                nearest_pair = None
                for i in range(len(ps)):
                    for j in range(i + 1, len(ps)):
                        point1 = np.array(ps[i])
                        point2 = np.array(ps[j])
                        distance = np.linalg.norm(point1 - point2)
                        if distance < min_distance:
                            min_distance = distance
                            nearest_pair = (ps[i], ps[j])
                return nearest_pair
            y_vector1 = [0,0,0]
            for i in range(len(points_in_set)):
                line1 = points_in_set[i]
                points_number = len(line1)
                if points_number == 1:
                    continue
                elif points_number == 2:
                    y_vector1 = list(np.array(line1[0])-np.array(line1[1]))
                    return y_vector1
                else: # more than 2 in a line
                    nearset_pair = find_nearest_two_points(line1)
                    y_vector1 = list(np.array(nearset_pair[0]) - np.array(nearset_pair[1]))
                    return y_vector1
            if y_vector1 == [0,0,0]:
                raise ValueError ('An error occurred when calculating the y which is smaller than cell parameters,two different judgement when judging the collinear')
            else:
                return y_vector1

        is_collinear = False
        line_direction_array = np.array(line_direction)
        for points_pair in combinations(points,2):
            vector = np.array(points_pair[0]) - np.array(points_pair[1])
            cross_value = np.cross(vector,line_direction_array)
            if np.linalg.norm(cross_value) < 1e-9:
                is_collinear = True
            else:
                continue

        if is_collinear is False:
            collinear_exist = False
            return collinear_exist, None,None
        else:
            collinear_exist = True
            classified_line_sets = PointsInLineSet(points,line_direction_array)
            y_vector = FindYVector(classified_line_sets)
            one_line_one_point = []
            for i in range(len(classified_line_sets)):
                line_i = classified_line_sets[i]
                one_line_one_point.append(line_i[0])
            return collinear_exist,y_vector,one_line_one_point

    def CollinearTwoLinesMethod(self,y_frac,points1,points2):
        '''
        if there are only a pair of chiral line in cell but more than one points in one line
        :param y_frac: y vector in frac, must be parallel with line direction but shorter
        :param points1: points in frac coordinates
        :param points2: points in frac coordinates
        :return: basic info for the folded-plane model
        '''
        line_direction_ortho = self.line_direction_ortho
        line_direction_frac = self.line_direction_frac
        cell_parameters = self.cell_parameters
        y_ortho = math_method.Orthogonalization(y_frac,cell_parameters['a'],cell_parameters['b'], cell_parameters['c'],
                                cell_parameters['alpha'], cell_parameters['beta'],cell_parameters['gamma'])
        y_value = np.linalg.norm(y_ortho)

        point1 = points1[0] #one points with line direction represent a line
        point2 = points2[0]
        enlarge_matrix = cell_method.EnlargeMatrix_2D(point1, point2, line_direction_frac)
        enlarged_points1 = cell_method.PointEnlarging(point1, enlarge_matrix)

        enlarged_points1_ortho = self.GetOrthoPointsCoordinate(enlarged_points1, cell_parameters['a'],cell_parameters['b'], cell_parameters['c'],
                                                               cell_parameters['alpha'], cell_parameters['beta'],cell_parameters['gamma'])
        point2_ortho = math_method.Orthogonalization(point2, cell_parameters['a'], cell_parameters['b'],cell_parameters['c'],
                                                     cell_parameters['alpha'], cell_parameters['beta'],cell_parameters['gamma'])

        print('')
        mn_output = self.GetMN(point2_ortho, enlarged_points1_ortho, line_direction_ortho)

        m = mn_output['m']
        n = mn_output['n']
        selected_two_lines = mn_output['selected_lines']

        point11 = selected_two_lines[0]
        point12 = selected_two_lines[1]
        x_vector = list(np.array(point11) - np.array(point12))
        x = np.linalg.norm(x_vector)

        x_vector_frac = math_method.Deorthogonalization(x_vector,cell_parameters['a'],cell_parameters['b'], cell_parameters['c'],cell_parameters['alpha'],
                                                  cell_parameters['beta'], cell_parameters['gamma'])

        cell_parameters = self.cell_parameters
        sigma_radius, sigma_label = self.GetSigma(x_vector, line_direction_ortho, cell_parameters['alpha'],
                                                  cell_parameters['beta'], cell_parameters['gamma'])
        if sigma_label != 'None':
            sigma_angle = cell_parameters[sigma_label]
        else:
            sigma_angle = (sigma_radius * 180) / np.pi


        sg_number = self.sg_number
        theta_radius = self.GetTheta(sg_number, x, y_value, m, n, sigma_radius)
        theta_angle = (theta_radius * 180) / np.pi
        theta_angle = round(theta_angle, 3)
        k = (m + n) / (x * np.sin(sigma_radius))
        output_dict = {'theta': theta_angle, 'sigma': sigma_angle, 'm': m, 'n': n, 'x': x, 'y': y_value,'x_vector':x_vector_frac,'k':k}
        return output_dict

    @staticmethod
    def RemoveDuplicateCollinearPoints(points,line_direction):
        line_direction_array = np.array(line_direction)
        lines = []
        remaining_points = points[:]  # Copy the original points list
        while remaining_points:
            point1 = remaining_points[0]
            one_line = [point1]
            points_on_line = []  # Track points to remove later
            for point2 in remaining_points[1:]:  # Iterate over remaining points (excluding point1)
                vector1 = np.array(point1) - np.array(point2)
                cross_value1 = np.cross(vector1, line_direction_array)
                if np.linalg.norm(cross_value1) < 1e-9:  # Use a small tolerance instead of == 0
                    one_line.append(point2)
                    points_on_line.append(point2)  # Mark point2 for removal
            remaining_points.remove(point1)
            for p in points_on_line:
                remaining_points.remove(p)
            lines.append(one_line)
        one_line_one_point = []
        for i in range(len(lines)):
            line_i = lines[i]
            one_line_one_point.append(line_i[0])
        return one_line_one_point


class overlapped_folded_plane:

    def __init__(self,sg_number,molecule_frac,ops,line_direction_frac,a,b,c,alpha,beta,gamma):

        self.sg_number = sg_number
        points_frac = folded_plane.GetFracPointsCoordinates(ops,molecule_frac)
        points_ortho = folded_plane.GetOrthoPointsCoordinate(points_frac,a,b,c,alpha,beta,gamma)
        self.points_frac = list(points_frac)
        self.points_ortho = list(points_ortho)
        self.line_direction_frac = list(line_direction_frac)
        line_direction_ortho = math_method.Orthogonalization(line_direction_frac,a,b,c,alpha,beta,gamma)
        self.line_direction_ortho = list(line_direction_ortho)
        self.cell_parameters = {'a': a, "b": b, 'c': c, 'alpha': alpha, 'beta': beta, 'gamma': gamma}

        points_number = len(ops)
        if points_number == 1:
            point = points_frac[0]
            if list(line_direction_frac) in [[1,0,0],[0,1,0],[0,0,1]]:
                y_value = np.linalg.norm(line_direction_ortho)
                info = self.NormalOnePointOneCell(y_value)
            else:
                info = self.CrossOnePointOneCell(point)
        elif points_number == 2:
            all_collinear,y_vector = folded_plane.WhetherAllCollinear(points_frac,line_direction_frac)
            if all_collinear is True:
                y_ortho = math_method.Orthogonalization(y_vector,a,b,c,alpha,beta,gamma)
                y_value = np.linalg.norm(y_ortho)
                info = self.NormalOnePointOneCell(y_value)
            else:
                y_value = np.linalg.norm(line_direction_ortho)
                info = self.TwoLineMethod(y_value,points_frac)
        else:
            collinear_exist,y_vector,one_line_one_point = folded_plane.WhetherCollinearExists(points_frac,line_direction_frac)

            if collinear_exist is False:
                y_value = np.linalg.norm(line_direction_ortho)
                info = self.GeneralMethod(y_value,points_frac)
            else:
                y_vector_ortho = math_method.Orthogonalization(y_vector,a,b,c,alpha,beta,gamma)
                y_value = np.linalg.norm(y_vector_ortho)
                info = self.GeneralMethod(y_value,one_line_one_point)
        theta_angle = info['theta']
        sigma_angle = info['sigma']
        m_value = info['m']
        n_value = info['n']
        x_value = info['x']
        y_value = info['y']
        x_vector_frac = info['x_vector']
        k_value = info['k']

        self.x_vector = x_vector_frac
        self.theta = theta_angle
        self.sigma = sigma_angle
        self.m = m_value
        self.n = n_value
        self.x = x_value
        self.y = y_value
        self.k = k_value
        return

    def GeneralMethod(self,y_value,points_frac):
        '''
        :param points: should be non-collinear.fractional
        :return:
        '''
        cell_parameters = self.cell_parameters
        points_ortho = folded_plane.GetOrthoPointsCoordinate(points_frac,cell_parameters['a'],cell_parameters['b'], cell_parameters['c'],
                                                                      cell_parameters['alpha'], cell_parameters['beta'],cell_parameters['gamma'])

        enlarge_matrix = [[1,0,0],[0,1,0],[0,0,1]]
        line_direction_frac = self.line_direction_frac
        line_direction_ortho = self.line_direction_ortho
        for i in range(3):
            e = enlarge_matrix[i]
            if e == line_direction_frac:
                enlarge_matrix[i] = [0,0,0]
            else:
                continue

        enlarged_points_frac = cell_method.PointEnlarging(points_frac,enlarge_matrix)

        enlarged_points_ortho = folded_plane.GetOrthoPointsCoordinate(enlarged_points_frac,cell_parameters['a'],cell_parameters['b'], cell_parameters['c'],
                                                                      cell_parameters['alpha'], cell_parameters['beta'],cell_parameters['gamma'])


        points_center_frac = math_method.GetCentroid(enlarged_points_frac,len(enlarged_points_frac))
        points_center_ortho = math_method.Orthogonalization(points_center_frac,cell_parameters['a'],cell_parameters['b'], cell_parameters['c'],
                                                            cell_parameters['alpha'], cell_parameters['beta'],cell_parameters['gamma'])

        distance_to_center = []
        for i in range(len(enlarged_points_ortho)):
            point_temp = enlarged_points_ortho[i]
            distance_temp = math_method.PPdistance(point_temp,points_center_ortho)
            distance_to_center.append(distance_temp)
        min_distance_to_center = min(distance_to_center)
        index_min_center = distance_to_center.index(min_distance_to_center)

        central_point = enlarged_points_ortho[index_min_center]


        # find the nearest two lines
        enlarged_points_ortho.remove(central_point)


        mn_output_info = folded_plane.GetMN(central_point,enlarged_points_ortho,line_direction_ortho)
        m = mn_output_info['m']
        n = mn_output_info['n']
        selected_two_lines = mn_output_info['selected_lines']
        x_pair_points_ortho,x_cell_direction = folded_plane.XPointsPairFix(cell_parameters,selected_two_lines,line_direction_frac)
        print('x_pair_points_ortho',x_pair_points_ortho)
        x_direction_ortho = list(np.array(x_pair_points_ortho[0]) - np.array(x_pair_points_ortho[1]))
        x = np.linalg.norm(x_direction_ortho)
        sigma_radius, sigma_label = folded_plane.GetSigma(x_direction_ortho,line_direction_ortho,cell_parameters['alpha'],
                                                          cell_parameters['beta'], cell_parameters['gamma'])
        if sigma_label != 'None':
            sigma_angle = cell_parameters[sigma_label]
        else:
            sigma_angle = (sigma_radius * 180) / np.pi
        sg_number = self.sg_number
        theta_radius = folded_plane.GetTheta(sg_number,x,y_value,m,n,sigma_radius)
        theta_angle = (theta_radius * 180) / np.pi
        theta_angle = round(theta_angle, 3)
        k = (m + n) / (x * np.sin(sigma_radius))
        output_dict = {'theta': theta_angle, 'sigma': sigma_angle, 'm': m, 'n': n, 'x': x, 'y': y_value,'k':k,'x_vector':x_cell_direction}
        return output_dict

    def TwoLineMethod(self,y_value,points_frac):
        def enlarging_cell(points,line_direction_frac):
            all_direction = [[1,0,0],[0,1,0],[0,0,1]]
            all_direction.remove(line_direction_frac)
            direction1_array = np.array(all_direction[0])
            direction2_array = np.array(all_direction[1])
            enlarged_points_frac = []
            for point in points:
                point_array = np.array(point)
                point11 = point
                enlarged_points_frac.append(point11)
                point12 = list(point_array+direction1_array)
                enlarged_points_frac.append(point12)
                point13 = list(point_array-direction1_array)
                enlarged_points_frac.append(point13)
                point14 = list(point_array+direction2_array)
                enlarged_points_frac.append(point14)
                point15 = list(point_array-direction2_array)
                enlarged_points_frac.append(point15)
            return enlarged_points_frac

        line_direction_frac = self.line_direction_frac
        line_direction_ortho = self.line_direction_ortho
        cell_parameters = self.cell_parameters
        point1 = points_frac[0]
        point1_ortho = math_method.Orthogonalization(point1,cell_parameters['a'],cell_parameters['b'], cell_parameters['c'],
                                                    cell_parameters['alpha'], cell_parameters['beta'],cell_parameters['gamma'])

        enlarged_points_frac = enlarging_cell(points_frac,line_direction_frac)
        enlarged_points_frac.remove(point1)
        enlarged_points_ortho = folded_plane.GetOrthoPointsCoordinate(enlarged_points_frac, cell_parameters['a'],
                                                                      cell_parameters['b'], cell_parameters['c'],
                                                                      cell_parameters['alpha'], cell_parameters['beta'],
                                                                      cell_parameters['gamma'])

        mn_output_info = folded_plane.GetMN(point1_ortho, enlarged_points_ortho, line_direction_ortho)
        m = mn_output_info['m']
        n = mn_output_info['n']
        selected_two_lines = mn_output_info['selected_lines']
        x_pair_points_ortho, x_cell_direction = folded_plane.XPointsPairFix(cell_parameters, selected_two_lines,
                                                                            line_direction_frac)

        x_direction_ortho = list(np.array(x_pair_points_ortho[0]) - np.array(x_pair_points_ortho[1]))
        x = np.linalg.norm(x_direction_ortho)
        sigma_radius, sigma_label = folded_plane.GetSigma(x_direction_ortho, line_direction_ortho,
                                                          cell_parameters['alpha'],
                                                          cell_parameters['beta'], cell_parameters['gamma'])

        if sigma_label != 'None':
            sigma_angle = cell_parameters[sigma_label]
        else:
            sigma_angle = (sigma_radius * 180) / np.pi
        sg_number = self.sg_number
        theta_radius = folded_plane.GetTheta(sg_number, x, y_value, m, n, sigma_radius)
        theta_angle = (theta_radius * 180) / np.pi
        theta_angle = round(theta_angle, 3)
        k = (m + n) / (x * np.sin(sigma_radius))
        output_dict = {'theta': theta_angle, 'sigma': sigma_angle, 'm': m, 'n': n, 'x': x, 'y': y_value, 'k': k,
                       'x_vector': x_cell_direction}
        return output_dict

    def NormalOnePointOneCell(self,y_value):
        '''
        for line direction which is parallel with [[1,0,0],[0,1,0],[0,0,1],[-1,0,0],[0,-1,0],[0,0,-1]]
        :return:
        '''
        points_frac = self.points_frac
        points_ortho = self.points_ortho
        point_frac = points_frac[0]
        point_ortho = points_ortho[0]
        line_direction_frac = self.line_direction_frac
        cell_parameters = self.cell_parameters
        a = cell_parameters['a']
        b = cell_parameters['b']
        c = cell_parameters['c']
        alpha = cell_parameters['alpha']
        beta = cell_parameters['beta']
        gamma = cell_parameters['gamma']
        line_direction_ortho = math_method.Orthogonalization(line_direction_frac, a, b, c, alpha, beta, gamma)

        enlarge_direction = [[1,0,0],[0,1,0],[0,0,1]]
        enlarge_direction.remove(line_direction_frac)
        direction1 = enlarge_direction[0]
        direction2 = enlarge_direction[1]
        direction3 = list(np.array(direction1) + np.array(direction2))
        direction4 = list(np.array(direction1) - np.array(direction2))
        potential_direction = [direction1,direction2,direction3,direction4]

        surrounding_points = []
        for i in range(len(potential_direction)):
            direction_i = potential_direction[i]
            x = direction_i[0] + point_frac[0]
            y = direction_i[1] + point_frac[1]
            z = direction_i[2] + point_frac[2]
            temp_point = [x, y, z]
            surrounding_points.append(temp_point)
        surrounding_points_ortho = []
        for i in range(len(surrounding_points)):
            surround_i = surrounding_points[i]
            surround_ortho = math_method.Orthogonalization(surround_i,a,b,c,alpha,beta,gamma)
            surrounding_points_ortho.append(surround_ortho)

        line_line_distance = []
        math_tools = math_method()
        for i in range(4):
            pointi = surrounding_points_ortho[i]
            ll_distance = math_tools.LineLineDistance(point_ortho,pointi,line_direction_ortho)
            line_line_distance.append(ll_distance)

        lengths2 = line_line_distance.copy()
        mn = sorted(lengths2)[:2]
        min1 = mn[0]
        min2 = mn[1]

        if min1 != min2: # flat
            print('m=!n',min1)
            m = n = min1
            m_lengths_index = line_line_distance.index(min1)
            m_vector = potential_direction[m_lengths_index]
            x_vector = list(2*np.array(m_vector))
            x = 2*m

            x_vector_ortho = math_method.Orthogonalization(x_vector,a,b,c,alpha,beta,gamma)
            sigma_radius,sigma_label = folded_plane.GetSigma(x_vector_ortho,line_direction_ortho,alpha,beta,gamma)

        else:
            indexes = [index for index, element in enumerate(line_line_distance) if element == min1]
            m_lengths_index = indexes[0]
            n_lengths_index = indexes[1]
            m = min1
            n = min2

            m_vector = potential_direction[m_lengths_index]
            n_vector = potential_direction[n_lengths_index]
            x_vector = list(np.array(n_vector) - np.array(m_vector))
            x_vector_ortho = math_method.Orthogonalization(x_vector, a, b, c, alpha, beta, gamma)
            x = np.linalg.norm(x_vector_ortho)
            sigma_radius, sigma_label = folded_plane.GetSigma(x_vector_ortho, line_direction_ortho, alpha, beta, gamma)
            print(sigma_radius,sigma_label)

        sg_number = self.sg_number
        theta_radius = folded_plane.GetTheta(sg_number,x,y_value,m,n,sigma_radius)

        theta_angle = (theta_radius * 180) / np.pi
        theta_angle = round(theta_angle, 3)

        k = (m + n) / (x * np.sin(sigma_radius))
        sigma_angle = sigma_radius*180/np.pi
        output_dict = {'theta': theta_angle, 'sigma': sigma_angle, 'm': m, 'n': n, 'x': x, 'y': y_value,'x_vector':x_vector,'k':k}
        return output_dict

    def CrossOnePointOneCell(self,point_frac):
        cell_parameters = self.cell_parameters
        line_direction_frac = self.line_direction_frac
        line_direction_ortho = self.line_direction_ortho
        y_value = np.linalg.norm(line_direction_ortho)

        all_directions = [[1,0,0],[0,1,0],[0,0,1],
                          [1,1,0],[1,-1,0],[1,0,1],[1,0,-1],[0,1,1],[0,1,-1],
                          [1,1,1],[-1,1,1],[1,-1,1],[1,1,-1]]
        all_directions.remove(line_direction_frac)

        one_line_one_point = folded_plane.RemoveDuplicateCollinearPoints(all_directions,line_direction_frac)
        enlarge_directions = []
        for i in range(len(one_line_one_point)):
            linei = one_line_one_point[i]
            enlarge_directions.append(linei[0])

        # enlarge the point to more than one
        enlarged_points = [point_frac]
        for i in range(len(enlarge_directions)):
            direction_i = enlarge_directions[i]
            x = direction_i[0]+point_frac[0]
            y = direction_i[1]+point_frac[1]
            z = direction_i[2]+point_frac[2]
            temp_point = [x,y,z]
            enlarged_points.append(temp_point)

        # do point orthogonal
        points_ortho = []
        for i in range(len(enlarged_points)):
            point_i = enlarged_points[i]
            point_i_ortho = math_method.Orthogonalization(point_i,cell_parameters['a'],cell_parameters['b'],cell_parameters['c'],
                                                          cell_parameters['alpha'],cell_parameters['beta'],cell_parameters['gamma'])
            points_ortho.append(point_i_ortho)

        center_point = math_method.GetCentroid(points_ortho,len(points_ortho))
        # find the nearest point to the points center
        pp_distances = []
        for i in range(len(points_ortho)):
            distance_i = math_method.PPdistance(center_point,points_ortho[i])
            pp_distances.append(distance_i)
        min_ppdistance = min(pp_distances)
        point_index = pp_distances.index(min_ppdistance)
        most_center_point = points_ortho(point_index)

        math_tools = math_method()
        index_most_center_point = points_ortho.index(most_center_point)
        points_ortho.remove(most_center_point)
        removed_frac_point = enlarged_points[index_most_center_point]
        enlarged_points.remove(removed_frac_point)

        line_line_distances = []
        for i in range(len(points_ortho)):
            ll_distance = math_tools.LineLineDistance(most_center_point,points_ortho[i],line_direction_ortho)
            line_line_distances.append(ll_distance)

        ll_distance_list_copy = line_line_distances.copy()
        mn = sorted(ll_distance_list_copy)[:2]
        min1 = mn[0]
        min2 = mn[1]

        if min1 != min2:  # flat
            m = n = min1
            m_lengths_index = lengths.index(min1)
            m_vector = list(np.array(enlarged_points[m_lengths_index])-np.array(most_center_point))
            x_vector = 2 * m_vector
            x = 2*m
            x_vector_ortho = math_method.Orthogonalization(x_vector,a,b,c,alpha,beta,gamma)
            sigma_radius,sigma_label = folded_plane.GetSigma(x_vector_ortho,line_direction_ortho,alpha,beta,gamma)

        else:
            indexes = [index for index, element in enumerate(line_line_distances) if element == min1]
            m_lengths_index = indexes[0]
            n_lengths_index = indexes[1]
            m = min1
            n = min2
            m_vector = list(np.array(enlarged_points[m_lengths_index])-np.array(most_center_point))
            n_vector = list(np.array(enlarged_points[n_lengths_index])-np.array(most_center_point))
            x_vector = list(np.array(n_vector) - np.array(m_vector))
            x_vector_ortho = math_method.Orthogonalization(x_vector, a, b, c, alpha, beta, gamma)
            sigma_radius, sigma_label = folded_plane.GetSigma(x_vector_ortho, line_direction_ortho, alpha, beta, gamma)

        sg_number = self.sg_number
        theta_radius = folded_plane.GetTheta(sg_number,x,y_value,m,n,sigma_radius)
        theta_angle = (theta_radius * 180) / np.pi
        theta_angle = round(theta_angle, 3)
        k = (m + n) / (x * np.sin(sigma_radius))
        output_dict = {'theta': theta_angle, 'sigma': sigma_angle, 'm': m, 'n': n, 'x': x, 'y': y_value,'x_vector':x_vector,'k':k}
        return output_dict
