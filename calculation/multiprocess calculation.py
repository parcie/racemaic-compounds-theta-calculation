import multiprocessing as mp
import pandas as pd
import json
import csv
from one_sample_calculation import a_racemic_sample
import openpyxl
import numpy as np

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

def get_ccdc_code_list(file_name):
    '''save the ccdc code list in a file and read it'''
    with open(file_name,mode='r') as file:
        csv_reader = csv.reader(file)
        code_list = list(csv_reader)
    return code_list

def output_one_sample_value(args):
    ccdc_code, sg_table = args  # Unpack the tuple
    try:
        racemic1 = a_racemic_sample(ccdc_code, sg_table)
        sg_name = racemic1.sg_name
        m = racemic1.m
        n = racemic1.n
        x = racemic1.x
        y = racemic1.y
        k = racemic1.k
        sigma = racemic1.sigma
        theta = racemic1.theta
        is_overlap = racemic1.is_overlap
        line_direction = racemic1.line_direction
        line_direction_string = str(line_direction)
        single_info = [ccdc_code, sg_name, is_overlap, line_direction_string, m, n, x, y, k, sigma,theta]
        return single_info
    except Exception as e:
        # Log the error and continue
        print(f"Error processing sample {args[0]}: {e}")
        return [ccdc_code,'error','error','error','error','error','error','error','error','error','error']  # Return None to skip the failed sample

def parallel_compute(ccdc_code_list,sg_file_name):

    sg_info_table = sg_table_read(sg_file_name)
    num_cores = mp.cpu_count()
    print('num_cores',num_cores)
    args = [(code, sg_info_table) for code in ccdc_code_list]

    with mp.Pool(num_cores) as pool:
        result_list = pool.map(output_one_sample_value,args)

    column_names = ['ccdc_code', 'spacegroup', 'is_overlapped', 'line_direction', 'm', 'n', 'x', 'y', 'k', 'sigma','theta']
    valid_results = [row for row in result_list if row and len(row) == len(column_names)]

    print(type(valid_results))
    print(len(valid_results))
    print('valid_results[0]',valid_results[0])
    print('result',valid_results)
    df1 = pd.DataFrame(valid_results,columns=column_names)
    return df1

def Generate_codelist_from_xlsx(file_path):
    xlsx_sg = file_path
    wb_sg = openpyxl.load_workbook(xlsx_sg)
    sheet_sg = wb_sg.active
    ccdc_codes = [cell.value for cell in sheet_sg['A']]
    ccdc_codes.remove('Refcode')
    return ccdc_codes

def Generate_codelist_from_csv(file_path):
    df = pd.read_csv(file_path)
    series = df['Refcode']
    code_list = list(series)
    return code_list


if __name__ == '__main__':
    sg_file = '/Users/chenfangyi/Documents/勉強する/theta_calculation/sg_test/1111.json'
    ccdc_code_list = Generate_codelist_from_csv('/Users/chenfangyi/Documents/勉強する/theta_calculation/chipi_running_result/crystal_classification/primary_classification_result/achiral_classification/sg_separation/group_56.csv')

    df = parallel_compute(ccdc_code_list, sg_file)
    df.to_csv('/Users/chenfangyi/Documents/勉強する/theta_calculation/chipi_running_result/crystal_classification/primary_classification_result/achiral_classification/sg_separation/sg43.csv', index=False)