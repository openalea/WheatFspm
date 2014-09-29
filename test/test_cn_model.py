# -*- coding: latin-1 -*-
'''
Test cnwheat.py. 

CSV files must contain only ASCII characters and ',' as separator.

Be sure to add the 'src' directory to your PYTHONPATH before running this script.
'''

import os

import numpy as np
import pandas as pd

from cnwheat import cnwheat
from cnwheat import organ

DATA_DIRPATH = 'data'
DESIRED_OUTPUT_FILENAME = 'desired_output.csv'
ACTUAL_OUTPUT_FILENAME = 'actual_output.csv'


RELATIVE_TOLERANCE = 10e-3
ABSOLUTE_TOLERANCE = 10e-3


def read_t_data(curr_data_dirpath, data_filename):
    data_filepath = os.path.join(curr_data_dirpath, data_filename)
    return pd.read_csv(data_filepath, index_col='t')


def compare_actual_to_desired(curr_data_dirpath, actual_output_df, save_actual_output=False):
    # read desired output
    desired_output_filepath = os.path.join(curr_data_dirpath, DESIRED_OUTPUT_FILENAME)
    desired_output_df = pd.read_csv(desired_output_filepath)
      
    # keep only the columns to test
    actual_output_df = actual_output_df[desired_output_df.columns]
      
    if save_actual_output:
        actual_output_filepath = os.path.join(curr_data_dirpath, ACTUAL_OUTPUT_FILENAME)
        actual_output_df.to_csv(actual_output_filepath, na_rep='NA', index=False)
      
    # compare to the desired output
    np.testing.assert_allclose(actual_output_df.values, desired_output_df.values, RELATIVE_TOLERANCE, ABSOLUTE_TOLERANCE)
    
    
def test_minimalist_system():
    curr_data_dirpath = os.path.join(DATA_DIRPATH, 'minimalist_system')   
    
    # create the chaff
    name='chaff'
    PAR_df = read_t_data(curr_data_dirpath, 'PAR_%s.csv' % name)
    chaff = organ.Chaff(Area=0.00075, Mstruct=0.21, PAR=PAR_df.PAR, 
                        STORAGE_0=0, SUCROSE_0=0, TRIOSESP_0=0, name=name)
     
    # create the internodes
    number_of_internodes = 3
    Areas = [(0.0012, 0.0003), 0.0004, 0.00025]
    Mstructs = [(0.148, 0.04), 0.18, 0.154]
    internodes = []
    for i in xrange(number_of_internodes):
        internode_index = i + 1
        name = 'internode%d' % internode_index
        if internode_index == 1:
            current_names = (name + '_enclosed', name + '_exposed')
            current_Areas = Areas[i]
            current_Mstructs = Mstructs[i]
        else:
            current_names = (name,)
            current_Areas = (Areas[i],)
            current_Mstructs = (Mstructs[i],)
        for j in xrange(len(current_names)):
            PAR_df = read_t_data(curr_data_dirpath, 'PAR_%s.csv' % current_names[j])
            internodes.append(organ.Internode(Area=current_Areas[j], Mstruct=current_Mstructs[j], 
                                              PAR=PAR_df.PAR, FRUCTAN_0=0, STORAGE_0=0, SUCROSE_0=0, 
                                              TRIOSESP_0=0, name=current_names[j]))
             
    # create the laminae
    number_of_laminae = 3
    Areas = [0.00346, 0.0034, 0.00228]
    Mstructs = [0.14, 0.09, 0.05]
    Rdarks = [1.5, 0.0, 0.0]
    laminae = []
    for i in xrange(number_of_laminae):
        name = 'lamina%d' % (i + 1)
        PAR_df = read_t_data(curr_data_dirpath, 'PAR_%s.csv' % name)
        laminae.append(organ.Lamina(Area=Areas[i], Mstruct=Mstructs[i], Rdark=Rdarks[i],
                                    PAR=PAR_df.PAR, 
                                    STORAGE_0=0, SUCROSE_0=0, TRIOSESP_0=0,
                                    name=name))
         
    # create the peduncles
    peduncles = []
    names = ['peduncle_enclosed', 'peduncle_exposed']
    Areas = [0.00155, 0.00085]
    Mstructs = [0.168, 0.089]
    for i in xrange(len(names)):
        PAR_df = read_t_data(curr_data_dirpath, 'PAR_%s.csv' % names[i])
        peduncles.append(organ.Peduncle(Area=Areas[i], Mstruct=Mstructs[i], PAR=PAR_df.PAR, 
                                        FRUCTAN_0=0, STORAGE_0=0, SUCROSE_0=0, TRIOSESP_0=0, 
                                        name=names[i]))
    
    # create the sheaths
    number_of_sheaths = 3
    sheaths = []
    Areas = [0.0006, 0.0005, 0.0004]
    Mstructs = [0.103, 0.069, 0.043]
    for i in xrange(number_of_sheaths):
        name = 'sheath%d' % (i + 1)
        PAR_df = read_t_data(curr_data_dirpath, 'PAR_%s.csv' % name)
        sheaths.append(organ.Sheath(Area=Areas[i], Mstruct=Mstructs[i], PAR=PAR_df.PAR, 
                                    FRUCTAN_0=0, STORAGE_0=0, SUCROSE_0=0, TRIOSESP_0=0, name=name))
         
    # create the grains
    grains = organ.Grains(STORAGE_0=0, STRUCTURE_0=10850, name='grains')
     
    # create the roots
    roots = organ.Roots(Mstruct=0.504, Sucrose_0=0, name='roots')
     
    # create the phloem
    phloem = organ.Phloem(SUCROSE_0=0, name='phloem')
    
    # get meteo data
    meteo_df = read_t_data(curr_data_dirpath, 'meteo.csv')
     
    # run the model
    actual_output_df = cnwheat.run(start_time=0, stop_time=48, number_of_output_steps=7,
                                    phloem=phloem,
                                    organs=[chaff]+internodes+laminae+peduncles+sheaths+[grains]+[roots], 
                                    meteo=meteo_df,
                                    photosynthesis_computation_interval=4)
    
    compare_actual_to_desired(curr_data_dirpath, actual_output_df)
    

if __name__ == '__main__':
    test_minimalist_system()


