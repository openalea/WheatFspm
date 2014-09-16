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


def read_assimilation_file(curr_data_dirpath, Assimilation_filename):
    Assimilation_filepath = os.path.join(curr_data_dirpath, Assimilation_filename)
    return pd.read_csv(Assimilation_filepath, index_col='t')

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
    

def test_Chaff():
    curr_data_dirpath = os.path.join(DATA_DIRPATH, 'chaff')
     
    # create the chaff
    Assimilation_df = read_assimilation_file(curr_data_dirpath, 'Assimilation.csv')
    chaff = organ.Chaff(Area=0.00075, Mstruct=0.21, Assimilation=Assimilation_df.An, 
                        STORAGE_0=0, SUCROSE_0=0, TRIOSESP_0=0)
     
    # run the model
    actual_output_df = cnwheat.run(start_time=0, stop_time=24, number_of_output_steps=7,
                                    organs=[chaff])
    
    compare_actual_to_desired(curr_data_dirpath, actual_output_df)
     
    
def test_Internode():
    curr_data_dirpath = os.path.join(DATA_DIRPATH, 'internode')
      
    # create the internode
    Assimilation_df = read_assimilation_file(curr_data_dirpath, 'Assimilation.csv')
    internode = organ.Internode(Area=0.0004, Mstruct=0.18, Assimilation=Assimilation_df.An, 
                                FRUCTAN_0=0, STORAGE_0=0, SUCROSE_0=0, TRIOSESP_0=0)
      
    # run the model
    actual_output_df = cnwheat.run(start_time=0, stop_time=24, number_of_output_steps=7,
                                    organs=[internode])
      
    compare_actual_to_desired(curr_data_dirpath, actual_output_df)
     
 
def test_Lamina():
    curr_data_dirpath = os.path.join(DATA_DIRPATH, 'lamina')
      
    # create the lamina
    Assimilation_df = read_assimilation_file(curr_data_dirpath, 'Assimilation.csv')
    lamina = organ.Lamina(Area=0.00346, Mstruct=0.14, Rdark=1.5, Assimilation=Assimilation_df.An, 
                          STORAGE_0=0, SUCROSE_0=0, TRIOSESP_0=0)
      
    # run the model
    actual_output_df = cnwheat.run(start_time=0, stop_time=24, number_of_output_steps=7,
                                    organs=[lamina])
     
    compare_actual_to_desired(curr_data_dirpath, actual_output_df)
     
 
def test_Peduncle():
    curr_data_dirpath = os.path.join(DATA_DIRPATH, 'peduncle')
      
    # create the peduncle
    Assimilation_df = read_assimilation_file(curr_data_dirpath, 'Assimilation.csv')
    peduncle = organ.Peduncle(Area=0.00155, Mstruct=0.168, Assimilation=Assimilation_df.An, 
                                FRUCTAN_0=0, STORAGE_0=0, SUCROSE_0=0, TRIOSESP_0=0)
      
    # run the model
    actual_output_df = cnwheat.run(start_time=0, stop_time=96, number_of_output_steps=13,
                                    organs=[peduncle])
      
    compare_actual_to_desired(curr_data_dirpath, actual_output_df)
     
     
def test_Sheath():
    curr_data_dirpath = os.path.join(DATA_DIRPATH, 'sheath')
      
    # create the sheath
    Assimilation_df = read_assimilation_file(curr_data_dirpath, 'Assimilation.csv')
    sheath = organ.Sheath(Area=0.0006, Mstruct=0.103, Assimilation=Assimilation_df.An, 
                          FRUCTAN_0=0, STORAGE_0=0, SUCROSE_0=0, TRIOSESP_0=0)
      
    # run the model
    actual_output_df = cnwheat.run(start_time=0, stop_time=24, number_of_output_steps=7,
                                    organs=[sheath])
      
    compare_actual_to_desired(curr_data_dirpath, actual_output_df)
     
     
def test_Grains():
    curr_data_dirpath = os.path.join(DATA_DIRPATH, 'grains')
      
    # create the grains
    grains = organ.Grains(STORAGE_0=0, STRUCTURE_0=10850)
      
    # run the model
    actual_output_df = cnwheat.run(start_time=0, stop_time=24, number_of_output_steps=7,
                                    organs=[grains])
      
    compare_actual_to_desired(curr_data_dirpath, actual_output_df)
     
     
def test_Roots():
    curr_data_dirpath = os.path.join(DATA_DIRPATH, 'roots')
      
    # create the roots
    roots = organ.Roots(Mstruct=0.504, Sucrose_0=0)
      
    # run the model
    actual_output_df = cnwheat.run(start_time=0, stop_time=24, number_of_output_steps=7,
                                    organs=[roots])
      
    compare_actual_to_desired(curr_data_dirpath, actual_output_df)
     
     
def test_Phloem():
    curr_data_dirpath = os.path.join(DATA_DIRPATH, 'phloem')
      
    # create the phloem
    phloem = organ.Phloem(SUCROSE_0=0, Respiration_0=0)
      
    # run the model
    actual_output_df = cnwheat.run(start_time=0, stop_time=24, number_of_output_steps=7,
                                    phloem=phloem)
      
    compare_actual_to_desired(curr_data_dirpath, actual_output_df)
 
 
def test_simplified_system():
    curr_data_dirpath = os.path.join(DATA_DIRPATH, 'simplified_system')
     
    # create the chaff
    name='chaff'
    Assimilation_df = read_assimilation_file(curr_data_dirpath, 'Assimilation_%s.csv' % name)
    chaff = organ.Chaff(Area=0.00075, Mstruct=0.21, Assimilation=Assimilation_df.An, 
                        STORAGE_0=0, SUCROSE_0=0, TRIOSESP_0=0, name=name)
     
    # create the internode
    name='internode'
    Assimilation_df = read_assimilation_file(curr_data_dirpath, 'Assimilation_%s.csv' % name)
    internode = organ.Internode(Area=0.0004, Mstruct=0.18, Assimilation=Assimilation_df.An, 
                                FRUCTAN_0=0, STORAGE_0=0, SUCROSE_0=0, TRIOSESP_0=0, name=name)
     
    # create the lamina
    name='lamina'
    Assimilation_df = read_assimilation_file(curr_data_dirpath, 'Assimilation_%s.csv' % name)
    lamina = organ.Lamina(Area=0.00346, Mstruct=0.14, Rdark=1.5, Assimilation=Assimilation_df.An, 
                          STORAGE_0=0, SUCROSE_0=0, TRIOSESP_0=0, name=name)
         
    # create the peduncle
    name='peduncle'
    Assimilation_df = read_assimilation_file(curr_data_dirpath, 'Assimilation_%s.csv' % name)
    peduncle = organ.Peduncle(Area=0.00155, Mstruct=0.168, Assimilation=Assimilation_df.An, 
                              FRUCTAN_0=0, STORAGE_0=0, SUCROSE_0=0, TRIOSESP_0=0, name=name)
     
    # create the sheath
    name='sheath'
    Assimilation_df = read_assimilation_file(curr_data_dirpath, 'Assimilation_%s.csv' % name)
    sheath = organ.Sheath(Area=0.0006, Mstruct=0.103, Assimilation=Assimilation_df.An, 
                          FRUCTAN_0=0, STORAGE_0=0, SUCROSE_0=0, TRIOSESP_0=0, name=name)
     
    # create the grains
    grains = organ.Grains(STORAGE_0=0, STRUCTURE_0=10850, name='grains')
     
    # create the roots
    roots = organ.Roots(Mstruct=0.504, Sucrose_0=0, name='roots')
     
    # create the phloem
    phloem = organ.Phloem(SUCROSE_0=0, Respiration_0=0, name='phloem')
     
    # run the model
    actual_output_df = cnwheat.run(start_time=0, stop_time=96, number_of_output_steps=13,
                                    organs=[chaff]+[internode]+[lamina]+[peduncle]+[sheath]+[grains]+[roots], 
                                    phloem=phloem)
     
    compare_actual_to_desired(curr_data_dirpath, actual_output_df)
    

if __name__ == '__main__':
    test_Chaff()
    test_Internode()
    test_Lamina()
    test_Peduncle()
    test_Sheath()
    test_Grains()
    test_Roots()
    test_Phloem()
    test_simplified_system()
    

