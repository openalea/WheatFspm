# -*- coding: latin-1 -*-
'''
Test cn_model.py. 

CSV files must contain only ASCII characters and ',' as separator.

Be sure to add the 'src' directory to your PYTHONPATH before running this script.
'''

import os

import numpy as np
import pandas as pd

import organ
import cn_model

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
    actual_output_df = cn_model.run(start_time=0, stop_time=24, number_of_output_steps=7,
                                    organs=[chaff])
    
    compare_actual_to_desired(curr_data_dirpath, actual_output_df)
     
    
def test_Internode():
    curr_data_dirpath = os.path.join(DATA_DIRPATH, 'internode')
     
    # create the internode
    Assimilation_df = read_assimilation_file(curr_data_dirpath, 'Assimilation.csv')
    internode = organ.Internode(Area=0.0004, Mstruct=0.18, Assimilation=Assimilation_df.An, 
                                FRUCTAN_0=0, STORAGE_0=0, SUCROSE_0=0, TRIOSESP_0=0)
     
    # run the model
    actual_output_df = cn_model.run(start_time=0, stop_time=24, number_of_output_steps=7,
                                    organs=[internode])
     
    compare_actual_to_desired(curr_data_dirpath, actual_output_df)
    

def test_Lamina():
    curr_data_dirpath = os.path.join(DATA_DIRPATH, 'lamina')
     
    # create the lamina
    Assimilation_df = read_assimilation_file(curr_data_dirpath, 'Assimilation.csv')
    lamina = organ.Lamina(Area=0.00346, Mstruct=0.14, Rdark=1.5, Assimilation=Assimilation_df.An, 
                          STORAGE_0=0, SUCROSE_0=0, TRIOSESP_0=0)
     
    # run the model
    actual_output_df = cn_model.run(start_time=0, stop_time=24, number_of_output_steps=7,
                                    organs=[lamina])
    
    compare_actual_to_desired(curr_data_dirpath, actual_output_df)
    

def test_Peduncle():
    curr_data_dirpath = os.path.join(DATA_DIRPATH, 'peduncle')
     
    # create the peduncle
    Assimilation_df = read_assimilation_file(curr_data_dirpath, 'Assimilation.csv')
    peduncle = organ.Peduncle(Area=0.00155, Mstruct=0.168, Assimilation=Assimilation_df.An, 
                                FRUCTAN_0=0, STORAGE_0=0, SUCROSE_0=0, TRIOSESP_0=0)
     
    # run the model
    actual_output_df = cn_model.run(start_time=0, stop_time=96, number_of_output_steps=13,
                                    organs=[peduncle])
     
    compare_actual_to_desired(curr_data_dirpath, actual_output_df)
    
    
def test_Sheath():
    curr_data_dirpath = os.path.join(DATA_DIRPATH, 'sheath')
     
    # create the sheath
    Assimilation_df = read_assimilation_file(curr_data_dirpath, 'Assimilation.csv')
    sheath = organ.Sheath(Area=0.0006, Mstruct=0.103, Assimilation=Assimilation_df.An, 
                          FRUCTAN_0=0, STORAGE_0=0, SUCROSE_0=0, TRIOSESP_0=0)
     
    # run the model
    actual_output_df = cn_model.run(start_time=0, stop_time=24, number_of_output_steps=7,
                                    organs=[sheath])
     
    compare_actual_to_desired(curr_data_dirpath, actual_output_df)
    
    
def test_Grains():
    curr_data_dirpath = os.path.join(DATA_DIRPATH, 'grains')
     
    # create the grains
    grains = organ.Grains(STORAGE_0=0, STRUCTURE_0=10850)
     
    # run the model
    actual_output_df = cn_model.run(start_time=0, stop_time=24, number_of_output_steps=7,
                                    organs=[grains])
     
    compare_actual_to_desired(curr_data_dirpath, actual_output_df)
    
    
def test_Roots():
    curr_data_dirpath = os.path.join(DATA_DIRPATH, 'roots')
     
    # create the roots
    roots = organ.Roots(Mstruct=0.504, Sucrose_0=0)
     
    # run the model
    actual_output_df = cn_model.run(start_time=0, stop_time=24, number_of_output_steps=7,
                                    organs=[roots])
     
    compare_actual_to_desired(curr_data_dirpath, actual_output_df)
    
    
def test_Phloem():
    curr_data_dirpath = os.path.join(DATA_DIRPATH, 'phloem')
     
    # create the phloem
    phloem = organ.Phloem(SUCROSE_0=0, Respiration_0=0)
     
    # run the model
    actual_output_df = cn_model.run(start_time=0, stop_time=24, number_of_output_steps=7,
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
    actual_output_df = cn_model.run(start_time=0, stop_time=96, number_of_output_steps=13,
                                    organs=[chaff]+[internode]+[lamina]+[peduncle]+[sheath]+[grains]+[roots], 
                                    phloem=phloem)
    
    compare_actual_to_desired(curr_data_dirpath, actual_output_df)

   
def test_complex_system():
    '''This system is equivalent to the system Distribution_CN_postflo_V2.mbk from ModelMaker'''
    curr_data_dirpath = os.path.join(DATA_DIRPATH, 'complex_system')
    
    # create the chaff
    name='chaff'
    curr_Assimilation_df = read_assimilation_file(curr_data_dirpath, 'Assimilation_%s.csv' % name)
    chaff = organ.Chaff(Area=0.00075, Mstruct=0.21, Assimilation=curr_Assimilation_df.An, 
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
            curr_Assimilation_df = read_assimilation_file(curr_data_dirpath, 'Assimilation_%s.csv' % current_names[j])
            internodes.append(organ.Internode(Area=current_Areas[j], Mstruct=current_Mstructs[j], 
                                              Assimilation=curr_Assimilation_df.An, FRUCTAN_0=0, STORAGE_0=0, SUCROSE_0=0, 
                                              TRIOSESP_0=0, name=current_names[j]))
            
    # create the laminae
    number_of_laminae = 3
    Areas = [0.00346, 0.0034, 0.00228]
    Mstructs = [0.14, 0.09, 0.05]
    Rdarks = [1.5, 0.0, 0.0]
    laminae = []
    for i in xrange(number_of_laminae):
        name = 'lamina%d' % (i + 1)
        curr_Assimilation_df = read_assimilation_file(curr_data_dirpath, 'Assimilation_%s.csv' % name)
        laminae.append(organ.Lamina(Area=Areas[i], Mstruct=Mstructs[i], Rdark=Rdarks[i],
                                    Assimilation=curr_Assimilation_df.An, 
                                    STORAGE_0=0, SUCROSE_0=0, TRIOSESP_0=0,
                                    name=name))
        
    # create the peduncles
    peduncles = []
    names = ['peduncle_enclosed', 'peduncle_exposed']
    Areas = [0.00155, 0.00085]
    Mstructs = [0.168, 0.089]
    for i in xrange(len(names)):
        curr_Assimilation_df = read_assimilation_file(curr_data_dirpath, 'Assimilation_%s.csv' % names[i])
        peduncles.append(organ.Peduncle(Area=Areas[i], Mstruct=Mstructs[i], Assimilation=curr_Assimilation_df.An, 
                                        FRUCTAN_0=0, STORAGE_0=0, SUCROSE_0=0, TRIOSESP_0=0, 
                                        name=names[i]))

    # create the sheaths
    number_of_sheaths = 3
    sheaths = []
    Areas = [0.0006, 0.0005, 0.0004]
    Mstructs = [0.103, 0.069, 0.043]
    for i in xrange(number_of_sheaths):
        name = 'sheath%d' % (i + 1)
        curr_Assimilation_df = read_assimilation_file(curr_data_dirpath, 'Assimilation_%s.csv' % name)
        sheaths.append(organ.Sheath(Area=Areas[i], Mstruct=Mstructs[i], Assimilation=curr_Assimilation_df.An, 
                                    FRUCTAN_0=0, STORAGE_0=0, SUCROSE_0=0, TRIOSESP_0=0, name=name))
        
    # create the grains
    grains = organ.Grains(STORAGE_0=0, STRUCTURE_0=10850, name='grains')
    
    # create the roots
    roots = organ.Roots(Mstruct=0.504, Sucrose_0=0, name='roots')
    
    # create the phloem
    phloem = organ.Phloem(SUCROSE_0=0, Respiration_0=0, name='phloem')
    
    # run the model
    actual_output_df = cn_model.run(start_time=0, stop_time=960, number_of_output_steps=241,
                                    organs=[chaff]+internodes+laminae+peduncles+sheaths+[grains]+[roots], 
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
    test_complex_system()
    

