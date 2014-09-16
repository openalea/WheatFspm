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
MODELMAKER_OUTPUT_FILENAME = 'ModelMaker_output.csv'
PYTHON_OUTPUT_FILENAME = 'python_output.csv'


relative_tolerance = 10e-3
absolute_tolerance = 10e-3


def read_assimilation_file(curr_data_dirpath, Assimilation_filename):
    Assimilation_filepath = os.path.join(curr_data_dirpath, Assimilation_filename)
    return pd.read_csv(Assimilation_filepath, index_col='t')


def test_1_lamina():
    # force parameters for backward compatibility
    organ.Organ.delta_t = 1
    
    curr_data_dirpath = os.path.join(DATA_DIRPATH, '1_lamina')
    
    # create the lamina
    Assimilation_df = read_assimilation_file(curr_data_dirpath, 'Assimilation.csv')
    lamina_1 = organ.Lamina(Area=0.00346, Mstruct=0.14, Rdark=0, Assimilation=Assimilation_df.An, 
                             STORAGE_0=0, SUCROSE_0=0, TRIOSESP_0=0, name='')

    # run the model
    python_output_df = cn_model.run(start_time=0, stop_time=960, number_of_output_steps=241, organs=[lamina_1])
    
    # read ModelMaker output
    ModelMaker_output_filepath = os.path.join(curr_data_dirpath, MODELMAKER_OUTPUT_FILENAME)
    ModelMaker_output_df = pd.read_csv(ModelMaker_output_filepath)
    
    # keep only the columns to test
    python_output_df = python_output_df[ModelMaker_output_df.columns]
    
    # save the results (optional)
    python_output_filepath = os.path.join(curr_data_dirpath, PYTHON_OUTPUT_FILENAME)
    python_output_df.to_csv(python_output_filepath, na_rep='NA', index=False)
    
    # compare to the output of ModelMaker
    np.testing.assert_allclose(python_output_df.values, ModelMaker_output_df.values, relative_tolerance, absolute_tolerance)
    
    
def test_3_laminae_phloem_ear():
    # force parameters for backward compatibility
    organ.Organ.delta_t = 1
    organ.Lamina.laminae_inflexion_points.clear()
    
    curr_data_dirpath = os.path.join(DATA_DIRPATH, '3_laminae_phloem_ear')
    
    # create the laminae
    number_of_laminae = 3
    Areas = [0.00346, 0.0034, 0.00228]
    Mstructs = [0.14, 0.09, 0.05]
    Rdarks = [0, 0, 0]
    laminae = []
    for i in xrange(number_of_laminae):
        name = 'lamina%d' % (i + 1)
        curr_Assimilation_df = read_assimilation_file(curr_data_dirpath, 'Assimilation_%s.csv' % name)
        laminae.append(organ.Lamina(Area=Areas[i], Mstruct=Mstructs[i], Rdark=Rdarks[i],
                                    Assimilation=curr_Assimilation_df.An, 
                                    STORAGE_0=0, SUCROSE_0=0, TRIOSESP_0=0,
                                    name=name))
        
    # create the ear
    ear = organ.Ear(value_0=0, name='ear')
    
    # create the phloem
    phloem = organ.Phloem(SUCROSE_0=0, Respiration_0=0, name='phloem')
    
    # run the model
    python_output_df = cn_model.run(start_time=0, stop_time=960, number_of_output_steps=241, organs=laminae+[ear], phloem=phloem)
    
    # read ModelMaker output
    ModelMaker_output_filepath = os.path.join(curr_data_dirpath, MODELMAKER_OUTPUT_FILENAME)
    ModelMaker_output_df = pd.read_csv(ModelMaker_output_filepath)
    
    # keep only the columns to test
    python_output_df = python_output_df[ModelMaker_output_df.columns]
    
    # save the results (optional)
    python_output_filepath = os.path.join(curr_data_dirpath, PYTHON_OUTPUT_FILENAME)
    python_output_df.to_csv(python_output_filepath, na_rep='NA', index=False)

    # compare to the output of ModelMaker
    np.testing.assert_allclose(python_output_df.values, ModelMaker_output_df.values, relative_tolerance, absolute_tolerance)
    
    
def test_Distribution_CN_postflo_V2():
    curr_data_dirpath = os.path.join(DATA_DIRPATH, 'Distribution_CN_postflo_V2')
    
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
    python_output_df = cn_model.run(start_time=0, stop_time=960, number_of_output_steps=241,
                                    organs=[chaff]+internodes+laminae+peduncles+sheaths+[grains]+[roots], 
                                    phloem=phloem)
    
    # read ModelMaker output
    ModelMaker_output_filepath = os.path.join(curr_data_dirpath, MODELMAKER_OUTPUT_FILENAME)
    ModelMaker_output_df = pd.read_csv(ModelMaker_output_filepath)
    
    # keep only the columns to test
    python_output_df = python_output_df[ModelMaker_output_df.columns]
    
    # save the results (optional)
    python_output_filepath = os.path.join(curr_data_dirpath, PYTHON_OUTPUT_FILENAME)
    python_output_df.to_csv(python_output_filepath, na_rep='NA', index=False)
    
    # compare to the output of ModelMaker
    np.testing.assert_allclose(python_output_df.values, ModelMaker_output_df.values, relative_tolerance, absolute_tolerance)
    

if __name__ == '__main__':
    test_1_lamina()
    test_3_laminae_phloem_ear()
    test_Distribution_CN_postflo_V2()

