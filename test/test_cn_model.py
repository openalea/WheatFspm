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

def test_1_lamina():
    curr_data_dirpath = os.path.join(DATA_DIRPATH, '1_lamina')
    
    # create the lamina
    Assimilation_filepath = os.path.join(curr_data_dirpath, 'Assimilation_lamina1.csv')
    Assimilation_df = pd.read_csv(Assimilation_filepath, index_col='t')
    lamina_1 = organ.Lamina(lamina_area=0.00346, Mstruct_lamina=0.14, Assimilation=Assimilation_df.An, 
                             STORAGE_0=0, SUCROSE_lamina_0=0, TRIOSESP_0=0)

    # run the model
    python_output_df = cn_model.run(start_time=0, stop_time=960, number_of_output_steps=241, organs=[lamina_1])
    
    # keep only the columns to test
    python_output_df = python_output_df[['t', 'Photosynthesis', 'D_storage', 'S_storage', 'S_sucrose', 'STORAGE', 'SUCROSE_lamina', 'TRIOSESP']]
    
    # save the results
    python_output_filepath = os.path.join(curr_data_dirpath, PYTHON_OUTPUT_FILENAME)
    python_output_df.to_csv(python_output_filepath, na_rep='NA', index=False)
    
    # compare to the output of ModelMaker
    ModelMaker_output_filepath = os.path.join(curr_data_dirpath, MODELMAKER_OUTPUT_FILENAME)
    ModelMaker_output_df = pd.read_csv(ModelMaker_output_filepath)
    np.testing.assert_allclose(python_output_df.values, ModelMaker_output_df.values, relative_tolerance, absolute_tolerance)
    
    
def test_3_laminae_phloem_ear():
    curr_data_dirpath = os.path.join(DATA_DIRPATH, '3_laminae_phloem_ear')
    
    # create the laminae
    number_of_leaves = 3
    lamina_areas = [0.00346, 0.0034, 0.00228]
    Mstruct_laminae = [0.14, 0.09, 0.05]
    laminae = []
    for i in range(number_of_leaves):
        lamina_name = str(i + 1)
        curr_Assimilation_filepath = os.path.join(curr_data_dirpath, 'Assimilation_lamina%s.csv' % lamina_name)
        curr_Assimilation_df = pd.read_csv(curr_Assimilation_filepath, index_col='t')
        laminae.append(organ.Lamina(lamina_area=lamina_areas[i], Mstruct_lamina=Mstruct_laminae[i], 
                                    Assimilation=curr_Assimilation_df.An, 
                                    STORAGE_0=0, SUCROSE_lamina_0=0, TRIOSESP_0=0,
                                    name=lamina_name))
        
    # create the ear
    ear = organ.Ear(Ear_value_0=0)
    
    # create the phloem
    phloem = organ.Phloem(SUCROSE_phloem_0=0)
    
    # run the model
    python_output_df = cn_model.run(start_time=0, stop_time=960, number_of_output_steps=241, organs=laminae+[ear], phloem=phloem)
    
    # save the results
    python_output_filepath = os.path.join(curr_data_dirpath, PYTHON_OUTPUT_FILENAME)
    python_output_df.to_csv(python_output_filepath, na_rep='NA', index=False)
    
    # compare to the output of ModelMaker
    ModelMaker_output_filepath = os.path.join(curr_data_dirpath, MODELMAKER_OUTPUT_FILENAME)
    ModelMaker_output_df = pd.read_csv(ModelMaker_output_filepath)
    np.testing.assert_allclose(python_output_df.values, ModelMaker_output_df.values, relative_tolerance, absolute_tolerance)
    
    
if __name__ == '__main__':
    test_1_lamina()
    test_3_laminae_phloem_ear()

