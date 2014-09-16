'''
Test cn_model.py. 

CSV files must contain only ASCII characters and ',' as separator.

Be sure to add the 'src' directory to your PYTHONPATH before running this script.
'''

import os

import numpy as np
import pandas as pd

import lamina
import cn_model

data_dirpath = 'data'
ModelMaker_output_filename = 'ModelMaker_output.csv'
python_output_filename = 'python_output.csv'

relative_tolerance = 10e-3
absolute_tolerance = 10e-3

def test_1_lamina():
    curr_data_dirpath = os.path.join(data_dirpath, '1_lamina')
    Assimilation_filepath = os.path.join(curr_data_dirpath, 'Assimilation_lamina1.csv')
    ModelMaker_output_filepath = os.path.join(curr_data_dirpath, ModelMaker_output_filename)
    python_output_filepath = os.path.join(curr_data_dirpath, python_output_filename)
    
    Assimilation_df = pd.read_csv(Assimilation_filepath, index_col='t')
    
    lamina_1 = lamina.Lamina(lamina_area=0.00346, Mstruct_lamina=0.14, Assimilation=Assimilation_df.An, 
                           STORAGE_0=0, SUCROSE_lamina_0=0, TRIOSESP_0=0)

    # run the model
    python_output_df = cn_model.run(start_time=0, stop_time=960, number_of_output_steps=241, lamina_1=lamina_1)
    python_output_df.to_csv(python_output_filepath, na_rep='NA', index=False)
    
    # compare to the output of ModelMaker
    ModelMaker_output_df = pd.read_csv(ModelMaker_output_filepath)
    np.testing.assert_allclose(python_output_df.values, ModelMaker_output_df.values, relative_tolerance, absolute_tolerance)
     
     
if __name__ == '__main__':
    test_1_lamina()

