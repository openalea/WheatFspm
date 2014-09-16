'''
Test cn_model.py. 

Be sure to add the 'src' directory to your PYTHONPATH before running this script.
'''

import os

import numpy as np
import pandas as pd

import cn_model

data_dirpath = 'data'
input_filepath = os.path.join(data_dirpath, 'input.csv') # ASCII characters only and ',' as separator
modelmaker_output_filepath = os.path.join(data_dirpath, 'modelmaker_output.csv')
python_output_filepath = os.path.join(data_dirpath, 'python_output.csv')

relative_tolerance = 10e-3
absolute_tolerance = 10e-3

def test_run():
    input_df = pd.read_csv(input_filepath, index_col='t')
    
    # initial conditions
    STORAGE_0 = 0
    SUCROSE_lamina_0 = 0
    TRIOSESP_0 = 0
    y0 = [STORAGE_0, SUCROSE_lamina_0, TRIOSESP_0]  
    
    # run the model
    python_output_df = cn_model.run(start_time=0, stop_time=960, number_of_output_steps=241,
                                    initial_conditions=y0, An=input_df.An)
    python_output_df.to_csv(python_output_filepath, na_rep='NA', index=False)
    
    # compare to the output of ModelMaker
    modelmaker_output_df = pd.read_csv(modelmaker_output_filepath)
    np.testing.assert_allclose(python_output_df.values, modelmaker_output_df.values, relative_tolerance, absolute_tolerance)
    
if __name__ == '__main__':
    test_run()

