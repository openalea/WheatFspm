"""
Use :func:`cnwheat.post_processing.plot_linear_regression` to compare:
    * the output computed using ModelMaker, 
    * and the output computed using :func:`cnwheat.cnwheat.run`.

The plot is saved in 'compare_outputs.png'.

The package cnwheat must be correctly installed before running this example. 

CSV files must contain only ASCII characters and ',' as separator.
"""

import os

import pandas as pd

from cnwheat import post_processing

DATA_DIRPATH = 'data'
MODELMAKER_OUTPUT_FILENAME = 'modelmaker_output.csv'
CNWHEAT_OUTPUT_FILENAME = 'cnwheat_output.csv'
PLOT_FILENAME = 'compare_outputs.png'
OUTPUT_NAME = 'SUCROSE_phloem'

modelmaker_output_filepath = os.path.join(DATA_DIRPATH, MODELMAKER_OUTPUT_FILENAME)
cnwheat_output_filepath = os.path.join(DATA_DIRPATH, CNWHEAT_OUTPUT_FILENAME)

ModelMaker_array = pd.read_csv(modelmaker_output_filepath)[OUTPUT_NAME].values
cnwheat_array = pd.read_csv(cnwheat_output_filepath)[OUTPUT_NAME].values

post_processing.plot_linear_regression(ModelMaker_array, 
                                       cnwheat_array, 
                                       x_label='modelmaker_{}'.format(OUTPUT_NAME), 
                                       y_label='cnwheat_{}'.format(OUTPUT_NAME), 
                                       plot_filepath=PLOT_FILENAME)

