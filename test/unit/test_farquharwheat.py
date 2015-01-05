# -*- coding: latin-1 -*-
"""
    test_farquhar_wheat
    ~~~~~~~~~~~~~~~~~~~

    Test the package farquharwheat.

    CSV files must contain only ASCII characters and ',' as separator.

    Be sure to add the 'farquharwheat' directory to your PYTHONPATH before running this script.

    :copyright: Copyright 2014 INRA-EGC, see AUTHORS.
    :license: TODO, see LICENSE for details.
"""

"""
    Information about this versioned file:
        $LastChangedBy$
        $LastChangedDate$
        $LastChangedRevision$
        $URL$
        $Id$
"""

import os

import numpy as np
import pandas as pd

from farquharwheat.model import PhotosynthesisModel

DATA_DIRPATH = 'data'
DESIRED_OUTPUT_FILENAME = 'desired_output.csv'
ACTUAL_OUTPUT_FILENAME = 'actual_output.csv'

PRECISION = 2
RELATIVE_TOLERANCE = 10**-PRECISION
ABSOLUTE_TOLERANCE = RELATIVE_TOLERANCE


def read_t_data(curr_data_dirpath, data_filename):
    data_filepath = os.path.join(curr_data_dirpath, data_filename)
    return pd.read_csv(data_filepath, sep=None, index_col='t', engine = 'python')


def compare_actual_to_desired(DATA_DIRPATH, actual_output_df, save_actual_output=False):
    # read desired output
    desired_output_filepath = os.path.join(DATA_DIRPATH, DESIRED_OUTPUT_FILENAME)
    desired_output_df = pd.read_csv(desired_output_filepath)

    # keep only the columns to test
    actual_output_df = actual_output_df[desired_output_df.columns]

    if save_actual_output:
        actual_output_filepath = os.path.join(DATA_DIRPATH, ACTUAL_OUTPUT_FILENAME)
        actual_output_df.to_csv(actual_output_filepath, na_rep='NA', index=False, float_format='%.{}f'.format(PRECISION))

    # compare to the desired output
    np.testing.assert_allclose(actual_output_df.values, desired_output_df.values, RELATIVE_TOLERANCE, ABSOLUTE_TOLERANCE)


def test_calculate_An():
    
    # organ dimensions
    organ_width = 0.018
    organ_height = 0.6
    
    # get meteo_df data
    meteo_df = read_t_data(DATA_DIRPATH, 'meteo.csv')
    
    # get the PAR
    PAR_series = read_t_data(DATA_DIRPATH, 'PAR_Lamina.csv').PAR
    
    time_grid = meteo_df.index
    
    # compute An and Tr for each t in the time grid
    An_Tr_list = []
    for t in time_grid:
        try:
            PAR_t = PAR_series[t]
        except:
            pass
        air_temperature_t = meteo_df['air_temperature'][t]
        humidity_t = meteo_df['humidity'][t]
        ambient_CO2_t = meteo_df['ambient_CO2'][t]
        Wind_top_canopy_t = meteo_df['Wind'][t]
        An_Tr = PhotosynthesisModel.calculate_An(t, organ_width, organ_height, PAR_t, 
                                                 air_temperature_t, ambient_CO2_t, 
                                                 humidity_t, Wind_top_canopy_t)
        An_Tr_list.append(An_Tr)
    
    An_Tr = np.array(An_Tr_list)
    actual_output_df = pd.DataFrame.from_items([('t', time_grid), ('An', An_Tr[:,0]), ('Tr', An_Tr[:,1])])
    
    compare_actual_to_desired(DATA_DIRPATH, actual_output_df, True)


if __name__ == '__main__':
    test_calculate_An()


