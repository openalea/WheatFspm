# -*- coding: latin-1 -*-
"""
    test_farquhar_wheat
    ~~~~~~~~~~~~~~~~~~~

    Test the Farquhar-Wheat model (standalone).

    You must first install :mod:`farquharwheat` and its dependencies
    before running this script with the command `python`.

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

from farquharwheat import model, simulation 

INPUTS_DIRPATH = 'inputs'

INPUTS_FILENAME = 'inputs.csv'

OUTPUTS_DIRPATH = 'outputs'

DESIRED_OUTPUTS_FILENAME = 'desired_outputs.csv'

ACTUAL_OUTPUTS_FILENAME = 'actual_outputs.csv'

PRECISION = 6
RELATIVE_TOLERANCE = 10**-PRECISION
ABSOLUTE_TOLERANCE = RELATIVE_TOLERANCE


def compare_actual_to_desired(outputs_dirpath, actual_output_df, desired_output_filename, actual_output_filename, save_actual_output=False):
    # read desired output
    desired_output_filepath = os.path.join(outputs_dirpath, desired_output_filename)
    desired_output_df = pd.read_csv(desired_output_filepath)

    if save_actual_output:
        actual_output_filepath = os.path.join(outputs_dirpath, actual_output_filename)
        actual_output_df.to_csv(actual_output_filepath, na_rep='NA', index=False, float_format='%.{}f'.format(PRECISION))

    # keep only numerical data
    for column in ('axis', 'organ', 'element'):
        if column in desired_output_df.columns:
            del desired_output_df[column]
            del actual_output_df[column]

    # compare to the desired output
    np.testing.assert_allclose(actual_output_df.values, desired_output_df.values, RELATIVE_TOLERANCE, ABSOLUTE_TOLERANCE)


def test_run():
    
    # create a simulation
    simulation_ = simulation.Simulation()
    # read inputs from Pandas dataframe
    inputs = pd.read_csv(os.path.join(INPUTS_DIRPATH, INPUTS_FILENAME))
    # initialize the simulation with the inputs
    simulation_.initialize(inputs)
    # run the simulation
    simulation_.run(Ta=18.8, ambient_CO2=360, RH=0.68, Ur=3.171, PARi=2262400)
    # format the outputs to Pandas dataframe
    outputs_df = simulation_.format_outputs()
    
    compare_actual_to_desired(OUTPUTS_DIRPATH, outputs_df, DESIRED_OUTPUTS_FILENAME, ACTUAL_OUTPUTS_FILENAME, save_actual_output=True)
    

if __name__ == '__main__':
    test_run()
