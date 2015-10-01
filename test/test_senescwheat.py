# -*- coding: latin-1 -*-
"""
    test_senescwheat
    ~~~~~~~~~~~~~~~

    Test the model Senesc-Wheat.

    You must first install :mod:`senescwheat` (and add it to your PYTHONPATH)
    before running this script with the command `python`.

    :copyright: Copyright 2014-2015 INRA-ECOSYS, see AUTHORS.
    :license: TODO, see LICENSE for details.
    
    .. seealso:: Barillot et al. 2015.
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

from senescwheat import model, simulation, converter

INPUTS_DIRPATH = 'inputs'
ROOTS_INPUTS_FILENAME = 'roots_inputs.csv'
ELEMENTS_INPUTS_FILENAME = 'elements_inputs.csv'

OUTPUTS_DIRPATH = 'outputs'
DESIRED_ROOTS_OUTPUTS_FILENAME = 'desired_roots_outputs.csv'
DESIRED_ELEMENTS_OUTPUTS_FILENAME = 'desired_elements_outputs.csv'
ACTUAL_ROOTS_OUTPUTS_FILENAME = 'actual_roots_outputs.csv'
ACTUAL_ELEMENTS_OUTPUTS_FILENAME = 'actual_elements_outputs.csv'


PRECISION = 6
RELATIVE_TOLERANCE = 10**-PRECISION
ABSOLUTE_TOLERANCE = RELATIVE_TOLERANCE


def compare_actual_to_desired(data_dirpath, actual_data_df, desired_data_filename, actual_data_filename=None):
    # read desired data
    desired_data_filepath = os.path.join(data_dirpath, desired_data_filename)
    desired_data_df = pd.read_csv(desired_data_filepath)

    if actual_data_filename is not None:
        actual_data_filepath = os.path.join(data_dirpath, actual_data_filename)
        actual_data_df.to_csv(actual_data_filepath, na_rep='NA', index=False)

    # keep only numerical data
    for column in ('axis', 'organ', 'element'):
        if column in desired_data_df.columns:
            del desired_data_df[column]
            del actual_data_df[column]

    # compare to the desired data
    np.testing.assert_allclose(actual_data_df.values, desired_data_df.values, RELATIVE_TOLERANCE, ABSOLUTE_TOLERANCE)


def test_run():

    # create a simulation
    simulation_ = simulation.Simulation()
    # read inputs from Pandas dataframe
    roots_inputs_df = pd.read_csv(os.path.join(INPUTS_DIRPATH, ROOTS_INPUTS_FILENAME))
    elements_inputs_df = pd.read_csv(os.path.join(INPUTS_DIRPATH, ELEMENTS_INPUTS_FILENAME))
    # convert the dataframe to simulation inputs format
    inputs = converter.from_dataframe(roots_inputs_df, elements_inputs_df)
    # initialize the simulation with the inputs
    simulation_.initialize(inputs)
    # convert the inputs to Pandas dataframe
    roots_inputs_reconverted_df, elements_inputs_reconverted_df = converter.to_dataframe(simulation_.inputs)
    # compare inputs
    compare_actual_to_desired('inputs', roots_inputs_reconverted_df, ROOTS_INPUTS_FILENAME)
    compare_actual_to_desired('inputs', elements_inputs_reconverted_df, ELEMENTS_INPUTS_FILENAME)
    # run the simulation
    simulation_.run()
    # convert the outputs to Pandas dataframe
    roots_outputs_df, elements_outputs_df = converter.to_dataframe(simulation_.outputs)
    # compare outputs
    compare_actual_to_desired('outputs', roots_outputs_df, DESIRED_ROOTS_OUTPUTS_FILENAME, ACTUAL_ROOTS_OUTPUTS_FILENAME)
    compare_actual_to_desired('outputs', elements_outputs_df, DESIRED_ELEMENTS_OUTPUTS_FILENAME, ACTUAL_ELEMENTS_OUTPUTS_FILENAME)


if __name__ == '__main__':
    test_run()
