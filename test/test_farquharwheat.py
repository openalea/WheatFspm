# -*- coding: latin-1 -*-
"""
    test_farquhar_wheat
    ~~~~~~~~~~~~~~~~~~~

    Test the Farquhar-Wheat model (standalone).

    You must first install :mod:`farquharwheat` and its dependencies
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

from farquharwheat import model, simulation, converter

ORGANS_INPUTS_FILENAME = 'organs_inputs.csv'
ELEMENTS_INPUTS_FILENAME = 'elements_inputs.csv'

DESIRED_ELEMENTS_OUTPUTS_FILENAME = 'desired_elements_outputs.csv'

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
    for column in ('axis', 'organ', 'element', 'label'):
        if column in desired_data_df.columns:
            del desired_data_df[column]
            del actual_data_df[column]

    # compare to the desired data
    np.testing.assert_allclose(actual_data_df.values, desired_data_df.values, RELATIVE_TOLERANCE, ABSOLUTE_TOLERANCE)


def test_run():

    # create a simulation and a converter
    simulation_ = simulation.Simulation()
    # read inputs from Pandas dataframe
    organs_inputs_df = pd.read_csv(ORGANS_INPUTS_FILENAME)
    elements_inputs_df = pd.read_csv(ELEMENTS_INPUTS_FILENAME)
    # convert the dataframe to simulation inputs format
    inputs = converter.from_dataframe(organs_inputs_df, elements_inputs_df)
    # initialize the simulation with the inputs
    simulation_.initialize(inputs)
    # convert the inputs to Pandas dataframe
    reconverted_organs_inputs, reconverted_elements_inputs = converter.to_dataframe(simulation_.inputs)
    # compare inputs
    compare_actual_to_desired('.', reconverted_organs_inputs, ORGANS_INPUTS_FILENAME)
    compare_actual_to_desired('.', reconverted_elements_inputs, ELEMENTS_INPUTS_FILENAME)
    # run the simulation
    simulation_.run(Ta=18.8, ambient_CO2=360, RH=0.530000, Ur=2.200000, PARi=3838000)
    # convert the outputs to Pandas dataframe
    _, elements_outputs_df = converter.to_dataframe(simulation_.outputs)
    # compare outputs
    compare_actual_to_desired('.', elements_outputs_df, DESIRED_ELEMENTS_OUTPUTS_FILENAME, ACTUAL_ELEMENTS_OUTPUTS_FILENAME)


if __name__ == '__main__':
    test_run()
