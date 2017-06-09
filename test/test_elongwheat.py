# -*- coding: latin-1 -*-
"""
    test_elongwheat
    ~~~~~~~~~~~~~~~~

    An example to show how to:

        * initialize and run the model Elong-Wheat,
        * format the outputs of Elong-Wheat.

    You must first install :mod:`elongwheat` and its dependencies
    before running this script with the command `python`.

    :copyright: Copyright 2014-2016 INRA-ECOSYS, see AUTHORS.
    :license: TODO, see LICENSE for details.

    .. seealso:: Barillot et al. 2016.

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

from elongwheat import simulation, converter

INPUTS_DIRPATH = 'inputs'
HIDDENZONE_INPUTS_FILENAME = 'hiddenzones_inputs.csv'
ORGAN_INPUTS_FILENAME = 'organs_inputs.csv'
SAM_INPUTS_FILENAME = 'SAM_inputs.csv'

OUTPUTS_DIRPATH = 'outputs'
DESIRED_HIDDENZONE_OUTPUTS_FILENAME = 'desired_hiddenzone_outputs.csv'
DESIRED_ORGAN_OUTPUTS_FILENAME = 'desired_organ_outputs.csv'
DESIRED_SAM_OUTPUTS_FILENAME = 'desired_SAM_outputs.csv'
ACTUAL_HIDDENZONE_OUTPUTS_FILENAME = 'actual_hiddenzone_outputs.csv'
ACTUAL_ORGAN_OUTPUTS_FILENAME = 'actual_organ_outputs.csv'
ACTUAL_SAM_OUTPUTS_FILENAME = 'actual_SAM_outputs.csv'

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
    for column in ('axis', 'organ', 'leaf_is_growing', 'internode_is_growing','leaf_is_emerged', 'is_growing','status'):
        if column in desired_data_df.columns:
            assert desired_data_df[column].equals(actual_data_df[column])
            del desired_data_df[column]
            del actual_data_df[column]

    # compare to the desired data
    np.testing.assert_allclose(actual_data_df.values, desired_data_df.values, RELATIVE_TOLERANCE, ABSOLUTE_TOLERANCE)


def test_run():

    # create a simulation
    simulation_ = simulation.Simulation(delta_t=3600)
    # read inputs from Pandas dataframe
    hiddenzone_inputs_df = pd.read_csv(os.path.join(INPUTS_DIRPATH, HIDDENZONE_INPUTS_FILENAME))
    organ_inputs_df = pd.read_csv(os.path.join(INPUTS_DIRPATH, ORGAN_INPUTS_FILENAME))
    SAM_inputs_df = pd.read_csv(os.path.join(INPUTS_DIRPATH, SAM_INPUTS_FILENAME))
    # convert the dataframe to simulation inputs format
    inputs = converter.from_dataframes(hiddenzone_inputs_df, organ_inputs_df, SAM_inputs_df)
    # initialize the simulation with the inputs
    simulation_.initialize(inputs)
##    # convert the inputs to Pandas dataframe
##    hiddenzone_inputs_reconverted_df, organ_inputs_reconverted_df = converter.to_dataframes(simulation_.inputs)
##    # compare inputs
##    compare_actual_to_desired('inputs', hiddenzone_inputs_reconverted_df, HIDDENZONE_INPUTS_FILENAME)
##    compare_actual_to_desired('inputs', organ_inputs_reconverted_df, ORGAN_INPUTS_FILENAME)
    # run the simulation
    simulation_.run(Ta=25,Ts=20)
    # convert the outputs to Pandas dataframe
    hiddenzone_outputs_df, organ_outputs_df , SAM_outputs_df = converter.to_dataframes(simulation_.outputs)
    # compare outputs
    compare_actual_to_desired('outputs', hiddenzone_outputs_df, DESIRED_HIDDENZONE_OUTPUTS_FILENAME, ACTUAL_HIDDENZONE_OUTPUTS_FILENAME)
    compare_actual_to_desired('outputs', organ_outputs_df, DESIRED_ORGAN_OUTPUTS_FILENAME, ACTUAL_ORGAN_OUTPUTS_FILENAME)
    compare_actual_to_desired('outputs', SAM_outputs_df, DESIRED_SAM_OUTPUTS_FILENAME, ACTUAL_SAM_OUTPUTS_FILENAME)


if __name__ == '__main__':
    test_run()
