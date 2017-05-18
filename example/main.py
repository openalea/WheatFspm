# -*- coding: latin-1 -*-
"""
    main
    ~~~~

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

from elongwheat import simulation as elongwheat_simulation, model as elongwheat_model, converter as elongwheat_converter

INPUTS_DIRPATH = 'inputs'

# elongwheat inputs at t0
HIDDENZONE_INPUTS_FILEPATH = os.path.join(INPUTS_DIRPATH, 'hiddenzones_inputs.csv')
ORGAN_INPUTS_FILEPATH = os.path.join(INPUTS_DIRPATH, 'organs_inputs.csv')
SAM_INPUTS_FILEPATH = os.path.join(INPUTS_DIRPATH, 'SAM_inputs.csv')

# elongwheat outputs
OUTPUTS_DIRPATH = 'outputs'
HIDDENZONE_OUTPUTS_FILENAME = 'hiddenzone_outputs.csv'
ORGAN_OUTPUTS_FILENAME = 'organ_outputs.csv'
SAM_OUTPUTS_FILENAME = 'SAM_outputs.csv'

# define the time step in hours for each elongwheat
elongwheat_ts = 1

# read elongwheat inputs at t0
elongwheat_hiddenzones_inputs_t0 = pd.read_csv(HIDDENZONE_INPUTS_FILEPATH)
elongwheat_organ_inputs_t0 = pd.read_csv(ORGAN_INPUTS_FILEPATH)
elongwheat_SAM_inputs_t0 = pd.read_csv(SAM_INPUTS_FILEPATH)

OUTPUTS_PRECISION = 8

if __name__ == '__main__':

    # Create population
    simulation_ = elongwheat_simulation.Simulation(delta_t=3600)
    # read inputs from Pandas dataframe
    hiddenzone_inputs_df = elongwheat_hiddenzones_inputs_t0
    organ_inputs_df = elongwheat_organ_inputs_t0
    SAM_inputs_df = elongwheat_SAM_inputs_t0
    # convert the dataframe to simulation inputs format
    inputs = elongwheat_converter.from_dataframes(hiddenzone_inputs_df, organ_inputs_df, SAM_inputs_df)
    # initialize the simulation with the inputs
    simulation_.initialize(inputs)
    # run the simulation
    simulation_.run()
    # convert the outputs to Pandas dataframe
    hiddenzone_outputs_df, organ_outputs_df, SAM_outputs_df = elongwheat_converter.to_dataframes(simulation_.outputs)
    # write the dataframe to CSV
    hiddenzone_outputs_df.to_csv(os.path.join(OUTPUTS_DIRPATH, HIDDENZONE_OUTPUTS_FILENAME), index=False, na_rep='NA', float_format='%.{}f'.format(OUTPUTS_PRECISION))
    organ_outputs_df.to_csv(os.path.join(OUTPUTS_DIRPATH, ORGAN_OUTPUTS_FILENAME), index=False, na_rep='NA', float_format='%.{}f'.format(OUTPUTS_PRECISION))
    SAM_outputs_df.to_csv(os.path.join(OUTPUTS_DIRPATH, SAM_OUTPUTS_FILENAME), index=False, na_rep='NA', float_format='%.{}f'.format(OUTPUTS_PRECISION))