# -*- coding: latin-1 -*-

'''
    main
    ~~~~

    An example to show how to:
        
        * initialize and run the model Senesc-Wheat,
        * format the outputs of Senesc-Wheat to Pandas dataframe.

    You must first install :mod:`senescwheat` and its dependencies 
    before running this script with the command `python`.

    :copyright: Copyright 2014-2016 INRA-ECOSYS, see AUTHORS.
    :license: TODO, see LICENSE for details.

    .. seealso:: Barillot et al. 2016.
'''

'''
    Information about this versioned file:
        $LastChangedBy$
        $LastChangedDate$
        $LastChangedRevision$
        $URL$
        $Id$
'''

import os

import pandas as pd

from senescwheat import model, simulation, converter

# inputs paths
INPUTS_DIRPATH = 'inputs'

ROOTS_INPUTS_FILENAME = 'roots_inputs.csv'
ELEMENTS_INPUTS_FILENAME = 'elements_inputs.csv'

# outputs paths
OUTPUTS_DIRPATH = 'outputs'

ROOTS_OUTPUTS_FILENAME = 'roots_outputs.csv'
ELEMENTS_OUTPUTS_FILENAME = 'elements_outputs.csv'

roots_outputs_filepath = os.path.join(OUTPUTS_DIRPATH, ROOTS_OUTPUTS_FILENAME)
elements_outputs_filepath = os.path.join(OUTPUTS_DIRPATH, ELEMENTS_OUTPUTS_FILENAME)


if __name__ == '__main__':
    
    # create a simulation
    simulation_ = simulation.Simulation(delta_t=3600)
    # read inputs from Pandas dataframes
    roots_inputs_df = pd.read_csv(os.path.join(INPUTS_DIRPATH, ROOTS_INPUTS_FILENAME))
    elements_inputs_df = pd.read_csv(os.path.join(INPUTS_DIRPATH, ELEMENTS_INPUTS_FILENAME))
    # convert the dataframes to simulation inputs format
    inputs = converter.from_dataframes(roots_inputs_df, elements_inputs_df)
    # initialize the simulation with the inputs
    simulation_.initialize(inputs)
    # run the simulation
    simulation_.run()
    # convert the outputs to Pandas dataframes
    roots_outputs_df, elements_outputs_df = converter.to_dataframes(simulation_.outputs)
    # write the dataframes to CSV
    roots_outputs_df.to_csv(os.path.join(OUTPUTS_DIRPATH, ROOTS_OUTPUTS_FILENAME), index=False, na_rep='NA')
    elements_outputs_df.to_csv(os.path.join(OUTPUTS_DIRPATH, ELEMENTS_OUTPUTS_FILENAME), index=False, na_rep='NA') 
    
