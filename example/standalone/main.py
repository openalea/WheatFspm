# -*- coding: latin-1 -*-

'''
    main
    ~~~~

    An example to show how to:
        
        * initialize and run the model Farquhar-Wheat,
        * format the outputs of Farquhar-Wheat.

    You must first install :mod:`farquharwheat` and its dependencies 
    before running this script with the command `python`.

    :copyright: Copyright 2014-2015 INRA-ECOSYS, see AUTHORS.
    :license: TODO, see LICENSE for details.

    .. seealso:: Barillot et al. 2015.
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

from farquharwheat import model, simulation, converter 

ORGANS_INPUTS_FILENAME = 'organs_inputs.csv'
ELEMENTS_INPUTS_FILENAME = 'elements_inputs.csv'

ELEMENTS_OUTPUTS_FILENAME = 'elements_outputs.csv'

OUTPUTS_PRECISION = 6

if __name__ == '__main__':
    
    # create a simulation and a converter
    simulation_ = simulation.Simulation()
    # read inputs from Pandas dataframe
    organs_inputs_df = pd.read_csv(ORGANS_INPUTS_FILENAME)
    elements_inputs_df = pd.read_csv(ELEMENTS_INPUTS_FILENAME)
    # convert the dataframe to simulation inputs format
    inputs = converter.from_dataframes(organs_inputs_df, elements_inputs_df)
    # initialize the simulation with the inputs
    simulation_.initialize(inputs)
    # run the simulation
    simulation_.run(Ta=18.8, ambient_CO2=360, RH=0.530000, Ur=2.200000, PARi=3838000)
    # convert the outputs to Pandas dataframes
    _, elements_outputs_df = converter.to_dataframes(simulation_.outputs)
    # write the dataframe to CSV
    elements_outputs_df.to_csv(ELEMENTS_OUTPUTS_FILENAME, index=False, na_rep='NA', float_format='%.{}f'.format(OUTPUTS_PRECISION)) 
    
