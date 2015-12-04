# -*- coding: latin-1 -*-

'''
    main
    ~~~~

    An example to show how to:
        
        * initialize and run the model Farquhar-Wheat from an Adel-Wheat MTG,
        * format the outputs of Farquhar-Wheat to both Adel-Wheat MTG and Pandas dataframes.

    You must first install :mod:`alinea.adel` and :mod:`farquharwheat` and its dependencies 
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
from alinea.adel import astk_interface

from farquharwheat import simulation, converter

INPUTS_DIRPATH = 'inputs'

# adelwheat inputs
ADELWHEAT_INPUTS_DIRPATH = os.path.join(INPUTS_DIRPATH, 'adelwheat') # the directory adelwheat must contain files 'adel0000.pckl' and 'scene0000.bgeom'

# farquharwheat inputs
FARQUHARWHEAT_INPUTS_DIRPATH = os.path.join(INPUTS_DIRPATH, 'farquharwheat')
FARQUHARWHEAT_ORGANS_INPUTS_FILEPATH = os.path.join(FARQUHARWHEAT_INPUTS_DIRPATH, 'organs_inputs.csv')
FARQUHARWHEAT_ELEMENTS_INPUTS_FILEPATH = os.path.join(FARQUHARWHEAT_INPUTS_DIRPATH, 'elements_inputs.csv')


OUTPUTS_DIRPATH = 'outputs'

# farquharwheat outputs
FARQUHARWHEAT_OUTPUTS_DIRPATH = os.path.join(OUTPUTS_DIRPATH, 'farquharwheat')
FARQUHARWHEAT_ELEMENTS_OUTPUTS_FILEPATH = os.path.join(FARQUHARWHEAT_OUTPUTS_DIRPATH, 'elements_outputs.csv')

OUTPUTS_PRECISION = 6


if __name__ == '__main__':
    
    # create a simulation and a converter
    simulation_ = simulation.Simulation()
    # read adelwheat inputs
    adel_wheat = astk_interface.AdelWheat(seed=1234)
    g = adel_wheat.load(dir=ADELWHEAT_INPUTS_DIRPATH)[0]
    # convert the MTG to Farquhar-Wheat inputs and initialize the simulation
    farquharwheat_organs_inputs_df = pd.read_csv(FARQUHARWHEAT_ORGANS_INPUTS_FILEPATH)
    farquharwheat_elements_inputs_df = pd.read_csv(FARQUHARWHEAT_ELEMENTS_INPUTS_FILEPATH)
    simulation_.initialize(converter.from_MTG(g, farquharwheat_organs_inputs_df, farquharwheat_elements_inputs_df))
    # run the simulation
    simulation_.run(Ta=18.8, ambient_CO2=360, RH=0.530000, Ur=2.200000, PARi=3838000)
    # update the MTG from Farquhar-Wheat outputs
    converter.update_MTG(simulation_.inputs, simulation_.outputs, g)
    # format Farquhar-Wheat outputs to Pandas dataframes
    _, farquharwheat_elements_outputs_df = converter.to_dataframe(simulation_.outputs)
    # write the dataframe to CSV
    farquharwheat_elements_outputs_df.to_csv(FARQUHARWHEAT_ELEMENTS_OUTPUTS_FILEPATH, index=False, na_rep='NA', float_format='%.{}f'.format(OUTPUTS_PRECISION)) 


