# -*- coding: latin-1 -*-

'''
    main
    ~~~~

    An example to show how to:
        
        * initialize and run the model Senesc-Wheat from an Adel-Wheat MTG,
        * format the outputs of Senesc-Wheat to both MTG and Pandas dataframes.
        
    You must first install :mod:`alinea.adel` and :mod:`senescwheat` and its dependencies 
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

from senescwheat import simulation, converter

INPUTS_DIRPATH = 'inputs'

# adelwheat inputs
ADELWHEAT_INPUTS_DIRPATH = os.path.join(INPUTS_DIRPATH, 'adelwheat') # the directory adelwheat must contain files 'adel0000.pckl' and 'scene0000.bgeom'

# senescwheat inputs
SENESCWHEAT_INPUTS_DIRPATH = os.path.join(INPUTS_DIRPATH, 'senescwheat')
SENESCWHEAT_ROOTS_INPUTS_FILEPATH = os.path.join(SENESCWHEAT_INPUTS_DIRPATH, 'roots_inputs.csv')
SENESCWHEAT_ELEMENTS_INPUTS_FILEPATH = os.path.join(SENESCWHEAT_INPUTS_DIRPATH, 'elements_inputs.csv')


OUTPUTS_DIRPATH = 'outputs'

# senescwheat outputs
SENESCWHEAT_OUTPUTS_DIRPATH = os.path.join(OUTPUTS_DIRPATH, 'senescwheat')
SENESCWHEAT_ROOTS_OUTPUTS_FILEPATH = os.path.join(SENESCWHEAT_OUTPUTS_DIRPATH, 'roots_outputs.csv')
SENESCWHEAT_ELEMENTS_OUTPUTS_FILEPATH = os.path.join(SENESCWHEAT_OUTPUTS_DIRPATH, 'elements_outputs.csv')


if __name__ == '__main__':
    
    # create a simulation
    simulation_ = simulation.Simulation(delta_t=3600)
    # read adelwheat inputs
    adel_wheat = astk_interface.AdelWheat(seed=1234)
    g = adel_wheat.load(dir=ADELWHEAT_INPUTS_DIRPATH)[0]
    # convert the MTG to Senesc-Wheat inputs and initialize the simulation
    senescwheat_roots_inputs_df = pd.read_csv(SENESCWHEAT_ROOTS_INPUTS_FILEPATH)
    senescwheat_elements_inputs_df = pd.read_csv(SENESCWHEAT_ELEMENTS_INPUTS_FILEPATH)
    simulation_.initialize(converter.from_MTG(g, senescwheat_roots_inputs_df, senescwheat_elements_inputs_df))
    # run the simulation
    simulation_.run()
    # update the MTG from Senesc-Wheat outputs
    converter.update_MTG(simulation_.inputs, simulation_.outputs, g)
    # format Senesc-Wheat outputs to Pandas dataframes
    roots_outputs_df, elements_outputs_df = converter.to_dataframes(simulation_.outputs)
    # write the dataframes to CSV
    roots_outputs_df.to_csv(SENESCWHEAT_ROOTS_OUTPUTS_FILEPATH, index=False, na_rep='NA')
    elements_outputs_df.to_csv(SENESCWHEAT_ELEMENTS_OUTPUTS_FILEPATH, index=False, na_rep='NA')


