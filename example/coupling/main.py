# -*- coding: latin-1 -*-

'''
    main
    ~~~~

    An example to show how to:
        
        * initialize and run the model Senesc-Wheat from an MTG,
        * format the outputs of Senesc-Wheat to both MTG and Pandas dataframes.
        
    You must first install :mod:`senescwheat` and its dependencies 
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

from senescwheat import simulation, converter

INPUTS_DIRPATH = 'inputs'
ROOTS_INPUTS_FILENAME = 'roots_inputs.csv'
ELEMENTS_INPUTS_FILENAME = 'elements_inputs.csv'

OUTPUTS_DIRPATH = 'outputs'
ROOTS_OUTPUTS_FILENAME = 'roots_outputs.csv'
ELEMENTS_OUTPUTS_FILENAME = 'elements_outputs.csv'


def setup_MTG():
    """Construct and fill a valid MTG.
    TODO: add roots to the MTG.
    TODO: load the fulfilled MTG from files. 
    """
    import random
    import numpy as np
    import pandas as pd
    from alinea.adel.astk_interface import initialise_stand
    
    random.seed(1234)
    np.random.seed(1234)
    
    elements_inputs_df = pd.read_csv(os.path.join(INPUTS_DIRPATH, ELEMENTS_INPUTS_FILENAME))
    
    g, wheat, domain_area, domain, convUnit = initialise_stand(1500)
    
    # add the properties which do not exist yet
    property_names = g.property_names()
    for input_name in converter.SENESCWHEAT_ELEMENTS_INPUTS:
        if input_name not in property_names:
            g.add_property(input_name)
    
    # traverse the MTG recursively from top
    for plant_vid in g.components_iter(g.root):
        plant_index = int(g.index(plant_vid))
        current_plant_group = elements_inputs_df[elements_inputs_df['plant'] == plant_index]
        if len(current_plant_group) == 0:
            # TODO: remove plant_vid and all its components
            continue
        for axis_vid in g.components_iter(plant_vid):
            axis_id = g.label(axis_vid)
            current_axis_group = current_plant_group[current_plant_group['axis'] == axis_id]
            if len(current_axis_group) == 0:
                # TODO: remove axis_vid and all its components
                continue
            for metamer_vid in g.components(axis_vid): 
                metamer_index = int(g.index(metamer_vid))
                current_metamer_group = current_axis_group[current_axis_group['metamer'] == metamer_index]
                if len(current_metamer_group) == 0:
                    # TODO: remove metamer_vid and all its components
                    continue
                for organ_vid in g.components_iter(metamer_vid):
                    organ_label = g.label(organ_vid)
                    current_organ_group = current_metamer_group[current_metamer_group['organ'] == organ_label]
                    if len(current_organ_group) == 0:
                        # TODO: remove organ_vid and all its components
                        continue
                    for element_vid in g.components_iter(organ_vid):
                        element_label = g.label(element_vid)
                        current_element_group = current_organ_group[current_organ_group['element'] == element_label]
                        if len(current_element_group) == 0:
                            # TODO: remove element_vid and all its components
                            continue
                        current_element_series = current_element_group.loc[current_element_group.first_valid_index()]
                        for property_ in converter.SENESCWHEAT_ELEMENTS_INPUTS:
                            g.property(property_)[element_vid] = current_element_series[property_]
    return g

if __name__ == '__main__':
    
    # create a simulation
    simulation_ = simulation.Simulation()
    # setup a MTG 
    g = setup_MTG()
    # convert the MTG to Senesc-Wheat inputs and initialize the simulation
    roots_inputs = pd.read_csv(os.path.join(INPUTS_DIRPATH, ROOTS_INPUTS_FILENAME)) #TODO: remove this code as soon as :mod:`openalea.mtg` permits to add/remove components.
    elements_inputs_df = pd.read_csv(os.path.join(INPUTS_DIRPATH, ELEMENTS_INPUTS_FILENAME)).astype(str)
    available_components = set(elements_inputs_df.groupby(['plant']).groups.keys() + \
                               elements_inputs_df.groupby(['plant', 'axis']).groups.keys() + \
                               elements_inputs_df.groupby(['plant', 'axis', 'metamer']).groups.keys() + \
                               elements_inputs_df.groupby(['plant', 'axis', 'metamer', 'organ']).groups.keys() + \
                               elements_inputs_df.groupby(converter.ELEMENTS_TOPOLOGY_COLUMNS).groups.keys())
    simulation_.initialize(converter.from_MTG(g, roots_inputs, available_components))
    # run the simulation
    simulation_.run()
    # format Senesc-Wheat outputs to Pandas dataframe
    roots_outputs_df, elements_outputs_df = converter.to_dataframes(simulation_.outputs)
    # update the MTG from Senesc-Wheat outputs
    converter.update_MTG(simulation_.outputs, g, roots_outputs_df, available_components)
    # write the dataframe to CSV
    roots_outputs_df.to_csv(os.path.join(OUTPUTS_DIRPATH, ROOTS_OUTPUTS_FILENAME), index=False, na_rep='NA')
    elements_outputs_df.to_csv(os.path.join(OUTPUTS_DIRPATH, ELEMENTS_OUTPUTS_FILENAME), index=False, na_rep='NA')  


