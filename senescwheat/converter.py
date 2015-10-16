# -*- coding: latin-1 -*-

from __future__ import division # use "//" to do integer division

"""
    senescwheat.converter
    ~~~~~~~~~~~~~~~~~~~~~~~

    The module :mod:`senescwheat.converter` defines functions to:
        
        * convert :class:`dataframes <pandas.DataFrame>` to/from SenescWheat inputs or outputs format.
        * convert a :class:`MTG <openalea.mtg.mtg.MTG>` to/from SenescWheat inputs or outputs format.
        
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

import numpy as np
import pandas as pd

import simulation


class ConverterError(Exception): pass


class NotModeledComponentError(ConverterError):
    '''MTG component not modeled by SenescWheat.
    TODO: make this warning an error'''
    def __init__(self, component_label, component_vertex_id):
        self.message = 'MTG component "{0}" is not modeled by SenescWheat (vertex {1}).'.format(component_label, component_vertex_id)
    def __str__(self):
        return repr(self.message)

class PropertyNotFoundError(ConverterError):
    '''Property not found in a vertex of a MTG.
    TODO: make this warning an error'''
    def __init__(self, vertex_property, vertex_id):
        self.message = 'Property "{0}" not found in vertex {1}.'.format(vertex_property, vertex_id)
    def __str__(self):
        return repr(self.message)
    
class MismatchedTopologiesError(ConverterError):
    '''Topologies mismatched between SenescWheat and MTG.
    TODO: make this warning an error'''
    def __init__(self, senescwheat_id, vertex_id):
        self.message = 'Topologies mismatched between SenescWheat and MTG: no mapping between SenescWheat object {0} and vertex {1}.'.format(senescwheat_id, vertex_id)
    def __str__(self):
        return repr(self.message)
    
    
#: the name of the photosynthetic organs modeled by SenescWheat 
SENESCWHEAT_PHOTOSYNTHETIC_ORGANS_NAMES = set(['internode', 'blade', 'sheath', 'peduncle', 'ear'])

#: the inputs needed by SenescWheat at roots scale
SENESCWHEAT_ROOTS_INPUTS = ['sucrose', 'amino_acids', 'mstruct', 'Nstruct', 'cytokinines']

#: the inputs needed by SenescWheat at element scale
SENESCWHEAT_ELEMENTS_INPUTS = ['green_area', 'proteins', 'mstruct', 'max_proteins', 'Nstruct', 'nitrates', 'amino_acids', 'starch', 'fructan', 'cytokinines', 'sucrose']

#: the inputs needed by SenescWheat
SENESCWHEAT_INPUTS = SENESCWHEAT_ROOTS_INPUTS + SENESCWHEAT_ELEMENTS_INPUTS

#: the outputs computed by SenescWheat at roots scale
SENESCWHEAT_ROOTS_OUTPUTS = ['mstruct_C_growth', 'Nstruct_N_growth', 'mstruct_death', 'mstruct', 'Nstruct', 'cytokinines']

#: the outputs computed by SenescWheat at elements scale
SENESCWHEAT_ELEMENTS_OUTPUTS = ['green_area', 'mstruct', 'Nstruct', 'starch', 'sucrose', 'fructan', 'proteins', 'amino_acids', 'cytokinines', 'surfacic_nitrogen']

#: the outputs computed by SenescWheat
SENESCWHEAT_OUTPUTS = SENESCWHEAT_ROOTS_OUTPUTS + SENESCWHEAT_ELEMENTS_OUTPUTS

#: the inputs and outputs of SenescWheat at roots scale
SENESCWHEAT_ROOTS_INPUTS_OUTPUTS = list(set(SENESCWHEAT_ROOTS_INPUTS + SENESCWHEAT_ROOTS_OUTPUTS))

#: the inputs and outputs of SenescWheat at elements scale
SENESCWHEAT_ELEMENTS_INPUTS_OUTPUTS = list(set(SENESCWHEAT_ELEMENTS_INPUTS + SENESCWHEAT_ELEMENTS_OUTPUTS))

#: the inputs and outputs of FarquharWheat. 
SENESCWHEAT_INPUTS_OUTPUTS = SENESCWHEAT_INPUTS + SENESCWHEAT_OUTPUTS

#: the columns which define the topology of a roots in the input/output dataframe
ROOTS_TOPOLOGY_COLUMNS = ['plant', 'axis']

#: the columns which define the topology of an element in the input/output dataframe
ELEMENTS_TOPOLOGY_COLUMNS = ['plant', 'axis', 'metamer', 'organ', 'element']


def from_dataframes(roots_dataframe, elements_dataframe):
    """
    Convert inputs/outputs from Pandas dataframes to Senesc-Wheat format.
    
    :Parameters:
        
        - `roots_dataframe` (:class:`pandas.DataFrame`) - Roots inputs dataframe to convert, with one line by roots.
        
        - `elements_dataframe` (:class:`pandas.DataFrame`) - Elements inputs dataframe to convert, with one line by element.
    
    :Returns:
        The inputs/outputs in a dictionary.
    
    :Returns Type:
        :class:`dict` of :class:`dict`
        
    .. seealso:: see :attr:`simulation.Simulation.inputs` and :attr:`simulation.Simulation.outputs` 
       for the structure of Senesc-Wheat inputs/outputs.
    
    """
    all_roots_dict = {}
    all_elements_dict = {}
    for (all_current_dict, current_dataframe, current_topology_columns) in ((all_roots_dict, roots_dataframe, ROOTS_TOPOLOGY_COLUMNS), 
                                                                            (all_elements_dict, elements_dataframe, ELEMENTS_TOPOLOGY_COLUMNS)):
        current_columns = current_dataframe.columns.difference(current_topology_columns)
        for current_id, current_group in current_dataframe.groupby(current_topology_columns):
            current_series = current_group.loc[current_group.first_valid_index()]
            current_dict = current_series[current_columns].to_dict()
            all_current_dict[current_id] = current_dict

    return {'roots': all_roots_dict, 'elements': all_elements_dict}


def to_dataframes(data_dict):
    """
    Convert inputs/outputs from Senesc-Wheat format to Pandas dataframe.
    
    :Parameters:
        
        - `data_dict` (:class:`dict`) - The inputs/outputs in Senesc-Wheat format.
    
    :Returns:
        One dataframe for roots inputs/outputs and one dataframe for elements inputs/outputs.
    
    :Returns Type:
        :class:`tuple` of :class:`pandas.DataFrame`
        
    .. seealso:: see :attr:`simulation.Simulation.inputs` and :attr:`simulation.Simulation.outputs` 
       for the structure of Senesc-Wheat inputs/outputs.
        
    """
    dataframes_dict = {}
    for (current_key, current_topology_columns, current_inputs_outputs_names) in (('roots', ROOTS_TOPOLOGY_COLUMNS, SENESCWHEAT_ROOTS_INPUTS_OUTPUTS), 
                                                                                  ('elements', ELEMENTS_TOPOLOGY_COLUMNS, SENESCWHEAT_ELEMENTS_INPUTS_OUTPUTS)):
        current_data_dict = data_dict[current_key]
        current_ids_df = pd.DataFrame(current_data_dict.keys(), columns=current_topology_columns)
        current_data_df = pd.DataFrame(current_data_dict.values())
        current_df = pd.concat([current_ids_df, current_data_df], axis=1)
        current_df.sort_index(by=current_topology_columns, inplace=True)
        current_columns_sorted = current_topology_columns + [input_output for input_output in current_inputs_outputs_names if input_output in current_df.columns]
        current_df = current_df.reindex_axis(current_columns_sorted, axis=1, copy=False)
        current_df.reset_index(drop=True, inplace=True)
        dataframes_dict[current_key] = current_df
        
    return dataframes_dict['roots'], dataframes_dict['elements']


def from_MTG(g, roots_inputs, available_components):
    """
    Convert a MTG to Senesc-Wheat inputs. 
    
    :Parameters:
        
            - g (:class:`openalea.mtg.mtg.MTG`) - A MTG which contains the inputs
              of Senesc-Wheat. These inputs are: :mod:`SENESCWHEAT_INPUTS`. 
              
            - `roots_dataframe` (:class:`pandas.DataFrame`) - Roots dataframe, with one line by roots.
            
            - `available_components` - TODO: remove this argument
              
    :Returns:
        The inputs of Senesc-Wheat.
        
    :Returns Type:
        :class:`dict` of :class:`dict`
        
    .. seealso:: see :attr:`simulation.Simulation.inputs` for the structure of Senesc-Wheat inputs.
    
    .. todo:: remove argument `roots_dataframe` as soon as :mod:`openalea.mtg` permits to add/remove components.
        
    """
    all_roots_inputs_dict = {}
    all_elements_inputs_dict = {}
    
    # traverse the MTG recursively from top ...
    for plant_vertex_id in g.components_iter(g.root):
        if g.index(plant_vertex_id) not in available_components: continue
        plant_index = int(g.index(plant_vertex_id))
        for axis_vertex_id in g.components_iter(plant_vertex_id):
            if (g.index(plant_vertex_id), g.label(axis_vertex_id)) not in available_components: continue
            axis_id = g.label(axis_vertex_id)
            for metamer_vertex_id in g.components_iter(axis_vertex_id):
                if (g.index(plant_vertex_id), g.label(axis_vertex_id), g.index(metamer_vertex_id)) not in available_components: continue
                metamer_index = int(g.index(metamer_vertex_id))
                for organ_vertex_id in g.components_iter(metamer_vertex_id):
                    if (g.index(plant_vertex_id), g.label(axis_vertex_id), g.index(metamer_vertex_id), g.label(organ_vertex_id)) not in available_components: continue
                    organ_label = g.label(organ_vertex_id)
                    if organ_label not in SENESCWHEAT_PHOTOSYNTHETIC_ORGANS_NAMES:
                        raise NotModeledComponentError(organ_label, organ_vertex_id)
                    for element_vertex_id in g.components_iter(organ_vertex_id):
                        if (g.index(plant_vertex_id), g.label(axis_vertex_id), g.index(metamer_vertex_id), g.label(organ_vertex_id), g.label(element_vertex_id)) not in available_components: continue
                        vertex_properties = g.get_vertex_property(element_vertex_id)
                        element_label = g.label(element_vertex_id)
                        current_element_inputs = {}
                        for element_input in SENESCWHEAT_ELEMENTS_INPUTS:
                            if element_input not in vertex_properties:
                                raise PropertyNotFoundError(element_input, element_vertex_id)
                            current_element_inputs[element_input] = vertex_properties[element_input]
                        # ... and set the inputs
                        all_elements_inputs_dict[(plant_index, axis_id, metamer_index, organ_label, element_label)] = current_element_inputs
    
    #TODO: remove the following code as soon as :mod:`openalea.mtg` permits to add/remove components.
    for current_roots_inputs_id, current_roots_inputs_group in roots_inputs.groupby(ROOTS_TOPOLOGY_COLUMNS):
        current_roots_inputs_series = current_roots_inputs_group.loc[current_roots_inputs_group.first_valid_index()]
        current_roots_inputs_dict = current_roots_inputs_series[SENESCWHEAT_ROOTS_INPUTS].to_dict()
        all_roots_inputs_dict[current_roots_inputs_id] = current_roots_inputs_dict
    
    return {'roots': all_roots_inputs_dict, 'elements': all_elements_inputs_dict}


def update_MTG(outputs, g, roots_outputs_dataframe, available_components):
    """
    Update a MTG from Senesc-Wheat outputs.
    
    :Parameters:
        
            - outputs (:class:`dict` of :class:`dict`) - Senesc-Wheat outputs. 
            These outputs are: :mod:`SENESCWHEAT_OUTPUTS`. 
        
            - `g` (:class:`openalea.mtg.mtg.MTG`) - The MTG to update from the outputs of SenescWheat.
            
            - `roots_outputs_dataframe` (:class:`pandas.DataFrame`) - Roots outputs dataframe, with one line by roots.
            
            - `available_components` - TODO: remove this argument
            
    .. seealso:: see :attr:`simulation.Simulation.outputs` for the structure of Senesc-Wheat outputs.
    
    .. todo:: remove argument `roots_outputs_dataframe` as soon as :mod:`openalea.mtg` permits to add/remove components.        
    
    .. todo:: update the MTG from SenescWheat roots
    
    """
    # add the properties if needed
    property_names = g.property_names()
    for senescwheat_output_name in SENESCWHEAT_OUTPUTS:
        if senescwheat_output_name not in property_names:
            g.add_property(senescwheat_output_name)
    
    elements_outputs_dict = outputs['elements']
    # traverse the MTG recursively from top ...
    for plant_vertex_id in g.components_iter(g.root):
        if g.index(plant_vertex_id) not in available_components: continue
        plant_index = int(g.index(plant_vertex_id))
        for axis_vertex_id in g.components_iter(plant_vertex_id):
            if (g.index(plant_vertex_id), g.label(axis_vertex_id)) not in available_components: continue
            axis_id = g.label(axis_vertex_id)
            for metamer_vertex_id in g.components(axis_vertex_id): 
                if (g.index(plant_vertex_id), g.label(axis_vertex_id), g.index(metamer_vertex_id)) not in available_components: continue
                metamer_index = int(g.index(metamer_vertex_id))
                for organ_vertex_id in g.components_iter(metamer_vertex_id):
                    if (g.index(plant_vertex_id), g.label(axis_vertex_id), g.index(metamer_vertex_id), g.label(organ_vertex_id)) not in available_components: continue
                    organ_label = g.label(organ_vertex_id)
                    if organ_label not in SENESCWHEAT_PHOTOSYNTHETIC_ORGANS_NAMES:
                        raise NotModeledComponentError(organ_label, organ_vertex_id)
                    for element_vertex_id in g.components_iter(organ_vertex_id):
                        if (g.index(plant_vertex_id), g.label(axis_vertex_id), g.index(metamer_vertex_id), g.label(organ_vertex_id), g.label(element_vertex_id)) not in available_components: continue
                        element_label = g.label(element_vertex_id)
                        element_id = (plant_index, axis_id, metamer_index, organ_label, element_label)
                        if element_id not in elements_outputs_dict:
                            raise MismatchedTopologiesError(element_id, element_vertex_id)
                        element_outputs = elements_outputs_dict[element_id]
                        for output_name, output_value in element_outputs.iteritems():
                            # ... and set the properties
                            g.property(output_name)[element_vertex_id] = output_value
                            
    #TODO: remove the following code as soon as :mod:`openalea.mtg` permits to add/remove components.
    roots_outputs_dict = outputs['roots']
    for current_roots_id, current_roots_inputs in roots_outputs_dict.iteritems():
        plant_index = current_roots_id[0]
        axis_id = current_roots_id[1]
        roots_outputs_dataframe.loc[(roots_outputs_dataframe['plant'] == plant_index) & (roots_outputs_dataframe['axis'] == axis_id), current_roots_inputs.keys()] = current_roots_inputs.values()

