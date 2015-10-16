# -*- coding: latin-1 -*-

from __future__ import division # use "//" to do integer division

"""
    farquharwheat.converter
    ~~~~~~~~~~~~~~~~~~~~~~~

    The module :mod:`farquharwheat.converter` defines functions to:
        
        * convert a :class:`dataframe <pandas.DataFrame>` to/from FarquharWheat inputs or outputs format.
        * convert a :class:`MTG <openalea.mtg.mtg.MTG>` to/from FarquharWheat inputs or outputs format.
        
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

import warnings

import pandas as pd


class ConverterWarning(UserWarning): pass


class NotModeledComponentWarning(ConverterWarning):
    '''MTG component not modeled by FarquharWheat.
    '''
    def __init__(self, component_label, component_vertex_id):
        self.message = 'MTG component "{0}" is not modeled by FarquharWheat (vertex {1}).'.format(component_label, component_vertex_id)
    def __str__(self):
        return repr(self.message)


warnings.simplefilter('always', ConverterWarning)


class ConverterError(Exception): pass


class PropertyNotFoundError(ConverterError):
    '''Property not found in a vertex of a MTG.
    '''
    def __init__(self, vertex_property, vertex_id):
        self.message = 'Property "{0}" not found in vertex {1}.'.format(vertex_property, vertex_id)
    def __str__(self):
        return repr(self.message)
    
    
class MismatchedTopologiesError(ConverterError):
    '''Topologies mismatched between FarquharWheat and MTG.
    '''
    def __init__(self, farquharwheat_id, vertex_id):
        self.message = 'Topologies mismatched between FarquharWheat and MTG: no mapping between FarquharWheat object {0} and vertex {1}.'.format(farquharwheat_id, vertex_id)
    def __str__(self):
        return repr(self.message)
    

#: the name of the organs modeled by FarquharWheat 
FARQUHARWHEAT_ORGANS_NAMES = set(['internode', 'blade', 'sheath', 'peduncle', 'ear'])

#: the inputs needed by FarquharWheat at organ scale
FARQUHARWHEAT_ORGANS_INPUTS = ['label']

#: the inputs needed by FarquharWheat at element scale
FARQUHARWHEAT_ELEMENTS_INPUTS = ['surfacic_nitrogen', 'width', 'height', 'STAR']

#: the inputs needed by FarquharWheat
FARQUHARWHEAT_INPUTS = FARQUHARWHEAT_ORGANS_INPUTS + FARQUHARWHEAT_ELEMENTS_INPUTS

#: the outputs computed by FarquharWheat at organ scale
FARQUHARWHEAT_ORGANS_OUTPUTS = []

#: the outputs computed by FarquharWheat at element scale
FARQUHARWHEAT_ELEMENTS_OUTPUTS = ['Ag', 'An', 'Rd', 'Tr', 'Ts', 'gs']

#: the outputs computed by FarquharWheat
FARQUHARWHEAT_OUTPUTS = FARQUHARWHEAT_ORGANS_OUTPUTS + FARQUHARWHEAT_ELEMENTS_OUTPUTS

#: the inputs and outputs of FarquharWheat at elements scale
FARQUHARWHEAT_ORGANS_INPUTS_OUTPUTS = FARQUHARWHEAT_ORGANS_INPUTS + FARQUHARWHEAT_ORGANS_OUTPUTS

#: the inputs and outputs of FarquharWheat at organs scale
FARQUHARWHEAT_ELEMENTS_INPUTS_OUTPUTS = FARQUHARWHEAT_ELEMENTS_INPUTS + FARQUHARWHEAT_ELEMENTS_OUTPUTS

#: the inputs and outputs of FarquharWheat. 
FARQUHARWHEAT_INPUTS_OUTPUTS = FARQUHARWHEAT_INPUTS + FARQUHARWHEAT_OUTPUTS

#: the columns which define the topology of an element in the input/output dataframe
ELEMENTS_TOPOLOGY_COLUMNS = ['plant', 'axis', 'metamer', 'organ', 'element']


def from_dataframe(data_df):
    """
    Convert inputs/outputs from Pandas dataframe to Farquhar-Wheat format.
    
    :Parameters:
        
        - `data_df` (:class:`pandas.DataFrame`) - The inputs/outputs to convert, with one row by element. 
    
    :Returns:
        The inputs/outputs in a dictionary.
    
    :Returns Type:
        :class:`dict` of :class:`dict`
        
    .. seealso:: see :attr:`simulation.Simulation.inputs` and :attr:`simulation.Simulation.outputs` 
       for the structure of Farquhar-Wheat inputs/outputs.
        
    """
    data_dict = {}
    columns_without_keys = data_df.columns.difference(ELEMENTS_TOPOLOGY_COLUMNS)
    for elements_id, element_group in data_df.groupby(ELEMENTS_TOPOLOGY_COLUMNS):
        data_dict[elements_id] = element_group.loc[element_group.first_valid_index(), columns_without_keys].to_dict()
    return data_dict


def to_dataframe(data_dict):
    """
    Convert inputs/outputs from Farquhar-Wheat format to Pandas dataframe.
    
    :Parameters:
        
        - `data_dict` (:class:`dict`) - The inputs/outputs in Farquhar-Wheat format.
    
    :Returns:
        The inputs/outputs in a dataframe, with one row by element.
    
    :Returns Type:
        :class:`pandas.DataFrame`
        
    .. seealso:: see :attr:`simulation.Simulation.inputs` and :attr:`simulation.Simulation.outputs` 
       for the structure of Farquhar-Wheat inputs/outputs.
        
    """
    elements_ids_df = pd.DataFrame(data_dict.keys(), columns=ELEMENTS_TOPOLOGY_COLUMNS)
    elements_data_df = pd.DataFrame(data_dict.values())
    data_df = pd.concat([elements_ids_df, elements_data_df], axis=1)
    data_df.sort_index(by=ELEMENTS_TOPOLOGY_COLUMNS, inplace=True)
    columns_sorted = ELEMENTS_TOPOLOGY_COLUMNS + [input_output for input_output in FARQUHARWHEAT_INPUTS_OUTPUTS if input_output in data_df.columns]
    data_df = data_df.reindex_axis(columns_sorted, axis=1, copy=False)
    data_df.reset_index(drop=True, inplace=True)
    return data_df


def from_MTG(g, available_components):
    """
    Convert a MTG to Farquhar-Wheat inputs. 
    
    :Parameters:
        
            - g (:class:`openalea.mtg.mtg.MTG`) - A MTG which contains the inputs
              of Farquhar-Wheat. These inputs are: :mod:`FARQUHARWHEAT_INPUTS`.
              
            - `available_components` - TODO: remove this argument 
              
    :Returns:
        The inputs of Farquhar-Wheat.
        
    :Returns Type:
        :class:`dict` of :class:`dict`
        
    .. seealso:: see :attr:`simulation.Simulation.inputs` for the structure of Farquhar-Wheat inputs.
        
    """
    inputs = {}
    
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
                    if organ_label not in FARQUHARWHEAT_ORGANS_NAMES:
                        warnings.warn(NotModeledComponentWarning(organ_label, organ_vertex_id))
                        continue
                    vertex_properties = g.get_vertex_property(organ_vertex_id)
                    organ_inputs = {}
                    for organ_input in FARQUHARWHEAT_ORGANS_INPUTS:
                        if organ_input not in vertex_properties:
                            raise PropertyNotFoundError(organ_input, organ_vertex_id)
                        organ_inputs[organ_input] = vertex_properties[organ_input]
                    for element_vertex_id in g.components_iter(organ_vertex_id):
                        if (g.index(plant_vertex_id), g.label(axis_vertex_id), g.index(metamer_vertex_id), g.label(organ_vertex_id), g.label(element_vertex_id)) not in available_components: continue
                        vertex_properties = g.get_vertex_property(element_vertex_id)
                        element_label = g.label(element_vertex_id)
                        current_inputs = organ_inputs.copy()
                        for element_input in FARQUHARWHEAT_ELEMENTS_INPUTS:
                            if element_input not in vertex_properties:
                                raise PropertyNotFoundError(element_input, element_vertex_id)
                            current_inputs[element_input] = vertex_properties[element_input]
                        # ... and set the inputs
                        inputs[(plant_index, axis_id, metamer_index, organ_label, element_label)] = current_inputs
    return inputs


def update_MTG(outputs, g, available_components):
    """
    Update a MTG from Farquhar-Wheat outputs.
    
    :Parameters:
        
            - outputs (:class:`dict` of :class:`dict`) - Farquhar-Wheat outputs. 
            These outputs are: :mod:`FARQUHARWHEAT_ELEMENTS_OUTPUTS`. 
        
            - `g` (:class:`openalea.mtg.mtg.MTG`) - The MTG to update from the outputs of FarquharWheat.
            
            - `available_components` - TODO: remove this argument
            
    .. seealso:: see :attr:`simulation.Simulation.outputs` for the structure of Farquhar-Wheat outputs.
            
    """
    # add the properties if needed
    property_names = g.property_names()
    for farquharwheat_output_name in FARQUHARWHEAT_OUTPUTS:
        if farquharwheat_output_name not in property_names:
            g.add_property(farquharwheat_output_name)
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
                    if organ_label not in FARQUHARWHEAT_ORGANS_NAMES:
                        warnings.warn(NotModeledComponentWarning(organ_label, organ_vertex_id))
                        continue
                    for element_vertex_id in g.components_iter(organ_vertex_id):
                        if (g.index(plant_vertex_id), g.label(axis_vertex_id), g.index(metamer_vertex_id), g.label(organ_vertex_id), g.label(element_vertex_id)) not in available_components: continue
                        element_label = g.label(element_vertex_id)
                        element_id = (plant_index, axis_id, metamer_index, organ_label, element_label)
                        if element_id not in outputs:
                            raise MismatchedTopologiesError(element_id, element_vertex_id)
                        element_outputs = outputs[element_id]
                        for output_name, output_value in element_outputs.iteritems():
                            # ... and set the properties
                            g.property(output_name)[element_vertex_id] = output_value

