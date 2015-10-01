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

class PropertyNotFoundWarning(ConverterWarning):
    '''Property not found in a vertex of a MTG.'''
    def __init__(self, property_, vertex_id):
        self.message = 'Property "{0}" not found in vertex {1}.'.format(property_, vertex_id)
    def __str__(self):
        return repr(self.message)
    
warnings.simplefilter('always', ConverterWarning)
    
    
class ConverterError(Exception): pass

class InvalidMTGError(ConverterError):
    '''The input MTG does not contain the required properties.'''
    def __init__(self, required_properties):
        self.message = 'The input MTG does not contain the required properties ({}).'.format(required_properties)
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

#: the outputs computed by FarquharWheat at element scale
FARQUHARWHEAT_ELEMENTS_OUTPUTS = ['Ag', 'An', 'Rd', 'Tr', 'Ts', 'gs']

#: the inputs and outputs of FarquharWheat. 
FARQUHARWHEAT_INPUTS_OUTPUTS = FARQUHARWHEAT_ORGANS_INPUTS + FARQUHARWHEAT_ELEMENTS_INPUTS + FARQUHARWHEAT_ELEMENTS_OUTPUTS

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
    return data_df


def from_MTG(g):
    """
    Convert a MTG to Farquhar-Wheat inputs. 
    
    :Parameters:
        
            - g (:class:`openalea.mtg.mtg.MTG`) - A MTG which contains the inputs
              of Farquhar-Wheat. These inputs are: :mod:`FARQUHARWHEAT_INPUTS`. 
              
    :Returns:
        The inputs of Farquhar-Wheat.
        
    :Returns Type:
        :class:`dict` of :class:`dict`
        
    .. seealso:: see :attr:`simulation.Simulation.inputs` for the structure of Farquhar-Wheat inputs.
        
    """
    inputs = {}
    
    # check needed properties
    if not set(FARQUHARWHEAT_ORGANS_INPUTS).issubset(g.properties()):
        raise InvalidMTGError(FARQUHARWHEAT_ORGANS_INPUTS)
    
    # traverse the MTG recursively from top ...
    for plant_vid in g.components_iter(g.root):
        plant_index = int(g.index(plant_vid))
        for axis_vid in g.components_iter(plant_vid):
            axis_id = g.label(axis_vid)
            for metamer_vid in g.components_iter(axis_vid):
                metamer_index = int(g.index(metamer_vid))
                for organ_vid in g.components_iter(metamer_vid):
                    organ_label = g.label(organ_vid)
                    if organ_label not in FARQUHARWHEAT_ORGANS_NAMES:
                        continue
                    vertex_properties = g.get_vertex_property(organ_vid)
                    organ_inputs = {}
                    is_missing = False
                    for mtg_property in FARQUHARWHEAT_ORGANS_INPUTS:
                        if mtg_property not in vertex_properties:
                            warnings.warn(PropertyNotFoundWarning(mtg_property, organ_vid))
                            is_missing = True
                            break
                        organ_inputs[mtg_property] = vertex_properties[mtg_property]
                    if is_missing:
                        continue
                    for element_vid in g.components_iter(organ_vid):
                        vertex_properties = g.get_vertex_property(element_vid)
                        element_label = g.label(element_vid)
                        current_inputs = organ_inputs.copy()
                        for mtg_property in FARQUHARWHEAT_ELEMENTS_INPUTS:
                            if mtg_property not in vertex_properties:
                                warnings.warn(PropertyNotFoundWarning(mtg_property, element_vid))
                                is_missing = True
                                break
                            current_inputs[mtg_property] = vertex_properties[mtg_property]
                        if is_missing:
                            continue
                        # ... and set the inputs
                        inputs[(plant_index, axis_id, metamer_index, organ_label, element_label)] = current_inputs
    return inputs


def update_MTG(outputs, g):
    """
    Update a MTG from Farquhar-Wheat outputs.
    
    :Parameters:
        
            - outputs (:class:`dict` of :class:`dict`) - Farquhar-Wheat outputs. 
            These outputs are: :mod:`FARQUHARWHEAT_ELEMENTS_OUTPUTS`. 
        
            - `g` (:class:`openalea.mtg.mtg.MTG`) - The MTG to update from the outputs of FarquharWheat.
            
    .. seealso:: see :attr:`simulation.Simulation.outputs` for the structure of Farquhar-Wheat outputs.
            
    """
    # add the properties if needed
    properties = g.properties()
    for farquharwheat_output_name in FARQUHARWHEAT_ELEMENTS_OUTPUTS:
        if farquharwheat_output_name not in properties:
            g.add_property(farquharwheat_output_name)
    # traverse the MTG recursively from top ...
    for plant_vid in g.components_iter(g.root):
        plant_index = int(g.index(plant_vid))
        for axis_vid in g.components_iter(plant_vid):
            axis_id = g.label(axis_vid)
            for metamer_vid in g.components(axis_vid): 
                metamer_index = int(g.index(metamer_vid))
                for organ_vid in g.components_iter(metamer_vid):
                    organ_label = g.label(organ_vid)
                    if organ_label not in FARQUHARWHEAT_ORGANS_NAMES:
                        continue
                    for element_vid in g.components_iter(organ_vid):
                        element_label = g.label(element_vid)
                        element_id = (plant_index, axis_id, metamer_index, organ_label, element_label)
                        if element_id not in outputs:
                            continue
                        element_outputs = outputs[element_id]
                        for output_name, output_value in element_outputs.iteritems():
                            # ... and set the properties
                            g.property(output_name)[element_vid] = output_value

