# -*- coding: latin-1 -*-

from __future__ import division # use "//" to do integer division

"""
    farquharwheat.converter
    ~~~~~~~~~~~~~~~~~~~~~~~

    The module :mod:`farquharwheat.converter` defines functions to:
        
        * convert :class:`dataframes <pandas.DataFrame>` to/from FarquharWheat :class:`inputs <simulation.Simulation.inputs>` or :class:`outputs <simulation.Simulation.outputs>`.
        * convert :class:`MTG <openalea.mtg.mtg.MTG>` to/from FarquharWheat :class:`inputs <simulation.Simulation.inputs>` or :class:`outputs <simulation.Simulation.outputs>`.
        
    :copyright: Copyright 2014 INRA-EGC, see AUTHORS.
    :license: TODO, see LICENSE for details.

    .. seealso:: Barillot et al. 2014.
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

import numpy as np
import pandas as pd
from alinea.astk import plantgl_utils

import simulation


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

#: the name, in the MTG, of the organs modeled by FarquharWheat
MTG_ORGANS_NAMES_SET = set(['internode', 'blade', 'sheath', 'peduncle', 'ear'])

#: the MTG properties, at organ scale, needed by FarquharWheat
MTG_ORGANS_PROPERTIES_NEEDED_BY_FARQUHARWHEAT = ['label']

#: the MTG properties, at element scale, needed by FarquharWheat
MTG_ELEMENTS_PROPERTIES_NEEDED_BY_FARQUHARWHEAT = ['surfacic_nitrogen', 'width', 'height', 'STAR']

#: the MTG properties needed by FarquharWheat
MTG_PROPERTIES_NEEDED_BY_FARQUHARWHEAT = MTG_ORGANS_PROPERTIES_NEEDED_BY_FARQUHARWHEAT + MTG_ELEMENTS_PROPERTIES_NEEDED_BY_FARQUHARWHEAT

#: the MTG properties computed by FarquharWheat
MTG_PROPERTIES_COMPUTED_BY_FARQUHARWHEAT = ['Ag', 'An', 'Rd', 'Tr', 'Ts', 'gs']

#: the inputs and outputs of FarquharWheat. 
FARQUHARWHEAT_INPUTS_OUTPUTS = MTG_ORGANS_PROPERTIES_NEEDED_BY_FARQUHARWHEAT + MTG_ELEMENTS_PROPERTIES_NEEDED_BY_FARQUHARWHEAT + MTG_PROPERTIES_COMPUTED_BY_FARQUHARWHEAT

#: the columns which define the topology of an element in the input/output dataframe
ELEMENTS_TOPOLOGY_COLUMNS = ['plant', 'axis', 'metamer', 'organ', 'element']


def from_dataframe(data_df):
    """
    Convert Pandas dataframe `data_df` to a dictionary.
    The dictionary has the same structure as :attr:`simulation.Simulation.inputs` 
    and :attr:`simulation.Simulation.outputs`.
    
    :Parameters:
        
        - `data_df` (:class:`pandas.DataFrame`) - The dataframe to convert, with one row by element. 
    
    :Returns:
        The data in a dictionary.
    
    :Returns Type:
        :class:`dict`
        
    """
    data_dict = {}
    columns_without_keys = data_df.columns.difference(ELEMENTS_TOPOLOGY_COLUMNS)
    for elements_id, element_group in data_df.groupby(ELEMENTS_TOPOLOGY_COLUMNS):
        data_dict[elements_id] = element_group.loc[element_group.first_valid_index(), columns_without_keys].to_dict()
    return data_dict


def to_dataframe(data_dict):
    """
    Convert the dictionary `data_dict` to Pandas dataframe.
    
    :Parameters:
        
        - `data_dict` (:class:`dict`) - The data to convert.
        The data has the same structure as :attr:`simulation.Simulation.inputs` 
        and :attr:`simulation.Simulation.outputs`.
    
    :Returns:
        The data in a dataframe, with one row by element.
    
    :Returns Type:
        :class:`pandas.DataFrame`
        
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
    Extract the inputs of Farquhar-Wheat from a MTG. 
    Farquhar-Wheat inputs are :mod:`MTG_PROPERTIES_COMPUTED_BY_FARQUHARWHEAT`.
    
    :Parameters:
        
            - g (:class:`openalea.mtg.mtg.MTG`) - A MTG which contains the inputs
              needed by Farquhar-Wheat to be run on each element.
              
    :Returns:
        The inputs of FarquharWheat by element.
        
        The inputs is a dictionary of dictionaries: {element1_id: element1_inputs, element2_id: element2_inputs, ..., elementN_id: elementN_inputs}, where:
         
            * elementi_id is a tuple: (plant_index, axis_id, metamer_index, organ_label, element_type),
            * and elementi_inputs is a dictionary: {'elementi_input1_name': elementi_input1_value, 'elementi_input2_name': elementi_input2_value, ..., 'elementi_inputN_name': elementi_inputN_value}.
         
        See :meth:`Model.calculate_An <farquharwheat.model.Model.calculate_An>` 
        for more information about the inputs. 
        
    :Returns Type:
        :class:`dict` of :class:`dict`
        
    """
    inputs = {}
    
    # check needed properties
    if not set(MTG_ORGANS_PROPERTIES_NEEDED_BY_FARQUHARWHEAT).issubset(g.properties()):
        raise InvalidMTGError(MTG_ORGANS_PROPERTIES_NEEDED_BY_FARQUHARWHEAT)
    
    # traverse the MTG recursively from top ...
    for plant_vid in g.components_iter(g.root):
        plant_index = int(g.index(plant_vid))
        for axis_vid in g.components_iter(plant_vid):
            axis_id = g.label(axis_vid)
            for metamer_vid in g.components_iter(axis_vid):
                metamer_index = int(g.index(metamer_vid))
                for organ_vid in g.components_iter(metamer_vid):
                    organ_label = g.label(organ_vid)
                    if organ_label not in MTG_ORGANS_NAMES_SET:
                        continue
                    vertex_properties = g.get_vertex_property(organ_vid)
                    organ_inputs = {}
                    is_missing = False
                    for mtg_property in MTG_ORGANS_PROPERTIES_NEEDED_BY_FARQUHARWHEAT:
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
                        for mtg_property in MTG_ELEMENTS_PROPERTIES_NEEDED_BY_FARQUHARWHEAT:
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
    Farquhar-Wheat outputs are :mod:`MTG_PROPERTIES_COMPUTED_BY_FARQUHARWHEAT`.
    
    :Parameters:
        
            - outputs (:class:`dict` of :class:`dict`) - The outputs by element. 
        
            `outputs` is a dictionary of dictionaries: {element1_id: element1_outputs, element2_id: element2_outputs, ..., elementN_id: elementN_outputs}, where:
             
                * elementi_id is a tuple: (plant_index, axis_id, metamer_index, organ_label, element_type),
                * and elementi_outputs is a dictionary: {'elementi_output1_name': elementi_output1_value, 'elementi_output2_name': elementi_output2_value, ..., 'elementi_outputN_name': elementi_outputN_value}.
        
            See :meth:`Model.calculate_An <farquharwheat.model.Model.calculate_An>` 
            for more information about the outputs.
        
            - `g` (:class:`openalea.mtg.mtg.MTG`) - The MTG to update from the `outputs` of FarquharWheat.
            
    """
    # add the properties if needed
    properties = g.properties()
    for farquharwheat_output_name in MTG_PROPERTIES_COMPUTED_BY_FARQUHARWHEAT:
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
                    if organ_label not in MTG_ORGANS_NAMES_SET:
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

