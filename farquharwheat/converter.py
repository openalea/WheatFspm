# -*- coding: latin-1 -*-

from __future__ import division # use "//" to do integer division

"""
    farquharwheat.converter
    ~~~~~~~~~~~~~~~~~~~~~~~

    The module :mod:`farquharwheat.converter` defines functions to:
        
        * convert a :class:`dataframe <pandas.DataFrame>` to/from FarquharWheat inputs or outputs format.
        * convert a :class:`MTG <openalea.mtg.mtg.MTG>` to/from FarquharWheat inputs or outputs format.
    
    Both dataframes and MTG follow AdelWheat naming convention.
        
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

import pandas as pd


#: the name of the organs modeled by FarquharWheat 
FARQUHARWHEAT_ORGANS_NAMES = set(['internode', 'blade', 'sheath', 'peduncle', 'ear'])

#: the inputs needed by FarquharWheat at organ scale
FARQUHARWHEAT_ORGANS_INPUTS = ['label']

#: the inputs needed by FarquharWheat at element scale
FARQUHARWHEAT_ELEMENTS_INPUTS = ['width', 'height', 'STAR', 'nitrates', 'amino_acids', 'proteins', 'Nstruct', 'green_area']

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

#: the columns which define the topology of an organ in the input/output dataframe
ORGANS_TOPOLOGY_COLUMNS = ['plant', 'axis', 'metamer', 'organ']

#: the columns which define the topology of an element in the input/output dataframe
ELEMENTS_TOPOLOGY_COLUMNS = ['plant', 'axis', 'metamer', 'organ', 'element']


def from_dataframes(organs_inputs, elements_inputs):
    """
    Convert inputs/outputs from Pandas dataframes to Farquhar-Wheat format.
    The column names of the dataframe respect the naming convention of AdelWheat.
    
    :Parameters:
        
        - `organs_inputs` (:class:`pandas.DataFrame`) - Organs inputs dataframe to convert, with one line by organ.
        
        - `elements_inputs` (:class:`pandas.DataFrame`) - Elements inputs dataframe to convert, with one line by element.
    
    :Returns:
        The inputs/outputs in a dictionary.
    
    :Returns Type:
        :class:`dict` of :class:`dict`
        
    .. seealso:: see :attr:`simulation.Simulation.inputs` and :attr:`simulation.Simulation.outputs` 
       for the structure of Farquhar-Wheat inputs/outputs.
    
    """
    all_organs_dict = {}
    all_elements_dict = {}
    for (all_current_dict, current_dataframe, current_topology_columns) in ((all_organs_dict, organs_inputs, ORGANS_TOPOLOGY_COLUMNS), 
                                                                            (all_elements_dict, elements_inputs, ELEMENTS_TOPOLOGY_COLUMNS)):
        current_columns = current_dataframe.columns.difference(current_topology_columns)
        for current_id, current_group in current_dataframe.groupby(current_topology_columns):
            current_series = current_group.loc[current_group.first_valid_index()]
            current_dict = current_series[current_columns].to_dict()
            all_current_dict[current_id] = current_dict

    return {'organs': all_organs_dict, 'elements': all_elements_dict}
    

def to_dataframes(data_dict):
    """
    Convert inputs/outputs from Farquhar-Wheat format to Pandas dataframe.
    The column names of the dataframe respect the naming convention of AdelWheat.
    
    :Parameters:
        
        - `data_dict` (:class:`dict`) - The inputs/outputs in Farquhar-Wheat format.
    
    :Returns:
        One dataframe for organs inputs/outputs and one dataframe for elements inputs/outputs.
    
    :Returns Type:
        :class:`tuple` of :class:`pandas.DataFrame`
        
    .. seealso:: see :attr:`simulation.Simulation.inputs` and :attr:`simulation.Simulation.outputs` 
       for the structure of Farquhar-Wheat inputs/outputs.
        
    """
    dataframes_dict = {}
    for (current_key, current_topology_columns, current_inputs_outputs_names) in (('organs', ORGANS_TOPOLOGY_COLUMNS, FARQUHARWHEAT_ORGANS_INPUTS_OUTPUTS), 
                                                                                  ('elements', ELEMENTS_TOPOLOGY_COLUMNS, FARQUHARWHEAT_ELEMENTS_INPUTS_OUTPUTS)):
        current_data_dict = data_dict[current_key]
        current_ids_df = pd.DataFrame(current_data_dict.keys(), columns=current_topology_columns)
        current_data_df = pd.DataFrame(current_data_dict.values())
        current_df = pd.concat([current_ids_df, current_data_df], axis=1)
        current_df.sort_index(by=current_topology_columns, inplace=True)
        current_columns_sorted = current_topology_columns + [input_output for input_output in current_inputs_outputs_names if input_output in current_df.columns]
        current_df = current_df.reindex_axis(current_columns_sorted, axis=1, copy=False)
        current_df.reset_index(drop=True, inplace=True)
        dataframes_dict[current_key] = current_df
        
    return dataframes_dict['organs'], dataframes_dict['elements']
    

def from_MTG(g, organs_inputs, elements_inputs):
    """
    Convert a MTG to Farquhar-Wheat inputs. 
    Use data in `organs_inputs` and `elements_inputs` if `g` is incomplete.
    The property names in the MTG respect the naming convention of AdelWheat.
    
    :Parameters:
        
            - g (:class:`openalea.mtg.mtg.MTG`) - A MTG which contains the inputs
              of Farquhar-Wheat. These inputs are: :mod:`FARQUHARWHEAT_INPUTS`.
              
            - `organs_inputs` (:class:`pandas.DataFrame`) - Organs dataframe, with one line by organ.
            
            - `elements_inputs` (:class:`pandas.DataFrame`) - Elements dataframe, with one line by element.
              
    :Returns:
        The inputs of Farquhar-Wheat.
        
    :Returns Type:
        :class:`dict` of :class:`dict`
        
    .. seealso:: see :attr:`simulation.Simulation.inputs` for the structure of Farquhar-Wheat inputs.
        
    """
    all_organs_inputs_dict = {}
    all_elements_inputs_dict = {}
    
    organs_inputs_grouped = organs_inputs.groupby(ORGANS_TOPOLOGY_COLUMNS)
    elements_inputs_grouped = elements_inputs.groupby(ELEMENTS_TOPOLOGY_COLUMNS)
    
    # traverse the MTG recursively from top ...
    for plant_vid in g.components_iter(g.root):
        plant_index = int(g.index(plant_vid))
        for axis_vid in g.components_iter(plant_vid):
            axis_label = g.label(axis_vid)
            for axis_component_vid in g.components_iter(axis_vid):
                if not g.label(axis_component_vid).startswith('metamer'): continue
                metamer_vid = axis_component_vid
                metamer_index = int(g.index(metamer_vid))
                for organ_vid in g.components_iter(metamer_vid):
                    vertex_properties = g.get_vertex_property(organ_vid)
                    organ_label = g.label(organ_vid)
                    if organ_label not in FARQUHARWHEAT_ORGANS_NAMES: continue
                    organ_id = (plant_index, axis_label, metamer_index, organ_label)
                    if organ_id in organs_inputs_grouped.groups:
                        organs_inputs_group = organs_inputs_grouped.get_group(organ_id)
                        organs_inputs_group_series = organs_inputs_group.loc[organs_inputs_group.first_valid_index(), organs_inputs_group.columns.intersection(FARQUHARWHEAT_ORGANS_INPUTS)]
                    else:
                        organs_inputs_group_series = pd.Series()
                    organ_inputs = {}
                    is_valid_organ = True
                    for organ_input_name in FARQUHARWHEAT_ORGANS_INPUTS:
                        if organ_input_name in vertex_properties:
                            # use the properties of the vertex
                            organ_inputs[organ_input_name] = vertex_properties[organ_input_name]
                        else:
                            # use the value in organs_inputs
                            if organ_input_name in organs_inputs_group_series:
                                organ_inputs[organ_input_name] = organs_inputs_group_series[organ_input_name]
                            else:
                                is_valid_organ = False
                                break
                    has_valid_element = False
                    for element_vid in g.components_iter(organ_vid):
                        vertex_properties = g.get_vertex_property(element_vid)
                        element_label = g.label(element_vid)
                        element_id = (plant_index, axis_label, metamer_index, organ_label, element_label)
                        if element_id in elements_inputs_grouped.groups:
                            elements_inputs_group = elements_inputs_grouped.get_group(element_id)
                            elements_inputs_group_series = elements_inputs_group.loc[elements_inputs_group.first_valid_index(), elements_inputs_group.columns.intersection(FARQUHARWHEAT_ELEMENTS_INPUTS)]
                        else:
                            elements_inputs_group_series = pd.Series()
                        element_inputs = {}
                        is_valid_element = True
                        for element_input_name in FARQUHARWHEAT_ELEMENTS_INPUTS:
                            if element_input_name in vertex_properties:
                                # use the properties of the vertex
                                element_inputs[element_input_name] = vertex_properties[element_input_name]
                            else:
                                # use the value in elements_inputs
                                if element_input_name in elements_inputs_group_series:
                                    element_inputs[element_input_name] = elements_inputs_group_series[element_input_name]
                                else:
                                    is_valid_element = False
                                    break
                        if is_valid_element:
                            has_valid_element = True
                            all_elements_inputs_dict[element_id] = element_inputs
                    if is_valid_organ and has_valid_element:
                        all_organs_inputs_dict[organ_id] = organ_inputs
    
    return {'organs': all_organs_inputs_dict, 'elements': all_elements_inputs_dict}


def update_MTG(inputs, outputs, g):
    """
    Update a MTG from Farquhar-Wheat inputs and outputs.
    The property names in the MTG respect the naming convention of AdelWheat.
    
    :Parameters:
    
            - inputs (:class:`dict` of :class:`dict`) - Farquhar-Wheat inputs. 
            These inputs are: :mod:`FARQUHARWHEAT_INPUTS`. 
        
            - outputs (:class:`dict` of :class:`dict`) - Farquhar-Wheat outputs. 
            These outputs are: :mod:`FARQUHARWHEAT_OUTPUTS`. 
        
            - `g` (:class:`openalea.mtg.mtg.MTG`) - The MTG to update from the inputs and outputs of FarquharWheat.
            
    .. seealso:: see :attr:`simulation.Simulation.inputs` and :attr:`simulation.Simulation.outputs` for the structure of Farquhar-Wheat inputs and outputs.
            
    """
    
    # add the properties if needed
    property_names = g.property_names()
    for farquharwheat_output_name in FARQUHARWHEAT_INPUTS_OUTPUTS:
        if farquharwheat_output_name not in property_names:
            g.add_property(farquharwheat_output_name)
    
    organs_inputs_dict = inputs['organs']
    elements_inputs_dict = inputs['elements']
    organs_outputs_dict = outputs['organs']
    elements_outputs_dict = outputs['elements']
    
    # traverse the MTG recursively from top ...
    for plant_vid in g.components_iter(g.root):
        plant_index = int(g.index(plant_vid))
        for axis_vid in g.components_iter(plant_vid):
            axis_label = g.label(axis_vid)
            for axis_component_vid in g.components_iter(axis_vid): 
                if not g.label(axis_component_vid).startswith('metamer'): continue
                metamer_vid = axis_component_vid
                metamer_index = int(g.index(metamer_vid))
                for organ_vid in g.components_iter(metamer_vid):
                    organ_label = g.label(organ_vid)
                    if organ_label not in FARQUHARWHEAT_ORGANS_NAMES: continue
                    organ_id = (plant_index, axis_label, metamer_index, organ_label)
                    if organ_id not in organs_outputs_dict: continue
                    # update the organ in the MTG
                    organ_inputs = organs_inputs_dict[organ_id]
                    for organ_input_name, organ_input_value in organ_inputs.iteritems():
                        g.property(organ_input_name)[organ_vid] = organ_input_value
                    organ_outputs = organs_outputs_dict[organ_id]
                    for organ_output_name, organ_output_value in organ_outputs.iteritems():
                        g.property(organ_output_name)[organ_vid] = organ_output_value
                    for element_vid in g.components_iter(organ_vid):
                        element_label = g.label(element_vid)
                        element_id = (plant_index, axis_label, metamer_index, organ_label, element_label)
                        if element_id not in elements_outputs_dict: continue
                        # update the element in the MTG
                        element_inputs = elements_inputs_dict[element_id]
                        for element_input_name, element_input_value in element_inputs.iteritems():
                            g.property(element_input_name)[element_vid] = element_input_value
                        element_outputs = elements_outputs_dict[element_id]
                        for element_output_name, element_output_value in element_outputs.iteritems():
                            g.property(element_output_name)[element_vid] = element_output_value
    
