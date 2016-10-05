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

import pandas as pd


#: the name of the organs modeled by FarquharWheat
FARQUHARWHEAT_ORGANS_NAMES = set(['internode', 'blade', 'sheath', 'peduncle', 'ear'])

#: the inputs needed by FarquharWheat at element scale
FARQUHARWHEAT_INPUTS = ['width', 'height', 'Eabsm2', 'nitrates', 'amino_acids', 'proteins', 'Nstruct', 'green_area']

#: the outputs computed by FarquharWheat
FARQUHARWHEAT_OUTPUTS = ['Ag', 'An', 'Rd', 'Tr', 'Ts', 'gs']

#: the inputs and outputs of FarquharWheat.
FARQUHARWHEAT_INPUTS_OUTPUTS = FARQUHARWHEAT_INPUTS + FARQUHARWHEAT_OUTPUTS

#: the columns which define the topology in the input/output dataframe
DATAFRAME_TOPOLOGY_COLUMNS = ['plant', 'axis', 'metamer', 'organ', 'element']


def from_dataframe(dataframe):
    """
    Convert inputs/outputs from Pandas dataframe to Farquhar-Wheat format.

    :Parameters:

        - `dataframe` (:class:`pandas.DataFrame`) - Inputs/outputs dataframe to convert, with one line by element.

    :Returns:
        The inputs/outputs in a dictionary.

    :Returns Type:
        :class:`dict` of :class:`dict`

    .. seealso:: see :attr:`simulation.Simulation.inputs` and :attr:`simulation.Simulation.outputs`
       for the structure of Farquhar-Wheat inputs/outputs.

    """
    data_dict = {}
    data_columns = dataframe.columns.difference(DATAFRAME_TOPOLOGY_COLUMNS)
    for current_id, current_group in dataframe.groupby(DATAFRAME_TOPOLOGY_COLUMNS):
        current_series = current_group.loc[current_group.first_valid_index()]
        current_dict = current_series[data_columns].to_dict()
        data_dict[current_id] = current_dict
    return data_dict


def to_dataframe(data_dict):
    """
    Convert inputs/outputs from Farquhar-Wheat format to Pandas dataframe.

    :Parameters:

        - `data_dict` (:class:`dict`) - The inputs/outputs in Farquhar-Wheat format.

    :Returns:
        The inputs/outputs in a dataframe.

    :Returns Type:
        :class:`pandas.DataFrame`

    .. seealso:: see :attr:`simulation.Simulation.inputs` and :attr:`simulation.Simulation.outputs`
       for the structure of Farquhar-Wheat inputs/outputs.

    """
    ids_df = pd.DataFrame(data_dict.keys(), columns=DATAFRAME_TOPOLOGY_COLUMNS)
    data_df = pd.DataFrame(data_dict.values())
    df = pd.concat([ids_df, data_df], axis=1)
    df.sort_values(by=DATAFRAME_TOPOLOGY_COLUMNS, inplace=True)
    columns_sorted = DATAFRAME_TOPOLOGY_COLUMNS + [column_name for column_name in FARQUHARWHEAT_INPUTS_OUTPUTS if column_name in df.columns]
    df = df.reindex_axis(columns_sorted, axis=1, copy=False)
    df.reset_index(drop=True, inplace=True)
    return df


def from_MTG(g, inputs):
    """
    Convert a MTG to Farquhar-Wheat inputs.
    Use data in `inputs` if `g` is incomplete.

    :Parameters:

            - g (:class:`openalea.mtg.mtg.MTG`) - A MTG which contains the inputs
              of Farquhar-Wheat. These inputs are: :mod:`FARQUHARWHEAT_INPUTS`.

            - `inputs` (:class:`pandas.DataFrame`) - Inputs dataframe, with one line by element.

    :Returns:
        The inputs of Farquhar-Wheat.

    :Returns Type:
        :class:`dict`

    .. seealso:: see :attr:`simulation.Simulation.inputs` for the structure of Farquhar-Wheat inputs.

    """
    all_inputs_dict = {}

    inputs_grouped = inputs.groupby(DATAFRAME_TOPOLOGY_COLUMNS)

    # traverse the MTG recursively from top ...
    for plant_vid in g.components_iter(g.root):
        plant_index = int(g.index(plant_vid))
        for axis_vid in g.components_iter(plant_vid):
            axis_label = g.label(axis_vid)
            for metamer_vid in g.components_iter(axis_vid):
                metamer_index = int(g.index(metamer_vid))
                for organ_vid in g.components_iter(metamer_vid):
                    vertex_properties = g.get_vertex_property(organ_vid)
                    organ_label = g.label(organ_vid)
                    if organ_label not in FARQUHARWHEAT_ORGANS_NAMES: continue
                    for element_vid in g.components_iter(organ_vid):
                        vertex_properties = g.get_vertex_property(element_vid)
                        element_label = g.label(element_vid)
                        element_id = (plant_index, axis_label, metamer_index, organ_label, element_label)
                        if element_id in inputs_grouped.groups:
                            elements_inputs_group = inputs_grouped.get_group(element_id)
                            elements_inputs_group_series = elements_inputs_group.loc[elements_inputs_group.first_valid_index(), elements_inputs_group.columns.intersection(FARQUHARWHEAT_INPUTS)]
                        else:
                            elements_inputs_group_series = pd.Series()
                        element_inputs = {}
                        is_valid_element = True

                        for element_input_name in FARQUHARWHEAT_INPUTS:
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
                            all_inputs_dict[element_id] = element_inputs

    return all_inputs_dict


def update_MTG(inputs, outputs, g):
    """
    Update a MTG from Farquhar-Wheat inputs and outputs.

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

    # traverse the MTG recursively from top ...
    for plant_vid in g.components_iter(g.root):
        plant_index = int(g.index(plant_vid))
        for axis_vid in g.components_iter(plant_vid):
            axis_label = g.label(axis_vid)
            for metamer_vid in g.components_iter(axis_vid):
                metamer_index = int(g.index(metamer_vid))
                for organ_vid in g.components_iter(metamer_vid):
                    organ_label = g.label(organ_vid)
                    if organ_label not in FARQUHARWHEAT_ORGANS_NAMES: continue
                    for element_vid in g.components_iter(organ_vid):
                        element_label = g.label(element_vid)
                        element_id = (plant_index, axis_label, metamer_index, organ_label, element_label)
                        if element_id not in outputs: continue
                        # update the element in the MTG
                        element_inputs = inputs[element_id]
                        for element_input_name, element_input_value in element_inputs.iteritems():
                            g.property(element_input_name)[element_vid] = element_input_value
                        element_outputs = outputs[element_id]
                        for element_output_name, element_output_value in element_outputs.iteritems():
                            g.property(element_output_name)[element_vid] = element_output_value

