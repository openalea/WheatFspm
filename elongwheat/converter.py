# -*- coding: latin-1 -*-

from __future__ import division # use "//" to do integer division

"""
    elongwheat.converter
    ~~~~~~~~~~~~~~~~~~~~~~~

    The module :mod:`elongwheat.converter` defines functions to convert
    :class:`dataframes <pandas.DataFrame>` to/from ElongWheat inputs or outputs format.

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

import simulation

#: the columns which define the topology in the input/output dataframe
HIDDENZONE_TOPOLOGY_COLUMNS = ['plant', 'axis', 'metamer']
ELEMENT_TOPOLOGY_COLUMNS = ['plant', 'axis', 'metamer', 'organ', 'element'] # Mature + emerging elements
SAM_TOPOLOGY_COLUMNS = ['plant', 'axis']


def from_dataframes(hiddenzone_inputs, element_inputs, SAM_inputs):
    """
    Convert inputs/outputs from Pandas dataframe to Elong-Wheat format.

    :Parameters:

        - `SAM_inputs` (:class:`pandas.DataFrame`) - Shoot Apical Meristem inputs dataframe to convert, with one line by SAM ie. one line per axis.
        - `hiddenzone_inputs` (:class:`pandas.DataFrame`) - Hidden zone inputs dataframe to convert, with one line by Hidden zone.
        - `element_inputs` (:class:`pandas.DataFrame`) - Emergeing and mature element inputs dataframe to convert, with one line by element.

    :Returns:
        The inputs in a dictionary.

    :Returns Type:
        :class:`dict` of :class:`dict`

    .. seealso:: see :attr:`simulation.Simulation.inputs` for the structure of Elong-Wheat inputs.

    """
    all_hiddenzone_dict = {}
    all_element_dict = {}
    all_SAM_dict = {}
    hiddenzone_L_calculation_dict = {}
    SAM_height_calculation_dict = {}

    hiddenzone_inputs_columns = hiddenzone_inputs.columns.difference(HIDDENZONE_TOPOLOGY_COLUMNS)
    emerging_element_inputs_columns = element_inputs.columns.difference(ELEMENT_TOPOLOGY_COLUMNS)
    SAM_inputs_columns = SAM_inputs.columns.difference(SAM_TOPOLOGY_COLUMNS)
    sheath_inputs_grouped = element_inputs[element_inputs.organ == 'sheath'].groupby(['plant', 'axis', 'metamer', 'organ'])

    for SAM_inputs_id, SAM_inputs_group in SAM_inputs.groupby(SAM_TOPOLOGY_COLUMNS):
        # SAM
        SAM_inputs_series = SAM_inputs_group.loc[SAM_inputs_group.first_valid_index()]
        SAM_inputs_dict = SAM_inputs_series[SAM_inputs_columns].to_dict()
        all_SAM_dict[SAM_inputs_id] = SAM_inputs_dict

        SAM_height_calculation_dict[SAM_inputs_id] = 0

    for element_inputs_id, element_inputs_group in element_inputs.groupby(ELEMENT_TOPOLOGY_COLUMNS):
        # Elements
        element_inputs_series = element_inputs_group.loc[element_inputs_group.first_valid_index()]
        element_inputs_dict = element_inputs_series[emerging_element_inputs_columns].to_dict()
        all_element_dict[element_inputs_id] = element_inputs_dict

    hiddenzone_inputs_grouped = hiddenzone_inputs.groupby(HIDDENZONE_TOPOLOGY_COLUMNS)
    for hiddenzone_inputs_id, hiddenzone_inputs_group in hiddenzone_inputs_grouped:
        # hiddenzone
        hiddenzone_inputs_series = hiddenzone_inputs_group.loc[hiddenzone_inputs_group.first_valid_index()]
        hiddenzone_inputs_dict = hiddenzone_inputs_series[hiddenzone_inputs_columns].to_dict()
        all_hiddenzone_dict[hiddenzone_inputs_id] = hiddenzone_inputs_dict

        # Get lengths required for the calculation of the distances before internode emergence and leaf emergence
        previous_hiddenzone_id = tuple(list(hiddenzone_inputs_id[:2]) + [hiddenzone_inputs_id[-1]-1])
        # previous distance to leaf emergence
        if hiddenzone_inputs_grouped.groups.has_key(previous_hiddenzone_id):
            previous_hiddenzone = hiddenzone_inputs_grouped.get_group(previous_hiddenzone_id)
            previous_leaf_dist_to_emerge = previous_hiddenzone.loc[previous_hiddenzone.first_valid_index(), 'leaf_dist_to_emerge']
        else:
            previous_leaf_dist_to_emerge = None

        # previous sheath length
        previous_sheath_id = tuple(list(hiddenzone_inputs_id[:2]) + [hiddenzone_inputs_id[-1]-1] + ['sheath'])
        if sheath_inputs_grouped.groups.has_key(previous_sheath_id):
            previous_sheath = sheath_inputs_grouped.get_group(previous_sheath_id)
            previous_sheath_visible_length =  previous_sheath[previous_sheath['element']=='StemElement']['length'].iloc[0]
            if not previous_leaf_dist_to_emerge: #: if no previous hiddenzone found, get the final hidden length of the previous sheath (assumes that no previous hiddenzone means a mature sheath)
                previous_sheath_hidden_length = previous_sheath[previous_sheath['element']=='HiddenElement']['length'].iloc[0]

        else:
            previous_sheath_visible_length = 0
            previous_sheath_hidden_length = 0

        # Growing internode length
        if hiddenzone_inputs_grouped.groups.has_key(hiddenzone_inputs_id):
            curr_hiddenzone = hiddenzone_inputs_grouped.get_group(hiddenzone_inputs_id)
            curr_internode_length = curr_hiddenzone.loc[curr_hiddenzone.first_valid_index(), 'internode_L']
            # For SAM height calculation
            SAM_id = tuple(curr_hiddenzone[SAM_TOPOLOGY_COLUMNS].iloc[0])
            SAM_height_calculation_dict[SAM_id] += curr_internode_length

        else:
            curr_internode_length = 0

        hiddenzone_L_calculation_dict[hiddenzone_inputs_id] = {'previous_leaf_dist_to_emerge': previous_leaf_dist_to_emerge,
                                                 'previous_sheath_visible_length': previous_sheath_visible_length,
                                                 'previous_sheath_hidden_length': previous_sheath_hidden_length,
                                                 'internode_length': curr_internode_length}

    return {'hiddenzone': all_hiddenzone_dict, 'elements': all_element_dict, 'SAM': all_SAM_dict , 'hiddenzone_L_calculation': hiddenzone_L_calculation_dict, 'SAM_height_calculation': SAM_height_calculation_dict}

def to_dataframes(data_dict):
    """
    Convert outputs from Elong-Wheat format to Pandas dataframe.

    :Parameters:

        - `data_dict` (:class:`dict`) - The outputs in Elong-Wheat format.

    :Returns:
        One dataframe for hiddenzone outputs, one dataframe for element outputs and one dataframe for SAM outputs.

    :Returns Type:
        :class:`tuple` of :class:`pandas.DataFrame`

    .. seealso:: see :attr:`simulation.Simulation.outputs` for the structure of Elong-Wheat outputs.

    """
    dataframes_dict = {}
    for (current_key, current_topology_columns, current_outputs_names) in (('hiddenzone', HIDDENZONE_TOPOLOGY_COLUMNS, simulation.HIDDENZONE_OUTPUTS),
                                                                           ('elements', ELEMENT_TOPOLOGY_COLUMNS, simulation.ELEMENT_OUTPUTS),
                                                                           ('SAM', SAM_TOPOLOGY_COLUMNS, simulation.SAM_OUTPUTS)):

        current_data_dict = data_dict[current_key]
        current_ids_df = pd.DataFrame(current_data_dict.keys(), columns=current_topology_columns)
        current_data_df = pd.DataFrame(current_data_dict.values())
        current_df = pd.concat([current_ids_df, current_data_df], axis=1)
        current_df.sort_values(by=current_topology_columns, inplace=True)
        current_columns_sorted = current_topology_columns + current_outputs_names
        current_df = current_df.reindex_axis(current_columns_sorted, axis=1, copy=False)
        current_df.reset_index(drop=True, inplace=True)
        dataframes_dict[current_key] = current_df

    return dataframes_dict['hiddenzone'], dataframes_dict['elements'], dataframes_dict['SAM']