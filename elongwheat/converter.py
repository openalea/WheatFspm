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
ORGAN_TOPOLOGY_COLUMNS = ['plant', 'axis', 'metamer', 'organ'] # exposed organs
SAM_TOPOLOGY_COLUMNS = ['plant', 'axis']


def from_dataframes(hiddenzone_inputs, organ_inputs, SAM_inputs):
    """
    Convert inputs/outputs from Pandas dataframe to Elong-Wheat format.

    :Parameters:

        - `SAM_inputs` (:class:`pandas.DataFrame`) - Shoot Apical Meristem inputs dataframe to convert, with one line by SAM ie. one line per axis.
        - `hiddenzone_inputs` (:class:`pandas.DataFrame`) - Hidden zone inputs dataframe to convert, with one line by Hidden zone.
        - `organ_inputs` (:class:`pandas.DataFrame`) - Exposed organ inputs dataframe to convert, with one line by organ.

    :Returns:
        The inputs in a dictionary.

    :Returns Type:
        :class:`dict` of :class:`dict`

    .. seealso:: see :attr:`simulation.Simulation.inputs` for the structure of Elong-Wheat inputs.

    """
    all_hiddenzone_dict = {}
    all_organ_dict = {}
    all_SAM_dict = {}
    all__dict = {}
    hiddenzone_L_calculation_dict = {}
    hiddenzone_inputs_columns = hiddenzone_inputs.columns.difference(HIDDENZONE_TOPOLOGY_COLUMNS)
    organ_inputs_columns = organ_inputs.columns.difference(ORGAN_TOPOLOGY_COLUMNS)
    SAM_inputs_columns = SAM_inputs.columns.difference(SAM_TOPOLOGY_COLUMNS)
    sheath_inputs_grouped = organ_inputs[organ_inputs.organ == 'sheath'].groupby(ORGAN_TOPOLOGY_COLUMNS)
    sheath_inputs_grouped_all_metamers = organ_inputs[organ_inputs.organ == 'sheath'].groupby(['plant', 'axis'])

    for SAM_inputs_id, SAM_inputs_group in SAM_inputs.groupby(SAM_TOPOLOGY_COLUMNS):
        # SAM
        SAM_inputs_series = SAM_inputs_group.loc[SAM_inputs_group.first_valid_index()]
        SAM_inputs_dict = SAM_inputs_series[SAM_inputs_columns].to_dict()
        all_SAM_dict[SAM_inputs_id] = SAM_inputs_dict

    for organ_inputs_id, organ_inputs_group in organ_inputs.groupby(ORGAN_TOPOLOGY_COLUMNS):
        # organ
        organ_inputs_series = organ_inputs_group.loc[organ_inputs_group.first_valid_index()]
        organ_inputs_dict = organ_inputs_series[organ_inputs_columns].to_dict()
        all_organ_dict[organ_inputs_id] = organ_inputs_dict

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
            previous_sheath_visible_length = previous_sheath.loc[previous_sheath.first_valid_index(), 'visible_length']
            if not previous_leaf_dist_to_emerge: #: if no previous hiddenzone found, get the final hidden length of the previous sheath (assumes that no previous hiddenzone means a mature sheath)
                previous_sheath_final_hidden_length = previous_sheath.loc[previous_sheath.first_valid_index(), 'final_hidden_length']

        else:
            previous_sheath_visible_length = 0
            previous_sheath_final_hidden_length = 0

        # internode length
        if hiddenzone_inputs_grouped.groups.has_key(hiddenzone_inputs_id):
            curr_hiddenzone = hiddenzone_inputs_grouped.get_group(hiddenzone_inputs_id)
            curr_internode_length = curr_hiddenzone.loc[curr_hiddenzone.first_valid_index(), 'internode_L']
        else:
            curr_internode_length = 0

        hiddenzone_L_calculation_dict[hiddenzone_inputs_id] = {'previous_leaf_dist_to_emerge': previous_leaf_dist_to_emerge,
                                                 'previous_sheath_visible_length': previous_sheath_visible_length,
                                                 'previous_sheath_final_hidden_length': previous_sheath_final_hidden_length,
                                                 'internode_length': curr_internode_length}

    return {'hiddenzone': all_hiddenzone_dict, 'organs': all_organ_dict, 'SAM': all_SAM_dict , 'hiddenzone_L_calculation': hiddenzone_L_calculation_dict}

def to_dataframes(data_dict):
    """
    Convert outputs from Elong-Wheat format to Pandas dataframe.

    :Parameters:

        - `data_dict` (:class:`dict`) - The outputs in Elong-Wheat format.

    :Returns:
        One dataframe for hiddenzone outputs and one dataframe for organ outputs.

    :Returns Type:
        :class:`tuple` of :class:`pandas.DataFrame`

    .. seealso:: see :attr:`simulation.Simulation.outputs` for the structure of Elong-Wheat outputs.

    """
    dataframes_dict = {}
    for (current_key, current_topology_columns, current_outputs_names) in (('hiddenzone', HIDDENZONE_TOPOLOGY_COLUMNS, simulation.HIDDENZONE_OUTPUTS),
                                                                           ('organs', ORGAN_TOPOLOGY_COLUMNS, simulation.ORGAN_OUTPUTS),
                                                                           ('SAM', SAM_TOPOLOGY_COLUMNS, simulation.SAM_OUTPUTS) ):
        current_data_dict = data_dict[current_key]
        current_ids_df = pd.DataFrame(current_data_dict.keys(), columns=current_topology_columns)
        current_data_df = pd.DataFrame(current_data_dict.values())
        current_df = pd.concat([current_ids_df, current_data_df], axis=1)
        current_df.sort_values(by=current_topology_columns, inplace=True)
        current_columns_sorted = current_topology_columns + current_outputs_names
        current_df = current_df.reindex_axis(current_columns_sorted, axis=1, copy=False)
        current_df.reset_index(drop=True, inplace=True)
        dataframes_dict[current_key] = current_df
    return dataframes_dict['hiddenzone'], dataframes_dict['organs'], dataframes_dict['SAM']
