# -*- coding: latin-1 -*-

from __future__ import division # use "//" to do integer division

"""
    farquharwheat.converter
    ~~~~~~~~~~~~~~~~~~~~~~~

    The module :mod:`farquharwheat.converter` defines functions to convert
    :class:`dataframes <pandas.DataFrame>` to/from FarquharWheat inputs or outputs format.

    :copyright: Copyright 2014-2016 INRA-ECOSYS, see AUTHORS.
    :license: TODO, see LICENSE for details.

    .. seealso:: Barillot et al. 2016.
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


#: the inputs needed by FarquharWheat at element scale
FARQUHARWHEAT_INPUTS = ['width', 'height', 'PARa', 'nitrates', 'amino_acids', 'proteins', 'Nstruct', 'green_area']

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

