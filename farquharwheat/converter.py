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
FARQUHARWHEAT_ELEMENTS_INPUTS = ['width', 'height', 'PARa', 'nitrates', 'amino_acids', 'proteins', 'Nstruct', 'green_area']
#: the inputs needed by FarquharWheat at axis scale
FARQUHARWHEAT_SAMS_INPUTS = ['T_SAM', 'H_Canopy']

#: the outputs computed by FarquharWheat
FARQUHARWHEAT_ELEMENTS_OUTPUTS = ['Ag', 'An', 'Rd', 'Tr', 'Ts', 'gs','width', 'height']

#: the inputs and outputs of FarquharWheat.
FARQUHARWHEAT_ELEMENTS_INPUTS_OUTPUTS = set(FARQUHARWHEAT_ELEMENTS_INPUTS + FARQUHARWHEAT_ELEMENTS_OUTPUTS)

#: the columns which define the topology in the input/output elements dataframe
ELEMENT_TOPOLOGY_COLUMNS = ['plant', 'axis', 'metamer', 'organ', 'element']
#: the columns which define the topology in the input/output elements dataframe
SAM_TOPOLOGY_COLUMNS = ['plant', 'axis']


def from_dataframe(element_inputs, SAM_inputs):
    """
    Convert inputs/outputs from Pandas dataframe to Farquhar-Wheat format.

    :Parameters:

        - `SAM_inputs` (:class:`pandas.DataFrame`) - Shoot Apical Meristem inputs dataframe to convert, with one line by SAM ie. one line per axis.
        - `element_inputs` (:class:`pandas.DataFrame`) - Emergeing and mature element inputs dataframe to convert, with one line by element.

    :Returns:
        The inputs/outputs in a dictionary.

    :Returns Type:
        :class:`dict` of :class:`dict`

    .. seealso:: see :attr:`simulation.Simulation.inputs` and :attr:`simulation.Simulation.outputs`
       for the structure of Farquhar-Wheat inputs/outputs.

    """
    all_elements_dict = {}
    data_columns = element_inputs.columns.difference( ELEMENT_TOPOLOGY_COLUMNS )
    for current_id, current_group in element_inputs.groupby( ELEMENT_TOPOLOGY_COLUMNS ):
        current_series = current_group.loc[current_group.first_valid_index()]
        current_dict = current_series[data_columns].to_dict()
        all_elements_dict[current_id] = current_dict

    all_SAMs_dict = {}
    data_columns = SAM_inputs.columns.difference(SAM_TOPOLOGY_COLUMNS)
    for current_id, current_group in SAM_inputs.groupby(SAM_TOPOLOGY_COLUMNS):
        current_series = current_group.loc[current_group.first_valid_index()]
        current_dict = current_series[data_columns].to_dict()
        all_SAMs_dict[current_id] = current_dict

    return {'elements': all_elements_dict, 'SAMs': all_SAMs_dict}


def to_dataframe(data_dict):
    """
    Convert inputs/outputs from Farquhar-Wheat format to Pandas dataframe.

    :Parameters:

        - `data_dict` (:class:`dict`) - The inputs/outputs in Farquhar-Wheat format.

    :Returns:
         - one dataframe for element outputs

    :Returns Type:
        :class:`pandas.DataFrame`

    .. seealso:: see :attr:`simulation.Simulation.inputs` and :attr:`simulation.Simulation.outputs`
       for the structure of Farquhar-Wheat inputs/outputs.

    """
    ids_df = pd.DataFrame(data_dict.keys(), columns=ELEMENT_TOPOLOGY_COLUMNS)
    data_df = pd.DataFrame(data_dict.values())
    df = pd.concat([ids_df, data_df], axis=1)
    df.sort_values(by=ELEMENT_TOPOLOGY_COLUMNS, inplace=True)
    columns_sorted = ELEMENT_TOPOLOGY_COLUMNS + [column_name for column_name in FARQUHARWHEAT_ELEMENTS_INPUTS_OUTPUTS if column_name in df.columns]
    df = df.reindex_axis(columns_sorted, axis=1, copy=False)
    df.reset_index(drop=True, inplace=True)
    return df

