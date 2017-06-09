# -*- coding: latin-1 -*-
"""
    fspmwheat.tools
    ~~~~~~~~~~~~~~~

    This module provides convenient tools needed by the facades.

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

def combine_dataframes_inplace(model_dataframe, shared_column_indexes, shared_dataframe_to_update):
    """Combine `model_dataframe` and `shared_dataframe_to_update` in-place:

           * re-index `model_dataframe` and `shared_dataframe_to_update` by `shared_column_indexes`,
           * use method pd.DataFrame.combine_first(),
           * reset to the right types in `shared_dataframe_to_update`,
           * reorder the columns: first columns in `shared_column_indexes`, then others columns alphabetically,
           * and reset the index in `shared_dataframe_to_update`.

    :Parameters:

        - `model_dataframe` (:class:`pandas.DataFrame`) - The dataframe to use for updating `shared_dataframe_to_update`.

        - `shared_column_indexes` (:class:`list`) - The indexes to re-index `model_dataframe` and `shared_dataframe_to_update` before combining them.

        - `shared_dataframe_to_update` (:class:`pandas.DataFrame`) - The dataframe to update.

    .. note:: `shared_dataframe_to_update` is updated in-place. Thus, `shared_dataframe_to_update` keeps the same object's memory address.

    """

    # re-index the dataframes to have common indexes
    if len(shared_dataframe_to_update) == 0:
        shared_dataframe_to_update_reindexed = shared_dataframe_to_update
    else:
        shared_dataframe_to_update.sort_values(shared_column_indexes, inplace=True) ##RB
        shared_dataframe_to_update_reindexed = pd.DataFrame(shared_dataframe_to_update.values,
                                                            index=sorted(shared_dataframe_to_update.groupby(shared_column_indexes).groups.keys()),
                                                            columns=shared_dataframe_to_update.columns)

    model_dataframe.sort_values(shared_column_indexes, inplace=True) ##RB
    model_dataframe_reindexed = pd.DataFrame(model_dataframe.values,
                                             index=sorted(model_dataframe.groupby(shared_column_indexes).groups.keys()),
                                             columns=model_dataframe.columns)

    # combine model and shared re-indexed dataframes
    if model_dataframe_reindexed.empty and shared_dataframe_to_update.empty: #TODO: check with Camille
        new_shared_dataframe = model_dataframe_reindexed.copy()
        for new_header in shared_dataframe_to_update_reindexed.columns.difference(model_dataframe_reindexed.columns):
            new_shared_dataframe[new_header] = ""
    else:
        new_shared_dataframe = model_dataframe_reindexed.combine_first(shared_dataframe_to_update_reindexed) # if there are booleans in shared_dataframe_to_update, then these booleans are converted to floats

    # reset to the right types in the combined dataframe
    dtypes = model_dataframe_reindexed.dtypes.combine_first(shared_dataframe_to_update_reindexed.dtypes)
    for column_name, data_type in dtypes.iteritems():
        new_shared_dataframe[column_name] = new_shared_dataframe[column_name].astype(data_type)

    # reorder the columns
    new_shared_dataframe = new_shared_dataframe.reindex_axis(shared_column_indexes + sorted(new_shared_dataframe.columns.difference(shared_column_indexes)), axis=1)

    # update the shared dataframe in-place
    shared_dataframe_to_update.drop(shared_dataframe_to_update.index, axis=0, inplace=True)
    shared_dataframe_to_update.drop(shared_dataframe_to_update.columns, axis=1, inplace=True)
    shared_dataframe_to_update['dataframe_to_update_index'] = new_shared_dataframe.index
    shared_dataframe_to_update.set_index('dataframe_to_update_index', inplace=True)
    for column in new_shared_dataframe.columns:
        shared_dataframe_to_update[column] = new_shared_dataframe[column]
    shared_dataframe_to_update.reset_index(0, drop=True, inplace=True)