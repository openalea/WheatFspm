import os
import pandas as pd

# Get the list of scenarii
titi = os.walk('.').next()[1]
scenarii = []
for i in titi:
    if i[:9] == 'Scenario_':
        scenarii.append(int(i[9:]))
scenarii_inputs = scenarii

INDEX_COLUMNS = ['t', 'plant', 'axis', 'metamer', 'organ', 'element']


def rearrange_postprocessing(postprocessing_tables, t=None, scenarii=None, merge_with_scenarii_list=True, scenarii_list_columns=None):
    """
    For each graph type, create a common directory with the graphs of all the scenarii.

    :param list postprocessing_tables: List of names of postprocessing tables that we want to compare for each scenario.
    :param int t: Time step of the model for which we want to compare the postprocessing of each scenario.
    :param list scenarii: List of scenarii (numbers) for which the graphs will be rearranged.
    :param bool merge_with_scenarii_list: if True, the output datafram will be merge with a subset of the scenarii_list table (only the columns in scenarii_list_columns)
    :param list scenarii_list_columns: List of columns name of the scenarii_list table that are of interest.
    """

    if scenarii is None:
        scenarii = scenarii_inputs

    # Create directory with general postprocessing
    general_res_dir = os.path.join('General_postprocessing')
    if not os.path.exists(general_res_dir):
        os.mkdir(general_res_dir)

    # read the description of the scenarii if necessary
    if merge_with_scenarii_list:
        if scenarii_list_columns is None:
            scenarii_list_columns = ['Scenario', 'Scenario_label', 'Plant_Density',
                                     'constant_Conc_Nitrates', 'fertilization_U_3500', 'fertilization_quantity', 'fertilization_interval', 'PAR']
        scenarii_df = pd.read_csv(os.path.join('inputs', 'scenarii_list.csv'), index_col='Scenario')
        scenarii_df['Scenario'] = scenarii_df.index

        scenarii_list_columns = [c for c in scenarii_list_columns if c in scenarii_df.columns] + [c for c in scenarii_df.columns if c.startswith('P_')]

    #: For each postprocessing table, extract the prostprocessings at t of the scenarii
    for pp_table in postprocessing_tables:
        res = []
        res_index_column = []
        res_file_name = None
        scenarii_with_res = []

        #: For each scenario, extract the prostprocessings at t
        for scenario in scenarii:

            scenario_name = 'Scenario_%.4d' % scenario # 'Scenario_' + str(scenario)
            scenario_postprocessing_dirpath = os.path.join(scenario_name, 'postprocessing')
            scenario_postprocessing_table_dirpath = os.path.join(scenario_postprocessing_dirpath, pp_table + '.csv')

            if not os.path.exists(scenario_postprocessing_table_dirpath):
                continue
            pp_res = pd.read_csv(scenario_postprocessing_table_dirpath)
            if 't' in pp_res.columns:
                if t not in pp_res.t:
                    continue
                scenarii_with_res.append(scenario)
                res.append(pp_res[pp_res.t == t])
                if not res_file_name:
                    res_file_name = pp_table + '_' + str(int(t))
            else:
                scenarii_with_res.append(scenario)
                res.append(pp_res)
                if not res_file_name:
                    res_file_name = pp_table

            # Find the index of the postprocessing table
            if not res_index_column:
                res_index_column = ['Scenario']
                res_index_column.extend([c for c in INDEX_COLUMNS if c in pp_res.columns])

        # convert list of postprocessing into datafram
        res_df = pd.concat(res, keys=scenarii_with_res, sort=False)
        res_df.reset_index(0, inplace=True)
        res_df.rename({'level_0': 'Scenario'}, axis=1, inplace=True)
        res_df = res_df.reindex(res_index_column + res_df.columns.difference(res_index_column).tolist(), axis=1, copy=False)

        # Add info about the scenario
        if merge_with_scenarii_list:
            scenarii_description = scenarii_df[scenarii_list_columns]
            res_df = scenarii_description.merge(res_df, how='outer', on='Scenario')

        # write the datafram
        res_dir = os.path.join(general_res_dir, res_file_name + '.csv')
        res_df.to_csv(res_dir, na_rep='NA', index=False)


if __name__ == '__main__':
    # rearrange_postprocessing(postprocessing_tables=['axes_postprocessing'], t=1999)
    # rearrange_postprocessing(postprocessing_tables=['axes_postprocessing'], t=3499)
    rearrange_postprocessing(postprocessing_tables=['organs_postprocessing'], t=3499)
    # rearrange_postprocessing(postprocessing_tables=['performance_indices'])
    # rearrange_postprocessing(postprocessing_tables=['leaf_traits'])
