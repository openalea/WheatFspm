# -*- coding: latin-1 -*-

import sys
import os
import pandas as pd
import main
import tools
import additional_graphs

my_dir = os.getcwd()

# Scenarii
scenarii_df = pd.read_csv( os.path.join('inputs','scenarii_list.csv'), index_col='Scenario')
scenarii_df['Scenario'] = scenarii_df.index
scenarii = scenarii_df.Scenario


def run_one_scenario(param):
    if type(param) != list:
        scenario = param
    elif len(param) > 1:
        scenario = int(param[1])
    else:
        sys.exit('Error : No Scenario value as input')

    scenario_name = 'Scenario_{}'.format(scenario)

    # Create the directory of the Scenario where results will be stored
    scenario_dirpath = scenario_name
    if not os.path.exists(scenario_dirpath):
        os.mkdir(scenario_dirpath)

    # Create directory paths for graphs, outputs and postprocessings of this scneario
    scenario_graphs_dirpath = os.path.join(scenario_dirpath, 'graphs')
    if not os.path.exists(scenario_graphs_dirpath):
        os.mkdir(scenario_graphs_dirpath)
    # Outputs
    scenario_outputs_dirpath = os.path.join(scenario_dirpath, 'outputs')
    if not os.path.exists(scenario_outputs_dirpath):
        os.mkdir(scenario_outputs_dirpath)
    # Postprocessings
    scenario_postprocessing_dirpath = os.path.join(scenario_dirpath, 'postprocessing')
    if not os.path.exists(scenario_postprocessing_dirpath):
        os.mkdir(scenario_postprocessing_dirpath)

    # Create dict of parameters for the scenario
    scenario_dict = tools.buildDic(scenarii_df.to_dict('index')[scenario])

    # Run main fspmwheat
    try:
        main.main(3000, forced_start_time=0, run_simu=True, run_postprocessing=True, generate_graphs=True, run_from_outputs=False,
                  option_static=False, show_3Dplant = False,
                  tillers_replications = None, heterogeneous_canopy = True,
                  N_fertilizations = {'constant_Conc_Nitrates': scenario_dict.get('constant_Conc_Nitrates') },
                  PLANT_DENSITY={1:scenario_dict.get('Plant_Density', 250.)},
                  INPUTS_DIRPATH='inputs',
                  METEO_FILENAME=scenario_dict.get('METEO_FILENAME'),
                  GRAPHS_DIRPATH = scenario_graphs_dirpath,
                  OUTPUTS_DIRPATH = scenario_outputs_dirpath,
                  POSTPROCESSING_DIRPATH=scenario_postprocessing_dirpath,
                  update_parameters_all_models = scenario_dict)
        additional_graphs.graph_RER(scenario)
        additional_graphs.graph_summary(scenario)
        additional_graphs.table_C_usages(scenario)

    except Exception as e:
        print e
    pass

if __name__ == '__main__':
    run_one_scenario(sys.argv)