from __future__ import print_function
import os
import sys
import getopt

import pandas as pd

import main
import tools
import additional_postprocessing


def run_fspmwheat(scenario_id, scenarii_file, simulation_length, outputs_dir_path):

    # Path of the directory which contains the inputs of the model
    INPUTS_DIRPATH = 'inputs'

    # Scenario to be run
    if scenariifile:  # Scenarii csv file provided by user
        scenarii_df = pd.read_csv(scenarii_file, index_col='Scenario')
    else:  # Use default scenarii csv file
        scenarii_df = pd.read_csv(os.path.join(INPUTS_DIRPATH, 'scenarii_list.csv'), index_col='Scenario')

    scenario = scenarii_df.loc[scenario_id].to_dict()
    scenario_name = 'Scenario_{}'.format(scenario['Scenario_label'])

    # Create dict of parameters for the scenario
    update_parameters = tools.buildDic(scenario)

    # Create the directory of the Scenario where results will be stored
    if outputs_dir_path:
        scenario_dirpath = os.path.join(outputs_dir_path, scenario_name)
    else:
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

    # -- SIMULATION PARAMETERS --

    # Do run the simulation?
    RUN_SIMU = True

    # Do run the postprocessing?
    RUN_POSTPROCESSING = False  #: TODO separate postprocessings coming from other models

    # Do generate the graphs?
    GENERATE_GRAPHS = False  #: TODO separate postprocessings coming from other models

    # Run main fspmwheat
    try:
        main.main(simulation_length, forced_start_time=0, run_simu=RUN_SIMU, run_postprocessing=RUN_POSTPROCESSING, generate_graphs=GENERATE_GRAPHS,
                  run_from_outputs=False, option_static=False,
                  tillers_replications=None, heterogeneous_canopy=True,
                  N_fertilizations={'constant_Conc_Nitrates': scenario.get('constant_Conc_Nitrates')},
                  PLANT_DENSITY={1: scenario.get('Plant_Density', 250.)},
                  INPUTS_DIRPATH=INPUTS_DIRPATH,
                  METEO_FILENAME=scenario.get('METEO_FILENAME'),
                  GRAPHS_DIRPATH=scenario_graphs_dirpath,
                  OUTPUTS_DIRPATH=scenario_outputs_dirpath,
                  POSTPROCESSING_DIRPATH=scenario_postprocessing_dirpath,
                  update_parameters_all_models=update_parameters)
        additional_postprocessing.table_C_usages(scenario_name)
        additional_postprocessing.calculate_performance_indices(scenario_name)

    except Exception as e:
        print(e)


if __name__ == '__main__':
    scenariifile = None
    simlength = 3000
    outputs = None

    try:
        opts, args = getopt.getopt(sys.argv[1:], "s:l:o:d", ["scenariifile=", "simlength=", "outputs="])
    except getopt.GetoptError as err:
        print(str(err))
        sys.exit(2)

    for opt, arg in opts:
        if opt in ("-s", "--scenariifile"):
            scenariifile = arg
        elif opt in ("-l", "--simlength"):
            simlength = int(arg)
        elif opt in ("-o", "--outputs"):
            outputs = arg
    # slurm = os.environ['SLURM_ARRAY_TASK_ID']
    run_fspmwheat(scenario_id=os.environ['SLURM_ARRAY_TASK_ID'], scenarii_file=scenariifile, simulation_length=simlength, outputs_dir_path=outputs)
