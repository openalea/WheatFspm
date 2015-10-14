# -*- coding: latin-1 -*-

'''
    main
    ~~~~

    An example to show how to couple models CN-Wheat, Farquhar-Wheat and Senesc-Wheat.
    This example uses the format MTG to exchange data between the models. 

    You must first install :mod:`cnwheat` and :mod:`farquharwheat` (and add them to your PYTHONPATH)
    before running this script with the command `python`.

    :copyright: Copyright 2014-2015 INRA-ECOSYS, see AUTHORS.
    :license: TODO, see LICENSE for details.

    .. seealso:: Barillot et al. 2015.
'''

'''
    Information about this versioned file:
        $LastChangedBy$
        $LastChangedDate$
        $LastChangedRevision$
        $URL$
        $Id$
'''

import os
import time, datetime
import profile, pstats

import logging

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from cnwheat import simulation as cnwheat_simulation, model as cnwheat_model, parameters as cnwheat_parameters, converter as cnwheat_converter, tools as cnwheat_tools
from farquharwheat import simulation as farquharwheat_simulation, model as farquharwheat_model, converter as farquharwheat_converter
from senescwheat import simulation as senescwheat_simulation, model as senescwheat_model, converter as senescwheat_converter

INPUTS_DIRPATH = 'inputs'
PLANTS_INPUTS_FILENAME = 'plants_inputs.csv' # TODO: load plants inputs from MTG file
AXES_INPUTS_FILENAME = 'axes_inputs.csv' # TODO: load axes inputs from MTG file 
PHYTOMERS_INPUTS_FILENAME = 'phytomers_inputs.csv' # TODO: load phytomers inputs from MTG file
ORGANS_INPUTS_OUTPUTS_FILENAME = 'organs_inputs_outputs.csv' # TODO: load organs inputs from MTG file ; TODO: store the state of system in a MTG
ELEMENTS_INPUTS_FILENAME = 'elements_inputs.csv' # TODO: load elements inputs from MTG file

PAR_FILENAME = 'PAR_Clermont_rebuild.csv'
METEO_FILENAME = 'meteo_Clermont_rebuild.csv'

OUTPUTS_DIRPATH = 'outputs' # the outputs dataframes also include inputs
PLANTS_OUTPUTS_FILENAME = 'plants_outputs.csv'
AXES_OUTPUTS_FILENAME = 'axes_outputs.csv'
PHYTOMERS_OUTPUTS_FILENAME = 'phytomers_outputs.csv'
ORGANS_OUTPUTS_FILENAME = 'organs_outputs.csv'
ELEMENTS_OUTPUTS_FILENAME = 'elements_outputs.csv'

GRAPHS_DIRPATH = 'graphs'

OUTPUTS_PRECISION = 6

LOGGING_CONFIG_FILEPATH = os.path.join('..', 'logging.json')

LOGGING_LEVEL = logging.INFO # set to logging.DEBUG for debugging or logging.INFO for INFO level

cnwheat_tools.setup_logging(LOGGING_CONFIG_FILEPATH, LOGGING_LEVEL, log_model=False, log_compartments=False, log_derivatives=False)


def setup_MTG():
    """Construct and fill a valid MTG.
    TODO: add roots, phloem, grains and soil to the MTG.
    TODO: load the fulfilled MTG from file. 
    """
    import random
    import numpy as np
    import pandas as pd
    from alinea.adel.astk_interface import initialise_stand
    
    random.seed(1234)
    np.random.seed(1234)
    
    elements_inputs_df = pd.read_csv(os.path.join(INPUTS_DIRPATH, ELEMENTS_INPUTS_FILENAME))
    
    g, wheat, domain_area, domain, convUnit = initialise_stand(1500)
    
    MTG_INPUTS_PROPERTIES = set(senescwheat_converter.SENESCWHEAT_INPUTS).union(set(farquharwheat_converter.FARQUHARWHEAT_INPUTS) - set(senescwheat_converter.SENESCWHEAT_INPUTS_OUTPUTS)).union(set(cnwheat_converter.CNWHEAT_INPUTS) - set(senescwheat_converter.SENESCWHEAT_INPUTS_OUTPUTS) - set(farquharwheat_converter.FARQUHARWHEAT_INPUTS_OUTPUTS))
    
    property_names = g.property_names()
    for cnwheat_variable_name in MTG_INPUTS_PROPERTIES:
        if cnwheat_variable_name not in property_names:
            g.add_property(cnwheat_variable_name)

    MTG_ELEMENTS_INPUTS_PROPERTIES = set(senescwheat_converter.SENESCWHEAT_ELEMENTS_INPUTS).union(set(farquharwheat_converter.FARQUHARWHEAT_ELEMENTS_INPUTS) - set(senescwheat_converter.SENESCWHEAT_ELEMENTS_INPUTS_OUTPUTS)).union(set(cnwheat_simulation.Simulation.ELEMENTS_STATE) - set(senescwheat_converter.SENESCWHEAT_ELEMENTS_INPUTS_OUTPUTS) - set(farquharwheat_converter.FARQUHARWHEAT_ELEMENTS_INPUTS_OUTPUTS)) 
    
    # traverse the MTG recursively from top
    for plant_vid in g.components_iter(g.root):
        plant_index = int(g.index(plant_vid))
        current_plant_group = elements_inputs_df[elements_inputs_df['plant'] == plant_index]
        if len(current_plant_group) == 0:
            # TODO: remove plant_vid and all its components
            continue
        for axis_vid in g.components_iter(plant_vid):
            axis_id = g.label(axis_vid)
            current_axis_group = current_plant_group[current_plant_group['axis'] == axis_id]
            if len(current_axis_group) == 0:
                # TODO: remove axis_vid and all its components
                continue
            for metamer_vid in g.components(axis_vid): 
                metamer_index = int(g.index(metamer_vid))
                current_metamer_group = current_axis_group[current_axis_group['phytomer'] == metamer_index]
                if len(current_metamer_group) == 0:
                    # TODO: remove metamer_vid and all its components
                    continue
                for organ_vid in g.components_iter(metamer_vid):
                    organ_label = g.label(organ_vid)
                    if organ_label not in cnwheat_converter.MTG_TO_CNWHEAT_PHYTOMERS_ORGANS_NAMES_MAPPING:
                        continue
                    current_organ_group = current_metamer_group[current_metamer_group['organ'] == organ_label]
                    if len(current_organ_group) == 0:
                        # TODO: remove organ_vid and all its components
                        continue
                    for element_vid in g.components_iter(organ_vid):
                        element_label = g.label(element_vid)
                        if element_label not in cnwheat_converter.MTG_TO_CNWHEAT_ELEMENTS_NAMES_MAPPING:
                            continue
                        current_element_group = current_organ_group[current_organ_group['element'] == element_label]
                        if len(current_element_group) == 0:
                            # TODO: remove element_vid and all its components
                            continue
                        current_element_series = current_element_group.loc[current_element_group.first_valid_index()]
                        for property_ in MTG_ELEMENTS_INPUTS_PROPERTIES:
                            g.property(property_)[element_vid] = current_element_series[property_]
    return g


def compute_CN_distrib(run_simu=True, make_graphs=True):

    plants_outputs_filepath = os.path.join(OUTPUTS_DIRPATH, PLANTS_OUTPUTS_FILENAME)
    axes_outputs_filepath = os.path.join(OUTPUTS_DIRPATH, AXES_OUTPUTS_FILENAME)
    phytomers_outputs_filepath = os.path.join(OUTPUTS_DIRPATH, PHYTOMERS_OUTPUTS_FILENAME)
    organs_outputs_filepath = os.path.join(OUTPUTS_DIRPATH, ORGANS_OUTPUTS_FILENAME)
    elements_outputs_filepath = os.path.join(OUTPUTS_DIRPATH, ELEMENTS_OUTPUTS_FILENAME)
    
    meteo_df = pd.read_csv(os.path.join(INPUTS_DIRPATH, METEO_FILENAME), index_col='t')
    PAR_df = pd.read_csv(os.path.join(INPUTS_DIRPATH, PAR_FILENAME), index_col=farquharwheat_converter.ELEMENTS_TOPOLOGY_COLUMNS)
    PARi_df = PAR_df.copy()
    # compute incident PAR ; TODO: why don't we read PARi directly ?
    PARi_df.PAR *= 0.9 * 0.95
    
    if run_simu:

        t0 = time.time()
        
        # define the time step of each simulator 
        senescwheat_ts = 2
        farquharwheat_ts = 2
        cnwheat_ts = 1
        
        g = setup_MTG()
        
        # create the simulators
        senescwheat_simulation_ = senescwheat_simulation.Simulation(time_step=senescwheat_ts)
        farquharwheat_simulation_ = farquharwheat_simulation.Simulation()
        cnwheat_simulation_ = cnwheat_simulation.Simulation()
        
        all_organs_inputs_ouputs_df = pd.read_csv(os.path.join(INPUTS_DIRPATH, ORGANS_INPUTS_OUTPUTS_FILENAME)) # TODO: remove this code as soon as :mod:`openalea.mtg` permits to add/remove components.

        PARi_grouped = PARi_df.groupby('t')

        # define the time grid of the simulators
        start_time = 0
        stop_time = 960

        all_plants_inputs_ouputs_df_list = []
        all_axes_inputs_ouputs_df_list = []
        all_phytomers_inputs_ouputs_df_list = []
        all_organs_inputs_ouputs_df_list = []
        all_elements_inputs_ouputs_df_list = []

        # TODO: this is to avoid converting error ; TODO: to remove as soon as :mod:`openalea.mtg` permits to add/remove components.
        elements_inputs_df = pd.read_csv(os.path.join(INPUTS_DIRPATH, ELEMENTS_INPUTS_FILENAME)).astype(str)
        available_components = set(elements_inputs_df.groupby(cnwheat_simulation.Simulation.PLANTS_INPUTS_INDEXES).groups.keys() + \
                                   elements_inputs_df.groupby(cnwheat_simulation.Simulation.AXES_INPUTS_INDEXES).groups.keys() + \
                                   elements_inputs_df.groupby(cnwheat_simulation.Simulation.PHYTOMERS_INPUTS_INDEXES).groups.keys() + \
                                   elements_inputs_df.groupby(cnwheat_simulation.Simulation.ORGANS_INPUTS_INDEXES).groups.keys() + \
                                   elements_inputs_df.groupby(cnwheat_simulation.Simulation.ELEMENTS_INPUTS_INDEXES).groups.keys())

        for t_senescwheat in xrange(start_time, stop_time, senescwheat_ts):

            senescwheat_roots_inputs_df = all_organs_inputs_ouputs_df.loc[all_organs_inputs_ouputs_df.organ == 'Roots'].dropna(axis=1).drop('organ', axis=1)
            senescwheat_simulation_.initialize(senescwheat_converter.from_MTG(g, senescwheat_roots_inputs_df, available_components))
            senescwheat_simulation_.run()
            senescwheat_roots_outputs_df, senescwheat_elements_outputs_df = senescwheat_converter.to_dataframes(senescwheat_simulation_.outputs)
            senescwheat_converter.update_MTG(senescwheat_simulation_.outputs, g, senescwheat_roots_outputs_df, available_components)
            all_organs_inputs_ouputs_df.loc[all_organs_inputs_ouputs_df.organ == 'Roots', senescwheat_roots_outputs_df.columns] = senescwheat_roots_outputs_df.loc[:, senescwheat_roots_outputs_df.columns].values
            _, senescwheat_elements_inputs_df = senescwheat_converter.to_dataframes(senescwheat_simulation_.inputs)
            senescwheat_elements_inputs_outputs_df = senescwheat_elements_outputs_df.combine_first(senescwheat_elements_inputs_df)
            
            for t_farquharwheat in xrange(t_senescwheat, t_senescwheat + senescwheat_ts, farquharwheat_ts):
                
                farquharwheat_simulation_.initialize(farquharwheat_converter.from_MTG(g, available_components))
                Ta, ambient_CO2, RH, Ur = meteo_df.loc[t_farquharwheat, ['air_temperature', 'ambient_CO2', 'humidity', 'Wind', ]]
                PARi_at_t_farquharwheat = PARi_grouped.get_group(t_farquharwheat).to_dict()['PAR']
                farquharwheat_simulation_.run(Ta, ambient_CO2, RH, Ur, PARi_at_t_farquharwheat)
                farquharwheat_converter.update_MTG(farquharwheat_simulation_.outputs, g, available_components)
                farquharwheat_elements_inputs_df = farquharwheat_converter.to_dataframe(farquharwheat_simulation_.inputs)
                farquharwheat_elements_outputs_df = farquharwheat_converter.to_dataframe(farquharwheat_simulation_.outputs)
                farquharwheat_elements_inputs_outputs_df = farquharwheat_elements_outputs_df.combine_first(farquharwheat_elements_inputs_df)
                senescwheat_farquharwheat_elements_inputs_outputs_combined_df = farquharwheat_elements_inputs_outputs_df.combine_first(senescwheat_elements_inputs_outputs_df) 
                
                for t_cnwheat in xrange(t_farquharwheat, t_farquharwheat + farquharwheat_ts, cnwheat_ts):
                    
                    cnwheat_simulation_.initialize(cnwheat_converter.from_MTG(g, all_organs_inputs_ouputs_df, available_components))
                    cnwheat_simulation_.run(start_time=t_cnwheat, stop_time=t_cnwheat+cnwheat_ts, number_of_output_steps=cnwheat_ts+1)
                    _, _, _, cnwheat_organs_inputs_outputs_df, _ = cnwheat_converter.to_dataframes(cnwheat_simulation_.population)
                    cnwheat_converter.update_MTG(cnwheat_simulation_.population, g, cnwheat_organs_inputs_outputs_df, available_components)
                    
                    (cnwheat_plants_postprocessing_df, 
                     cnwheat_axes_postprocessing_df,
                     cnwheat_phytomers_postprocessing_df,
                     cnwheat_organs_postprocessing_df,
                     cnwheat_elements_postprocessing_df) = cnwheat_simulation_.postprocessings()
                    
                    all_plants_inputs_ouputs_df = cnwheat_plants_postprocessing_df
                    all_axes_inputs_ouputs_df = cnwheat_axes_postprocessing_df
                    all_phytomers_inputs_ouputs_df = cnwheat_phytomers_postprocessing_df
                    
                    cnwheat_organs_postprocessing_df = cnwheat_organs_postprocessing_df.loc[cnwheat_organs_postprocessing_df.t == t_cnwheat, :].reset_index(drop=True)
                    all_organs_inputs_ouputs_df = cnwheat_organs_postprocessing_df.combine_first(cnwheat_organs_inputs_outputs_df).combine_first(all_organs_inputs_ouputs_df)
                    
                    cnwheat_elements_postprocessing_df = cnwheat_elements_postprocessing_df.loc[cnwheat_elements_postprocessing_df.t == t_cnwheat, :]
                    senescwheat_farquharwheat_elements_inputs_outputs_combined_df.rename_axis({'metamer': 'phytomer'}, axis=1, inplace=True) # TODO: should we standardize the names of the variables? 
                    senescwheat_farquharwheat_elements_inputs_outputs_combined_df.replace({'HiddenElement': 'enclosed', 'StemElement': 'exposed', 'LeafElement1': 'exposed'}, inplace=True) # TODO: should we standardize the names of the variables?
                    senescwheat_farquharwheat_elements_inputs_outputs_combined_df.replace({'internode': 'Internode', 'blade': 'Lamina', 'sheath': 'Sheath', 'ear': 'Chaff', 'peduncle': 'Peduncle'}, inplace=True) # TODO: should we standardize the names of the variables?
                    senescwheat_farquharwheat_elements_inputs_outputs_combined_df.sort_index(by=cnwheat_simulation.Simulation.ELEMENTS_INPUTS_INDEXES, inplace=True)
                    senescwheat_farquharwheat_elements_inputs_outputs_combined_df.reset_index(drop=True, inplace=True)
                    all_elements_inputs_ouputs_df = cnwheat_elements_postprocessing_df.combine_first(senescwheat_farquharwheat_elements_inputs_outputs_combined_df)
                    
                    all_plants_inputs_ouputs_df_list.append(all_plants_inputs_ouputs_df)
                    all_axes_inputs_ouputs_df_list.append(all_axes_inputs_ouputs_df)
                    all_phytomers_inputs_ouputs_df_list.append(all_phytomers_inputs_ouputs_df)
                    all_organs_inputs_ouputs_df_list.append(all_organs_inputs_ouputs_df)
                    all_elements_inputs_ouputs_df_list.append(all_elements_inputs_ouputs_df)
                    
        
        global_plants_df = pd.concat(all_plants_inputs_ouputs_df_list, ignore_index=True)
        global_plants_df.drop_duplicates(subset=cnwheat_simulation.Simulation.PLANTS_OUTPUTS_INDEXES, inplace=True)
        global_plants_df.to_csv(plants_outputs_filepath, na_rep='NA', index=False, float_format='%.{}f'.format(OUTPUTS_PRECISION))

        global_axes_df = pd.concat(all_axes_inputs_ouputs_df_list, ignore_index=True)
        global_axes_df.drop_duplicates(subset=cnwheat_simulation.Simulation.AXES_OUTPUTS_INDEXES, inplace=True)
        global_axes_df.to_csv(axes_outputs_filepath, na_rep='NA', index=False, float_format='%.{}f'.format(OUTPUTS_PRECISION))

        global_phytomers_df = pd.concat(all_phytomers_inputs_ouputs_df_list, ignore_index=True)
        global_phytomers_df.drop_duplicates(subset=cnwheat_simulation.Simulation.PHYTOMERS_OUTPUTS_INDEXES, inplace=True)
        global_phytomers_df.to_csv(phytomers_outputs_filepath, na_rep='NA', index=False, float_format='%.{}f'.format(OUTPUTS_PRECISION))

        global_organs_df = pd.concat(all_organs_inputs_ouputs_df_list, ignore_index=True)
        global_organs_df.drop_duplicates(subset=cnwheat_simulation.Simulation.ORGANS_OUTPUTS_INDEXES, inplace=True)
        global_organs_df.to_csv(organs_outputs_filepath, na_rep='NA', index=False, float_format='%.{}f'.format(OUTPUTS_PRECISION))

        global_elements_df = pd.concat(all_elements_inputs_ouputs_df_list, ignore_index=True)
        global_elements_df.drop_duplicates(subset=cnwheat_simulation.Simulation.ELEMENTS_OUTPUTS_INDEXES, inplace=True)
        all_elements_inputs_outputs_variables = [column for column in global_elements_df.columns if column not in cnwheat_simulation.Simulation.ELEMENTS_OUTPUTS_INDEXES]
        global_elements_df = global_elements_df.reindex_axis(cnwheat_simulation.Simulation.ELEMENTS_OUTPUTS_INDEXES + all_elements_inputs_outputs_variables, axis=1, copy=False)
        global_elements_df.sort_index(by=cnwheat_simulation.Simulation.ELEMENTS_OUTPUTS_INDEXES, inplace=True)
        global_elements_df.reset_index(drop=True, inplace=True)
#         global_elements_df = global_elements_df.convert_objects(copy=False)
#         global_elements_df[['plant', 'phytomer']] = global_elements_df[['plant', 'phytomer']].astype(int)
        global_elements_df.to_csv(elements_outputs_filepath, na_rep='NA', index=False, float_format='%.{}f'.format(OUTPUTS_PRECISION))
        
        execution_time = int(time.time()-t0)
        print '\n', 'Model executed in ', str(datetime.timedelta(seconds=execution_time))

    ######POST-PROCESSING##
    if make_graphs:
        
        if not os.path.isdir(GRAPHS_DIRPATH):
            os.mkdir(GRAPHS_DIRPATH)
        
        x_name = 't'
        x_label='Time (Hour)'
        
        # 1) Photosynthetic organs
        ph_elements_output_df = pd.read_csv(elements_outputs_filepath)

        graph_variables_ph_elements = {'Ag': u'Gross photosynthesis (µmol m$^{-2}$ s$^{-1}$)','An': u'Net photosynthesis (µmol m$^{-2}$ s$^{-1}$)', 'Tr':u'Organ surfacic transpiration rate (mmol H$_{2}$0 m$^{-2}$ s$^{-1}$)', 'Transpiration':u'Organ transpiration rate (mmol H$_{2}$0 s$^{-1}$)', 'Rd': u'Mitochondrial respiration rate of organ in light (µmol C h$^{-1}$)', 'Ts': u'Temperature surface (°C)', 'gs': u'Conductance stomatique (mol m$^{-2}$ s$^{-1}$)',
                           'Conc_TriosesP': u'[TriosesP] (µmol g$^{-1}$ mstruct)', 'Conc_Starch':u'[Starch] (µmol g$^{-1}$ mstruct)', 'Conc_Sucrose':u'[Sucrose] (µmol g$^{-1}$ mstruct)', 'Conc_Fructan':u'[Fructan] (µmol g$^{-1}$ mstruct)',
                           'Conc_Nitrates': u'[Nitrates] (µmol g$^{-1}$ mstruct)', 'Conc_Amino_Acids': u'[Amino_Acids] (µmol g$^{-1}$ mstruct)', 'Conc_Proteins': u'[Proteins] (g g$^{-1}$ mstruct)', 'SLN': u'Surfacic nitrogen content (g m$^{-2}$)',
                           'Nitrates_import': u'Total nitrates imported (µmol h$^{-1}$)', 'Amino_Acids_import': u'Total amino acids imported (µmol N h$^{-1}$)',
                           'S_Amino_Acids': u'[Rate of amino acids synthesis] (µmol N g$^{-1}$ mstruct h$^{-1}$)', 'S_Proteins': u'Rate of protein synthesis (µmol N g$^{-1}$ mstruct h$^{-1}$)', 'D_Proteins': u'Rate of protein degradation (µmol N g$^{-1}$ mstruct h$^{-1}$)', 'k_proteins': u'Relative rate of protein degradation (µmol N g$^{-1}$ mstruct h$^{-1}$)',
                           'Loading_Sucrose': u'Loading Sucrose (µmol C sucrose h$^{-1}$)', 'Loading_Amino_Acids': u'Loading Amino acids (µmol N amino acids h$^{-1}$)',
                           'green_area': u'Green area (m$^{2}$)', 'R_phloem_loading': u'Respiration phloem loading (µmol C h$^{-1}$)', 'R_Nnit_red': u'Respiration nitrate reduction (µmol C h$^{-1}$)', 'R_residual': u'Respiration residual (µmol C h$^{-1}$)', 'R_maintenance': u'Respiration residual (µmol C h$^{-1}$)',
                           'mstruct': u'Structural mass (g)', 'Nstruct': u'Structural N mass (g)',
                           'Conc_cytokinines':u'[cytokininess] (UA g$^{-1}$ mstruct)', 'D_cytokinines':u'cytokinines degradation (UA g$^{-1}$ mstruct)', 'cytokinines_import':u'cytokinines import (UA)'}


        for org_ph in (['Lamina'], ['Sheath'], ['Internode'], ['Peduncle', 'Chaff']):
            for variable_name, variable_label in graph_variables_ph_elements.iteritems():
                graph_name = variable_name + '_' + '_'.join(org_ph) + '.PNG'
                cnwheat_tools.plot_cnwheat_ouputs(ph_elements_output_df,
                              x_name = x_name,
                              y_name = variable_name,
                              x_label=x_label,
                              y_label=variable_label,
                              filters={'organ': org_ph},
                              plot_filepath=os.path.join(GRAPHS_DIRPATH, graph_name),
                              explicit_label=False)

        # 2) Roots, grains and phloem
        organs_output_df = pd.read_csv(organs_outputs_filepath)

        graph_variables_organs = {'Conc_Sucrose':u'[Sucrose] (µmol g$^{-1}$ mstruct)', 'Dry_Mass':'Dry mass (g)',
                            'Conc_Nitrates': u'[Nitrates] (µmol g$^{-1}$ mstruct)', 'Conc_Amino_Acids':u'[Amino Acids] (µmol g$^{-1}$ mstruct)', 'Proteins_N_Mass': u'[N Proteins] (g)',
                            'Uptake_Nitrates':u'Nitrates uptake (µmol h$^{-1}$)', 'Unloading_Sucrose':u'Unloaded sucrose (µmol C g$^{-1}$ mstruct h$^{-1}$)', 'Unloading_Amino_Acids':u'Unloaded Amino Acids (µmol N AA g$^{-1}$ mstruct h$^{-1}$)',
                            'S_Amino_Acids': u'Rate of amino acids synthesis (µmol N g$^{-1}$ mstruct h$^{-1}$)', 'S_Proteins': u'Rate of protein synthesis (µmol N h$^{-1}$)', 'Export_Nitrates': u'Total export of nitrates (µmol N h$^{-1}$)', 'Export_Amino_Acids': u'Total export of Amino acids (µmol N h$^{-1}$)',
                            'R_Nnit_upt': u'Respiration nitrates uptake (µmol C h$^{-1}$)', 'R_Nnit_red': u'Respiration nitrate reduction (µmol C h$^{-1}$)', 'R_residual': u'Respiration residual (µmol C h$^{-1}$)', 'R_maintenance': u'Respiration residual (µmol C h$^{-1}$)',
                            'R_grain_growth_struct': u'Respiration grain structural growth (µmol C h$^{-1}$)', 'R_grain_growth_starch': u'Respiration grain starch growth (µmol C h$^{-1}$)',
                            'R_growth': u'Growth respiration of roots (µmol C h$^{-1}$)', 'Nstruct_N_growth': u'Growth of Nstruct (µmol N h$^{-1}$)', 'mstruct_C_growth': u'Growth of structural dry mass (µmol C g$^{-1}$ mstruct h$^{-1}$)', 'mstruct_death': u'Death of root structural dry mass (g)', 'mstruct': u'Structural mass (g)',
                            'C_exudation': u'Carbon lost by root exudation (µmol C g$^{-1}$ mstruct h$^{-1}$', 'N_exudation': u'Nitrogen lost by root exudation (µmol N g$^{-1}$ mstruct h$^{-1}$',
                            'Conc_cytokinines':u'[cytokinines] (UA g$^{-1}$ mstruct)', 'S_cytokinines':u'Rate of cytokinines synthesis (UA g$^{-1}$ mstruct s$^{-1}$)', 'Export_cytokinines': 'Export of cytokinines from roots (UA h$^{-1}$)',
                            'HATS_LATS': u'Potential uptake (µmol h$^{-1}$)' , 'regul_transpiration':'Regulating transpiration function'}

        for org in (['Roots'], ['Grains'], ['Phloem']):
            for variable_name, variable_label in graph_variables_organs.iteritems():
                graph_name = variable_name + '_' + '_'.join(org) + '.PNG'
                cnwheat_tools.plot_cnwheat_ouputs(organs_output_df,
                              x_name = x_name,
                              y_name = variable_name,
                              x_label=x_label,
                              y_label=variable_label,
                              filters={'organ': org},
                              plot_filepath=os.path.join(GRAPHS_DIRPATH, graph_name),
                              explicit_label=False)

        # 3) Soil
        graph_name = 'Conc_Nitrates_Soil.PNG'
        cnwheat_tools.plot_cnwheat_ouputs(organs_output_df,
                              x_name = x_name,
                              y_name = 'Conc_Nitrates',
                              x_label=x_label,
                              y_label=u'[Nitrates] (µmol m$^{-3}$)',
                              filters={'organ': 'Soil'},
                              plot_filepath=os.path.join(GRAPHS_DIRPATH, graph_name),
                              explicit_label=False)

        ## Data NEMA
        t_NEMA = [0, 408, 648, 864, 1200]
        dry_mass_ph_NEMA = [1.92, 2.26, 1.97, 1.69, 1.56] # Photosynthetic organs i.e. laminae + stems + chaff (g)
        dry_mass_grains_NEMA_H3 = [0.12, 0.53, 1.22, 1.46, 1.54] # Grains (g)
        dry_mass_grains_NEMA_H15 = [0.12, 0.58, 1.29, 1.41, 1.50] # Grains (g)

        green_area_lamina1_H3 = [34.57712963, 36.62500501, 36.94493193, 18.09]
        green_area_lamina2_H3 = [33.96263889, 36.53136532, 28.85906228, 18.1]
        green_area_lamina3_H3 = [22.75143519, 23.95125661, 15.75888889]
        green_area_lamina1_H15 = [34.57712963, 34.66222222, 35.11382576, 22.26702485]
        green_area_lamina2_H15 = [33.96263889, 33.51185185, 31.36950758, 21.36666667]
        green_area_lamina3_H15 = [22.75143519, 23.01336941, 18.76296839]

        N_mass_ph_NEMA_H3 = [29.40, 26.80, 17.25, 8.84, 6.98] # Photosynthetic organs i.e. laminae + stems + chaff (mg)
        N_mass_grains_NEMA_H3 = [2.4, 11.10306667, 25.81933333, 32.171, 36.33444] # # Grains (mg)
        N_mass_ph_NEMA_H15 = [29.40, 29.69, 21.58, 11.35, 8.49] # Photosynthetic organs i.e. laminae + stems + chaff (mg)
        N_mass_grains_NEMA_H15 = [2.4, 10.01116667, 22.2441, 30.8332, 32.83882] # Grains (mg)


        # 4) Total dry mass accumulation
        fig, (ax1) = plt.subplots(1)

        TRIOSESP_MOLAR_MASS_C_RATIO = 0.21
        SUCROSE_MOLAR_MASS_C_RATIO = 0.42
        HEXOSE_MOLAR_MASS_C_RATIO = 0.4
        NITRATES_MOLAR_MASS_N_RATIO = 0.23
        AMINO_ACIDS_MOLAR_MASS_N_RATIO = 0.145


        ## Photosynthetic organs
        org_ph = ph_elements_output_df.groupby('t').sum()
        sum_dry_mass_org_ph =  ((org_ph['triosesP'] * 1E-6*cnwheat_parameters.OrganParameters.C_MOLAR_MASS)/TRIOSESP_MOLAR_MASS_C_RATIO +
                                (org_ph['sucrose'] * 1E-6*cnwheat_parameters.OrganParameters.C_MOLAR_MASS)/SUCROSE_MOLAR_MASS_C_RATIO +
                                (org_ph['starch'] * 1E-6*cnwheat_parameters.OrganParameters.C_MOLAR_MASS)/HEXOSE_MOLAR_MASS_C_RATIO +
                                (org_ph['fructan'] * 1E-6*cnwheat_parameters.OrganParameters.C_MOLAR_MASS)/HEXOSE_MOLAR_MASS_C_RATIO +
                                (org_ph['nitrates'] * 1E-6*cnwheat_parameters.OrganParameters.N_MOLAR_MASS)/NITRATES_MOLAR_MASS_N_RATIO +
                                (org_ph['amino_acids'] * 1E-6*cnwheat_parameters.OrganParameters.N_MOLAR_MASS)/AMINO_ACIDS_MOLAR_MASS_N_RATIO +
                                (org_ph['proteins'] * 1E-6*cnwheat_parameters.OrganParameters.N_MOLAR_MASS)/AMINO_ACIDS_MOLAR_MASS_N_RATIO +
                                [org_ph['mstruct'][0]]*len(org_ph.index))

        ax1.plot(sum_dry_mass_org_ph.index, sum_dry_mass_org_ph, label = r'$\sum$ (tp,i)', linestyle = '-', color = 'g')
        ax1.plot(t_NEMA, dry_mass_ph_NEMA, marker = 'o', color = 'g', linestyle = '')

        ## Roots
        roots = organs_output_df[organs_output_df['organ']=='Roots'].groupby('t').sum()
        sum_dry_mass_roots =    ((roots['sucrose'] * 1E-6*cnwheat_parameters.OrganParameters.C_MOLAR_MASS)/SUCROSE_MOLAR_MASS_C_RATIO +
                                (roots['nitrates'] * 1E-6*cnwheat_parameters.OrganParameters.N_MOLAR_MASS)/NITRATES_MOLAR_MASS_N_RATIO +
                                (roots['amino_acids'] * 1E-6*cnwheat_parameters.OrganParameters.N_MOLAR_MASS)/AMINO_ACIDS_MOLAR_MASS_N_RATIO +
                                roots['mstruct'])

        ax1.plot(sum_dry_mass_roots.index, sum_dry_mass_roots, label = 'Roots', linestyle = '-', color = 'k')

        ## Grains
        grains = organs_output_df[organs_output_df['organ']=='Grains'].groupby('t').sum()
        sum_dry_mass_grains = grains['Dry_Mass']

        ax1.plot(sum_dry_mass_grains.index, sum_dry_mass_grains, label = 'Grains', linestyle = '-', color = 'y')
        ax1.plot(t_NEMA, dry_mass_grains_NEMA_H3, label= 'H3', marker = 'v', color = 'y', linestyle = '')
        ax1.plot(t_NEMA, dry_mass_grains_NEMA_H15, label= 'H15', marker = 's', color = 'y', linestyle = '')

        ## Phloem
        phloem = organs_output_df[organs_output_df['organ']=='Phloem'].groupby('t').sum()
        sum_dry_mass_phloem =    ((phloem['sucrose'] * 1E-6*cnwheat_parameters.OrganParameters.C_MOLAR_MASS)/SUCROSE_MOLAR_MASS_C_RATIO +
                                 (phloem['amino_acids'] * 1E-6*cnwheat_parameters.OrganParameters.N_MOLAR_MASS)/AMINO_ACIDS_MOLAR_MASS_N_RATIO)

        ax1.plot(sum_dry_mass_phloem.index, sum_dry_mass_phloem, label = 'Phloem', linestyle = '-', color = 'b')
        ax1.plot(sum_dry_mass_org_ph.index, (sum_dry_mass_org_ph + sum_dry_mass_phloem), label = r'$\sum$ (tp,i) + phloem', linestyle = '--', color = 'g')

        ## Total
        total_dry_mass = sum_dry_mass_org_ph + sum_dry_mass_roots + sum_dry_mass_grains + sum_dry_mass_phloem
        ax1.plot(total_dry_mass.index, total_dry_mass, label = 'Total', linestyle = '-', color = 'r')

        ## Formatting
        ax1.set_ylabel('Dry mass (g)')
        ax1.legend(prop={'size':12}, bbox_to_anchor=(0.05, .6, 0.9, .5), loc='upper center', ncol=4, mode="expand", borderaxespad=0.)
        ax1.axvline(cnwheat_parameters.GrainsParameters.FILLING_INIT, color='k', linestyle='--')
        ax1.set_xlabel('Time from flowering (hour)')
        plt.savefig(os.path.join(GRAPHS_DIRPATH, 'Total_dry_mass.PNG'), format='PNG', bbox_inches='tight')
        plt.close()

        ## Dry mass production
        fig, (ax1) = plt.subplots(1)

        day, dry_mass_production = [], []
        prec_dry_mass = total_dry_mass[0]

        for t, dry_mass in total_dry_mass.iteritems():
            if t!=0 and t%23==0:
                day.append(int(t//24))
                dry_mass_production.append(dry_mass-prec_dry_mass)
                prec_dry_mass = dry_mass


        ax1.plot(day, dry_mass_production, linestyle = '-')
        ax1.set_ylabel('Dry mass production (g day$^{-1}$)')
        ax1.axvline(cnwheat_parameters.GrainsParameters.FILLING_INIT//24, color='k', linestyle='--')
        ax1.set_xlabel('Time from flowering (day)')
        plt.savefig(os.path.join(GRAPHS_DIRPATH, 'Dry_mass_production.PNG'), format='PNG', bbox_inches='tight')
        plt.close()

        # 5) Total N accumulation
        fig, (ax1) = plt.subplots(1)

        ## Photosynthetic organs
        org_ph = ph_elements_output_df.groupby('t').sum()
        sum_N_org_ph =  ((org_ph['nitrates'] * 1E-3*cnwheat_parameters.OrganParameters.N_MOLAR_MASS) +
                         (org_ph['amino_acids'] * 1E-3*cnwheat_parameters.OrganParameters.N_MOLAR_MASS) +
                         (org_ph['proteins'] * 1E-3*cnwheat_parameters.OrganParameters.N_MOLAR_MASS) +
                         [org_ph['Nstruct'][0]*1E3]*len(org_ph.index))

        ax1.plot(sum_N_org_ph.index, sum_N_org_ph, label = r'$\sum$ (tp,i)', linestyle = '-', color = 'g')
        ax1.plot(t_NEMA, N_mass_ph_NEMA_H3, marker = 'v', color = 'g', linestyle = '')
        ax1.plot(t_NEMA, N_mass_ph_NEMA_H15, marker = 's', color = 'g', linestyle = '')

        ## Roots
        roots = organs_output_df[organs_output_df['organ']=='Roots'].groupby('t').sum()
        sum_N_roots =    ((roots['nitrates'] * 1E-3*cnwheat_parameters.OrganParameters.N_MOLAR_MASS) +
                          (roots['amino_acids'] * 1E-3*cnwheat_parameters.OrganParameters.N_MOLAR_MASS) +
                          (roots['Nstruct'] * 1E3))

        ax1.plot(sum_N_roots.index, sum_N_roots, label = 'Roots', linestyle = '-', color = 'k')

        ## Grains
        grains = organs_output_df[organs_output_df['organ']=='Grains'].groupby('t').sum()
        sum_N_grains = grains['proteins'] * 1E-3*cnwheat_parameters.OrganParameters.N_MOLAR_MASS

        ax1.plot(sum_N_grains.index, sum_N_grains, label = 'Grains', linestyle = '-', color = 'y')
        ax1.plot(t_NEMA, N_mass_grains_NEMA_H3, label= 'H3', marker = 'v', color = 'y', linestyle = '')
        ax1.plot(t_NEMA, N_mass_grains_NEMA_H15, label= 'H15', marker = 's', color = 'y', linestyle = '')

        ## Phloem
        phloem = organs_output_df[organs_output_df['organ']=='Phloem'].groupby('t').sum()
        sum_N_phloem = phloem['amino_acids'] * 1E-3*cnwheat_parameters.OrganParameters.N_MOLAR_MASS

        ax1.plot(sum_N_phloem.index, sum_N_phloem, label = 'Phloem', linestyle = '-', color = 'b')

        ## Total
        total_N_mass = sum_N_org_ph + sum_N_roots + sum_N_grains + sum_N_phloem
        ax1.plot(total_N_mass.index, total_N_mass, label = 'Total', linestyle = '-', color = 'r')

        ## Formatting
        ax1.set_ylabel('N mass (mg)')
        ax1.legend(prop={'size':12}, bbox_to_anchor=(0.05, .6, 0.9, .5), loc='upper center', ncol=4, mode="expand", borderaxespad=0.)
        ax1.axvline(cnwheat_parameters.GrainsParameters.FILLING_INIT, color='k', linestyle='--')
        ax1.set_xlabel('Time from flowering (hour)')
        plt.savefig(os.path.join(GRAPHS_DIRPATH, 'Total_N_mass.PNG'), format='PNG', bbox_inches='tight')
        plt.close()

        # 6) Respiration plots
        ph_elements_output_df['day'] = ph_elements_output_df['t']//24
        organs_output_df['day'] = organs_output_df['t']//24
        days = ph_elements_output_df['day'].unique()
        green_area = ph_elements_output_df[ph_elements_output_df['organ']=='Lamina'].groupby('day')['green_area'].sum()/24

        R_phloem_loading = ph_elements_output_df.groupby('day')['R_phloem_loading'].sum()*12E-9/green_area
        R_Nnit_red = ph_elements_output_df.groupby('day')['R_Nnit_red'].sum()*12E-9/green_area + organs_output_df.groupby('day')['R_Nnit_red'].sum()*12E-9/green_area
        R_residual = ph_elements_output_df.groupby('day')['R_residual'].sum()*12E-9/green_area + organs_output_df.groupby('day')['R_residual'].sum()*12E-9/green_area
        R_Nnit_upt = organs_output_df.groupby('day')['R_Nnit_upt'].sum()*12E-9/green_area
        R_grain_growth_struct =  organs_output_df.groupby('day')['R_grain_growth_struct'].sum()*12E-9/green_area
        R_grain_growth_starch =  organs_output_df.groupby('day')['R_grain_growth_starch'].sum()*12E-9/green_area
        fig, ax1 = plt.subplots(1, sharex=True)
        ax1.plot(days, R_phloem_loading, label='R_phloem_loading')
        ax1.plot(days, R_Nnit_red, label='R_Nnit_red')
        ax1.plot(days, R_residual, label='R_residual')
        ax1.plot(days, R_Nnit_upt, label='R_Nnit_upt')
        ax1.plot(days, R_grain_growth_struct, label='R_grain_growth_struct')
        ax1.plot(days, R_grain_growth_starch, label='R_grain_growth_starch')

        ## Formatting
        ax1.legend(prop={'size':10}, framealpha=0.5)
        ax1.set_ylabel(u'Total tiller respiration (kg C m$^{-2}$ d$^{-1}$)')
        ax1.set_xlabel('Day')
        plt.tight_layout()
        plt.savefig(os.path.join(GRAPHS_DIRPATH, 'Respiration_total.PNG'), dpi=200, format='PNG')
        plt.close()

        ## 2nd plot
        R_phloem_loading = ph_elements_output_df.groupby('day')['R_phloem_loading'].mean()
        R_Nnit_red = ph_elements_output_df.groupby('day')['R_Nnit_red'].mean() + organs_output_df.groupby('day')['R_Nnit_red'].mean()
        R_residual = ph_elements_output_df.groupby('day')['R_residual'].mean() + organs_output_df.groupby('day')['R_residual'].mean()
        R_Nnit_upt = organs_output_df.groupby('day')['R_Nnit_upt'].mean()
        R_grain_growth_struct =  organs_output_df.groupby('day')['R_grain_growth_struct'].mean()
        R_grain_growth_starch =  organs_output_df.groupby('day')['R_grain_growth_starch'].mean()
        fig, ax1 = plt.subplots(1, sharex=True)
        ax1.plot(days, R_phloem_loading, label='R_phloem_loading')
        ax1.plot(days, R_Nnit_red, label='R_Nnit_red')
        ax1.plot(days, R_residual, label='R_residual')
        ax1.plot(days, R_Nnit_upt, label='R_Nnit_upt')
        ax1.plot(days, R_grain_growth_struct, label='R_grain_growth_struct')
        ax1.plot(days, R_grain_growth_starch, label='R_grain_growth_starch')

        ## Formatting
        ax1.legend(prop={'size':10}, framealpha=0.5)
        ax1.set_ylabel(u'Mean hourly tiller respiration (µmol C h$^{-1}$)')
        ax1.set_xlabel('Day')
        plt.tight_layout()
        plt.savefig(os.path.join(GRAPHS_DIRPATH, 'Respiration_total2.PNG'), dpi=200, format='PNG')
        plt.close()

        # 7) Daily N uptake plot
        organs_output_df['day'] = organs_output_df['t']//24+1
        fig, ax1 = plt.subplots(1, sharex=True)
        days = organs_output_df['day'].unique()
        daily_N_uptk = organs_output_df.groupby('day')['Uptake_Nitrates'].aggregate(np.sum)
        ax1.plot(days, daily_N_uptk)

        ## Formatting
        ax1.legend(prop={'size':10}, framealpha=0.5)
        ax1.set_ylabel(u'Daily nitrates uptake (µmol N day$^{-1}$)')
        ax1.set_xlabel('Day')
        plt.tight_layout()
        plt.savefig(os.path.join(GRAPHS_DIRPATH, 'Uptake_Nitrates_Roots_daily.PNG'), dpi=200, format='PNG')
        plt.close()

        # 7) PAR interception
        PARi_df_t = PARi_df.loc[PARi_df.t.isin(ph_elements_output_df.t), :]
        meteo_df_t = meteo_df.loc[ph_elements_output_df.t, :]
        ph_elements_output_df['day'] = ph_elements_output_df['t']//24
        days = ph_elements_output_df['day'].unique()
        ph_elements_output_df['integrated_PAR'] = PARi_df_t['PAR'].values * ph_elements_output_df['green_area'].values * 3600
        PAR_absorbed_cumul = ph_elements_output_df.groupby('day')['integrated_PAR'].sum()*410
        meteo_df_t.loc[:, 'day'] = meteo_df_t.index // 24
        PAR_incident_cumul = meteo_df_t.groupby('day')['PAR_incident'].sum()
        ratio_PAR_absorbed = PAR_absorbed_cumul/PAR_incident_cumul
        fig, (ax1) = plt.subplots(1)
        ax1.plot(meteo_df_t['day'].unique(), ratio_PAR_absorbed, label = 'Absorbed PAR', linestyle = '-', color = 'k', marker='o')

        ## Formatting
        ax1.set_ylabel(u'Fraction of incident PAR absorbed by the tiller')
        ax1.set_xlabel('Day')
        plt.tight_layout()
        plt.savefig(os.path.join(GRAPHS_DIRPATH, 'PAR_absorption.PNG'), dpi=200, format='PNG')
        plt.close()

        # Lamina green area
        fig, (ax1) = plt.subplots(1)

        green_area_lamina1_model = ph_elements_output_df[(ph_elements_output_df['organ']=='Lamina') & (ph_elements_output_df['phytomer']==1)]
        green_area_lamina2_model = ph_elements_output_df[(ph_elements_output_df['organ']=='Lamina') & (ph_elements_output_df['phytomer']==2)]
        green_area_lamina3_model = ph_elements_output_df[(ph_elements_output_df['organ']=='Lamina') & (ph_elements_output_df['phytomer']==3)]

        ax1.plot(green_area_lamina1_model['t'], green_area_lamina1_model['green_area']*10000, label = 'Lamina 1', linestyle = '-', color = 'b')
        ax1.plot(green_area_lamina2_model['t'], green_area_lamina2_model['green_area']*10000, label = 'Lamina 2', linestyle = '-', color = 'g')
        ax1.plot(green_area_lamina3_model['t'], green_area_lamina3_model['green_area']*10000, label = 'Lamina 3', linestyle = '-', color = 'r')

        # NEMA
        ax1.plot(t_NEMA[:-1], green_area_lamina1_H3, label = 'La1 H3', marker = 'v', color = 'b', linestyle = '')
        ax1.plot(t_NEMA[:-1], green_area_lamina2_H3, label = 'La2 H3', marker = 'v', color = 'g', linestyle = '')
        ax1.plot(t_NEMA[:-2], green_area_lamina3_H3, label = 'La3 H3', marker = 'v', color = 'r', linestyle = '')

        ax1.plot(t_NEMA[:-1], green_area_lamina1_H15, label = 'La1 H15', marker = 's', color = 'b', linestyle = '')
        ax1.plot(t_NEMA[:-1], green_area_lamina2_H15, label = 'La2 H15', marker = 's', color = 'g', linestyle = '')
        ax1.plot(t_NEMA[:-2], green_area_lamina3_H15, label = 'La3 H15', marker = 's', color = 'r', linestyle = '')

        ax1.set_ylabel(u'Photosynthetic area (cm$^{2}$)')
        ax1.set_xlabel('Time from flowering (hour)')
        ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        ax1.axvline(360, color='k', linestyle='--')
        # Shrink current axis by 20%
        box = ax1.get_position()
        ax1.set_position([box.x0, box.y0, box.width * 0.85, box.height])
        plt.savefig(os.path.join(GRAPHS_DIRPATH, 'green_area_Lamina_comparison.PNG'), dpi=200, format='PNG')
        plt.close()
        
        
if __name__ == '__main__':
    compute_CN_distrib(run_simu=True, make_graphs=True)
#     compute_CN_distrib(run_simu=False, make_graphs=True)

##    # Profiling
##    filename = 'profile.pstats'
##    profile.run('compute_CN_distrib(make_graphs=True)', filename)
##    stats = pstats.Stats(filename)
##    stats.sort_stats('time')
