# -*- coding: latin-1 -*-

'''
    main
    ~~~~

    Use of the modelling framework for the simulations presented in the paper Barillot et al. (2016).

    You must first install :mod:`alinea.adel`, :mod:`cnwheat`, :mod:`farquharwheat` and :mod:`senescwheat` (and add them to your PYTHONPATH)
    before running this script with the command `python`.

    :copyright: Copyright 2014-2015 INRA-ECOSYS, see AUTHORS.
    :license: TODO, see LICENSE for details.

    .. seealso:: Barillot et al. 2016.
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

from alinea.adel import astk_interface

from cnwheat import simulation as cnwheat_simulation, model as cnwheat_model, parameters as cnwheat_parameters, converter as cnwheat_converter, tools as cnwheat_tools
from farquharwheat import simulation as farquharwheat_simulation, model as farquharwheat_model, converter as farquharwheat_converter
from senescwheat import simulation as senescwheat_simulation, model as senescwheat_model, converter as senescwheat_converter



def force_MTG_data(g, df):
    """
    Force MTG data from dataframe
    """

    groups_df = df.groupby(['plant', 'axis', 'metamer', 'organ', 'element'])

    for vid in g.components_at_scale(g.root, scale=5):
        pid = int(g.index(g.complex_at_scale(vid, scale =1)))
        axid = g.property('label')[g.complex_at_scale(vid, scale =2)]
        mid = int(g.index(g.complex_at_scale(vid, scale =3)))
        org = g.property('label')[g.complex_at_scale(vid, scale =4)]
        elid = g.property('label')[vid]
        id_map = (pid, axid, mid, org, elid)
        if groups_df.groups.has_key(id_map):
            g.property('green_area')[vid] = (groups_df.get_group(id_map)['green_area'].iloc[0])*10000
            g.property('area')[vid] = (groups_df.get_group(id_map)['green_area'].iloc[0])*10000
        else:
            pass
    return g


INPUTS_DIRPATH = 'inputs'
GRAPHS_DIRPATH = 'graphs'
SCREENSHOT_DIRPATH = 'screenshots'

# adelwheat inputs at t0
ADELWHEAT_INPUTS_DIRPATH = os.path.join(INPUTS_DIRPATH, 'adelwheat') # the directory adelwheat must contain files 'adel0000.pckl' and 'scene0000.bgeom'

# cnwheat inputs at t0
CNWHEAT_INPUTS_DIRPATH = os.path.join(INPUTS_DIRPATH, 'cnwheat')
CNWHEAT_ORGANS_INPUTS_FILEPATH = os.path.join(CNWHEAT_INPUTS_DIRPATH, 'organs_inputs.csv')
CNWHEAT_ELEMENTS_INPUTS_FILEPATH = os.path.join(CNWHEAT_INPUTS_DIRPATH, 'elements_inputs.csv')
CNWHEAT_SOILS_INPUTS_FILEPATH = os.path.join(CNWHEAT_INPUTS_DIRPATH, 'soils_inputs.csv')

# farquharwheat inputs at t0
FARQUHARWHEAT_INPUTS_DIRPATH = os.path.join(INPUTS_DIRPATH, 'farquharwheat')
FARQUHARWHEAT_INPUTS_FILEPATH = os.path.join(FARQUHARWHEAT_INPUTS_DIRPATH, 'inputs.csv')
METEO_FILEPATH = os.path.join(FARQUHARWHEAT_INPUTS_DIRPATH, 'meteo_Clermont_rebuild.csv')

# senescwheat inputs at t0
SENESCWHEAT_INPUTS_DIRPATH = os.path.join(INPUTS_DIRPATH, 'senescwheat')
SENESCWHEAT_ROOTS_INPUTS_FILEPATH = os.path.join(SENESCWHEAT_INPUTS_DIRPATH, 'roots_inputs.csv')
SENESCWHEAT_ELEMENTS_INPUTS_FILEPATH = os.path.join(SENESCWHEAT_INPUTS_DIRPATH, 'elements_inputs.csv')

# inputs and outputs of all models for each step
OUTPUTS_DIRPATH = 'outputs'
AXES_INPUTS_OUTPUTS_FILEPATH = os.path.join(OUTPUTS_DIRPATH, 'axes_inputs_outputs.csv')
ORGANS_INPUTS_OUTPUTS_FILEPATH = os.path.join(OUTPUTS_DIRPATH, 'organs_inputs_outputs.csv')
ELEMENTS_INPUTS_OUTPUTS_FILEPATH = os.path.join(OUTPUTS_DIRPATH, 'elements_inputs_outputs.csv')
SOILS_INPUTS_OUTPUTS_FILEPATH = os.path.join(OUTPUTS_DIRPATH, 'soils_inputs_outputs.csv')

INPUTS_OUTPUTS_PRECISION = 6

LOGGING_CONFIG_FILEPATH = os.path.join('..', 'logging.json')

LOGGING_LEVEL = logging.WARNING # can be one of: DEBUG, INFO, WARNING, ERROR, CRITICAL

cnwheat_tools.setup_logging(LOGGING_CONFIG_FILEPATH, LOGGING_LEVEL, log_model=False, log_compartments=False, log_derivatives=False)


meteo_df = pd.read_csv(METEO_FILEPATH, index_col='t')

current_time_of_the_system = time.time()

def compute_CN_distrib(run_simu=True, make_graphs=True):

    # define the time step in hours for each simulator
    senescwheat_ts = 2
    farquharwheat_ts = 2
    cnwheat_ts = 1

    hour_to_second_conversion_factor = 3600
    # create the simulators
    senescwheat_simulation_ = senescwheat_simulation.Simulation(senescwheat_ts * hour_to_second_conversion_factor)
    farquharwheat_simulation_ = farquharwheat_simulation.Simulation()
    cnwheat_simulation_ = cnwheat_simulation.Simulation(cnwheat_ts * hour_to_second_conversion_factor)

    # read adelwheat inputs at t0
    adel_wheat = astk_interface.AdelWheat(seed=1234)
    g = adel_wheat.load(dir=ADELWHEAT_INPUTS_DIRPATH)[0]

    # read cnwheat inputs at t0
    cnwheat_organs_inputs_df = pd.read_csv(CNWHEAT_ORGANS_INPUTS_FILEPATH)
    cnwheat_elements_inputs_df = pd.read_csv(CNWHEAT_ELEMENTS_INPUTS_FILEPATH)
    cnwheat_soils_inputs_df = pd.read_csv(CNWHEAT_SOILS_INPUTS_FILEPATH)

    # read farquharwheat inputs at t0
    farquharwheat_inputs_df = pd.read_csv(FARQUHARWHEAT_INPUTS_FILEPATH)

    # read senescwheat inputs at t0
    senescwheat_roots_inputs_df = pd.read_csv(SENESCWHEAT_ROOTS_INPUTS_FILEPATH)
    senescwheat_elements_inputs_df = pd.read_csv(SENESCWHEAT_ELEMENTS_INPUTS_FILEPATH)

    # define the start and the end of the whole simulation (in hours)
    start_time = 0
    stop_time = 1200

    # Redefine MTG data
    g = force_MTG_data(g, cnwheat_elements_inputs_df)

    # define organs for which the variable 'max_proteins' is fixed
    forced_max_protein_elements = set(((1,'MS',9,'blade', 'LeafElement1'), (1,'MS',10,'blade', 'LeafElement1'), (1,'MS',11,'blade', 'LeafElement1')))

    # define lists of dataframes to store the inputs and outputs of cnwheat, farquharwheat and senescwheat at each step, and use them for postprocessings
    axes_inputs_ouputs_df_list = []
    organs_inputs_ouputs_df_list = []
    elements_inputs_ouputs_df_list = []
    soils_inputs_ouputs_df_list = []

    # Initialize dataframes to store the inputs and outputs of cnwheat, farquharwheat and senescwheat at a given step.
    # Since the topology is static, these dataframes keep the same shape during all the simulation.
    organs_inputs_ouputs_df = pd.DataFrame(index=range(len(cnwheat_organs_inputs_df.index)),
                                           columns=cnwheat_simulation.Simulation.ORGANS_OUTPUTS_INDEXES + sorted(set(cnwheat_simulation.Simulation.ORGANS_STATE + cnwheat_simulation.Simulation.ORGANS_INTERMEDIATE_VARIABLES + cnwheat_simulation.Simulation.ORGANS_FLUXES + cnwheat_simulation.Simulation.ORGANS_INTEGRATIVE_VARIABLES + cnwheat_simulation.Simulation.ORGANS_POSTPROCESSING_VARIABLES + senescwheat_converter.SENESCWHEAT_ROOTS_INPUTS_OUTPUTS)))
    organs_inputs_indexes = cnwheat_organs_inputs_df.loc[:, cnwheat_simulation.Simulation.ORGANS_INPUTS_INDEXES]
    organs_inputs_indexes.sort_index(by=cnwheat_simulation.Simulation.ORGANS_INPUTS_INDEXES, inplace=True)
    organs_inputs_indexes.reset_index(drop=True, inplace=True)
    organs_inputs_ouputs_df.loc[:, cnwheat_simulation.Simulation.ORGANS_INPUTS_INDEXES] = organs_inputs_indexes.loc[:, cnwheat_simulation.Simulation.ORGANS_INPUTS_INDEXES]

    elements_inputs_ouputs_df = pd.DataFrame(index=cnwheat_elements_inputs_df.index,
                                             columns=cnwheat_simulation.Simulation.ELEMENTS_OUTPUTS_INDEXES + sorted(set(cnwheat_simulation.Simulation.ELEMENTS_STATE + cnwheat_simulation.Simulation.ELEMENTS_INTERMEDIATE_VARIABLES + cnwheat_simulation.Simulation.ELEMENTS_FLUXES + cnwheat_simulation.Simulation.ELEMENTS_INTEGRATIVE_VARIABLES + cnwheat_simulation.Simulation.ELEMENTS_POSTPROCESSING_VARIABLES + farquharwheat_converter.FARQUHARWHEAT_INPUTS_OUTPUTS + senescwheat_converter.SENESCWHEAT_ELEMENTS_INPUTS_OUTPUTS)))
    elements_inputs_ouputs_df.loc[:, cnwheat_simulation.Simulation.ELEMENTS_INPUTS_INDEXES] = cnwheat_elements_inputs_df.loc[:, cnwheat_simulation.Simulation.ELEMENTS_INPUTS_INDEXES]

    if run_simu:
        # run the simulators
        for t_senescwheat in xrange(start_time, stop_time, senescwheat_ts):
            # initialize and run senescwheat, and update the global mtg
            senescwheat_simulation_.initialize(senescwheat_converter.from_MTG(g, senescwheat_roots_inputs_df, senescwheat_elements_inputs_df))
            senescwheat_simulation_.run()
            senescwheat_converter.update_MTG(senescwheat_simulation_.inputs, senescwheat_simulation_.outputs, g)
            # fill the global dataframes for post-processings
            senescwheat_roots_inputs_df, senescwheat_elements_inputs_df = senescwheat_converter.to_dataframes(senescwheat_simulation_.inputs)
            senescwheat_roots_outputs_df, senescwheat_elements_outputs_df = senescwheat_converter.to_dataframes(senescwheat_simulation_.outputs)
            senescwheat_roots_inputs_outputs_df = senescwheat_roots_outputs_df.combine_first(senescwheat_roots_inputs_df)
            senescwheat_elements_inputs_outputs_df = senescwheat_elements_outputs_df.combine_first(senescwheat_elements_inputs_df)
            organs_inputs_ouputs_df.loc[organs_inputs_ouputs_df.organ == 'roots', senescwheat_converter.SENESCWHEAT_ROOTS_INPUTS_OUTPUTS] = senescwheat_roots_inputs_outputs_df.loc[:, senescwheat_converter.SENESCWHEAT_ROOTS_INPUTS_OUTPUTS].values
            elements_inputs_ouputs_df.loc[:, senescwheat_converter.SENESCWHEAT_ELEMENTS_INPUTS_OUTPUTS] = senescwheat_elements_inputs_outputs_df.loc[:, senescwheat_converter.SENESCWHEAT_ELEMENTS_INPUTS_OUTPUTS].values

            for t_farquharwheat in xrange(t_senescwheat, t_senescwheat + senescwheat_ts, farquharwheat_ts):
                # get the meteo of the current step
                Ta, ambient_CO2, RH, Ur, PARi = meteo_df.loc[t_farquharwheat, ['air_temperature', 'ambient_CO2', 'humidity', 'Wind', 'PARi']]
                # initialize and run farquharwheat, and update the global mtg
                farquharwheat_simulation_.initialize(farquharwheat_converter.from_MTG(g, farquharwheat_inputs_df))
                farquharwheat_simulation_.run(Ta, ambient_CO2, RH, Ur, PARi)
                farquharwheat_converter.update_MTG(farquharwheat_simulation_.inputs, farquharwheat_simulation_.outputs, g)
                # fill the global dataframes for post-processings
                farquharwheat_inputs_df = farquharwheat_converter.to_dataframe(farquharwheat_simulation_.inputs)
                farquharwheat_outputs_df = farquharwheat_converter.to_dataframe(farquharwheat_simulation_.outputs)
                farquharwheat_inputs_outputs_df = farquharwheat_outputs_df.combine_first(farquharwheat_inputs_df)
                elements_inputs_ouputs_df.loc[:, farquharwheat_converter.FARQUHARWHEAT_INPUTS_OUTPUTS] = farquharwheat_inputs_outputs_df.loc[:, farquharwheat_converter.FARQUHARWHEAT_INPUTS_OUTPUTS].values

                for t_cnwheat in xrange(t_farquharwheat, t_farquharwheat + farquharwheat_ts, cnwheat_ts):
                    # initialize and run cnwheat, and update the global mtg
                    cnwheat_simulation_.initialize(cnwheat_converter.from_MTG(g, cnwheat_organs_inputs_df, cnwheat_elements_inputs_df),
                                                   cnwheat_converter.from_dataframes(soils_inputs=cnwheat_soils_inputs_df))
                    cnwheat_simulation_.run(start_time=t_cnwheat, stop_time=t_cnwheat+cnwheat_ts, number_of_output_steps=cnwheat_ts+1)
                    cnwheat_converter.update_MTG(cnwheat_simulation_.population, g)
                    cnwheat_soils_inputs_df = cnwheat_converter.to_dataframes(soils=cnwheat_simulation_.soils)
                    # run cnwheat post-processings
                    (_,
                     cnwheat_axes_postprocessing_df,
                     _,
                     cnwheat_organs_postprocessing_df,
                     cnwheat_elements_postprocessing_df,
                     cnwheat_soils_postprocessing_df) = cnwheat_simulation_.postprocessings()
                    # fill the global dataframes for post-processings
                    axes_inputs_ouputs_df = cnwheat_axes_postprocessing_df.loc[cnwheat_axes_postprocessing_df.t == t_cnwheat, :].reset_index(drop=True)
                    cnwheat_organs_postprocessing_df = cnwheat_organs_postprocessing_df.loc[cnwheat_organs_postprocessing_df.t == t_cnwheat, :].reset_index(drop=True)
                    organs_inputs_ouputs_df.loc[organs_inputs_ouputs_df.organ.isin(cnwheat_converter.MTG_TO_CNWHEAT_AXES_ORGANS_MAPPING), cnwheat_simulation.Simulation.ORGANS_ALL_VARIABLES] = cnwheat_organs_postprocessing_df.loc[:, cnwheat_simulation.Simulation.ORGANS_ALL_VARIABLES].values
                    organs_inputs_ouputs_df.loc[:, cnwheat_simulation.Simulation.T_INDEX] = t_cnwheat
                    cnwheat_elements_postprocessing_df = cnwheat_elements_postprocessing_df.loc[cnwheat_elements_postprocessing_df.t == t_cnwheat, :].reset_index(drop=True)
                    elements_inputs_ouputs_df = cnwheat_elements_postprocessing_df.combine_first(elements_inputs_ouputs_df)
                    soils_inputs_ouputs_df = cnwheat_soils_postprocessing_df.loc[cnwheat_soils_postprocessing_df.t == t_cnwheat, :].reset_index(drop=True)

                    axes_inputs_ouputs_df_list.append(axes_inputs_ouputs_df)
                    organs_inputs_ouputs_df_list.append(organs_inputs_ouputs_df)
                    elements_inputs_ouputs_df_list.append(elements_inputs_ouputs_df)
                    soils_inputs_ouputs_df_list.append(soils_inputs_ouputs_df)

                    organs_inputs_ouputs_df = organs_inputs_ouputs_df.copy()
                    elements_inputs_ouputs_df = elements_inputs_ouputs_df.copy()

                    # Screenshots
                    if t_cnwheat%24 == 0:
                        cnwheat_tools.color_MTG_Nitrogen(g, cnwheat_elements_inputs_df, t_cnwheat, SCREENSHOT_DIRPATH)

        global_axes_df = pd.concat(axes_inputs_ouputs_df_list, ignore_index=True)
        axes_inputs_outputs_variable_names = sorted([column for column in global_axes_df.columns if column not in cnwheat_simulation.Simulation.AXES_OUTPUTS_INDEXES])
        global_axes_df = global_axes_df.reindex_axis(cnwheat_simulation.Simulation.AXES_OUTPUTS_INDEXES + axes_inputs_outputs_variable_names, axis=1, copy=False)
        global_axes_df.to_csv(AXES_INPUTS_OUTPUTS_FILEPATH, na_rep='NA', index=False, float_format='%.{}f'.format(INPUTS_OUTPUTS_PRECISION))

        global_organs_df = pd.concat(organs_inputs_ouputs_df_list, ignore_index=True)
        organs_inputs_outputs_variable_names = sorted([column for column in global_organs_df.columns if column not in cnwheat_simulation.Simulation.ORGANS_OUTPUTS_INDEXES])
        global_organs_df = global_organs_df.reindex_axis(cnwheat_simulation.Simulation.ORGANS_OUTPUTS_INDEXES + organs_inputs_outputs_variable_names, axis=1, copy=False)
        global_organs_df.to_csv(ORGANS_INPUTS_OUTPUTS_FILEPATH, na_rep='NA', index=False, float_format='%.{}f'.format(INPUTS_OUTPUTS_PRECISION))

        global_elements_df = pd.concat(elements_inputs_ouputs_df_list, ignore_index=True)
        elements_inputs_outputs_variable_names = [column for column in global_elements_df.columns if column not in cnwheat_simulation.Simulation.ELEMENTS_OUTPUTS_INDEXES]
        global_elements_df = global_elements_df.reindex_axis(cnwheat_simulation.Simulation.ELEMENTS_OUTPUTS_INDEXES + elements_inputs_outputs_variable_names, axis=1, copy=False)
        global_elements_df.to_csv(ELEMENTS_INPUTS_OUTPUTS_FILEPATH, na_rep='NA', index=False, float_format='%.{}f'.format(INPUTS_OUTPUTS_PRECISION))

        global_soils_df = pd.concat(soils_inputs_ouputs_df_list, ignore_index=True)
        soils_inputs_outputs_variable_names = [column for column in global_soils_df.columns if column not in cnwheat_simulation.Simulation.SOILS_OUTPUTS_INDEXES]
        global_soils_df = global_soils_df.reindex_axis(cnwheat_simulation.Simulation.SOILS_OUTPUTS_INDEXES + soils_inputs_outputs_variable_names, axis=1, copy=False)
        global_soils_df.to_csv(SOILS_INPUTS_OUTPUTS_FILEPATH, na_rep='NA', index=False, float_format='%.{}f'.format(INPUTS_OUTPUTS_PRECISION))

        execution_time = int(time.time() - current_time_of_the_system)
        print '\n', 'Model executed in ', str(datetime.timedelta(seconds=execution_time))

    ########POST-PROCESSING##
    if make_graphs:
        from cnwheat import parameters
        from cnwheat import tools
        x_name = 't'
        x_label='Time (Hour)'

        # 1) Photosynthetic organs
        ph_elements_output_df = pd.read_csv(ELEMENTS_INPUTS_OUTPUTS_FILEPATH)

        graph_variables_ph_elements = {'Ag': u'Gross photosynthesis (µmol m$^{-2}$ s$^{-1}$)','An': u'Net photosynthesis (µmol m$^{-2}$ s$^{-1}$)', 'Tr':u'Organ surfacic transpiration rate (mmol H$_{2}$0 m$^{-2}$ s$^{-1}$)', 'Transpiration':u'Organ transpiration rate (mmol H$_{2}$0 s$^{-1}$)', 'Rd': u'Mitochondrial respiration rate of organ in light (µmol C h$^{-1}$)', 'Ts': u'Temperature surface (°C)', 'gs': u'Conductance stomatique (mol m$^{-2}$ s$^{-1}$)',
                           'Conc_TriosesP': u'[TriosesP] (µmol g$^{-1}$ mstruct)', 'Conc_Starch':u'[Starch] (µmol g$^{-1}$ mstruct)', 'Conc_Sucrose':u'[Sucrose] (µmol g$^{-1}$ mstruct)', 'Conc_Fructan':u'[Fructan] (µmol g$^{-1}$ mstruct)',
                           'Conc_Nitrates': u'[Nitrates] (µmol g$^{-1}$ mstruct)', 'Conc_Amino_Acids': u'[Amino_Acids] (µmol g$^{-1}$ mstruct)', 'Conc_Proteins': u'[Proteins] (g g$^{-1}$ mstruct)',
                           'Nitrates_import': u'Total nitrates imported (µmol h$^{-1}$)', 'Amino_Acids_import': u'Total amino acids imported (µmol N h$^{-1}$)',
                           'S_Amino_Acids': u'[Rate of amino acids synthesis] (µmol N g$^{-1}$ mstruct h$^{-1}$)', 'S_Proteins': u'Rate of protein synthesis (µmol N g$^{-1}$ mstruct h$^{-1}$)', 'D_Proteins': u'Rate of protein degradation (µmol N g$^{-1}$ mstruct h$^{-1}$)', 'k_proteins': u'Relative rate of protein degradation (s$^{-1}$)',
                           'Loading_Sucrose': u'Loading Sucrose (µmol C sucrose h$^{-1}$)', 'Loading_Amino_Acids': u'Loading Amino acids (µmol N amino acids h$^{-1}$)',
                           'green_area': u'Green area (m$^{2}$)', 'R_phloem_loading': u'Respiration phloem loading (µmol C h$^{-1}$)', 'R_Nnit_red': u'Respiration nitrate reduction (µmol C h$^{-1}$)', 'R_residual': u'Respiration residual (µmol C h$^{-1}$)', 'R_maintenance': u'Respiration residual (µmol C h$^{-1}$)',
                           'mstruct': u'Structural mass (g)', 'Nstruct': u'Structural N mass (g)',
                           'Conc_cytokinins':u'[cytokinins] (UA g$^{-1}$ mstruct)', 'D_cytokinins':u'Cytokinin degradation (UA g$^{-1}$ mstruct)', 'cytokinins_import':u'Cytokinin import (UA)'}


        for org_ph in (['blade'], ['sheath'], ['internode'], ['peduncle', 'ear']):
            for variable_name, variable_label in graph_variables_ph_elements.iteritems():
                graph_name = variable_name + '_' + '_'.join(org_ph) + '.PNG'
                tools.plot_cnwheat_ouputs(ph_elements_output_df,
                              x_name = x_name,
                              y_name = variable_name,
                              x_label=x_label,
                              y_label=variable_label,
                              filters={'organ': org_ph},
                              plot_filepath=os.path.join(GRAPHS_DIRPATH, graph_name),
                              explicit_label=False)

        # 2) Roots, grains and phloem
        organs_output_df = pd.read_csv(ORGANS_INPUTS_OUTPUTS_FILEPATH)

        graph_variables_organs = {'Conc_Sucrose':u'[Sucrose] (µmol g$^{-1}$ mstruct)', 'Dry_Mass':'Dry mass (g)',
                            'Conc_Nitrates': u'[Nitrates] (µmol g$^{-1}$ mstruct)', 'Conc_Amino_Acids':u'[Amino Acids] (µmol g$^{-1}$ mstruct)', 'Proteins_N_Mass': u'[N Proteins] (g)',
                            'Uptake_Nitrates':u'Nitrates uptake (µmol h$^{-1}$)', 'Unloading_Sucrose':u'Unloaded sucrose (µmol C g$^{-1}$ mstruct h$^{-1}$)', 'Unloading_Amino_Acids':u'Unloaded Amino Acids (µmol N AA g$^{-1}$ mstruct h$^{-1}$)',
                            'S_Amino_Acids': u'Rate of amino acids synthesis (µmol N g$^{-1}$ mstruct h$^{-1}$)', 'S_Proteins': u'Rate of protein synthesis (µmol N h$^{-1}$)', 'Export_Nitrates': u'Total export of nitrates (µmol N h$^{-1}$)', 'Export_Amino_Acids': u'Total export of Amino acids (µmol N h$^{-1}$)',
                            'R_Nnit_upt': u'Respiration nitrates uptake (µmol C h$^{-1}$)', 'R_Nnit_red': u'Respiration nitrate reduction (µmol C h$^{-1}$)', 'R_residual': u'Respiration residual (µmol C h$^{-1}$)', 'R_maintenance': u'Respiration residual (µmol C h$^{-1}$)',
                            'R_grain_growth_struct': u'Respiration grain structural growth (µmol C h$^{-1}$)', 'R_grain_growth_starch': u'Respiration grain starch growth (µmol C h$^{-1}$)',
                            'R_growth': u'Growth respiration of roots (µmol C h$^{-1}$)', 'Nstruct_N_growth': u'Growth of Nstruct (µmol N h$^{-1}$)', 'mstruct_C_growth': u'Growth of structural dry mass (µmol C g$^{-1}$ mstruct h$^{-1}$)', 'mstruct_death': u'Death of root structural dry mass (g)', 'mstruct': u'Structural mass (g)',
                            'C_exudation': u'Carbon lost by root exudation (µmol C g$^{-1}$ mstruct h$^{-1}$', 'N_exudation': u'Nitrogen lost by root exudation (µmol N g$^{-1}$ mstruct h$^{-1}$',
                            'Conc_cytokinins':u'[cytokinins] (UA g$^{-1}$ mstruct)', 'S_cytokinins':u'Rate of cytokinins synthesis (UA g$^{-1}$ mstruct)', 'Export_cytokinins': 'Export of cytokinins from roots (UA h$^{-1}$)',
                            'HATS_LATS': u'Potential uptake (µmol h$^{-1}$)' , 'regul_transpiration':'Regulating transpiration function'}

        for org in (['roots'], ['grains'], ['phloem']):
            for variable_name, variable_label in graph_variables_organs.iteritems():
                graph_name = variable_name + '_' + '_'.join(org) + '.PNG'
                tools.plot_cnwheat_ouputs(organs_output_df,
                              x_name = x_name,
                              y_name = variable_name,
                              x_label=x_label,
                              y_label=variable_label,
                              filters={'organ': org},
                              plot_filepath=os.path.join(GRAPHS_DIRPATH, graph_name),
                              explicit_label=False)

        # 3) Soil
        soil_output_df = pd.read_csv(SOILS_INPUTS_OUTPUTS_FILEPATH)

        fig, (ax1) = plt.subplots(1)
        conc_nitrates_soil = soil_output_df['Conc_Nitrates_Soil']*14E-6
        ax1.plot(soil_output_df['t'], conc_nitrates_soil)
        ax1.set_ylabel(u'[Nitrates] (g m$^{-3}$)')
        ax1.set_xlabel('Time from flowering (hour)')
        ax1.set_title = 'Conc Nitrates Soil'
        plt.savefig(os.path.join(GRAPHS_DIRPATH, 'Conc_Nitrates_Soil.PNG'), format='PNG', bbox_inches='tight')
        plt.close()

        TRIOSESP_MOLAR_MASS_C_RATIO = 0.21
        SUCROSE_MOLAR_MASS_C_RATIO = 0.42
        HEXOSE_MOLAR_MASS_C_RATIO = 0.4
        NITRATES_MOLAR_MASS_N_RATIO = 0.23
        AMINO_ACIDS_MOLAR_MASS_N_RATIO = 0.145
        phloem_shoot_root = 0.75

        # Data NEMA
        t_NEMA = [0, 408, 648, 864, 1200]
        dry_mass_ph_NEMA_H0 = [1.95, 2.34, 1.97, 1.73, 1.67] # Photosynthetic organs i.e. laminae + stems + chaff (g)
        dry_mass_ph_NEMA_H0_SD = [0.11, 0.02, 0.08, 0.06, 0.14]
        dry_mass_grains_NEMA_H0 = [0.58, 1.14, 1.54, 1.48] # Grains (g)
        dry_mass_grains_NEMA_H0_SD = [0.08, 0.09, 0.10, 0.04]
        dry_mass_tot_NEMA_H0 = [1.95, 2.92, 3.11, 3.27, 3.15]
        dry_mass_tot_NEMA_H0_SD = [0.11, 0.06, 0.17, 0.16, 0.14]

        green_area_lamina1_H0 = [34.6, 32.6, 27.2, 13.3]
        green_area_lamina1_H0_SD = [3.8, 2.7, 3.7, 0]
        green_area_lamina2_H0 = [34, 32.8, 22.3]
        green_area_lamina2_H0_SD = [2.7, 3.2, 2.1]
        green_area_lamina3_H0 = [22.8, 22.6, 20.5]
        green_area_lamina3_H0_SD = [1.8, 1.3, 0]
        green_area_lamina4_H0 = [16, 17.9]
        green_area_lamina4_H0_SD = [2.1, 5.1]

        N_tot_lamina1_H0 = [6.84, 4.55, 2.29, 0.88, 0.83]
        N_tot_lamina1_H0_SD = [0.81, 0.68, 0.51, 0.12, 0.04]
        N_tot_lamina2_H0 = [4.08, 2.80, 1.29, 0.62, 0.66]
        N_tot_lamina2_H0_SD = [0.34, 0.39, 0.15, 0.03, 0.08]
        N_tot_lamina3_H0 = [1.85, 1.18, 0.40, 0.36, 0.41]
        N_tot_lamina3_H0_SD = [0.15, 0.17, 0.07, 0.06, 0.05]
        N_tot_lamina4_H0 = [0.51, 0.50, 0.40, 0.29, 0.25]
        N_tot_lamina4_H0_SD = [0.22, 0.07, 0.08, 0.03, 0.10]
        N_tot_chaff_H0 = [5.17, 3.29, 1.81, 1.53, 1.75]
        N_tot_chaff_H0_SD = [0.36, 0.26, 0.16, 0.13, 0.73]
        N_tot_stem_H0 = [11.47, 8.78, 6.15, 3.44, 2.95]
        N_tot_stem_H0_SD = [0.98, 0.60, 0.85, 0.26, 0.29]

        N_mass_ph_NEMA_H0 = [29.91, 21.09, 12.34, 7.11, 6.84] # Photosynthetic organs i.e. laminae + stems + chaff (mg)
        N_mass_ph_NEMA_H0_SD = [2.14, 0.87, 1.27, 0.21, 1.17]
        N_mass_grains_NEMA_H0 = [9.15, 16.39, 24.68, 26.12] # Grains (mg)
        N_mass_grains_NEMA_H0_SD = [1.18, 1.58, 1.51, 0.33]
        N_tot_NEMA_H0 = [29.91, 30.24, 28.73, 31.79, 32.97]
        N_tot_NEMA_H0_SD = [2.14, 1.66, 2.63, 1.71, 1.32]

        DM_tot_lamina1_H0 = [0.18, 0.18, 0.17, 0.13, 0.13]
        DM_tot_lamina1_H0_SD = [0.02, 0.01, 0.01, 0.00, 0.00]
        DM_tot_lamina2_H0 = [0.13, 0.13, 0.12, 0.09, 0.08]
        DM_tot_lamina2_H0_SD = [0.01, 0.01, 0.00, 0.00, 0.00]
        DM_tot_lamina3_H0 = [0.08, 0.08, 0.06, 0.06, 0.05]
        DM_tot_lamina3_H0_SD = [0, 0.01, 0.00, 0.00, 0.00]
        DM_tot_lamina4_H0 = [0.03, 0.06, 0.05, 0.04, 0.03]
        DM_tot_lamina4_H0_SD = [0.01, 0.01, 0.01, 0.00, 0.01]
        DM_tot_stem_H0 = [1.27, 1.59, 1.30, 1.11, 1.02]
        DM_tot_stem_H0_SD = [0.07, 0.06, 0.04, 0.04, 0.01]
        DM_tot_chaff_H0 = [0.26, 0.31, 0.27, 0.29, 0.36]
        DM_tot_chaff_H0_SD = [0.02, 0.03, 0.03, 0.02, 0.12]

        # 4) Total dry mass accumulation
        fig, (ax1) = plt.subplots(1)

        ## Photosynthetic organs
        org_ph = ph_elements_output_df.groupby('t').sum()
        sum_dry_mass_org_ph =  ((org_ph['triosesP'] * 1E-6*parameters.OrganParameters.C_MOLAR_MASS)/TRIOSESP_MOLAR_MASS_C_RATIO +
                                (org_ph['sucrose'] * 1E-6*parameters.OrganParameters.C_MOLAR_MASS)/SUCROSE_MOLAR_MASS_C_RATIO +
                                (org_ph['starch'] * 1E-6*parameters.OrganParameters.C_MOLAR_MASS)/HEXOSE_MOLAR_MASS_C_RATIO +
                                (org_ph['fructan'] * 1E-6*parameters.OrganParameters.C_MOLAR_MASS)/HEXOSE_MOLAR_MASS_C_RATIO +
                                (org_ph['nitrates'] * 1E-6*parameters.OrganParameters.N_MOLAR_MASS)/NITRATES_MOLAR_MASS_N_RATIO +
                                (org_ph['amino_acids'] * 1E-6*parameters.OrganParameters.N_MOLAR_MASS)/AMINO_ACIDS_MOLAR_MASS_N_RATIO +
                                (org_ph['proteins'] * 1E-6*parameters.OrganParameters.N_MOLAR_MASS)/AMINO_ACIDS_MOLAR_MASS_N_RATIO +
                                [org_ph['mstruct'][0]]*len(org_ph.index))

        ax1.plot(sum_dry_mass_org_ph.index, sum_dry_mass_org_ph, label = r'$\sum$ (tp,i)', linestyle = '-', color = 'g')
        ax1.errorbar(t_NEMA, dry_mass_ph_NEMA_H0, yerr=dry_mass_ph_NEMA_H0_SD, marker = 'o', color = 'g', linestyle = '')

        ## Roots
        roots = organs_output_df[organs_output_df['organ']=='roots'].groupby('t').sum()
        sum_dry_mass_roots =    ((roots['sucrose'] * 1E-6*parameters.OrganParameters.C_MOLAR_MASS)/SUCROSE_MOLAR_MASS_C_RATIO +
                                (roots['nitrates'] * 1E-6*parameters.OrganParameters.N_MOLAR_MASS)/NITRATES_MOLAR_MASS_N_RATIO +
                                (roots['amino_acids'] * 1E-6*parameters.OrganParameters.N_MOLAR_MASS)/AMINO_ACIDS_MOLAR_MASS_N_RATIO +
                                roots['mstruct'])

        ax1.plot(sum_dry_mass_roots.index, sum_dry_mass_roots, label = 'Roots', linestyle = '-', color = 'k')

        ## Grains
        grains = organs_output_df[organs_output_df['organ']=='grains'].groupby('t').sum()
        sum_dry_mass_grains = grains['Dry_Mass']

        ax1.plot(sum_dry_mass_grains.index, sum_dry_mass_grains, label = 'Grains', linestyle = '-', color = 'y')
        ax1.errorbar(t_NEMA[1:], dry_mass_grains_NEMA_H0, yerr=dry_mass_grains_NEMA_H0_SD, label= 'H0', marker = 's', color = 'y', linestyle = '')

        ## Phloem
        phloem = organs_output_df[organs_output_df['organ']=='phloem'].groupby('t').sum()
        sum_dry_mass_phloem =    ((phloem['sucrose'] * 1E-6*parameters.OrganParameters.C_MOLAR_MASS)/SUCROSE_MOLAR_MASS_C_RATIO +
                                 (phloem['amino_acids'] * 1E-6*parameters.OrganParameters.N_MOLAR_MASS)/AMINO_ACIDS_MOLAR_MASS_N_RATIO)

        ax1.plot(sum_dry_mass_phloem.index, sum_dry_mass_phloem, label = 'Phloem', linestyle = '-', color = 'b')
        sum_dry_mass_ph_phloem = sum_dry_mass_org_ph.add(sum_dry_mass_phloem*phloem_shoot_root, fill_value=0)
        ax1.plot(sum_dry_mass_ph_phloem.index, sum_dry_mass_ph_phloem, label = r'$\sum$ (tp,i) + phloem', linestyle = '--', color = 'g')
        sum_dry_mass_roots_phloem = sum_dry_mass_roots.add(sum_dry_mass_phloem*(1-phloem_shoot_root), fill_value=0)
        ax1.plot(sum_dry_mass_roots_phloem.index, sum_dry_mass_roots_phloem, label = r'$\sum$ roots + phloem', linestyle = '--', color = 'k')

        ## Total aerial
        total_dry_mass = sum_dry_mass_org_ph + sum_dry_mass_grains + sum_dry_mass_phloem*phloem_shoot_root
        ax1.plot(total_dry_mass.index, total_dry_mass, label = 'Total aerial', linestyle = '-', color = 'r')
        ax1.errorbar(t_NEMA[1:], dry_mass_tot_NEMA_H0[1:], yerr=dry_mass_tot_NEMA_H0_SD[1:], label= 'H0', marker = 's', color = 'r', linestyle = '')

        ## Formatting
        ax1.set_ylabel('Dry mass (g)')
        ax1.legend(prop={'size':12}, bbox_to_anchor=(0.05, .6, 0.9, .5), loc='upper center', ncol=4, mode="expand", borderaxespad=0.)
        ax1.axvline(parameters.GrainsParameters.FILLING_INIT, color='k', linestyle='--')
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
        ax1.axvline(parameters.GrainsParameters.FILLING_INIT//24, color='k', linestyle='--')
        ax1.set_xlabel('Time from flowering (day)')
        plt.savefig(os.path.join(GRAPHS_DIRPATH, 'Dry_mass_production.PNG'), format='PNG', bbox_inches='tight')
        plt.close()

        # 5) Total N accumulation
        fig, (ax1) = plt.subplots(1)

        ## Photosynthetic organs
        org_ph = ph_elements_output_df.groupby('t').sum()
        sum_N_org_ph =  ((org_ph['nitrates'] * 1E-3*parameters.OrganParameters.N_MOLAR_MASS) +
                         (org_ph['amino_acids'] * 1E-3*parameters.OrganParameters.N_MOLAR_MASS) +
                         (org_ph['proteins'] * 1E-3*parameters.OrganParameters.N_MOLAR_MASS) +
                         [org_ph['Nstruct'][0]*1E3]*len(org_ph.index))

        ax1.plot(sum_N_org_ph.index, sum_N_org_ph, label = r'$\sum$ (tp,i)', linestyle = '-', color = 'g')
        ax1.errorbar(t_NEMA, N_mass_ph_NEMA_H0, yerr=N_mass_ph_NEMA_H0_SD,  marker = 's', color = 'g', linestyle = '')

        ## Roots
        roots = organs_output_df[organs_output_df['organ']=='roots'].groupby('t').sum()
        sum_N_roots =    ((roots['nitrates'] * 1E-3*parameters.OrganParameters.N_MOLAR_MASS) +
                          (roots['amino_acids'] * 1E-3*parameters.OrganParameters.N_MOLAR_MASS) +
                          (roots['Nstruct'] * 1E3))

        ax1.plot(sum_N_roots.index, sum_N_roots, label = 'Roots', linestyle = '-', color = 'k')

        ## Grains
        grains = organs_output_df[organs_output_df['organ']=='grains'].groupby('t').sum()
        sum_N_grains = grains['proteins'] * 1E-3*parameters.OrganParameters.N_MOLAR_MASS

        ax1.plot(sum_N_grains.index, sum_N_grains, label = 'Grains', linestyle = '-', color = 'y')
        ax1.errorbar(t_NEMA[1:], N_mass_grains_NEMA_H0, yerr=N_mass_grains_NEMA_H0_SD, label= 'H0', marker = 's', color = 'y', linestyle = '')

        ## Phloem
        phloem = organs_output_df[organs_output_df['organ']=='phloem'].groupby('t').sum()
        sum_N_phloem = phloem['amino_acids'] * 1E-3*parameters.OrganParameters.N_MOLAR_MASS
        sum_N_ph_phloem = sum_N_org_ph.add(sum_N_phloem*phloem_shoot_root, fill_value=0)
        sum_N_roots_phloem = sum_N_roots.add(sum_N_phloem*(1-phloem_shoot_root), fill_value=0)
        ax1.plot(sum_N_ph_phloem.index, sum_N_ph_phloem, label = r'$\sum$ (tp,i) + phloem', linestyle = '--', color = 'g')
        ax1.plot(sum_N_phloem.index, sum_N_phloem, label = 'Phloem', linestyle = '-', color = 'b')
        ax1.plot(sum_N_roots_phloem.index, sum_N_roots_phloem, label = r'$\sum$ roots + phloem', linestyle = '--', color = 'k')

        ## Total aerial
        total_N_mass = sum_N_org_ph + sum_N_grains + sum_N_phloem*phloem_shoot_root
        ax1.plot(total_N_mass.index, total_N_mass, label = 'Total aerial', linestyle = '-', color = 'r')
        ax1.errorbar(t_NEMA[1:], N_tot_NEMA_H0[1:], yerr=N_tot_NEMA_H0_SD[1:], label = 'H0', marker = 's', color = 'r', linestyle = '')

        ## Formatting
        ax1.set_ylabel('N mass (mg)')
        ax1.legend(prop={'size':12}, bbox_to_anchor=(0.05, .6, 0.9, .5), loc='upper center', ncol=4, mode="expand", borderaxespad=0.)
        ax1.axvline(parameters.GrainsParameters.FILLING_INIT, color='k', linestyle='--')
        ax1.set_xlabel('Time from flowering (hour)')
        plt.savefig(os.path.join(GRAPHS_DIRPATH, 'Total_N_mass.PNG'), format='PNG', bbox_inches='tight')
        plt.close()

        # 6) Respiration plots
        ph_elements_output_df['day'] = ph_elements_output_df['t']//24
        organs_output_df['day'] = organs_output_df['t']//24
        green_area = ph_elements_output_df[ph_elements_output_df['organ']=='blade'].groupby('day')['green_area'].sum()/24

        R_phloem_loading = ph_elements_output_df.groupby('day')['R_phloem_loading'].sum()*12E-9/green_area
        R_Nnit_red_ph = ph_elements_output_df.groupby('day')['R_Nnit_red'].sum()*12E-9/green_area
        R_Nnit_red_org = organs_output_df.groupby('day')['R_Nnit_red'].sum()*12E-9/green_area
        R_Nnit_red = R_Nnit_red_ph.add(R_Nnit_red_org, fill_value=0)
        R_residual_ph = ph_elements_output_df.groupby('day')['R_residual'].sum()*12E-9/green_area
        R_residual_org = organs_output_df.groupby('day')['R_residual'].sum()*12E-9/green_area
        R_residual = R_residual_ph.add(R_residual_org, fill_value=0)
        R_Nnit_upt = organs_output_df.groupby('day')['R_Nnit_upt'].sum()*12E-9/green_area
        R_grain_growth_struct =  organs_output_df.groupby('day')['R_grain_growth_struct'].sum()*12E-9/green_area
        R_grain_growth_starch =  organs_output_df.groupby('day')['R_grain_growth_starch'].sum()*12E-9/green_area
        fig, ax1 = plt.subplots(1, sharex=True)
        ax1.plot(R_phloem_loading.index, R_phloem_loading, label='R_phloem_loading')
        ax1.plot(R_Nnit_red.index, R_Nnit_red, label='R_Nnit_red')
        ax1.plot(R_residual.index, R_residual, label='R_residual')
        ax1.plot(R_Nnit_upt.index, R_Nnit_upt, label='R_Nnit_upt')
        ax1.plot(R_grain_growth_struct.index, R_grain_growth_struct, label='R_grain_growth_struct')
        ax1.plot(R_grain_growth_starch.index, R_grain_growth_starch, label='R_grain_growth_starch')

        ## Formatting
        ax1.legend(prop={'size':10}, framealpha=0.5)
        ax1.set_ylabel(u'Total tiller respiration (kg C m$^{-2}$ d$^{-1}$)')
        ax1.set_xlabel('Day')
        plt.tight_layout()
        plt.savefig(os.path.join(GRAPHS_DIRPATH, 'Respiration_total.PNG'), dpi=200, format='PNG')
        plt.close()

        ## 2nd plot
        R_phloem_loading = ph_elements_output_df.groupby('day')['R_phloem_loading'].mean()
        R_Nnit_red_ph = ph_elements_output_df.groupby('day')['R_Nnit_red'].mean()
        R_Nnit_red_org = organs_output_df.groupby('day')['R_Nnit_red'].mean()
        R_Nnit_red = R_Nnit_red_ph.add(R_Nnit_red_org, fill_value=0)
        R_residual_ph = ph_elements_output_df.groupby('day')['R_residual'].mean()
        R_residual_org = organs_output_df.groupby('day')['R_residual'].mean()
        R_residual = R_residual_ph.add(R_residual_org, fill_value=0)
        R_Nnit_upt = organs_output_df.groupby('day')['R_Nnit_upt'].mean()
        R_grain_growth_struct =  organs_output_df.groupby('day')['R_grain_growth_struct'].mean()
        R_grain_growth_starch =  organs_output_df.groupby('day')['R_grain_growth_starch'].mean()
        fig, ax1 = plt.subplots(1, sharex=True)
        ax1.plot(R_phloem_loading.index, R_phloem_loading, label='R_phloem_loading')
        ax1.plot(R_Nnit_red.index, R_Nnit_red, label='R_Nnit_red')
        ax1.plot(R_residual.index, R_residual, label='R_residual')
        ax1.plot(R_Nnit_upt.index, R_Nnit_upt, label='R_Nnit_upt')
        ax1.plot(R_grain_growth_struct.index, R_grain_growth_struct, label='R_grain_growth_struct')
        ax1.plot(R_grain_growth_starch.index, R_grain_growth_starch, label='R_grain_growth_starch')

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

        # 8) PAR interception
        ph_elements_output_df['day'] = ph_elements_output_df.t//24
        ph_elements_output_df['projected_area'] = ph_elements_output_df['Eabsm2'] * ph_elements_output_df['green_area'] # m²
        tiller_projected_area = ph_elements_output_df.groupby('day')['projected_area'].sum() / 24 # moyenne jour (m²)
        projected_area_square_meter = tiller_projected_area * parameters.SoilParameters.CULM_DENSITY # moyenne jour sur 1 m² sol (m²)

        days = ph_elements_output_df['day'].unique()
        fig, (ax1) = plt.subplots(1)
        ax1.plot(days, projected_area_square_meter, label = 'Absorbed PAR', linestyle = '-', color = 'k', marker='o')

        ## Formatting
        ax1.set_ylabel(u'Fraction of incident PAR absorbed by the tiller')
        ax1.set_xlabel('Day')
        plt.tight_layout()
        plt.savefig(os.path.join(GRAPHS_DIRPATH, 'PAR_absorption.PNG'), dpi=200, format='PNG')
        plt.close()

        # 9) Laminae
        lamina1_model = ph_elements_output_df[(ph_elements_output_df['organ']=='blade') & (ph_elements_output_df['metamer']==12)].reset_index()
        lamina2_model = ph_elements_output_df[(ph_elements_output_df['organ']=='blade') & (ph_elements_output_df['metamer']==11)].reset_index()
        lamina3_model = ph_elements_output_df[(ph_elements_output_df['organ']=='blade') & (ph_elements_output_df['metamer']==10)].reset_index()
        lamina4_model = ph_elements_output_df[(ph_elements_output_df['organ']=='blade') & (ph_elements_output_df['metamer']==9)].reset_index()

        # Green area
        fig, (ax1) = plt.subplots(1)
        ax1.plot(lamina1_model['t'], lamina1_model['green_area']*10000, label = 'Lamina 1', linestyle = '-', color = 'b')
        ax1.plot(lamina2_model['t'], lamina2_model['green_area']*10000, label = 'Lamina 2', linestyle = '-', color = 'g')
        ax1.plot(lamina3_model['t'], lamina3_model['green_area']*10000, label = 'Lamina 3', linestyle = '-', color = 'r')
        ax1.plot(lamina4_model['t'], lamina4_model['green_area']*10000, label = 'Lamina 4', linestyle = '-', color = 'c')

        ## NEMA
        ax1.errorbar(t_NEMA[:-1], green_area_lamina1_H0, yerr=green_area_lamina1_H0_SD, label = 'La1 H0', marker = 's', color = 'b', linestyle = '')
        ax1.errorbar(t_NEMA[:-2], green_area_lamina2_H0, yerr=green_area_lamina2_H0_SD, label = 'La2 H0', marker = 's', color = 'g', linestyle = '')
        ax1.errorbar(t_NEMA[:-2], green_area_lamina3_H0, yerr=green_area_lamina3_H0_SD, label = 'La3 H0', marker = 's', color = 'r', linestyle = '')
        ax1.errorbar(t_NEMA[:-3], green_area_lamina4_H0, yerr=green_area_lamina4_H0_SD, label = 'La4 H0', marker = 's', color = 'c', linestyle = '')

        ax1.set_ylabel(u'Photosynthetic area (cm$^{2}$)')
        ax1.set_xlabel('Time from flowering (hour)')
        ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        ax1.axvline(360, color='k', linestyle='--')
        ## Shrink current axis by 20%
        box = ax1.get_position()
        ax1.set_position([box.x0, box.y0, box.width * 0.85, box.height])
        plt.savefig(os.path.join(GRAPHS_DIRPATH, 'green_area_blade.PNG'), dpi=200, format='PNG')
        plt.close()

        # N tot
        phloem_model = organs_output_df[organs_output_df['organ']== 'phloem'].reset_index()
        area_tot_ph = ph_elements_output_df[ph_elements_output_df['t']==0]['green_area'].sum()
        contrib_area_lamina1 = (lamina1_model['green_area'][0]/area_tot_ph)*phloem_shoot_root
        contrib_area_lamina2 = (lamina2_model['green_area'][0]/area_tot_ph)*phloem_shoot_root
        contrib_area_lamina3 = (lamina3_model['green_area'][0]/area_tot_ph)*phloem_shoot_root
        contrib_area_lamina4 = (lamina4_model['green_area'][0]/area_tot_ph)*phloem_shoot_root


        sum_N_lamina1 =  ((lamina1_model['nitrates'] * 1E-3*parameters.OrganParameters.N_MOLAR_MASS) +
                         (lamina1_model['amino_acids'] * 1E-3*parameters.OrganParameters.N_MOLAR_MASS) +
                         (lamina1_model['proteins'] * 1E-3*parameters.OrganParameters.N_MOLAR_MASS))
        sum_N_lamina1 = sum_N_lamina1.add(phloem_model['amino_acids'] * 1E-3*parameters.OrganParameters.N_MOLAR_MASS*contrib_area_lamina1, fill_value=0)
        sum_N_lamina1 = sum_N_lamina1.add([lamina1_model['Nstruct'].iloc[0]*1E3]*len(phloem_model.index), fill_value=0)


        sum_N_lamina2 =  ((lamina2_model['nitrates'] * 1E-3*parameters.OrganParameters.N_MOLAR_MASS) +
                         (lamina2_model['amino_acids'] * 1E-3*parameters.OrganParameters.N_MOLAR_MASS) +
                         (lamina2_model['proteins'] * 1E-3*parameters.OrganParameters.N_MOLAR_MASS))
        sum_N_lamina2 = sum_N_lamina2.add(phloem_model['amino_acids'] * 1E-3*parameters.OrganParameters.N_MOLAR_MASS*contrib_area_lamina2, fill_value=0)
        sum_N_lamina2 = sum_N_lamina2.add([lamina2_model['Nstruct'].iloc[0]*1E3]*len(phloem_model.index), fill_value=0)


        sum_N_lamina3 =  ((lamina3_model['nitrates'] * 1E-3*parameters.OrganParameters.N_MOLAR_MASS) +
                         (lamina3_model['amino_acids'] * 1E-3*parameters.OrganParameters.N_MOLAR_MASS) +
                         (lamina3_model['proteins'] * 1E-3*parameters.OrganParameters.N_MOLAR_MASS))
        sum_N_lamina3 = sum_N_lamina3.add(phloem_model['amino_acids'] * 1E-3*parameters.OrganParameters.N_MOLAR_MASS*contrib_area_lamina3, fill_value=0)
        sum_N_lamina3 = sum_N_lamina3.add([lamina3_model['Nstruct'].iloc[0]*1E3]*len(phloem_model.index), fill_value=0)

        sum_N_lamina4 =  ((lamina4_model['nitrates'] * 1E-3*parameters.OrganParameters.N_MOLAR_MASS) +
                         (lamina4_model['amino_acids'] * 1E-3*parameters.OrganParameters.N_MOLAR_MASS) +
                         (lamina4_model['proteins'] * 1E-3*parameters.OrganParameters.N_MOLAR_MASS))
        sum_N_lamina4 = sum_N_lamina4.add(phloem_model['amino_acids'] * 1E-3*parameters.OrganParameters.N_MOLAR_MASS*contrib_area_lamina4, fill_value=0)
        sum_N_lamina4 = sum_N_lamina4.add([lamina4_model['Nstruct'].iloc[0]*1E3]*len(phloem_model.index), fill_value=0)

        fig, (ax1) = plt.subplots(1)
        ax1.plot(sum_N_lamina1.index, sum_N_lamina1, label = 'Lamina 1', linestyle = '-', color = 'b')
        ax1.plot(sum_N_lamina2.index, sum_N_lamina2, label = 'Lamina 2', linestyle = '-', color = 'g')
        ax1.plot(sum_N_lamina3.index, sum_N_lamina3, label = 'Lamina 3', linestyle = '-', color = 'r')
        ax1.plot(sum_N_lamina4.index, sum_N_lamina4, label = 'Lamina 4', linestyle = '-', color = 'c')

        ## NEMA
        ax1.errorbar(t_NEMA, N_tot_lamina1_H0, yerr=N_tot_lamina1_H0_SD, label = 'La1 H0', marker = 's', color = 'b', linestyle = '')
        ax1.errorbar(t_NEMA, N_tot_lamina2_H0, yerr=N_tot_lamina2_H0_SD, label = 'La2 H0', marker = 's', color = 'g', linestyle = '')
        ax1.errorbar(t_NEMA, N_tot_lamina3_H0, yerr=N_tot_lamina3_H0_SD, label = 'La3 H0', marker = 's', color = 'r', linestyle = '')
        ax1.errorbar(t_NEMA, N_tot_lamina4_H0, yerr=N_tot_lamina4_H0_SD, label = 'La4 H0', marker = 's', color = 'c', linestyle = '')

        ax1.set_ylabel(u'Lamina N mass (mg N)')
        ax1.set_xlabel('Time from flowering (hour)')
        ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        ax1.axvline(360, color='k', linestyle='--')
        ## Shrink current axis by 20%
        box = ax1.get_position()
        ax1.set_position([box.x0, box.y0, box.width * 0.85, box.height])
        plt.savefig(os.path.join(GRAPHS_DIRPATH, 'N_tot_laminae.PNG'), dpi=200, format='PNG')
        plt.close()

        # 10) Stem & chaff
        stem_model = ph_elements_output_df[(ph_elements_output_df['organ']!= 'blade') & (ph_elements_output_df['organ']!= 'ear')].groupby('t').sum().reset_index()
        contrib_area_stem = (stem_model['green_area'][0]/area_tot_ph)*phloem_shoot_root
        chaff_model = ph_elements_output_df[(ph_elements_output_df['organ']== 'ear')].reset_index()
        contrib_area_chaff = (chaff_model['green_area'][0]/area_tot_ph)*phloem_shoot_root

        # N tot
        sum_N_stem =  ((stem_model['nitrates'] * 1E-3*parameters.OrganParameters.N_MOLAR_MASS) +
                         (stem_model['amino_acids'] * 1E-3*parameters.OrganParameters.N_MOLAR_MASS) +
                         (stem_model['proteins'] * 1E-3*parameters.OrganParameters.N_MOLAR_MASS))
        sum_N_stem = sum_N_stem.add(phloem_model['amino_acids'] * 1E-3*parameters.OrganParameters.N_MOLAR_MASS*contrib_area_stem, fill_value=0)
        sum_N_stem = sum_N_stem.add([stem_model['Nstruct'].iloc[0]*1E3]*len(phloem_model.index), fill_value=0)

        sum_N_chaff =  ((chaff_model['nitrates'] * 1E-3*parameters.OrganParameters.N_MOLAR_MASS) +
                         (chaff_model['amino_acids'] * 1E-3*parameters.OrganParameters.N_MOLAR_MASS) +
                         (chaff_model['proteins'] * 1E-3*parameters.OrganParameters.N_MOLAR_MASS))
        sum_N_chaff = sum_N_chaff.add(phloem_model['amino_acids'] * 1E-3*parameters.OrganParameters.N_MOLAR_MASS*contrib_area_chaff, fill_value=0)
        sum_N_chaff = sum_N_chaff.add([chaff_model['Nstruct'].iloc[0]*1E3]*len(phloem_model.index), fill_value=0)

        fig, (ax1) = plt.subplots(1)
        ax1.plot(sum_N_stem.index, sum_N_stem, label = 'Stem', linestyle = '-', color = 'b')
        ax1.plot(sum_N_chaff.index, sum_N_chaff, label = 'Chaff', linestyle = '-', color = 'g')

        ## NEMA
        ax1.errorbar(t_NEMA, N_tot_stem_H0, yerr=N_tot_stem_H0_SD, label = 'Stem H0', marker = 's', color = 'b', linestyle = '')
        ax1.errorbar(t_NEMA, N_tot_chaff_H0, yerr=N_tot_chaff_H0_SD, label = 'Chaff H0', marker = 's', color = 'g', linestyle = '')

        ax1.set_ylabel(u'Stem N mass (mg N)')
        ax1.set_xlabel('Time from flowering (hour)')
        ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        ax1.axvline(360, color='k', linestyle='--')
        ## Shrink current axis by 20%
        box = ax1.get_position()
        ax1.set_position([box.x0, box.y0, box.width * 0.85, box.height])
        plt.savefig(os.path.join(GRAPHS_DIRPATH, 'N_tot_stem.PNG'), dpi=200, format='PNG')
        plt.close()

        # Dry mass
        ## Laminae
        sum_DM_lamina1 = ((lamina1_model['triosesP'] * 1E-6*parameters.OrganParameters.C_MOLAR_MASS)/TRIOSESP_MOLAR_MASS_C_RATIO +
                          (lamina1_model['sucrose'] * 1E-6*parameters.OrganParameters.C_MOLAR_MASS)/SUCROSE_MOLAR_MASS_C_RATIO +
                          (lamina1_model['starch'] * 1E-6*parameters.OrganParameters.C_MOLAR_MASS)/HEXOSE_MOLAR_MASS_C_RATIO +
                          (lamina1_model['fructan'] * 1E-6*parameters.OrganParameters.C_MOLAR_MASS)/HEXOSE_MOLAR_MASS_C_RATIO +
                          (lamina1_model['nitrates'] * 1E-6*parameters.OrganParameters.N_MOLAR_MASS)/NITRATES_MOLAR_MASS_N_RATIO +
                          (lamina1_model['amino_acids'] * 1E-6*parameters.OrganParameters.N_MOLAR_MASS)/AMINO_ACIDS_MOLAR_MASS_N_RATIO +
                          (lamina1_model['proteins'] * 1E-6*parameters.OrganParameters.N_MOLAR_MASS)/AMINO_ACIDS_MOLAR_MASS_N_RATIO).reset_index(drop=True)
        sum_DM_lamina1 = sum_DM_lamina1.add(sum_dry_mass_phloem*contrib_area_lamina1, fill_value=0)
        sum_DM_lamina1 = sum_DM_lamina1.add([lamina1_model['mstruct'].iloc[0]]*len(phloem_model.index), fill_value=0)

        sum_DM_lamina2 =  ((lamina2_model['triosesP'] * 1E-6*parameters.OrganParameters.C_MOLAR_MASS)/TRIOSESP_MOLAR_MASS_C_RATIO +
                          (lamina2_model['sucrose'] * 1E-6*parameters.OrganParameters.C_MOLAR_MASS)/SUCROSE_MOLAR_MASS_C_RATIO +
                          (lamina2_model['starch'] * 1E-6*parameters.OrganParameters.C_MOLAR_MASS)/HEXOSE_MOLAR_MASS_C_RATIO +
                          (lamina2_model['fructan'] * 1E-6*parameters.OrganParameters.C_MOLAR_MASS)/HEXOSE_MOLAR_MASS_C_RATIO +
                          (lamina2_model['nitrates'] * 1E-6*parameters.OrganParameters.N_MOLAR_MASS)/NITRATES_MOLAR_MASS_N_RATIO +
                          (lamina2_model['amino_acids'] * 1E-6*parameters.OrganParameters.N_MOLAR_MASS)/AMINO_ACIDS_MOLAR_MASS_N_RATIO +
                          (lamina2_model['proteins'] * 1E-6*parameters.OrganParameters.N_MOLAR_MASS)/AMINO_ACIDS_MOLAR_MASS_N_RATIO).reset_index(drop=True)
        sum_DM_lamina2 = sum_DM_lamina2.add(sum_dry_mass_phloem*contrib_area_lamina2, fill_value=0)
        sum_DM_lamina2 = sum_DM_lamina2.add([lamina2_model['mstruct'].iloc[0]]*len(phloem_model.index), fill_value=0)

        sum_DM_lamina3 =  ((lamina3_model['triosesP'] * 1E-6*parameters.OrganParameters.C_MOLAR_MASS)/TRIOSESP_MOLAR_MASS_C_RATIO +
                          (lamina3_model['sucrose'] * 1E-6*parameters.OrganParameters.C_MOLAR_MASS)/SUCROSE_MOLAR_MASS_C_RATIO +
                          (lamina3_model['starch'] * 1E-6*parameters.OrganParameters.C_MOLAR_MASS)/HEXOSE_MOLAR_MASS_C_RATIO +
                          (lamina3_model['fructan'] * 1E-6*parameters.OrganParameters.C_MOLAR_MASS)/HEXOSE_MOLAR_MASS_C_RATIO +
                          (lamina3_model['nitrates'] * 1E-6*parameters.OrganParameters.N_MOLAR_MASS)/NITRATES_MOLAR_MASS_N_RATIO +
                          (lamina3_model['amino_acids'] * 1E-6*parameters.OrganParameters.N_MOLAR_MASS)/AMINO_ACIDS_MOLAR_MASS_N_RATIO +
                          (lamina3_model['proteins'] * 1E-6*parameters.OrganParameters.N_MOLAR_MASS)/AMINO_ACIDS_MOLAR_MASS_N_RATIO).reset_index(drop=True)
        sum_DM_lamina3 = sum_DM_lamina3.add(sum_dry_mass_phloem*contrib_area_lamina3, fill_value=0)
        sum_DM_lamina3 = sum_DM_lamina3.add([lamina3_model['mstruct'].iloc[0]]*len(phloem_model.index), fill_value=0)

        sum_DM_lamina4 =  ((lamina4_model['triosesP'] * 1E-6*parameters.OrganParameters.C_MOLAR_MASS)/TRIOSESP_MOLAR_MASS_C_RATIO +
                          (lamina4_model['sucrose'] * 1E-6*parameters.OrganParameters.C_MOLAR_MASS)/SUCROSE_MOLAR_MASS_C_RATIO +
                          (lamina4_model['starch'] * 1E-6*parameters.OrganParameters.C_MOLAR_MASS)/HEXOSE_MOLAR_MASS_C_RATIO +
                          (lamina4_model['fructan'] * 1E-6*parameters.OrganParameters.C_MOLAR_MASS)/HEXOSE_MOLAR_MASS_C_RATIO +
                          (lamina4_model['nitrates'] * 1E-6*parameters.OrganParameters.N_MOLAR_MASS)/NITRATES_MOLAR_MASS_N_RATIO +
                          (lamina4_model['amino_acids'] * 1E-6*parameters.OrganParameters.N_MOLAR_MASS)/AMINO_ACIDS_MOLAR_MASS_N_RATIO +
                          (lamina4_model['proteins'] * 1E-6*parameters.OrganParameters.N_MOLAR_MASS)/AMINO_ACIDS_MOLAR_MASS_N_RATIO).reset_index(drop=True)
        sum_DM_lamina4 = sum_DM_lamina4.add(sum_dry_mass_phloem*contrib_area_lamina4, fill_value=0)
        sum_DM_lamina4 = sum_DM_lamina4.add([lamina4_model['mstruct'].iloc[0]]*len(phloem_model.index), fill_value=0)

        fig, (ax1) = plt.subplots(1)
        ax1.plot(sum_DM_lamina1.index, sum_DM_lamina1, label = 'Lamina 1', linestyle = '-', color = 'b')
        ax1.plot(sum_DM_lamina2.index, sum_DM_lamina2, label = 'Lamina 2', linestyle = '-', color = 'g')
        ax1.plot(sum_DM_lamina3.index, sum_DM_lamina3, label = 'Lamina 3', linestyle = '-', color = 'r')
        ax1.plot(sum_DM_lamina4.index, sum_DM_lamina4, label = 'Lamina 4', linestyle = '-', color = 'c')

        ## NEMA
        ax1.errorbar(t_NEMA, DM_tot_lamina1_H0, yerr=DM_tot_lamina1_H0_SD, label = 'La1 H0', marker = 's', color = 'b', linestyle = '')
        ax1.errorbar(t_NEMA, DM_tot_lamina2_H0, yerr=DM_tot_lamina2_H0_SD, label = 'La2 H0', marker = 's', color = 'g', linestyle = '')
        ax1.errorbar(t_NEMA, DM_tot_lamina3_H0, yerr=DM_tot_lamina3_H0_SD, label = 'La3 H0', marker = 's', color = 'r', linestyle = '')
        ax1.errorbar(t_NEMA, DM_tot_lamina4_H0, yerr=DM_tot_lamina4_H0_SD, label = 'La4 H0', marker = 's', color = 'c', linestyle = '')

        ax1.set_ylabel(u'Lamina dry mass mass (g)')
        ax1.set_xlabel('Time from flowering (hour)')
        ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        ax1.axvline(360, color='k', linestyle='--')
        ## Shrink current axis by 20%
        box = ax1.get_position()
        ax1.set_position([box.x0, box.y0, box.width * 0.85, box.height])
        plt.savefig(os.path.join(GRAPHS_DIRPATH, 'Dry_mass_laminae.PNG'), dpi=200, format='PNG')
        plt.close()

        ## Stem & chaff
        stem_model = ph_elements_output_df[(ph_elements_output_df['organ']!= 'blade') & (ph_elements_output_df['organ']!= 'ear')].groupby('t').sum()
        chaff_model = ph_elements_output_df[(ph_elements_output_df['organ']== 'ear')]

        sum_DM_stem =  ((stem_model['triosesP'] * 1E-6*parameters.OrganParameters.C_MOLAR_MASS)/TRIOSESP_MOLAR_MASS_C_RATIO +
                          (stem_model['sucrose'] * 1E-6*parameters.OrganParameters.C_MOLAR_MASS)/SUCROSE_MOLAR_MASS_C_RATIO +
                          (stem_model['starch'] * 1E-6*parameters.OrganParameters.C_MOLAR_MASS)/HEXOSE_MOLAR_MASS_C_RATIO +
                          (stem_model['fructan'] * 1E-6*parameters.OrganParameters.C_MOLAR_MASS)/HEXOSE_MOLAR_MASS_C_RATIO +
                          (stem_model['nitrates'] * 1E-6*parameters.OrganParameters.N_MOLAR_MASS)/NITRATES_MOLAR_MASS_N_RATIO +
                          (stem_model['amino_acids'] * 1E-6*parameters.OrganParameters.N_MOLAR_MASS)/AMINO_ACIDS_MOLAR_MASS_N_RATIO +
                          (stem_model['proteins'] * 1E-6*parameters.OrganParameters.N_MOLAR_MASS)/AMINO_ACIDS_MOLAR_MASS_N_RATIO).reset_index(drop=True)
        sum_DM_stem = sum_DM_stem.add(sum_dry_mass_phloem*contrib_area_stem, fill_value=0)
        sum_DM_stem = sum_DM_stem.add([stem_model['mstruct'].iloc[0]]*len(phloem_model.index), fill_value=0)

        sum_DM_chaff =  ((chaff_model['triosesP'] * 1E-6*parameters.OrganParameters.C_MOLAR_MASS)/TRIOSESP_MOLAR_MASS_C_RATIO +
                          (chaff_model['sucrose'] * 1E-6*parameters.OrganParameters.C_MOLAR_MASS)/SUCROSE_MOLAR_MASS_C_RATIO +
                          (chaff_model['starch'] * 1E-6*parameters.OrganParameters.C_MOLAR_MASS)/HEXOSE_MOLAR_MASS_C_RATIO +
                          (chaff_model['fructan'] * 1E-6*parameters.OrganParameters.C_MOLAR_MASS)/HEXOSE_MOLAR_MASS_C_RATIO +
                          (chaff_model['nitrates'] * 1E-6*parameters.OrganParameters.N_MOLAR_MASS)/NITRATES_MOLAR_MASS_N_RATIO +
                          (chaff_model['amino_acids'] * 1E-6*parameters.OrganParameters.N_MOLAR_MASS)/AMINO_ACIDS_MOLAR_MASS_N_RATIO +
                          (chaff_model['proteins'] * 1E-6*parameters.OrganParameters.N_MOLAR_MASS)/AMINO_ACIDS_MOLAR_MASS_N_RATIO).reset_index(drop=True)
        sum_DM_chaff = sum_DM_chaff.add(sum_dry_mass_phloem*contrib_area_chaff, fill_value=0)
        sum_DM_chaff = sum_DM_chaff.add([chaff_model['mstruct'].iloc[0]]*len(phloem_model.index), fill_value=0)

        fig, (ax1) = plt.subplots(1)
        ax1.plot(sum_DM_stem.index, sum_DM_stem, label = 'Stem', linestyle = '-', color = 'b')
        ax1.plot(sum_DM_chaff.index, sum_DM_chaff, label = 'Chaff', linestyle = '-', color = 'b')

        ## NEMA
        ax1.errorbar(t_NEMA, DM_tot_stem_H0, yerr=DM_tot_stem_H0_SD, label = 'Stem H0', marker = 's', color = 'b', linestyle = '')
        ax1.errorbar(t_NEMA, DM_tot_chaff_H0, yerr=DM_tot_chaff_H0_SD, label = 'Chaff H0', marker = 'o', color = 'b', linestyle = '')

        ax1.set_ylabel(u'Dry mass (g)')
        ax1.set_xlabel('Time from flowering (hour)')
        ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        ax1.axvline(360, color='k', linestyle='--')
        ## Shrink current axis by 20%
        box = ax1.get_position()
        ax1.set_position([box.x0, box.y0, box.width * 0.85, box.height])
        plt.savefig(os.path.join(GRAPHS_DIRPATH, 'Dry_mass_stem_chaff.PNG'), dpi=200, format='PNG')
        plt.close()


if __name__ == '__main__':
    compute_CN_distrib(run_simu=True, make_graphs=False)
##    # Profiling
##    filename = 'profile.pstats'
##    profile.run('compute_CN_distrib(make_graphs=True)', filename)
##    stats = pstats.Stats(filename)
##    stats.sort_stats('time')
