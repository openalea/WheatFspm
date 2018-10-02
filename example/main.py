# -*- coding: latin-1 -*-
"""
    main
    ~~~~

    A standalone script to:

        * run several time Elong-Wheat alone

    You must first install :mod:`elongwheat` and its dependencies
    before running this script with the command `python`.

    :copyright: Copyright 2014-2016 INRA-ECOSYS, see AUTHORS.
    :license: TODO, see LICENSE for details.

    .. seealso:: Barillot et al. 2016.

"""

# --- PREAMBLE

OPTION_SHOW_ADEL = False


import os

import numpy as np
import pandas as pd

import time

from tempwheat import simulation as tempwheat_simulation, model as tempwheat_model, converter as tempwheat_converter
from elongwheat import simulation as elongwheat_simulation, model as elongwheat_model, converter as elongwheat_converter



INPUTS_DIRPATH = 'inputs'
GRAPHS_DIRPATH = 'graphs'

if OPTION_SHOW_ADEL:
   from fspmwheat import elongwheat_facade
   from alinea.adel.adel_dynamic import AdelWheatDyn
   # adelwheat inputs at t0
   ADELWHEAT_INPUTS_DIRPATH = os.path.join(INPUTS_DIRPATH, 'adelwheat') # the directory adelwheat must contain files 'adel0000.pckl' and 'scene0000.bgeom'
   adel_wheat = AdelWheatDyn(seed=1234, convUnit=1)
   g = adel_wheat.load(dir=ADELWHEAT_INPUTS_DIRPATH)
   adel_wheat.domain = g.get_vertex_property(0)['meta']['domain'] # temp (until Christian's commit)
   adel_wheat.nplants = g.get_vertex_property(0)['meta']['nplants'] # temp (until Christian's commit)


# elongwheat inputs
HIDDENZONE_INPUTS_FILEPATH = os.path.join(INPUTS_DIRPATH, 'hiddenzones_inputs.csv')
ELEMENT_INPUTS_FILEPATH = os.path.join(INPUTS_DIRPATH, 'elements_inputs.csv')
SAM_INPUTS_FILEPATH = os.path.join(INPUTS_DIRPATH, 'SAM_inputs.csv')
METEO_FILEPATH = os.path.join(INPUTS_DIRPATH, 'meteo_test.csv')

# elongwheat outputs
OUTPUTS_DIRPATH = 'outputs'
HIDDENZONE_OUTPUTS_FILEPATH = os.path.join(OUTPUTS_DIRPATH, 'all_hiddenzone_outputs.csv')
ELEMENT_OUTPUTS_FILEPATH = os.path.join(OUTPUTS_DIRPATH, 'all_element_outputs.csv')
SAM_OUTPUTS_FILEPATH = os.path.join(OUTPUTS_DIRPATH, 'all_SAM_outputs.csv')

# general output dataframes
all_hiddenzone_outputs_df = pd.DataFrame()
all_element_outputs_df = pd.DataFrame()
all_SAM_outputs_df = pd.DataFrame()

# read elongwheat inputs from a given timestep
run_from_outputs = False

if run_from_outputs:
    hiddenzones_inputs_desired_t = pd.read_csv(HIDDENZONE_OUTPUTS_FILEPATH)
    element_inputs_desired_t = pd.read_csv(ELEMENT_OUTPUTS_FILEPATH)
    SAM_inputs_desired_t = pd.read_csv(SAM_OUTPUTS_FILEPATH)

    assert 't_step' in hiddenzones_inputs_desired_t.columns
    desired_t_step = 487 #max(hiddenzones_inputs_desired_t['t_step'])
    hiddenzone_inputs_df = hiddenzones_inputs_desired_t[hiddenzones_inputs_desired_t.t_step == desired_t_step].drop( ['t_step'], axis = 1)
    element_inputs_df = element_inputs_desired_t[element_inputs_desired_t.t_step == desired_t_step].drop( ['t_step'], axis = 1)
    SAM_inputs_df = SAM_inputs_desired_t[SAM_inputs_desired_t.t_step == desired_t_step].drop( ['t_step'], axis = 1)
    inputs = elongwheat_converter.from_dataframes(hiddenzone_inputs_df, element_inputs_df, SAM_inputs_df)

else:
    hiddenzone_inputs_df = pd.read_csv(HIDDENZONE_INPUTS_FILEPATH)
    element_inputs_df = pd.read_csv(ELEMENT_INPUTS_FILEPATH)
    SAM_inputs_df = pd.read_csv(SAM_INPUTS_FILEPATH)
    inputs = elongwheat_converter.from_dataframes(hiddenzone_inputs_df, element_inputs_df, SAM_inputs_df)
    desired_t_step = 0

meteo = pd.read_csv(METEO_FILEPATH, index_col='t')

## Write first line of dfs with input values
hiddenzone_outputs_df, element_outputs_df, SAM_outputs_df = elongwheat_converter.to_dataframes(inputs)
hiddenzone_outputs_df['t_step'] = desired_t_step
element_outputs_df['t_step'] = desired_t_step
SAM_outputs_df['t_step'] = desired_t_step

## increment the general output dataframes
all_hiddenzone_outputs_df = all_hiddenzone_outputs_df.append(hiddenzone_outputs_df)
all_element_outputs_df = all_element_outputs_df.append(element_outputs_df)
all_SAM_outputs_df = all_SAM_outputs_df.append(SAM_outputs_df)


# define the time step in hours for each elongwheat
elongwheat_ts = 1


# --- SETUP RUN

# setup outup precision
OUTPUTS_PRECISION = 8

# delta_t
delta_t = 3600

# end
loop_end = 800

# --- MAIN

# Initialization simulation
## Create population
tempwheat_simulation_ = tempwheat_simulation.Simulation(delta_t=delta_t)
simulation_ = elongwheat_simulation.Simulation(delta_t=delta_t)

start_time = time.time()


# Loop for several runs
for t_step in range(desired_t_step+1, loop_end+1, elongwheat_ts):

    print(t_step)

    ## -- tempwheat
    Ta, Tsol = meteo.loc[t_step, ['air_temperature', 'soil_temperature']]  # TODO: Add soil temperature in the weather input file
    inputs_tempwheat = tempwheat_converter.from_dataframes( SAM_inputs_df)
    tempwheat_simulation_.initialize(inputs_tempwheat)
    tempwheat_simulation_.run(Ta, Tsol)
    outputs_tempwheat = tempwheat_converter.to_dataframes( tempwheat_simulation_.outputs)

    ## -- elongwheat

    ## convert the dataframe to simulation inputs format
    new_SAM_inputs_df = pd.merge(outputs_tempwheat, SAM_inputs_df[['plant','axis']+[i for i in SAM_inputs_df.columns if i not in outputs_tempwheat.columns]], on =['plant','axis'])
    inputs = elongwheat_converter.from_dataframes(hiddenzone_inputs_df, element_inputs_df, new_SAM_inputs_df)

    ## initialize the simulation with the inputs
    simulation_.initialize(inputs)

    ## run the simulation
    simulation_.run()

    ## convert the outputs to Pandas dataframe
    hiddenzone_outputs_df, element_outputs_df, SAM_outputs_df = elongwheat_converter.to_dataframes(simulation_.outputs)

    new_SAM_outputs_df = pd.merge(SAM_outputs_df, outputs_tempwheat[['plant','axis']+[i for i in outputs_tempwheat.columns if i not in SAM_outputs_df.columns]], on =['plant','axis'])

    ## update MTG
    if OPTION_SHOW_ADEL:
       elongwheat_facade_._update_shared_MTG(simulation_.outputs['hiddenzone'], simulation_.outputs['elements'], simulation_.outputs['SAM'])

    ## use output as input for the next step
    hiddenzone_inputs_df, element_inputs_df, SAM_inputs_df = hiddenzone_outputs_df, element_outputs_df, new_SAM_outputs_df

    ## add column with time step
    hiddenzone_outputs_df['t_step'] = t_step
    element_outputs_df['t_step'] = t_step
    new_SAM_outputs_df['t_step'] = t_step

    ## increment the general output dataframes
    ## increment the general output dataframes, only for current output dataframes which are not empty
    if len(hiddenzone_outputs_df) != 0:
        all_hiddenzone_outputs_df = all_hiddenzone_outputs_df.append(hiddenzone_outputs_df)
    if len(element_outputs_df) != 0:
        all_element_outputs_df = all_element_outputs_df.append(element_outputs_df)
    if len(new_SAM_outputs_df) != 0:
        all_SAM_outputs_df = all_SAM_outputs_df.append(new_SAM_outputs_df)

# --- ADEL
if OPTION_SHOW_ADEL:
   adel_wheat.update_geometry(g)
   adel_wheat.plot(g)

# --- RESUTS

# Save write the dataframe to CSV
all_hiddenzone_outputs_df.to_csv(HIDDENZONE_OUTPUTS_FILEPATH, index=False, na_rep='NA', float_format='%.{}f'.format(OUTPUTS_PRECISION))
all_element_outputs_df.to_csv(ELEMENT_OUTPUTS_FILEPATH, index=False, na_rep='NA', float_format='%.{}f'.format(OUTPUTS_PRECISION))
all_SAM_outputs_df.to_csv(SAM_OUTPUTS_FILEPATH, index=False, na_rep='NA', float_format='%.{}f'.format(OUTPUTS_PRECISION))

print("--- %s seconds ---" % (time.time() - start_time))

# --- GRAPHS
execfile("graphs.py")