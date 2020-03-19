# -*- coding: latin-1 -*-

import pandas as pd
import numpy as np
import os

from cnwheat import model as cnwheat_model

def table_C_usages(scenario_postprocessing_dirpath):
    """ Calculate C usage from postprocessings and save it to a CSV file

    :param str scenario_postprocessing_dirpath: the path to the CSV file describing all scenarii

    """
    # --- Import simulations prostprocessings
    df_axe = pd.read_csv(os.path.join(scenario_postprocessing_dirpath, 'axes_postprocessing.csv'))
    df_elt = pd.read_csv(os.path.join(scenario_postprocessing_dirpath, 'elements_postprocessing.csv'))
    df_org = pd.read_csv(os.path.join(scenario_postprocessing_dirpath, 'organs_postprocessing.csv'))
    df_hz = pd.read_csv(os.path.join(scenario_postprocessing_dirpath, 'hiddenzones_postprocessing.csv'))

    df_roots = df_org[df_org['organ'] == 'roots'].copy()
    df_phloem = df_org[df_org['organ'] == 'phloem'].copy()

    # --- C usages relatif to Net Photosynthesis
    AMINO_ACIDS_C_RATIO = cnwheat_model.EcophysiologicalConstants.AMINO_ACIDS_C_RATIO  #: Mean number of mol of C in 1 mol of the major amino acids of plants (Glu, Gln, Ser, Asp, Ala, Gly)
    AMINO_ACIDS_N_RATIO = cnwheat_model.EcophysiologicalConstants.AMINO_ACIDS_N_RATIO  #: Mean number of mol of N in 1 mol of the major amino acids of plants (Glu, Gln, Ser, Asp, Ala, Gly)

    # Photosynthesis
    df_elt['Photosynthesis_tillers'] = df_elt['Photosynthesis'].fillna(0) * df_elt['nb_replications'].fillna(1.)
    Tillers_Photosynthesis_Ag = df_elt.groupby(['t'], as_index=False).agg({'Photosynthesis_tillers': 'sum'})
    C_usages = pd.DataFrame({'t': Tillers_Photosynthesis_Ag['t']})
    C_usages['C_produced'] = np.cumsum(Tillers_Photosynthesis_Ag.Photosynthesis_tillers)

    # Respiration
    C_usages['Respi_roots'] = np.cumsum(df_axe.C_respired_roots)
    C_usages['Respi_shoot'] = np.cumsum(df_axe.C_respired_shoot)

    # Exudation
    C_usages['exudation'] = np.cumsum(df_axe.C_exudated.fillna(0))

    # Structural growth
    C_consumption_mstruct_roots = df_roots.sucrose_consumption_mstruct.fillna(0) + df_roots.AA_consumption_mstruct.fillna(0) * AMINO_ACIDS_C_RATIO / AMINO_ACIDS_N_RATIO
    C_usages['Structure_roots'] = np.cumsum(C_consumption_mstruct_roots.reset_index(drop=True))

    df_hz['C_consumption_mstruct'] = df_hz.sucrose_consumption_mstruct.fillna(0) + df_hz.AA_consumption_mstruct.fillna(0) * AMINO_ACIDS_C_RATIO / AMINO_ACIDS_N_RATIO
    df_hz['C_consumption_mstruct_tillers'] = df_hz['C_consumption_mstruct'] * df_hz['nb_replications']
    C_consumption_mstruct_shoot = df_hz.groupby(['t'])['C_consumption_mstruct_tillers'].sum()
    C_usages['Structure_shoot'] = np.cumsum(C_consumption_mstruct_shoot.reset_index(drop=True))

    # Non structural C
    df_phloem['C_NS'] = df_phloem.sucrose.fillna(0) + df_phloem.amino_acids.fillna(0) * AMINO_ACIDS_C_RATIO / AMINO_ACIDS_N_RATIO
    C_NS_phloem_init = df_phloem.C_NS - df_phloem.C_NS.reset_index(drop=True)[0]
    C_usages['NS_phloem'] = C_NS_phloem_init.reset_index(drop=True)

    df_elt['C_NS'] = df_elt.sucrose.fillna(0) + df_elt.fructan.fillna(0) + df_elt.starch.fillna(0) + (
            df_elt.amino_acids.fillna(0) + df_elt.proteins.fillna(0)) * AMINO_ACIDS_C_RATIO / AMINO_ACIDS_N_RATIO
    df_elt['C_NS_tillers'] = df_elt['C_NS'] * df_elt['nb_replications'].fillna(1.)
    C_elt = df_elt.groupby(['t']).agg({'C_NS_tillers': 'sum'})

    df_hz['C_NS'] = df_hz.sucrose.fillna(0) + df_hz.fructan.fillna(0) + (df_hz.amino_acids.fillna(0) + df_hz.proteins.fillna(0)) * AMINO_ACIDS_C_RATIO / AMINO_ACIDS_N_RATIO
    df_hz['C_NS_tillers'] = df_hz['C_NS'] * df_hz['nb_replications'].fillna(1.)
    C_hz = df_hz.groupby(['t']).agg({'C_NS_tillers': 'sum'})

    df_roots['C_NS'] = df_roots.sucrose.fillna(0) + df_roots.amino_acids.fillna(0) * AMINO_ACIDS_C_RATIO / AMINO_ACIDS_N_RATIO

    C_NS_autre = df_roots.C_NS.reset_index(drop=True) + C_elt.C_NS_tillers.reset_index(drop=True) + C_hz.C_NS_tillers.reset_index(drop=True)
    C_NS_autre_init = C_NS_autre - C_NS_autre.reset_index(drop=True)[0]
    C_usages['NS_other'] = C_NS_autre_init.reset_index(drop=True)

    # Total
    C_usages['C_budget'] = (C_usages.Respi_roots + C_usages.Respi_shoot + C_usages.exudation + C_usages.Structure_roots + C_usages.Structure_shoot + C_usages.NS_phloem + C_usages.NS_other) /\
                           C_usages.C_produced

    C_usages.to_csv(os.path.join(scenario_postprocessing_dirpath, 'C_usages.csv'), index=False)


def calculate_performance_indices(scenario_postprocessing_dirpath, meteo_dirpath, plant_density):
    """
    Average RUE and photosynthetic yield for the whole cycle.

    :param str scenario_postprocessing_dirpath: the path to the CSV file describing all scenarii
    :param str meteo_dirpath: the path to the CSV meteo file
    :param int plant_density: the plant density (plant m-2)
    """

    # --- Import simulations prostprocessings
    df_elt = pd.read_csv(os.path.join(scenario_postprocessing_dirpath, 'elements_postprocessing.csv'))
    df_axe = pd.read_csv(os.path.join(scenario_postprocessing_dirpath, 'axes_postprocessing.csv'))

    # --- Import meteo file for incident PAR
    df_meteo = pd.read_csv(meteo_dirpath)

    # --- RUE (g DM. MJ-1 PARa)
    df_elt['PARa_MJ'] = df_elt['PARa'] * df_elt['green_area'] * df_elt['nb_replications'].fillna(
        1.) * 3600 / 4.6 * 10 ** -6  # Si tallage, il faut alors utiliser les calculcs green_area et PARa des talles.
    df_elt['RGa_MJ'] = df_elt['PARa'] * df_elt['green_area'] * df_elt['nb_replications'].fillna(
        1.) * 3600 / 2.02 * 10 ** -6  # Si tallage, il faut alors utiliser les calculcs green_area et PARa des talles.
    PARa = df_elt.groupby(['t'])['PARa_MJ'].agg('sum')
    PARa_cum = np.cumsum(PARa)

    RUE_shoot = np.polyfit(PARa_cum, df_axe.sum_dry_mass_shoot, 1)[0]
    RUE_plant = np.polyfit(PARa_cum, df_axe.sum_dry_mass, 1)[0]

    # --- RUE (g DM. MJ-1 RGint estimated from LAI using Beer-Lambert's law with extinction coefficient of 0.4)

    # Beer-Lambert
    df_LAI = df_elt[(df_elt.element == 'LeafElement1')].groupby(['t']).agg({'green_area': 'sum'})
    df_LAI['LAI'] = df_LAI.green_area * plant_density
    df_LAI['t'] = df_LAI.index

    toto = df_meteo[['t', 'PARi_MA4']].merge(df_LAI[['t', 'LAI']], on='t', how='inner')
    toto['PARi_caribu'] = toto.PARi_MA4
    ts_caribu = range(0, toto.shape[0], 4)
    save = toto.at[0, 'PARi_MA4']
    for i in range(0, toto.shape[0]):
        if i in ts_caribu:
            save = toto.at[i, 'PARi_MA4']
        toto.at[i, 'PARi_caribu'] = save

    toto['PARint_BL'] = toto.PARi_caribu * (1 - np.exp(-0.4 * toto.LAI))
    toto['RGint_BL_MJ'] = toto['PARint_BL'] * 3600 / 2.02 * 10 ** -6
    RGint_BL_cum = np.cumsum(toto.RGint_BL_MJ)

    df_axe['sum_dry_mass_shoot_couvert'] = df_axe.sum_dry_mass_shoot * plant_density
    df_axe['sum_dry_mass_couvert'] = df_axe.sum_dry_mass * plant_density

    RUE_shoot_couvert = np.polyfit(RGint_BL_cum, df_axe.sum_dry_mass_shoot_couvert, 1)[0]
    RUE_plant_couvert = np.polyfit(RGint_BL_cum, df_axe.sum_dry_mass_couvert, 1)[0]

    # ---  Photosynthetic efficiency of the plant
    df_elt['Photosynthesis_tillers'] = df_elt.Ag * df_elt.green_area * df_elt.nb_replications.fillna(1.)
    df_elt['PARa_tot_tillers'] = df_elt.PARa * df_elt.green_area * df_elt.nb_replications.fillna(1.)
    # df_elt['green_area_tillers'] = df_elt.green_area * df_elt.nb_replications.fillna(1.)
    # photo_y = df_elt.groupby(['t'],as_index=False).agg({'Photosynthesis_tillers':'sum', 'PARa_tot_tillers':'sum', 'green_area_tillers':'sum'})
    # photo_y['Photosynthetic_yield_plante'] = photo_y.Photosynthesis_tillers / photo_y.PARa_tot_tillers

    PARa2 = df_elt.groupby(['t'])['PARa_tot_tillers'].agg('sum')
    PARa2_cum = np.cumsum(PARa2)
    Photosynthesis = df_elt.groupby(['t'])['Photosynthesis_tillers'].agg('sum')
    Photosynthesis_cum = np.cumsum(Photosynthesis)

    avg_photo_y = np.polyfit(PARa2_cum, Photosynthesis_cum, 1)[0]

    # --- Photosynthetic C allocated to Respiration and to Exudation
    C_usages_path = os.path.join(scenario_postprocessing_dirpath, 'C_usages.csv')
    C_usages = pd.read_csv(C_usages_path)
    C_usages_div = C_usages.div(C_usages.C_produced, axis=0)

    # ---  Write results into a table
    res_df = pd.DataFrame.from_dict({'RUE_plant_MJ_PAR': [RUE_plant],
                                     'RUE_shoot_MJ_PAR': [RUE_shoot],
                                     'RUE_plant_MJ_RGint': [RUE_plant_couvert],
                                     'RUE_shoot_MJ_RGint': [RUE_shoot_couvert],
                                     'Photosynthetic_efficiency': [avg_photo_y],
                                     'C_usages_Respi_roots': C_usages_div.loc[max(C_usages_div.index), 'Respi_roots'],
                                     'C_usages_Respi_shoot': C_usages_div.loc[max(C_usages_div.index), 'Respi_shoot'],
                                     'C_usages_Respi': C_usages_div.loc[max(C_usages_div.index), 'Respi_shoot'] + C_usages_div.loc[max(C_usages_div.index), 'Respi_roots'],
                                     'C_usages_exudation': C_usages_div.loc[max(C_usages_div.index), 'exudation']})

    res_df.to_csv(os.path.join(scenario_postprocessing_dirpath, 'performance_indices.csv'), index=False)


# def all_scenraii_postprocessings(scenarii_list_dirpath):
#     # ------- Run the above functions for all the scenarii
#     # Import scenarii list and description
#     scenarii_df = pd.read_csv(scenarii_list_dirpath, index_col='Scenario')
#     scenarii_df['Scenario'] = scenarii_df.index
#
#     if 'Scenario_label' not in scenarii_df.keys():
#         scenarii_df['Scenario_label'] = ''
#     else:
#         scenarii_df['Scenario_label'] = scenarii_df['Scenario_label'].fillna('')
#     scenarii = scenarii_df.Scenario
#
#
#     for scenario in scenarii:
#         table_C_usages(int(scenario))
#         calculate_performance_indices(int(scenario))
