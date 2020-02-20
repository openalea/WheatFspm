# -*- coding: latin-1 -*-

import pandas as pd
import numpy as np
import os

import statsmodels.api as sm

import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import matplotlib.image as mpimg

from elongwheat import parameters as elongwheat_parameters
import tools

## Import scenarii list and description
scenarii_df = pd.read_csv( os.path.join('inputs','scenarii_list.csv'), index_col='Scenario')
scenarii_df['Scenario'] = scenarii_df.index
if 'Scenario_label' not in scenarii_df.keys():
    scenarii_df['Scenario_label'] =  ''
else:
    scenarii_df['Scenario_label'] = scenarii_df['Scenario_label'].fillna('')
scenarii = scenarii_df.Scenario

## ------- Def functions to make graphs
def graph_RER(scenario):

    scenario_name = 'Scenario_{}'.format(scenario)
    scenario_graphs_dirpath = os.path.join(scenario_name, 'graphs')

    ## ------ RERmax from parameters : from parameters file or from scenarii_list.csv
    scenario_dict = tools.buildDic(scenarii_df.to_dict('index')[scenario])

    rer_param = {}
    if scenario_dict.get('elongwheat'):
        rer_param = scenario_dict['elongwheat'].get('RERmax', {} )
    if rer_param == {}:
        rer_param = dict((k, v) for k, v in elongwheat_parameters.RERmax.iteritems())

    ## ------ Simulated RER

    # import simulation outputs
    data_RER = pd.read_csv(os.path.join(scenario_name, 'outputs', 'hiddenzones_states.csv'))
    data_RER = data_RER[(data_RER.axis == 'MS') & (data_RER.metamer >= 4)].copy()
    data_RER.sort_values(['t','metamer'], inplace=True)
    data_teq = pd.read_csv(os.path.join(scenario_name, 'outputs', 'SAM_states.csv'))
    data_teq = data_teq[data_teq.axis == 'MS'].copy()

    ## Time previous leaf emergence
    tmp = data_RER[data_RER.leaf_is_emerged]
    leaf_em = tmp.groupby('metamer', as_index = False)['t'].min()
    leaf_em['t_em'] = leaf_em.t
    prev_leaf_em = leaf_em
    prev_leaf_em.metamer = leaf_em.metamer + 1

    data_RER2 = pd.merge(data_RER, prev_leaf_em[['metamer','t_em']], on = 'metamer')
    data_RER2= data_RER2[data_RER2.t <= data_RER2.t_em]

    ## SumTimeEq
    data_teq['SumTimeEq'] = np.cumsum(data_teq.delta_teq)
    data_RER3 = pd.merge(data_RER2, data_teq[['t', 'SumTimeEq']], on = 't')

    ## logL
    data_RER3['logL'] = np.log(data_RER3.leaf_L)

    ## Estimate RER
    RER_sim = {}
    for l in data_RER3.metamer.drop_duplicates():
        Y = data_RER3.logL[data_RER3.metamer == l]
        X = data_RER3.SumTimeEq[data_RER3.metamer == l]
        X = sm.add_constant(X)
        mod = sm.OLS(Y, X)
        fit_RER = mod.fit()
        RER_sim[ int(l) ] = fit_RER.params['SumTimeEq']

    ## ------ Graph
    fig, (ax1) = plt.subplots(1, figsize=(4, 3))
    ax1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))

    x, y = zip(*sorted(RER_sim.items()))
    ax1.plot(x,y , label=r'Simulated RER', linestyle='-', color='g')

    ## obs
    # ax1.errorbar(rer_obs.leaf, rer_obs.RER, yerr=rer_obs.RER_confint, marker='o', color='g', linestyle='', label="Observed RER", markersize=2)

    ## Parameters
    ax1.plot(rer_param.keys(), rer_param.values(), marker='*', color='k', linestyle='', label="Model parameters")

    ## Formatting
    ax1.set_ylabel(u'Relative Elongation Rate at 12°C (s$^{-1}$)')
    ax1.legend(prop={'size': 12}, bbox_to_anchor=(0.05, .6, 0.9, .5), loc='upper center', ncol=3, mode="expand", borderaxespad=0.)
    ax1.legend(loc='upper left')
    ax1.set_xlabel('Phytomer rank')
    ax1.set_xlim(left=4)
    ax1.set_ylim(bottom=0., top=6e-6)
    plt.savefig(os.path.join(scenario_graphs_dirpath, 'RER_comparison.PNG'), format='PNG', bbox_inches='tight', dpi=600)
    plt.close()

def graph_C_usages(scenario):

    scenario_name = 'Scenario_{}'.format(scenario)
    scenario_graphs_dirpath = os.path.join(scenario_name, 'graphs')

    # --- Import simulations prostprocessings
    df_elt = pd.read_csv(os.path.join(scenario_name, 'postprocessing', 'elements_postprocessing.csv'))
    df_elt['day'] = df_elt['t'] // 24 + 1
    df_org = pd.read_csv(os.path.join(scenario_name, 'postprocessing', 'organs_postprocessing.csv'))
    df_org['day'] = df_org['t'] // 24 + 1
    df_axe = pd.read_csv(os.path.join(scenario_name, 'postprocessing', 'axes_postprocessing.csv'))
    df_axe['day'] = df_axe['t'] // 24 + 1
    df_hz = pd.read_csv(os.path.join(scenario_name, 'postprocessing', 'hiddenzones_postprocessing.csv'))
    df_hz['day'] = df_hz['t'] // 24 + 1

    # --- Import simulations prostprocessings
    df_elt = pd.read_csv(os.path.join(scenario_name, 'postprocessing', 'elements_postprocessing.csv'))
    df_org = pd.read_csv(os.path.join(scenario_name, 'postprocessing', 'organs_postprocessing.csv'))
    df_axe = pd.read_csv(os.path.join(scenario_name, 'postprocessing', 'axes_postprocessing.csv'))
    df_hz = pd.read_csv(os.path.join(scenario_name, 'postprocessing', 'hiddenzones_postprocessing.csv'))
    df_roots = df_org[df_org['organ'] == 'roots'].copy()
    df_phloem = df_org[df_org['organ'] == 'phloem'].copy()

    # --- C usages relatif to Net Photosynthesis
    AMINO_ACIDS_C_RATIO = 4.15              #: Mean number of mol of C in 1 mol of the major amino acids of plants (Glu, Gln, Ser, Asp, Ala, Gly)
    AMINO_ACIDS_N_RATIO = 1.25              #: Mean number of mol of N in 1 mol of the major amino acids of plants (Glu, Gln, Ser, Asp, Ala, Gly)

    # Photosynthesis
    df_elt['Photosynthesis_tillers'] = df_elt['Photosynthesis'].fillna(0) * df_elt['nb_replications'].fillna(1.)
    Tillers_Photosynthesis_Ag = df_elt.groupby(['t'], as_index=False).agg({'Photosynthesis_tillers': 'sum'})
    C_usages = pd.DataFrame( {'t' : Tillers_Photosynthesis_Ag['t']})
    C_usages['C_produced'] = np.cumsum(Tillers_Photosynthesis_Ag.Photosynthesis_tillers)

    # Respiration
    C_usages['Respi_roots'] = np.cumsum(df_axe.C_respired_roots)
    C_usages['Respi_shoot'] = np.cumsum(df_axe.C_respired_shoot)

    # Exudation
    C_usages['exudation'] = np.cumsum(df_axe.C_exudated.fillna(0))

    # Structural growth
    C_consumption_mstruct_roots = df_roots.sucrose_consumption_mstruct.fillna(0) + df_roots.AA_consumption_mstruct.fillna(0)* AMINO_ACIDS_C_RATIO/AMINO_ACIDS_N_RATIO
    C_usages['Structure_roots'] = np.cumsum(C_consumption_mstruct_roots.reset_index(drop=True))

    df_hz['C_consumption_mstruct'] = df_hz.sucrose_consumption_mstruct.fillna(0) + df_hz.AA_consumption_mstruct.fillna(0) * AMINO_ACIDS_C_RATIO / AMINO_ACIDS_N_RATIO
    df_hz['C_consumption_mstruct_tillers'] = df_hz['C_consumption_mstruct'] * df_hz['nb_replications']
    C_consumption_mstruct_shoot = df_hz.groupby(['t'])['C_consumption_mstruct_tillers'].sum()
    C_usages['Structure_shoot'] = np.cumsum(C_consumption_mstruct_shoot.reset_index(drop=True))

    # Non structural C
    df_phloem['C_NS'] = df_phloem.sucrose.fillna(0) + df_phloem.amino_acids.fillna(0) * AMINO_ACIDS_C_RATIO/AMINO_ACIDS_N_RATIO
    C_NS_phloem_init = df_phloem.C_NS - df_phloem.C_NS[0]
    C_usages['NS_phloem'] = C_NS_phloem_init.reset_index(drop=True)

    df_elt['C_NS'] = df_elt.sucrose.fillna(0) + df_elt.fructan.fillna(0) + df_elt.starch.fillna(0) + (df_elt.amino_acids.fillna(0) + df_elt.proteins.fillna(0))* AMINO_ACIDS_C_RATIO/AMINO_ACIDS_N_RATIO
    df_elt['C_NS_tillers'] = df_elt['C_NS'] * df_elt['nb_replications'].fillna(1.)
    C_elt = df_elt.groupby(['t']).agg({'C_NS_tillers': 'sum'})

    df_hz['C_NS'] = df_hz.sucrose.fillna(0) + df_hz.fructan.fillna(0) + (df_hz.amino_acids.fillna(0) + df_hz.proteins.fillna(0))* AMINO_ACIDS_C_RATIO/AMINO_ACIDS_N_RATIO
    df_hz['C_NS_tillers'] = df_hz['C_NS'] * df_hz['nb_replications'].fillna(1.)
    C_hz = df_hz.groupby(['t']).agg({'C_NS_tillers': 'sum'})

    df_roots['C_NS'] = df_roots.sucrose.fillna(0) + df_roots.amino_acids.fillna(0) * AMINO_ACIDS_C_RATIO/AMINO_ACIDS_N_RATIO

    C_NS_autre = df_roots.C_NS.reset_index(drop=True) + C_elt.C_NS_tillers + C_hz.C_NS_tillers
    C_NS_autre_init = C_NS_autre - C_NS_autre[0]
    C_usages['NS_other'] = C_NS_autre_init.reset_index(drop=True)

    # Total
    C_usages['C_budget'] = (C_usages.Respi_roots + C_usages.Respi_shoot + C_usages.exudation + C_usages.Structure_roots + C_usages.Structure_shoot + C_usages.NS_phloem + C_usages.NS_other) / \
                           C_usages.C_produced

    # ----- Graph
    fig, ax = plt.subplots()
    ax.plot(C_usages.t,  C_usages.Structure_shoot / C_usages.C_produced * 100,
            label=u'Structural mass - Shoot', color='g')
    ax.plot(C_usages.t, C_usages.Structure_roots / C_usages.C_produced * 100,
            label=u'Structural mass - Roots', color='r')
    ax.plot(C_usages.t, (C_usages.NS_phloem + C_usages.NS_other) / C_usages.C_produced * 100, label=u'Non-structural C', color='darkorange')
    ax.plot(C_usages.t, (C_usages.Respi_roots + C_usages.Respi_shoot) / C_usages.C_produced * 100, label=u'C loss by respiration', color='b')
    ax.plot(C_usages.t, C_usages.exudation / C_usages.C_produced * 100, label=u'C loss by exudation', color='c')
    # ax.plot(C_usages.t, C_usages.C_budget * 100, label=u'C consumption vs. production', color='k')

    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax.set_xlabel('Time (h)')
    ax.set_ylabel(u'Carbon usages : Photosynthesis (%)')
    ax.set_ylim(bottom=0, top=100.)

    fig.suptitle(u'Total cumulated usages are ' + str(round(C_usages.C_budget.tail(1)*100,0)) + u' % of Photosynthesis')

    plt.savefig(os.path.join(scenario_graphs_dirpath, 'C_usages_cumulated.PNG'), format='PNG', bbox_inches='tight')
    plt.close()

def graph_summary(scenario, graph_list=None):
    if graph_list is None:
        graph_list = ['LAI', 'sum_dry_mass_axis', 'shoot_roots_ratio_axis', 'N_content_shoot_axis', 'Conc_Amino_acids_phloem', 'Conc_Sucrose_phloem', 'leaf_Lmax',
                      'green_area_blade']
    scenario_name = 'Scenario_{}'.format(scenario)
    scenario_label = scenarii_df['Scenario_label'][scenario]
    scenario_graphs_dirpath = os.path.join(scenario_name, 'graphs')

    nb_graphs = len(graph_list)
    if nb_graphs<=4:
        nrow = 2
        ncol =2
    elif nb_graphs <= 6:
        nrow = 2
        ncol = 3
    elif  nb_graphs <= 8:
        nrow = 2
        ncol = 4
    elif  nb_graphs == 9:
        nrow = 3
        ncol = 3
    else:
        raise AttributeError('Too many graphs to diplay.')

    fig, axs = plt.subplots(nrow,ncol, figsize = (ncol*4,nrow*3) )
    plt.suptitle('Scenario ' + str(scenario) + ' ' + scenario_label , y = 1)
    i = 0
    for (x,y), ax in np.ndenumerate(axs):
        image = mpimg.imread(os.path.join(scenario_graphs_dirpath, graph_list[i] + '.PNG'))
        ax.imshow(image)
        ax.axis('off')
        i += 1
    plt.tight_layout(pad=0.05)
    plt.subplots_adjust(wspace=0, hspace=0)
    plt.subplots_adjust(top=0.95)
    # plt.show()
    plt.savefig( os.path.join(scenario_graphs_dirpath,'Summary.PNG'), format='PNG', bbox_inches='tight', dpi=200)
    plt.close()



if __name__ == '__main__':

    ## ------- Make the graphs and tables for all the scenarii

    for scenario in scenarii:
        graph_RER(int(scenario))
        graph_C_usages(int(scenario))
        graph_summary(int(scenario), graph_list=['LAI','sum_dry_mass_axis','shoot_roots_ratio_axis','N_content_shoot_axis','Conc_Amino_acids_phloem','Conc_Sucrose_phloem', 'leaf_Lmax',
                      'green_area_blade'] )
