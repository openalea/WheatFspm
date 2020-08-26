# -*- coding: latin-1 -*-


import os

import numpy as np
import pandas as pd

from cnwheat import tools
import matplotlib.pyplot as plt

GRAPHS_DIRPATH = 'graphs'

# elongwheat outputs
OUTPUTS_DIRPATH = 'outputs'
HIDDENZONE_OUTPUTS_FILENAME = 'all_hiddenzone_outputs.csv'
ELEMENT_OUTPUTS_FILENAME = 'all_element_outputs.csv'
AXIS_OUTPUTS_FILENAME = 'all_axes_outputs.csv'

# Charts

x_name = 't_step'
x_label = 'Time (h)'

# SAM outputs for SumTT
all_axis_outputs_df = pd.read_csv(os.path.join(OUTPUTS_DIRPATH, AXIS_OUTPUTS_FILENAME))
all_axis_outputs_df = all_axis_outputs_df[all_axis_outputs_df['axis'] == 'MS']

# 4) Hidden zones
all_hiddenzone_outputs_df = pd.read_csv(os.path.join(OUTPUTS_DIRPATH, HIDDENZONE_OUTPUTS_FILENAME))
all_hiddenzone_outputs_df = all_hiddenzone_outputs_df[all_hiddenzone_outputs_df['axis'] == 'MS']
graph_variables_hiddenzones = {'leaf_pseudostem_length': u'Length for leaf emergence (m)', 'leaf_L': u'Leaf length (m)', 'delta_leaf_L': u'Delta leaf length (m)',
                               'internode_distance_to_emerge': u'Length for internode emergence (m)', 'internode_L': u'Internode length (m)', 'delta_internode_L': u'Delta internode length (m)'}

for variable_name, variable_label in graph_variables_hiddenzones.items():
    graph_name = variable_name + '_hz' + '.PNG'
    tools.plot_cnwheat_ouputs(all_hiddenzone_outputs_df,
                              x_name=x_name,
                              y_name=variable_name,
                              x_label=x_label,
                              y_label=variable_label,
                              filters={'plant': 1, 'axis': 'MS'},
                              plot_filepath=os.path.join(GRAPHS_DIRPATH, graph_name),
                              explicit_label=False)

# 5) elements
all_element_outputs_df = pd.read_csv(os.path.join(OUTPUTS_DIRPATH, ELEMENT_OUTPUTS_FILENAME))
all_element_outputs_df = all_element_outputs_df[all_element_outputs_df['axis'] == 'MS']
graph_variables_elements = {'length': u'Length (m)'}

for organ_label in list(all_element_outputs_df['organ'].unique()):
    subdata = all_element_outputs_df[all_element_outputs_df['organ'] == organ_label]

    for element_label in list(all_element_outputs_df['element'].unique()):
        subsubdata = subdata[subdata['element'] == element_label]

        variable_name = 'length'
        if element_label in ('LeafElement1', 'StemElement'):
            variable_label = u'Visible Length (m)'
        else:
            variable_label = u'Hidden Length (m)'

        graph_name = organ_label + '_' + element_label + '.PNG'
        tools.plot_cnwheat_ouputs(subsubdata,
                                  x_name=x_name,
                                  y_name=variable_name,
                                  x_label=x_label,
                                  y_label=variable_label,
                                  filters={'plant': 1, 'axis': 'MS', 'element': element_label},
                                  plot_filepath=os.path.join(GRAPHS_DIRPATH, graph_name),
                                  explicit_label=False)

# --- CALCULATION OUTPUTS
tmp = all_hiddenzone_outputs_df[all_hiddenzone_outputs_df.leaf_is_emerged == 1]

res = tmp['t_step'].groupby(tmp.metamer).min()
res = pd.DataFrame(res)
res['metamer'] = res.index
res = res.merge(all_axis_outputs_df[['t_step', 'sum_TT']], on='t_step', how='left')
res = res.rename(columns={'sum_TT': 'sum_TT_em'})

tmp = res.copy()
tmp.metamer = res.metamer + 1
tmp = tmp.rename(columns={'sum_TT_em': 'sum_TT_em_prec'})
res = res.merge(tmp, on='metamer', how='left')
res['phyllochrone'] = res.sum_TT_em - res.sum_TT_em_prec

res = res[res.phyllochrone > 10]

tmp = all_hiddenzone_outputs_df[['leaf_Lmax', 'internode_Lmax', 'sheath_Lmax', 'lamina_Lmax']].groupby(all_hiddenzone_outputs_df.metamer, as_index=False).min()
tmp['metamer'] = tmp.index
res = res.merge(tmp, on='metamer', how='left')

# Plot Phyllochrone
plt.plot(res.metamer, res.phyllochrone, marker='o')
plt.axis([int(min(res.metamer) - 1), int(max(res.metamer) + 1), 0, 200])
plt.xlabel('Phytomer')
plt.ylabel('Phyllochrone (TT entre emergence feuilles successives)')
plt.title('Phyllochone')
for index, row in res[['metamer', 'phyllochrone']].iterrows():
    plt.text(row['metamer'], row['phyllochrone'], row['phyllochrone'].astype(int).astype(str))
plt.savefig(os.path.join(GRAPHS_DIRPATH, 'phyllo.png'))
plt.close()

# Comparison Ljutovac 2002

bchmk = pd.read_csv('Ljutovac2002.csv')
res = all_hiddenzone_outputs_df
res = res[(res['axis'] == 'MS') & (res['plant'] == 1) & ~np.isnan(res.leaf_Lmax)].copy()
res_IN = res[~ np.isnan(res.internode_Lmax)]
last_value_idx = res.groupby(['metamer'])['t_step'].transform(max) == res['t_step']
res = res[last_value_idx].copy()
bchmk = bchmk[bchmk.metamer >= min(res.metamer)]
last_value_idx = res_IN.groupby(['metamer'])['t_step'].transform(max) == res_IN['t_step']
res_IN = res_IN[last_value_idx].copy()
res = res[['metamer', 'leaf_Lmax', 'lamina_Lmax', 'sheath_Lmax']].merge(res_IN[['metamer', 'internode_Lmax']], left_on='metamer',
                                                                        right_on='metamer', how='outer').copy()

var_list = ['leaf_Lmax', 'lamina_Lmax', 'sheath_Lmax', 'internode_Lmax']
for var in list(var_list):
    fig, ax = plt.subplots()
    plt.xlim((int(min(res.metamer) - 1), int(max(res.metamer) + 1)))
    plt.ylim(ymin=0, ymax=np.nanmax(list(res[var] * 100 * 1.05) + list(bchmk[var] * 1.05)))

    tmp = res[['metamer', var]].drop_duplicates()

    line1 = ax.plot(tmp.metamer, tmp[var] * 100, color='c', marker='o')
    line2 = ax.plot(bchmk.metamer, bchmk[var], color='orange', marker='o')

    ax.set_ylabel(var + ' (cm)')
    ax.set_title(var)
    ax.legend((line1[0], line2[0]), ('Simulation', 'Ljutovac 2002'), loc=2)
    plt.savefig(os.path.join(GRAPHS_DIRPATH, var + '.PNG'))
    plt.close()
