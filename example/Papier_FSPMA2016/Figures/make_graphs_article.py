# -*- coding: latin-1 -*-

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import stats
from ast import literal_eval

from alinea.adel.astk_interface import AdelWheat

GRAPHS_DIRPATH = r'D:\Documents\PostDoc_Grignon\Modeles\FSPM_Wheat\example\Papier_FSPMA2016_R1\Figures'

C_MOLAR_MASS = 12
N_MOLAR_MASS = 14

TRIOSESP_MOLAR_MASS_C_RATIO = 0.21
SUCROSE_MOLAR_MASS_C_RATIO = 0.42
HEXOSE_MOLAR_MASS_C_RATIO = 0.4
NITRATES_MOLAR_MASS_N_RATIO = 0.23
AMINO_ACIDS_MOLAR_MASS_N_RATIO = 0.145


def cm2inch(*tupl):
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)


fontsize = 10
thickness = 1.25
markersize=6

t_end = 1201
phloem_shoot_root = 0.75

DATA_DIRPATH = r'D:\Documents\PostDoc_Grignon\Modeles\FSPM_Wheat\example\Papier_FSPMA2016_R1'
DENSITIES = [200, 410, 600, 800]
STANDS = ['Plano_diffus', 'Erect_diffus', 'Asso_diffus']

AXES_MAP = {200:0, 410:1, 600:2, 800:3}
markers = {'Soissons': None, 'Caphorn': None, 'Soissons_Mixed': None, 'Caphorn_Mixed': None}
colors = {'Soissons': 'r', 'Caphorn': 'b', 'Soissons_Mixed': 'r', 'Caphorn_Mixed': 'b'}
lines = {'Soissons': '-', 'Caphorn': '-', 'Soissons_Mixed': '--', 'Caphorn_Mixed': '--'}
label_mapping = {'Caphorn_Mixed':"Mixed erectophile", 'Soissons_Mixed':"Mixed planophile", 'Caphorn': "Erectophile", 'Soissons': "Planophile"}

label_mapping = {'Asso_mixte_1':"Mixed erectophile", 'Asso_mixte_2':"Mixed planophile", 'Erect_mixte': "Erectophile", 'Plano_mixte': "Planophile",
                 'Asso_diffus_1':"Mixed erectophile", 'Asso_diffus_2':"Mixed planophile", 'Erect_diffus': "Erectophile", 'Plano_diffus': "Planophile",
                 'Asso_direct_1':"Mixed erectophile", 'Asso_direct_2':"Mixed planophile", 'Erect_direct': "Erectophile", 'Plano_direct': "Planophile"}

color_mapping = {'Asso_diffus_1':"b--", 'Asso_diffus_2':"r--", 'Erect_diffus': "b-", 'Plano_diffus': "r-",}

labelx = -0.15 # axes coords


### Inclination
##inclination_data = r'D:\Documents\PostDoc_Grignon\Modeles\FSPM_Wheat\example\Papier_FSPMA2016\Reconstruction_interception\S2V\Analyse_inclinaison.csv'
##inclination_df = pd.read_csv(inclination_data)
##inclination_df_grouped = inclination_df.groupby(['Cultivar'])
##
##fig = plt.figure(figsize=cm2inch(30,15))
##gs = gridspec.GridSpec(1,3, wspace=0)
##ax1 = plt.subplot(gs[0,0])
##ax2 = plt.subplot(gs[0,1])
##ax3 = plt.subplot(gs[0,2])
##
##markers = {'Soissons': 's', 'Caphorn': '^', 'Mixture': 'o'}
##colors = {'Soissons': 'r', 'Caphorn': 'b', 'Mixture': 'g'}
##labels = {'Soissons': 'Planophile', 'Caphorn': 'Erectophile', 'Mixture': 'Mixture'}
##
##for cv, group in inclination_df_grouped:
##    ax1.plot(group['LAD'], group['Height'], marker = markers[cv], color = colors[cv], linewidth = thickness, label=labels[cv])
##    ax2.plot(group['Leaf inclination'], group['Height'], marker = markers[cv],  color = colors[cv], linewidth = thickness, label=labels[cv])
##    ax3.plot(group['Plant inclination'], group['Height'], marker = markers[cv],  color = colors[cv], linewidth = thickness, label=labels[cv])
##
##ax1.set_xlabel('Leaf Area Density (m$^{2}$ m$^{-3}$)', fontsize=fontsize+5)
##ax1.set_ylabel('Height from soil (m)', fontsize=fontsize+5)
##[i.set_linewidth(thickness) for i in ax1.spines.itervalues()]
##ax1.tick_params(width = thickness, length = 10, labelsize=fontsize+2)
##ax1.set_xticks(np.arange(0, 0.25, 0.05))
##ax1.set_xticks(ax1.get_xticks()[:-1])
##ax1.set_yticks(np.arange(0, 1., 0.2))
##ax1.text(0.1, 0.95, 'A', transform=ax1.transAxes, fontsize=fontsize+5, verticalalignment='top')
##ax1.legend(prop={'size':11}, bbox_to_anchor=(1, 0.95), ncol=1, frameon=False, numpoints = 1)
##
##
##ax2.set_xlabel('Lamina inclination ($^{o}$)', fontsize=fontsize+5)
##[i.set_linewidth(thickness) for i in ax2.spines.itervalues()]
##ax2.tick_params(width = thickness, length = 10, labelsize=fontsize+2, labelleft='off')
##ax2.set_xticks(np.arange(0, 105, 15))
##ax2.set_xticks(ax2.get_xticks()[:-1])
##ax2.set_yticks(np.arange(0, 1., 0.2))
##ax2.text(0.1, 0.95, 'B', transform=ax2.transAxes, fontsize=fontsize+5, verticalalignment='top')
##
##ax3.set_xlabel('Whole culm inclination ($^{o}$)', fontsize=fontsize+5)
##[i.set_linewidth(thickness) for i in ax3.spines.itervalues()]
##ax3.tick_params(width = thickness, length = 10, labelsize=fontsize+2, labelleft='off')
##ax3.set_xticks(np.arange(0, 105, 15))
##ax3.set_yticks(np.arange(0, 1., 0.2))
##ax3.text(0.1, 0.95, 'C', transform=ax3.transAxes, fontsize=fontsize+5, verticalalignment='top')
##
##plt.savefig(r'D:\Documents\PostDoc_Grignon\Modeles\FSPM_Wheat\example\Papier_FSPMA2016\Figures\Figure_02_LAD_Inclination.TIFF', dpi=300, format='TIFF', bbox_inches='tight')
##plt.close()



# # PARa direct vs diffus
#
# PARi_theoritical_df = pd.read_csv('PARi_theorique.csv')
#
# DATA_DIRPATH = r'D:\Documents\PostDoc_Grignon\Modeles\FSPM_Wheat\example\Papier_FSPMA2016\Reconstruction_interception\outputs'
# DENSITIES = [200, 410, 600, 800]
# STANDS = ['Caphorn', 'Soissons', 'Mixture']
# ORGAN_SELECT = [('blade', 4), ('blade', 3), ('internode', 6)]
#
# AXES_MAP = {200:0, 410:1, 600:2, 800:3}
# markers = {'Soissons': None, 'Caphorn': None, 'Soissons_Mixed': None, 'Caphorn_Mixed': None}
# colors = {'Soissons': 'r', 'Caphorn': 'b', 'Soissons_Mixed': 'r', 'Caphorn_Mixed': 'b'}
# lines = {'Soissons': '-', 'Caphorn': '-', 'Soissons_Mixed': '--', 'Caphorn_Mixed': '--'}
# label_mapping = {'Caphorn_Mixed':"Mixed erectophile", 'Soissons_Mixed':"Mixed planophile", 'Caphorn': "Erectophile", 'Soissons': "Planophile"}
# alpha= 0.2
#
# # Direct
# fig = plt.figure(figsize=cm2inch(30,25))
# gs = gridspec.GridSpec(3, 4, hspace=0, wspace=0)
#
# axes=[]
# IC95 = {'density':[], 't':[], 'species':[], 'organ':[], 'metamer':[], 'IC95_up':[], 'IC95_bot':[]}
#
# for density in DENSITIES:
#     for stand in STANDS:
#         # CSV files
#         PARa_direct_csv = os.path.join(DATA_DIRPATH, stand, 'sun_DOY[150, 199]_H[5, 19]_1000_plants_Density_{}.csv'.format(density))
#
#         # dfs
#         PARa_direct_df_all_data = pd.read_csv(PARa_direct_csv)
#
#         # Select data for 5<=t<=19 and listed organs
#         PARa_direct_df = PARa_direct_df_all_data[(5 <= PARa_direct_df_all_data['t']) & (PARa_direct_df_all_data['t'] <= 19)]
#
#         # Computation of theoritical PARa
#         ids = []
#         PARa_direct_theoritical_list = []
#         IC95_list_bot = []
#         IC95_list_up = []
#         for group_name, group_data in PARa_direct_df.groupby(['t', 'species', 'organ', 'metamer']):
#             if (group_name[2], group_name[3]) not in ORGAN_SELECT: continue
#             ids.append(group_name)
#             t = group_name[0]
#             PARi_theoritical = PARi_theoritical_df['PARi theorique'][PARi_theoritical_df['hour']==t].iloc[0]
#             PARa_direct_theoritical_list.append(group_data['Eabs'].iloc[0] * PARi_theoritical * 3600 * 1E-6)
#             if group_data['IC95'].iloc[0] == '(nan, nan)':
#                 IC95_list_bot.append(group_data['Eabs'].iloc[0] * PARi_theoritical * 3600 * 1E-6)
#                 IC95_list_up.append(group_data['Eabs'].iloc[0] * PARi_theoritical * 3600 * 1E-6)
#             else:
#                 IC95_list_bot.append(literal_eval(group_data['IC95'].iloc[0])[0] * PARi_theoritical * 3600 * 1E-6)
#                 IC95_list_up.append(literal_eval(group_data['IC95'].iloc[0])[1] * PARi_theoritical * 3600 * 1E-6)
#
#             # # if isinstance(group_data['IC95'].iloc[0], float):
#             # IC95['density'].append(density)
#             # IC95['t'].append(t)
#             # IC95['species'].append(group_name[1])
#             # IC95['organ'].append(group_name[2])
#             # IC95['metamer'].append(group_name[3])
#             # IC95['IC95_bot'].append(literal_eval(group_data['IC95'].iloc[0])[0] * PARi_theoritical * 3600 * 1E-6)
#             # IC95['IC95_up'].append(literal_eval(group_data['IC95'].iloc[0])[1] * PARi_theoritical * 3600 * 1E-6)
#
#         ids_df = pd.DataFrame(ids, columns=['t', 'species', 'organ', 'metamer'])
#         data_df = pd.DataFrame({'PARa theorique': PARa_direct_theoritical_list, 'IC95_bot': IC95_list_bot, 'IC95_up': IC95_list_up})
#         PARa_direct_df = pd.concat([ids_df, data_df], axis=1)
#         PARa_direct_df.sort_values(['t', 'species', 'organ', 'metamer'], inplace=True)
#
#         #Graphs
#         for species in PARa_direct_df['species'].unique():
#             if stand == 'Mixture':
#                 cv_map = species + '_Mixed'
#             else:
#                 cv_map = species
#
#             ax1 = plt.subplot(gs[0, AXES_MAP[density]])
#             ax2 = plt.subplot(gs[1, AXES_MAP[density]])
#             ax3 = plt.subplot(gs[2, AXES_MAP[density]])
#             ## Blade 4
#             blade4_direct = PARa_direct_df[(PARa_direct_df['species']==species) & (PARa_direct_df['metamer']==4) & (PARa_direct_df['organ']=='blade')]
#             PARa_blade4_direct = blade4_direct['PARa theorique']
#             ax1.plot(PARa_direct_df['t'].unique(), PARa_blade4_direct, label=label_mapping[cv_map], marker=markers[cv_map], color=colors[cv_map], linestyle=lines[cv_map], linewidth = thickness)
#             ax1.fill_between(PARa_direct_df['t'].unique(), blade4_direct['IC95_bot'], blade4_direct['IC95_up'], alpha=alpha, color=colors[cv_map], linewidth=0.0)
#             ## Blade 3
#             blade3_direct = PARa_direct_df[(PARa_direct_df['species']==species) & (PARa_direct_df['metamer']==3) & (PARa_direct_df['organ']=='blade')]
#             PARa_blade3_direct  = blade3_direct['PARa theorique']
#             ax2.plot(PARa_direct_df['t'].unique(), PARa_blade3_direct, label=label_mapping[cv_map], marker=markers[cv_map], color=colors[cv_map], linestyle=lines[cv_map], linewidth = thickness)
#             ax2.fill_between(PARa_direct_df['t'].unique(), blade3_direct['IC95_bot'], blade3_direct['IC95_up'], alpha=alpha, color=colors[cv_map], linewidth=0.0)
#             ## Chaff
#             chaff_direct = PARa_direct_df[(PARa_direct_df['species']==species) & (PARa_direct_df['metamer']==6) & (PARa_direct_df['organ']=='internode')]
#             PARa_chaff_direct = chaff_direct['PARa theorique']
#             ax3.plot(PARa_direct_df['t'].unique(), PARa_chaff_direct, label=label_mapping[cv_map], marker=markers[cv_map], color=colors[cv_map], linestyle=lines[cv_map], linewidth = thickness)
#             ax3.fill_between(PARa_direct_df['t'].unique(), chaff_direct['IC95_bot'], chaff_direct['IC95_up'], alpha=alpha, color=colors[cv_map], linewidth=0.0)
#
#             # Store axes
#             for ax in (ax1, ax2, ax3):
#                 axes.append(ax)
#
#
# IC95_df = pd.DataFrame(IC95)
# mean_IC95 = {}
# for density in DENSITIES:
#     mean_IC95[density] = {}
#     for organ_id in ORGAN_SELECT:
#         mean_IC95[density][organ_id] = {'t':[], 'IC95':[]}
#
# for group_name, group_data in IC95_df.groupby(['density', 't', 'organ', 'metamer']):
#     density = group_name[0]
#     organ_id = (group_name[2], group_name[3])
#     mean_IC95[density][organ_id]['t'].append(group_name[1])
#     mean_IC95[density][organ_id]['IC95'].append(group_data['IC95'].mean())
#
# # for density in mean_IC95.keys():
# #     ax = plt.subplot(gs[0, AXES_MAP[density]])
# #     t = mean_IC95[density][('blade', 4)]['t']
# #     ax.bar(t, mean_IC95[density][('blade', 4)]['IC95'])
# #     ax.set_yticks(np.arange(0, 4, 1))
#
# for ax in axes:
#     [i.set_linewidth(thickness) for i in ax.spines.itervalues()]
#     ax.tick_params(width = thickness, length = 10, labelbottom='off', labelleft='off', top='off', right='off')
#     ax.set_yticks(np.arange(0, 4, 1))
#     ax.set_xticks([5,9,13,17])
#
# # Labels
# ax_label_left = (plt.subplot(gs[0, 0]), plt.subplot(gs[1, 0]), plt.subplot(gs[2, 0]))
# for ax in ax_label_left:
#     ax.tick_params(labelsize=fontsize+2, labelleft='on')
#     if ax !=plt.subplot(gs[0, 0]):
#         ax.set_yticks(ax.get_yticks()[:-1])
# ax_label_bottom = (plt.subplot(gs[2, 0]), plt.subplot(gs[2, 1]), plt.subplot(gs[2, 2]), plt.subplot(gs[2, 3]))
# for ax in ax_label_bottom:
#     ax.tick_params(labelsize=fontsize+2, labelbottom='on')
#
# plt.subplot(gs[1, 0]).set_ylabel('PAR absorbed (mol m$^{-2}$ h$^{-1}$)', fontsize=fontsize+5)
# plt.text(30, -3.5, 'Time (hour UTC)', fontsize=fontsize+5)
#
# plt.text(8, 6.2, 'Density 200', fontsize=fontsize+5)
# plt.text(24, 6.2, 'Density 410', fontsize=fontsize+5)
# plt.text(40, 6.2, 'Density 600', fontsize=fontsize+5)
# plt.text(56, 6.2, 'Density 800', fontsize=fontsize+5)
#
# plt.subplot(gs[0, 1]).text(17, 2.5, 'Lamina n', fontsize=fontsize+5, backgroundcolor='w', verticalalignment='top')
# plt.subplot(gs[0, 1]).set_zorder(1)
# plt.subplot(gs[1, 1]).text(17, 2.5, 'Lamina n-1', fontsize=fontsize+5, backgroundcolor='w', verticalalignment='top')
# plt.subplot(gs[1, 1]).set_zorder(1)
# plt.subplot(gs[2, 1]).text(18, 2.5, 'Chaff', fontsize=fontsize+5, backgroundcolor='w', verticalalignment='top')
# plt.subplot(gs[2, 1]).set_zorder(1)
#
# plt.savefig(r'D:\Documents\PostDoc_Grignon\Modeles\FSPM_Wheat\example\Papier_FSPMA2016\Figures\Figure_03_PAR_direct.TIFF', dpi=300, format='TIFF', bbox_inches='tight')
# plt.close()
#
#
# # Diffus
# fig = plt.figure(figsize=cm2inch(30,25))
# gs = gridspec.GridSpec(3, 4, hspace=0, wspace=0)
#
# axes=[]
# IC95 = {'density':[], 't':[], 'species':[], 'organ':[], 'metamer':[], 'IC95_up':[], 'IC95_bot':[]}
#
# for density in DENSITIES:
#     for stand in STANDS:
#         # CSV files
#         PARa_diffus_csv = os.path.join(DATA_DIRPATH, stand, 'soc_a4z5_1000_plants_Density_{}.csv'.format(density))
#
#         # dfs
#         PARa_diffus_df_all_data = pd.read_csv(PARa_diffus_csv)
#
#         # Select data for 5<=t<=19 and listed organs
#         PARa_diffus_df = PARa_diffus_df_all_data[(5 <= PARa_diffus_df_all_data['t']) & (PARa_diffus_df_all_data['t'] <= 19)]
#
#         # Computation of theoritical PARa
#         ids = []
#         PARa_diffus_theoritical_list = []
#         IC95_list_bot = []
#         IC95_list_up = []
#         for group_name, group_data in PARa_diffus_df.groupby(['t', 'species', 'organ', 'metamer']):
#             if (group_name[2], group_name[3]) not in ORGAN_SELECT: continue
#             ids.append(group_name)
#             t = group_name[0]
#             PARi_theoritical = PARi_theoritical_df['PARi theorique'][PARi_theoritical_df['hour']==t].iloc[0]
#             PARa_diffus_theoritical_list.append(group_data['Eabs'].iloc[0] * PARi_theoritical * 3600 * 1E-6)
#             if group_data['IC95'].iloc[0] == '(nan, nan)':
#                 IC95_list_bot.append(group_data['Eabs'].iloc[0] * PARi_theoritical * 3600 * 1E-6)
#                 IC95_list_up.append(group_data['Eabs'].iloc[0] * PARi_theoritical * 3600 * 1E-6)
#             else:
#                 IC95_list_bot.append(literal_eval(group_data['IC95'].iloc[0])[0] * PARi_theoritical * 3600 * 1E-6)
#                 IC95_list_up.append(literal_eval(group_data['IC95'].iloc[0])[1] * PARi_theoritical * 3600 * 1E-6)
#
#             # # if isinstance(group_data['IC95'].iloc[0], float):
#             # IC95['density'].append(density)
#             # IC95['t'].append(t)
#             # IC95['species'].append(group_name[1])
#             # IC95['organ'].append(group_name[2])
#             # IC95['metamer'].append(group_name[3])
#             # IC95['IC95_bot'].append(literal_eval(group_data['IC95'].iloc[0])[0] * PARi_theoritical * 3600 * 1E-6)
#             # IC95['IC95_up'].append(literal_eval(group_data['IC95'].iloc[0])[1] * PARi_theoritical * 3600 * 1E-6)
#
#         ids_df = pd.DataFrame(ids, columns=['t', 'species', 'organ', 'metamer'])
#         data_df = pd.DataFrame({'PARa theorique': PARa_diffus_theoritical_list, 'IC95_bot': IC95_list_bot, 'IC95_up': IC95_list_up})
#         PARa_diffus_df = pd.concat([ids_df, data_df], axis=1)
#         PARa_diffus_df.sort_values(['t', 'species', 'organ', 'metamer'], inplace=True)
#
#         #Graphs
#         for species in PARa_diffus_df['species'].unique():
#             if stand == 'Mixture':
#                 cv_map = species + '_Mixed'
#             else:
#                 cv_map = species
#
#             ax1 = plt.subplot(gs[0, AXES_MAP[density]])
#             ax2 = plt.subplot(gs[1, AXES_MAP[density]])
#             ax3 = plt.subplot(gs[2, AXES_MAP[density]])
#             ## Blade 4
#             blade4_diffus = PARa_diffus_df[(PARa_diffus_df['species']==species) & (PARa_diffus_df['metamer']==4) & (PARa_diffus_df['organ']=='blade')]
#             PARa_blade4_diffus = blade4_diffus['PARa theorique']
#             ax1.plot(PARa_diffus_df['t'].unique(), PARa_blade4_diffus, label=label_mapping[cv_map], marker=markers[cv_map], color=colors[cv_map], linestyle=lines[cv_map], linewidth = thickness)
#             ax1.fill_between(PARa_diffus_df['t'].unique(), blade4_diffus['IC95_bot'], blade4_diffus['IC95_up'], alpha=alpha, color=colors[cv_map], linewidth=0.0)
#             ## Blade 3
#             blade3_diffus = PARa_diffus_df[(PARa_diffus_df['species']==species) & (PARa_diffus_df['metamer']==3) & (PARa_diffus_df['organ']=='blade')]
#             PARa_blade3_diffus  = blade3_diffus['PARa theorique']
#             ax2.plot(PARa_diffus_df['t'].unique(), PARa_blade3_diffus, label=label_mapping[cv_map], marker=markers[cv_map], color=colors[cv_map], linestyle=lines[cv_map], linewidth = thickness)
#             ax2.fill_between(PARa_diffus_df['t'].unique(), blade3_diffus['IC95_bot'], blade3_diffus['IC95_up'], alpha=alpha, color=colors[cv_map], linewidth=0.0)
#             ## Chaff
#             chaff_diffus = PARa_diffus_df[(PARa_diffus_df['species']==species) & (PARa_diffus_df['metamer']==6) & (PARa_diffus_df['organ']=='internode')]
#             PARa_chaff_diffus = chaff_diffus['PARa theorique']
#             ax3.plot(PARa_diffus_df['t'].unique(), PARa_chaff_diffus, label=label_mapping[cv_map], marker=markers[cv_map], color=colors[cv_map], linestyle=lines[cv_map], linewidth = thickness)
#             ax3.fill_between(PARa_diffus_df['t'].unique(), chaff_diffus['IC95_bot'], chaff_diffus['IC95_up'], alpha=alpha, color=colors[cv_map], linewidth=0.0)
#
#             # Store axes
#             for ax in (ax1, ax2, ax3):
#                 axes.append(ax)
#
#
# IC95_df = pd.DataFrame(IC95)
# mean_IC95 = {}
# for density in DENSITIES:
#     mean_IC95[density] = {}
#     for organ_id in ORGAN_SELECT:
#         mean_IC95[density][organ_id] = {'t':[], 'IC95':[]}
#
# for group_name, group_data in IC95_df.groupby(['density', 't', 'organ', 'metamer']):
#     density = group_name[0]
#     organ_id = (group_name[2], group_name[3])
#     mean_IC95[density][organ_id]['t'].append(group_name[1])
#     mean_IC95[density][organ_id]['IC95'].append(group_data['IC95'].mean())
#
# # for density in mean_IC95.keys():
# #     ax = plt.subplot(gs[0, AXES_MAP[density]])
# #     t = mean_IC95[density][('blade', 4)]['t']
# #     ax.bar(t, mean_IC95[density][('blade', 4)]['IC95'])
# #     ax.set_yticks(np.arange(0, 4, 1))
#
# for ax in axes:
#     [i.set_linewidth(thickness) for i in ax.spines.itervalues()]
#     ax.tick_params(width = thickness, length = 10, labelbottom='off', labelleft='off', top='off', right='off')
#     ax.set_yticks(np.arange(0, 4, 1))
#     ax.set_xticks([5,9,13,17])
#
# # Labels
# ax_label_left = (plt.subplot(gs[0, 0]), plt.subplot(gs[1, 0]), plt.subplot(gs[2, 0]))
# for ax in ax_label_left:
#     ax.tick_params(labelsize=fontsize+2, labelleft='on')
#     if ax !=plt.subplot(gs[0, 0]):
#         ax.set_yticks(ax.get_yticks()[:-1])
# ax_label_bottom = (plt.subplot(gs[2, 0]), plt.subplot(gs[2, 1]), plt.subplot(gs[2, 2]), plt.subplot(gs[2, 3]))
# for ax in ax_label_bottom:
#     ax.tick_params(labelsize=fontsize+2, labelbottom='on')
#
# plt.subplot(gs[1, 0]).set_ylabel('PAR absorbed (mol m$^{-2}$ h$^{-1}$)', fontsize=fontsize+5)
# plt.text(30, -3.5, 'Time (hour UTC)', fontsize=fontsize+5)
#
# plt.text(8, 6.2, 'Density 200', fontsize=fontsize+5)
# plt.text(24, 6.2, 'Density 410', fontsize=fontsize+5)
# plt.text(40, 6.2, 'Density 600', fontsize=fontsize+5)
# plt.text(56, 6.2, 'Density 800', fontsize=fontsize+5)
#
# plt.subplot(gs[0, 1]).text(17, 2.5, 'Lamina n', fontsize=fontsize+5, backgroundcolor='w', verticalalignment='top')
# plt.subplot(gs[0, 1]).set_zorder(1)
# plt.subplot(gs[1, 1]).text(17, 2.5, 'Lamina n-1', fontsize=fontsize+5, backgroundcolor='w', verticalalignment='top')
# plt.subplot(gs[1, 1]).set_zorder(1)
# plt.subplot(gs[2, 1]).text(18, 2.5, 'Chaff', fontsize=fontsize+5, backgroundcolor='w', verticalalignment='top')
# plt.subplot(gs[2, 1]).set_zorder(1)
#
# plt.savefig(r'D:\Documents\PostDoc_Grignon\Modeles\FSPM_Wheat\example\Papier_FSPMA2016\Figures\Figure_04_PAR_diffus.TIFF', dpi=300, format='TIFF', bbox_inches='tight')
# plt.close()








##fig = plt.figure(figsize=cm2inch(30,15))
##gs = gridspec.GridSpec(2, 4, hspace=0, wspace=0)
##
##for density in DENSITIES:
##    for stand in STANDS:
##        # mtg
##        adel_wheat = AdelWheat(seed=1234, convUnit=1)
##        g = adel_wheat.load(dir=os.path.join(DATA_DIRPATH, 'Densite_{}'.format(density), stand, 'inputs', 'adelwheat'))
##        # Photosynthetic elements
##        ph_elements_path = os.path.join(DATA_DIRPATH, 'Densite_{}'.format(density), stand, 'outputs\elements_states.csv')
##        ph_elements = pd.read_csv(ph_elements_path)
##        # Organs
##        organs_path = os.path.join(DATA_DIRPATH, 'Densite_{}'.format(density), stand, 'outputs\organs_states.csv')
##        organs = pd.read_csv(organs_path)
##        # Soil
##        soil_path = os.path.join(DATA_DIRPATH, 'Densite_{}'.format(density), stand, 'outputs\soils_states.csv')
##        soil = pd.read_csv(soil_path)
##
##        # 1.1) Lamina area
##        ax1 = plt.subplot(gs[0,AXES_MAP[density]])
##
##        if stand == 'Asso_diffus':
##            plant1_vid, plant2_vid = g.components(g.root)[0], g.components(g.root)[1]
##            plant1_species = g.property('species')[plant1_vid]
##            plant2_species = g.property('species')[plant2_vid]
##
##            if plant1_species=='Caphorn' and plant2_species=='Soissons':
##                green_area_lamina1 = ph_elements[(ph_elements['organ']=='blade') & (ph_elements['plant']==1)].groupby('t')['green_area'].aggregate(np.sum)*10000
##                green_area_lamina2 = ph_elements[(ph_elements['organ']=='blade') & (ph_elements['plant']==2)].groupby('t')['green_area'].aggregate(np.sum)*10000
##
##            else:
##                green_area_lamina1 = ph_elements[(ph_elements['organ']=='blade') & (ph_elements['plant']==2)].groupby('t')['green_area'].aggregate(np.sum)*10000
##                green_area_lamina2 = ph_elements[(ph_elements['organ']=='blade') & (ph_elements['plant']==1)].groupby('t')['green_area'].aggregate(np.sum)*10000
##
##            ax1.plot(green_area_lamina1.index, green_area_lamina1, color_mapping[stand + '_1'], label = label_mapping[stand + '_1'], linewidth = thickness)
##            ax1.plot(green_area_lamina2.index, green_area_lamina2, color_mapping[stand + '_2'], label = label_mapping[stand + '_2'], linewidth = thickness)
##
##
##        else:
##            green_area_lamina = ph_elements[(ph_elements['organ']=='blade') & (ph_elements['plant']==1)].groupby('t')['green_area'].aggregate(np.sum)*10000
##            ax1.plot(green_area_lamina.index, green_area_lamina, color_mapping[stand], label = label_mapping[stand], linewidth = thickness)
##
##        # Formatting
##        ax1.set_xticks(np.arange(0, 1300, 200))
##        ax1.set_yticks(np.arange(0, 150, 25))
##        if AXES_MAP[density]==0:
##            ax1.tick_params(labelbottom= 'off', width = thickness, labelsize=fontsize)
##            ax1.set_ylabel(u'Lamina green area (cm$^{2}$)', fontsize=fontsize+2)
##        else:
##            ax1.tick_params(labelleft= 'off', labelbottom= 'off', width = thickness, labelsize=fontsize)
##        ax1.yaxis.set_label_coords(labelx, 0.5)
##        [i.set_linewidth(thickness) for i in ax1.spines.itervalues()]
##        ax1.axvline(360, color='k', linestyle='--')
##
##
##        # 1.2) Grain dry mass
##        ax2 = plt.subplot(gs[1,AXES_MAP[density]])
##
##        if stand == 'Asso_diffus':
##            if plant1_species=='Caphorn' and plant2_species=='Soissons':
##                grains1 = organs[(organs['organ']=='grains') & (organs['plant']==1)].groupby('t').sum()
##                sum_dry_mass_grains1 = grains1['Dry_Mass']
##                grains2 = organs[(organs['organ']=='grains') & (organs['plant']==2)].groupby('t').sum()
##                sum_dry_mass_grains2 = grains2['Dry_Mass']
##            else:
##                grains1 = organs[(organs['organ']=='grains') & (organs['plant']==2)].groupby('t').sum()
##                sum_dry_mass_grains1 = grains1['Dry_Mass']
##                grains2 = organs[(organs['organ']=='grains') & (organs['plant']==1)].groupby('t').sum()
##                sum_dry_mass_grains2 = grains2['Dry_Mass']
##
##            ax2.plot(sum_dry_mass_grains1.index, sum_dry_mass_grains1, color_mapping[stand + '_1'], label = label_mapping[stand + '_1'], linewidth = thickness)
##            ax2.plot(sum_dry_mass_grains2.index, sum_dry_mass_grains2, color_mapping[stand + '_2'], label = label_mapping[stand + '_2'], linewidth = thickness)
##
##        else:
##            grains = organs[(organs['organ']=='grains')].groupby('t')['Dry_Mass']
##            sum_dry_mass_grains = grains.mean()
##            ax2.plot(sum_dry_mass_grains.index, sum_dry_mass_grains, color_mapping[stand], label = label_mapping[stand], linewidth = thickness)
##
##        # Formatting
##        ax2.set_xticks(np.arange(0, 1300, 200))
##        ax2.set_xlabel('Time from flowering (hour)', fontsize=fontsize)
##        ax2.set_yticks(np.arange(0, 4.5, 1))
##        ax2.set_yticks(ax2.get_yticks()[:-1])
##        if AXES_MAP[density]==0:
##            ax2.tick_params(labelbottom='on', width = thickness, labelsize=fontsize)
##            ax2.set_ylabel('Grain dry mass (g)', fontsize=fontsize+2)
##            ax2.set_xticks(ax2.get_xticks()[:-1])
##        else:
##            ax2.tick_params(labelleft= 'off', labelbottom= 'on', width = thickness, labelsize=fontsize)
##            if AXES_MAP[density]!=3:
##                ax2.set_xticks(ax2.get_xticks()[:-1])
##        ax2.yaxis.set_label_coords(labelx, 0.5)
##        [i.set_linewidth(thickness) for i in ax2.spines.itervalues()]
##        ax2.axvline(360, color='k', linestyle='--')
##
##
##plt.subplot(gs[0, 0]).text(300, 140, 'Density 200', fontsize=fontsize+2, backgroundcolor='w', verticalalignment='top')
##plt.subplot(gs[0, 0]).set_zorder(1)
##plt.subplot(gs[0, 1]).text(300, 140, 'Density 410', fontsize=fontsize+2, backgroundcolor='w', verticalalignment='top')
##plt.subplot(gs[0, 1]).set_zorder(1)
##plt.subplot(gs[0, 2]).text(300, 140, 'Density 600', fontsize=fontsize+2, backgroundcolor='w', verticalalignment='top')
##plt.subplot(gs[0, 2]).set_zorder(1)
##plt.subplot(gs[0, 3]).text(300, 140, 'Density 800', fontsize=fontsize+2, backgroundcolor='w', verticalalignment='top')
##plt.subplot(gs[0, 3]).set_zorder(1)
##
####plt.legend(prop={'size':11}, bbox_to_anchor=(0.5, 1.2), ncol=4, frameon=True)
##plt.savefig(os.path.join(GRAPHS_DIRPATH, 'Figure_06_Culm_scale.TIFF'), dpi=300, format='TIFF', bbox_inches='tight')
##plt.close()
##
### C non struct
##phloem_shoot_root = 0.75
##C_MOLAR_MASS = 12
##fig = plt.figure(figsize=cm2inch(30,20))
##gs = gridspec.GridSpec(3, 4, hspace=0, wspace=0)
##
##for density in DENSITIES:
##    for stand in STANDS:
##        axes = []
##        # mtg
##        adel_wheat = AdelWheat(seed=1234, convUnit=1)
##        g = adel_wheat.load(dir=os.path.join(DATA_DIRPATH, 'Densite_{}'.format(density), stand, 'inputs', 'adelwheat'))
##        # Photosynthetic elements
##        ph_elements_path = os.path.join(DATA_DIRPATH, 'Densite_{}'.format(density), stand, 'outputs\elements_states.csv')
##        ph_elements = pd.read_csv(ph_elements_path)
##        # Organs
##        organs_path = os.path.join(DATA_DIRPATH, 'Densite_{}'.format(density), stand, 'outputs\organs_states.csv')
##        organs = pd.read_csv(organs_path)
##        # Soil
##        soil_path = os.path.join(DATA_DIRPATH, 'Densite_{}'.format(density), stand, 'outputs\soils_states.csv')
##        soil = pd.read_csv(soil_path)
##
##        ax1 = plt.subplot(gs[0,AXES_MAP[density]])
##        ax2 = plt.subplot(gs[1,AXES_MAP[density]])
##        ax3 = plt.subplot(gs[2,AXES_MAP[density]])
##        axes.append(ax1)
##        axes.append(ax2)
##        axes.append(ax3)
##
##        if stand == 'Asso_diffus':
##            plant1_vid, plant2_vid = g.components(g.root)[0], g.components(g.root)[1]
##            plant1_species = g.property('species')[plant1_vid]
##            plant2_species = g.property('species')[plant2_vid]
##
##            if plant1_species=='Caphorn' and plant2_species=='Soissons':
##                ## Stem
##                stem_1 = ph_elements[ph_elements.organ.isin(['internode', 'sheath', 'peduncle']) & (ph_elements['plant']==1)]
##                stem_grouped_1 = stem_1[['t', 'triosesP', 'starch', 'fructan', 'sucrose', 'green_area']].groupby('t').aggregate(np.sum)
##                area_tot_ph_1 = stem_1.groupby('t')['green_area'].sum()
##                contrib_area_stem_1 = (stem_grouped_1['green_area'].reset_index(drop=True) / area_tot_ph_1) * phloem_shoot_root
##                sum_C_1 = (stem_grouped_1['triosesP'] + stem_grouped_1['starch'] + stem_grouped_1['fructan'] + stem_grouped_1['sucrose']) * (1E-3*C_MOLAR_MASS)
##                phloem_1 = organs[(organs['organ']== 'phloem') & (organs['plant']==1)].reset_index()
##                sum_C_1 = sum_C_1.add(phloem_1['sucrose'] * (1E-3*C_MOLAR_MASS) * contrib_area_stem_1, fill_value=0)
##                ax1.plot(sum_C_1.index, sum_C_1, color_mapping[stand + '_1'], label = label_mapping[stand + '_1'], linewidth = thickness)
##                ## Roots
##                roots_1 = organs[(organs['organ']== 'phloem') & (organs['plant']==1)].reset_index()
##                sum_C_1 = (roots_1['sucrose'] * (1E-3*C_MOLAR_MASS)).reset_index(drop=True)
##                sum_C_1 = sum_C_1.add(phloem_1['sucrose'] * (1E-3*C_MOLAR_MASS)*(1 - phloem_shoot_root), fill_value=0)
##                ax3.plot(sum_C_1.index, sum_C_1, color_mapping[stand + '_1'], label = label_mapping[stand + '_1'], linewidth = thickness)
##                ## Phloem
##                sum_C_1 = (phloem_1['sucrose'] * (1E-3*C_MOLAR_MASS)).reset_index(drop=True)
##                ax2.plot(sum_C_1.index, sum_C_1, color_mapping[stand + '_1'], label = label_mapping[stand + '_1'], linewidth = thickness)
##
##
##                ## Stem
##                stem_2 = ph_elements[ph_elements.organ.isin(['internode', 'sheath', 'peduncle']) & (ph_elements['plant']==2)]
##                stem_grouped_2 = stem_2[['t', 'triosesP', 'starch', 'fructan', 'sucrose', 'green_area']].groupby('t').aggregate(np.sum)
##                area_tot_ph_2 = stem_2.groupby('t')['green_area'].sum()
##                contrib_area_stem_2 = (stem_grouped_2['green_area'].reset_index(drop=True) / area_tot_ph_2) * phloem_shoot_root
##                sum_C_2 = (stem_grouped_2['triosesP'] + stem_grouped_2['starch'] + stem_grouped_2['fructan'] + stem_grouped_2['sucrose']) * (1E-3*C_MOLAR_MASS)
##                phloem_2 = organs[(organs['organ']== 'phloem') & (organs['plant']==2)].reset_index()
##                sum_C_2 = sum_C_2.add(phloem_2['sucrose'] * (1E-3*C_MOLAR_MASS) * contrib_area_stem_2, fill_value=0)
##                ax1.plot(sum_C_2.index, sum_C_2, color_mapping[stand + '_2'], label = label_mapping[stand + '_2'], linewidth = thickness)
##                ## Roots
##                roots_2 = organs[(organs['organ']== 'phloem') & (organs['plant']==2)].reset_index()
##                sum_C_2 = (roots_2['sucrose'] * (1E-3*C_MOLAR_MASS)).reset_index(drop=True)
##                sum_C_2 = sum_C_2.add(phloem_2['sucrose'] * (1E-3*C_MOLAR_MASS)*(1 - phloem_shoot_root), fill_value=0)
##                ax3.plot(sum_C_2.index, sum_C_2, color_mapping[stand + '_2'], label = label_mapping[stand + '_2'], linewidth = thickness)
##                ## Phloem
##                sum_C_2 = (phloem_2['sucrose'] * (1E-3*C_MOLAR_MASS)).reset_index(drop=True)
##                ax2.plot(sum_C_2.index, sum_C_2, color_mapping[stand + '_2'], label = label_mapping[stand + '_2'], linewidth = thickness)
##
##            else:
##                ## Stem
##                stem_1 = ph_elements[ph_elements.organ.isin(['internode', 'sheath', 'peduncle']) & (ph_elements['plant']==2)]
##                stem_grouped_1 = stem_1[['t', 'triosesP', 'starch', 'fructan', 'sucrose', 'green_area']].groupby('t').aggregate(np.sum)
##                area_tot_ph_1 = stem_1.groupby('t')['green_area'].sum()
##                contrib_area_stem_1 = (stem_grouped_1['green_area'].reset_index(drop=True) / area_tot_ph_1) * phloem_shoot_root
##                sum_C_1 = (stem_grouped_1['triosesP'] + stem_grouped_1['starch'] + stem_grouped_1['fructan'] + stem_grouped_1['sucrose']) * (1E-3*C_MOLAR_MASS)
##                phloem_1 = organs[(organs['organ']== 'phloem') & (organs['plant']==2)].reset_index()
##                sum_C_1 = sum_C_1.add(phloem_1['sucrose'] * (1E-3*C_MOLAR_MASS) * contrib_area_stem_1, fill_value=0)
##                ax1.plot(sum_C_1.index, sum_C_1, color_mapping[stand + '_1'], label = label_mapping[stand + '_1'], linewidth = thickness)
##                ## Roots
##                roots_1 = organs[(organs['organ']== 'phloem') & (organs['plant']==2)].reset_index()
##                sum_C_1 = (roots_1['sucrose'] * (1E-3*C_MOLAR_MASS)).reset_index(drop=True)
##                sum_C_1 = sum_C_1.add(phloem_1['sucrose'] * (1E-3*C_MOLAR_MASS)*(1 - phloem_shoot_root), fill_value=0)
##                ax3.plot(sum_C_1.index, sum_C_1, color_mapping[stand + '_1'], label = label_mapping[stand + '_1'], linewidth = thickness)
##                ## Phloem
##                sum_C_1 = (phloem_1['sucrose'] * (1E-3*C_MOLAR_MASS)).reset_index(drop=True)
##                ax2.plot(sum_C_1.index, sum_C_1, color_mapping[stand + '_1'], label = label_mapping[stand + '_1'], linewidth = thickness)
##
##
##                ## Stem
##                stem_2 = ph_elements[ph_elements.organ.isin(['internode', 'sheath', 'peduncle']) & (ph_elements['plant']==1)]
##                stem_grouped_2 = stem_2[['t', 'triosesP', 'starch', 'fructan', 'sucrose', 'green_area']].groupby('t').aggregate(np.sum)
##                area_tot_ph_2 = stem_2.groupby('t')['green_area'].sum()
##                contrib_area_stem_2 = (stem_grouped_2['green_area'].reset_index(drop=True) / area_tot_ph_2) * phloem_shoot_root
##                sum_C_2 = (stem_grouped_2['triosesP'] + stem_grouped_2['starch'] + stem_grouped_2['fructan'] + stem_grouped_2['sucrose']) * (1E-3*C_MOLAR_MASS)
##                phloem_2 = organs[(organs['organ']== 'phloem') & (organs['plant']==1)].reset_index()
##                sum_C_2 = sum_C_2.add(phloem_2['sucrose'] * (1E-3*C_MOLAR_MASS) * contrib_area_stem_2, fill_value=0)
##                ax1.plot(sum_C_2.index, sum_C_2, color_mapping[stand + '_2'], label = label_mapping[stand + '_2'], linewidth = thickness)
##                ## Roots
##                roots_2 = organs[(organs['organ']== 'phloem') & (organs['plant']==1)].reset_index()
##                sum_C_2 = (roots_2['sucrose'] * (1E-3*C_MOLAR_MASS)).reset_index(drop=True)
##                sum_C_2 = sum_C_2.add(phloem_2['sucrose'] * (1E-3*C_MOLAR_MASS)*(1 - phloem_shoot_root), fill_value=0)
##                ax3.plot(sum_C_2.index, sum_C_2, color_mapping[stand + '_2'], label = label_mapping[stand + '_2'], linewidth = thickness)
##                ## Phloem
##                sum_C_2 = (phloem_2['sucrose'] * (1E-3*C_MOLAR_MASS)).reset_index(drop=True)
##                ax2.plot(sum_C_2.index, sum_C_2, color_mapping[stand + '_2'], label = label_mapping[stand + '_2'], linewidth = thickness)
##
##        else:
##            ## Stem
##            stem = ph_elements[ph_elements.organ.isin(['internode', 'sheath', 'peduncle']) & (ph_elements['plant']==1)]
##            stem_grouped = stem[['t', 'triosesP', 'starch', 'fructan', 'sucrose', 'green_area']].groupby('t').aggregate(np.sum)
##            area_tot_ph = stem.groupby('t')['green_area'].sum()
##            contrib_area_stem = (stem_grouped['green_area'].reset_index(drop=True) / area_tot_ph) * phloem_shoot_root
##            sum_C = (stem_grouped['triosesP'] + stem_grouped['starch'] + stem_grouped['fructan'] + stem_grouped['sucrose']) * (1E-3*C_MOLAR_MASS)
##            phloem = organs[(organs['organ']== 'phloem')].reset_index()
##            sum_C = sum_C.add(phloem['sucrose'] * (1E-3*C_MOLAR_MASS) * contrib_area_stem, fill_value=0)
##            ax1.plot(sum_C.index, sum_C, color_mapping[stand + ''], label = label_mapping[stand + ''], linewidth = thickness)
##            ## Roots
##            roots = organs[(organs['organ']== 'phloem') & (organs['plant']==1)].reset_index()
##            sum_C = (roots['sucrose'] * (1E-3*C_MOLAR_MASS)).reset_index(drop=True)
##            sum_C = sum_C.add(phloem['sucrose'] * (1E-3*C_MOLAR_MASS)*(1 - phloem_shoot_root), fill_value=0)
##            ax3.plot(sum_C.index, sum_C, color_mapping[stand + ''], label = label_mapping[stand + ''], linewidth = thickness)
##            ## Phloem
##            sum_C = (phloem['sucrose'] * (1E-3*C_MOLAR_MASS)).reset_index(drop=True)
##            ax2.plot(sum_C.index, sum_C, color_mapping[stand + ''], label = label_mapping[stand + ''], linewidth = thickness)
##
##
##        for ax in axes:
##            [i.set_linewidth(thickness) for i in ax.spines.itervalues()]
##            ax.tick_params(width = thickness, length = 10, labelbottom='off', labelleft='off', top='off', right='off', labelsize=fontsize)
##            ax.set_yticks(np.arange(0, 800, 200))
##            ax.set_xticks(np.arange(0, 1300, 200))
##            ax.axvline(360, color='k', linestyle='--')
##
##
##ax_label_left = (plt.subplot(gs[0, 0]), plt.subplot(gs[1, 0]), plt.subplot(gs[2, 0]))
##for ax in ax_label_left:
##    ax.tick_params(labelleft='on')
##    if ax == plt.subplot(gs[1, 0]):
##        ax.set_ylabel('Non structural C mass (mg C)', fontsize=fontsize+2)
##    if ax != plt.subplot(gs[0, 0]):
##        ax.set_yticks(ax.get_yticks()[:-1])
##
##
##ax_label_bottom = (plt.subplot(gs[2, 0]), plt.subplot(gs[2, 1]), plt.subplot(gs[2, 2]), plt.subplot(gs[2, 3]))
##for ax in ax_label_bottom:
##    ax.tick_params(labelbottom='on')
##    ax.set_xlabel('Time from flowering (hour)', fontsize=fontsize)
##    if ax != plt.subplot(gs[2, 3]):
##        ax.set_xticks(ax.get_xticks()[:-1])
##
##xtext = 300
##ytext = 675
##plt.subplot(gs[0, 0]).text(xtext, ytext, 'Density 200', fontsize=fontsize+3, backgroundcolor='w', verticalalignment='top')
##plt.subplot(gs[0, 0]).set_zorder(1)
##plt.subplot(gs[0, 1]).text(xtext, ytext, 'Density 410', fontsize=fontsize+3, backgroundcolor='w', verticalalignment='top')
##plt.subplot(gs[0, 1]).set_zorder(1)
##plt.subplot(gs[0, 2]).text(xtext, ytext, 'Density 600', fontsize=fontsize+3, backgroundcolor='w', verticalalignment='top')
##plt.subplot(gs[0, 2]).set_zorder(1)
##plt.subplot(gs[0, 3]).text(xtext, ytext, 'Density 800', fontsize=fontsize+3, backgroundcolor='w', verticalalignment='top')
##plt.subplot(gs[0, 3]).set_zorder(1)
##
##xtext = 50
##ytext = 575
##plt.subplot(gs[0, 0]).text(xtext, ytext, 'Stem', fontsize=fontsize+2, backgroundcolor='w', verticalalignment='top')
##plt.subplot(gs[0, 0]).set_zorder(1)
##plt.subplot(gs[1, 0]).text(xtext, ytext, 'Phloem', fontsize=fontsize+2, backgroundcolor='w', verticalalignment='top')
##plt.subplot(gs[1, 0]).set_zorder(1)
##plt.subplot(gs[2, 0]).text(xtext, ytext, 'Roots', fontsize=fontsize+2, backgroundcolor='w', verticalalignment='top')
##plt.subplot(gs[2, 0]).set_zorder(1)
##
###plt.legend(prop={'size':11}, bbox_to_anchor=(0.5, 1.2), ncol=4, frameon=True)
##plt.savefig(os.path.join(GRAPHS_DIRPATH, 'Figure_08_C_dynamics.TIFF'), dpi=300, format='TIFF', bbox_inches='tight')
##plt.close()
##
##
##




### PARa vs Cass
##data = r'D:\Documents\PostDoc_Grignon\Modeles\FSPM_Wheat\example\Papier_FSPMA2016_R1\Figures\Relation PAR-Cass.csv'
##data_df = pd.read_csv(data)
##data_df_grouped = data_df.groupby(['Density', 'Leaf posture'])
##
##fig, ax = plt.subplots()
##
##markers = {'Erect':'^', 'Plano':'s'}
##colors = {'D200':'b', 'D410':'g', 'D600':'y', 'D800':'r'}
##
##for group_name, group in data_df_grouped:
##    density = group_name[0]
##    leaf_inclination = group_name[1]
##    ax.plot(group['Total_PARa'], group['Photosynthesis'], marker = markers[leaf_inclination], color = colors[density], label=density)
##
##slope, intercept, r_value, p_value, std_err = stats.linregress(data_df['Total_PARa'], data_df['Photosynthesis'])
##line = slope * data_df['Total_PARa'] + intercept
##ax.plot(data_df['Total_PARa'], line, 'k--')
##ax.text(0.75, 0.8, 'y = {}x + {}'.format(round(slope, 3), round(intercept, 3)), transform=ax.transAxes, fontsize=fontsize, verticalalignment='top')
##
##[i.set_linewidth(thickness) for i in ax.spines.itervalues()]
##ax.set_xticks(np.arange(0, 12, 2))
##ax.set_yticks(np.arange(0, 300, 50))
##ax.set_xlabel('Total absorbed PAR (mol)', fontsize=fontsize+5)
##ax.set_ylabel('Total assimilated C (mmol)', fontsize=fontsize+5)
##
##plt.savefig(r'D:\Documents\PostDoc_Grignon\Modeles\FSPM_Wheat\example\Papier_FSPMA2016_R1\Figures\Figure_05_PARavsCass.TIFF', dpi=300, format='TIFF', bbox_inches='tight')
##plt.close()
##
# Cass vs grain BM
# data = r'D:\Documents\PostDoc_Grignon\Modeles\FSPM_Wheat\example\Papier_FSPMA2016_R1\Figures\Cass_vs_grain_BM.csv'
# data_df = pd.read_csv(data)
# data_df_grouped = data_df.groupby(['Density (culm m-2)', 'posture'])
#
# fig, ax = plt.subplots()
#
# markers = {'Erect':'^', 'Plano':'s'}
# colors = {200:'b', 410:'g', 600:'y', 800:'r'}
#
# for group_name, group in data_df_grouped:
#     density = group_name[0]
#     leaf_inclination = group_name[1]
#     ax.plot(group['Sum C (g)'], group['grain DM (g per culm)'], marker = markers[leaf_inclination], color = colors[density], label=density)
#
# z = np.polyfit(data_df['Sum C (g)'], data_df['grain DM (g per culm)'], 2)
# p = np.poly1d(z)
# xp = np.linspace(min(data_df['Sum C (g)']), max(data_df['Sum C (g)']), 100)
# ax.plot(xp, p(xp), 'k--')
# ##ax.text(0.4, 0.9, 'y = {}x$^2$ + {}x + {}'.format(round(z[0], 3), round(z[1], 3), round(z[2], 3)), transform=ax.transAxes, fontsize=fontsize, verticalalignment='top')
#
# [i.set_linewidth(thickness) for i in ax.spines.itervalues()]
# ax.set_xticks(np.arange(0, 4, 1))
# ax.set_yticks(np.arange(0, 4, 1))
# ax.set_xlabel('Total assimilated C (g)', fontsize=fontsize+5)
# ax.set_ylabel('Grain final dry mass (g)', fontsize=fontsize+5)
#
# plt.savefig(r'D:\Documents\PostDoc_Grignon\Modeles\FSPM_Wheat\example\Papier_FSPMA2016_R1\Figures\Figure_07_Cass_vs_grain_mass.TIFF', dpi=300, format='TIFF', bbox_inches='tight')
# plt.close()
