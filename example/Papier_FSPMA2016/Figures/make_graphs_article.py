# -*- coding: latin-1 -*-

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from ast import literal_eval
from scipy import stats
import statsmodels.stats.api as sms
import gc

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
markersize = 6

t_end = 1201
phloem_shoot_root = 0.75

DATA_DIRPATH = r'D:\Documents\PostDoc_Grignon\Modeles\FSPM_Wheat\example\Papier_FSPMA2016'
GRAPHS_DIRPATH = r'D:\Documents\PostDoc_Grignon\Modeles\FSPM_Wheat\example\Papier_FSPMA2016\Figures'

# DENSITIES = ['200', '410', '600', '800']
DENSITIES = ['200_10_plantes', '410_10_plantes', '600_10_plantes', '800_10_plantes']
# DENSITIES = ['200', '410', '600', '800', '200_10_plantes', '410_10_plantes', '600_10_plantes', '800_10_plantes']

TREATMENTS = ['Plano_diffus', 'Plano_direct', 'Plano_mixte', 'Erect_diffus', 'Erect_direct', 'Erect_mixte', 'Asso_diffus', 'Asso_direct', 'Asso_mixte']

markers = {'Soissons': 's', 'Caphorn': '^', 'Mixture': 'o', 'Plano_diffus': 's', 'Erect_diffus': '^',
                      'Plano_mixte': 's', 'Erect_mixte': '^', 'Plano_direct': 's', 'Erect_direct': '^', 'Plano': 's', 'Erect': '^'}
colors = {'Soissons': 'r', 'Caphorn': 'b', 'Soissons_Mixed': 'r', 'Caphorn_Mixed': 'b', 'Mixture': 'g',
          '200': 'b', '410': 'g', '600': 'y', '800': 'r', 'Erect_diffus': 'b', 'Plano_diffus': 'r',
          '200_10_plantes': 'b', '410_10_plantes': 'g', '600_10_plantes': 'y', '800_10_plantes': 'r'}
lines = {'Soissons': '-', 'Caphorn': '-', 'Soissons_Mixed': '--', 'Caphorn_Mixed': '--', 'Plano_diffus': '-', 'Erect_diffus': '-'}
label_mapping = {'Asso_mixte_1':"Mixed erectophile", 'Asso_mixte_2':"Mixed planophile", 'Erect_mixte': "Erectophile", 'Plano_mixte': "Planophile",
                 'Asso_diffus_1':"Mixed erectophile", 'Asso_diffus_2':"Mixed planophile", 'Erect_diffus': "Erectophile", 'Plano_diffus': "Planophile",
                 'Asso_direct_1':"Mixed erectophile", 'Asso_direct_2':"Mixed planophile", 'Erect_direct': "Erectophile", 'Plano_direct': "Planophile",
                 'Caphorn_Mixed': "Mixed erectophile", 'Soissons_Mixed': "Mixed planophile", 'Soissons': 'Planophile', 'Caphorn': 'Erectophile', 'Mixture': 'Mixture'}

labelx = -0.2  # axes coords
AXES_MAP = {DENSITIES[i]: i for i in range(len(DENSITIES))}
alpha = 0.15


def figure_inclination(num):
    inclination_data = os.path.join(DATA_DIRPATH, 'Reconstruction_interception\S2V\Analyse_inclinaison.csv')
    inclination_df = pd.read_csv(inclination_data)
    inclination_df_grouped = inclination_df.groupby(['Cultivar'])

    plt.figure(figsize=cm2inch(30, 15))
    gs = gridspec.GridSpec(1, 3, wspace=0)
    ax1 = plt.subplot(gs[0, 0])
    ax2 = plt.subplot(gs[0, 1])
    ax3 = plt.subplot(gs[0, 2])

    for cv, group in inclination_df_grouped:
        ax1.plot(group['LAD'], group['Height'], marker=markers[cv], color=colors[cv], linewidth=thickness, label=label_mapping[cv])
        ax2.plot(group['Leaf inclination'], group['Height'], marker=markers[cv],  color=colors[cv], linewidth=thickness, label=label_mapping[cv])
        ax3.plot(group['Plant inclination'], group['Height'], marker=markers[cv],  color=colors[cv], linewidth=thickness, label=label_mapping[cv])

    ax1.set_xlabel('Leaf Area Density (m$^{2}$ m$^{-3}$)', fontsize=fontsize+5)
    ax1.set_ylabel('Height from soil (m)', fontsize=fontsize+5)
    [i.set_linewidth(thickness) for i in ax1.spines.itervalues()]
    ax1.tick_params(width=thickness, length=10, labelsize=fontsize+2)
    ax1.set_xticks(np.arange(0, 0.25, 0.05))
    ax1.set_xticks(ax1.get_xticks()[:-1])
    ax1.set_yticks(np.arange(0, 1., 0.2))
    ax1.text(0.1, 0.95, 'A', transform=ax1.transAxes, fontsize=fontsize+5, verticalalignment='top')
    ax1.legend(prop={'size': 11}, bbox_to_anchor=(1, 0.95), ncol=1, frameon=False, numpoints=1)

    ax2.set_xlabel('Lamina inclination ($^{o}$)', fontsize=fontsize+5)
    [i.set_linewidth(thickness) for i in ax2.spines.itervalues()]
    ax2.tick_params(width=thickness, length=10, labelsize=fontsize+2, labelleft='off')
    ax2.set_xticks(np.arange(0, 105, 15))
    ax2.set_xticks(ax2.get_xticks()[:-1])
    ax2.set_yticks(np.arange(0, 1., 0.2))
    ax2.text(0.1, 0.95, 'B', transform=ax2.transAxes, fontsize=fontsize+5, verticalalignment='top')

    ax3.set_xlabel('Whole culm inclination ($^{o}$)', fontsize=fontsize+5)
    [i.set_linewidth(thickness) for i in ax3.spines.itervalues()]
    ax3.tick_params(width=thickness, length=10, labelsize=fontsize+2, labelleft='off')
    ax3.set_xticks(np.arange(0, 105, 15))
    ax3.set_yticks(np.arange(0, 1., 0.2))
    ax3.text(0.1, 0.95, 'C', transform=ax3.transAxes, fontsize=fontsize+5, verticalalignment='top')

    plt.savefig(os.path.join(GRAPHS_DIRPATH,'Figure_{}_LAD_Inclination.TIFF'.format(num)), dpi=300, format='TIFF', bbox_inches='tight')
    plt.close()


def figure_PARa_dynamic_direct(num):

    PARi_theoritical_df = pd.read_csv('PARi_theorique.csv')

    stands = ['Caphorn', 'Soissons', 'Mixture']
    organ_select = [('blade', 4), ('blade', 3), ('internode', 6)]

    plt.figure(figsize=cm2inch(30,25))
    gs = gridspec.GridSpec(3, 4, hspace=0, wspace=0)

    axes = []
    IC95 = {'density': [], 't': [], 'species': [], 'organ': [], 'metamer': [], 'IC95_up': [], 'IC95_bot': []}

    for density in DENSITIES:
        for stand in stands:
            # CSV files
            PARa_direct_csv = os.path.join(DATA_DIRPATH, 'Reconstruction_interception\outputs', stand, 'sun_DOY[150, 199]_H[5, 19]_1000_plants_Density_{}.csv'.format(density))

            # dfs
            PARa_direct_df_all_data = pd.read_csv(PARa_direct_csv)

            # Select data for 5<=t<=19 and listed organs
            PARa_direct_df = PARa_direct_df_all_data[(5 <= PARa_direct_df_all_data['t']) & (PARa_direct_df_all_data['t'] <= 19)]

            # Computation of theoritical PARa
            ids = []
            PARa_direct_theoritical_list = []
            IC95_list_bot = []
            IC95_list_up = []
            for group_name, group_data in PARa_direct_df.groupby(['t', 'species', 'organ', 'metamer']):
                if (group_name[2], group_name[3]) not in organ_select: continue
                ids.append(group_name)
                t = group_name[0]
                PARi_theoritical = PARi_theoritical_df['PARi theorique'][PARi_theoritical_df['hour'] ==t ].iloc[0]
                PARa_direct_theoritical_list.append(group_data['Eabs'].iloc[0] * PARi_theoritical * 3600 * 1E-6)
                if group_data['IC95'].iloc[0] == '(nan, nan)':
                    IC95_list_bot.append(group_data['Eabs'].iloc[0] * PARi_theoritical * 3600 * 1E-6)
                    IC95_list_up.append(group_data['Eabs'].iloc[0] * PARi_theoritical * 3600 * 1E-6)
                else:
                    IC95_list_bot.append(literal_eval(group_data['IC95'].iloc[0])[0] * PARi_theoritical * 3600 * 1E-6)
                    IC95_list_up.append(literal_eval(group_data['IC95'].iloc[0])[1] * PARi_theoritical * 3600 * 1E-6)

            ids_df = pd.DataFrame(ids, columns=['t', 'species', 'organ', 'metamer'])
            data_df = pd.DataFrame({'PARa theorique': PARa_direct_theoritical_list, 'IC95_bot': IC95_list_bot, 'IC95_up': IC95_list_up})
            PARa_direct_df = pd.concat([ids_df, data_df], axis=1)
            PARa_direct_df.sort_values(['t', 'species', 'organ', 'metamer'], inplace=True)

            # Graphs
            for species in PARa_direct_df['species'].unique():
                if stand == 'Mixture':
                    cv_map = species + '_Mixed'
                else:
                    cv_map = species

                ax1 = plt.subplot(gs[0, AXES_MAP[density]])
                ax2 = plt.subplot(gs[1, AXES_MAP[density]])
                ax3 = plt.subplot(gs[2, AXES_MAP[density]])
                # Blade 4
                blade4_direct = PARa_direct_df[(PARa_direct_df['species'] == species) & (PARa_direct_df['metamer'] == 4) & (PARa_direct_df['organ'] == 'blade')]
                PARa_blade4_direct = blade4_direct['PARa theorique']
                ax1.plot(PARa_direct_df['t'].unique(), PARa_blade4_direct, label=label_mapping[cv_map], marker=None, color=colors[cv_map], linestyle=lines[cv_map], linewidth=thickness)
                ax1.fill_between(PARa_direct_df['t'].unique(), blade4_direct['IC95_bot'], blade4_direct['IC95_up'], alpha=alpha, color=colors[cv_map], linewidth=0.0)
                # Blade 3
                blade3_direct = PARa_direct_df[(PARa_direct_df['species'] == species) & (PARa_direct_df['metamer'] == 3) & (PARa_direct_df['organ'] == 'blade')]
                PARa_blade3_direct = blade3_direct['PARa theorique']
                ax2.plot(PARa_direct_df['t'].unique(), PARa_blade3_direct, label=label_mapping[cv_map], marker=None, color=colors[cv_map], linestyle=lines[cv_map], linewidth=thickness)
                ax2.fill_between(PARa_direct_df['t'].unique(), blade3_direct['IC95_bot'], blade3_direct['IC95_up'], alpha=alpha, color=colors[cv_map], linewidth=0.0)
                # Chaff
                chaff_direct = PARa_direct_df[(PARa_direct_df['species'] == species) & (PARa_direct_df['metamer']==6) & (PARa_direct_df['organ'] == 'internode')]
                PARa_chaff_direct = chaff_direct['PARa theorique']
                ax3.plot(PARa_direct_df['t'].unique(), PARa_chaff_direct, label=label_mapping[cv_map], marker=None, color=colors[cv_map], linestyle=lines[cv_map], linewidth=thickness)
                ax3.fill_between(PARa_direct_df['t'].unique(), chaff_direct['IC95_bot'], chaff_direct['IC95_up'], alpha=alpha, color=colors[cv_map], linewidth=0.0)

                # Store axes
                for ax in (ax1, ax2, ax3):
                    axes.append(ax)


    IC95_df = pd.DataFrame(IC95)
    mean_IC95 = {}
    for density in DENSITIES:
        mean_IC95[density] = {}
        for organ_id in organ_select:
            mean_IC95[density][organ_id] = {'t': [], 'IC95': []}

    for group_name, group_data in IC95_df.groupby(['density', 't', 'organ', 'metamer']):
        density = group_name[0]
        organ_id = (group_name[2], group_name[3])
        mean_IC95[density][organ_id]['t'].append(group_name[1])
        mean_IC95[density][organ_id]['IC95'].append(group_data['IC95'].mean())

    for ax in axes:
        [i.set_linewidth(thickness) for i in ax.spines.itervalues()]
        ax.tick_params(width=thickness, length=10, labelbottom='off', labelleft='off', top='off', right='off')
        ax.set_yticks(np.arange(0, 4, 1))
        ax.set_xticks([5, 9, 13, 17])

    # Labels
    ax_label_left = (plt.subplot(gs[0, 0]), plt.subplot(gs[1, 0]), plt.subplot(gs[2, 0]))
    for ax in ax_label_left:
        ax.tick_params(labelsize=fontsize+2, labelleft='on')
        if ax != plt.subplot(gs[0, 0]):
            ax.set_yticks(ax.get_yticks()[:-1])
    ax_label_bottom = (plt.subplot(gs[2, 0]), plt.subplot(gs[2, 1]), plt.subplot(gs[2, 2]), plt.subplot(gs[2, 3]))
    for ax in ax_label_bottom:
        ax.tick_params(labelsize=fontsize+2, labelbottom='on')

    plt.subplot(gs[1, 0]).set_ylabel('PAR absorbed (mol m$^{-2}$ h$^{-1}$)', fontsize=fontsize+5)
    plt.text(30, -3.5, 'Time (hour UTC)', fontsize=fontsize+5)

    plt.text(8, 6.2, 'Density 200', fontsize=fontsize+5)
    plt.text(24, 6.2, 'Density 410', fontsize=fontsize+5)
    plt.text(40, 6.2, 'Density 600', fontsize=fontsize+5)
    plt.text(56, 6.2, 'Density 800', fontsize=fontsize+5)

    plt.subplot(gs[0, 1]).text(17, 2.5, 'Lamina n', fontsize=fontsize+5, backgroundcolor='w', verticalalignment='top')
    plt.subplot(gs[0, 1]).set_zorder(1)
    plt.subplot(gs[1, 1]).text(17, 2.5, 'Lamina n-1', fontsize=fontsize+5, backgroundcolor='w', verticalalignment='top')
    plt.subplot(gs[1, 1]).set_zorder(1)
    plt.subplot(gs[2, 1]).text(18, 2.5, 'Chaff', fontsize=fontsize+5, backgroundcolor='w', verticalalignment='top')
    plt.subplot(gs[2, 1]).set_zorder(1)

    plt.savefig(os.path.join(GRAPHS_DIRPATH, 'Figure_{}_PAR_direct.TIFF'.format(num)), dpi=300, format='TIFF', bbox_inches='tight')
    plt.close()


def figure_PARa_dynamic_diffuse(num):
    PARi_theoritical_df = pd.read_csv('PARi_theorique.csv')

    stands = ['Caphorn', 'Soissons', 'Mixture']
    organ_select = [('blade', 4), ('blade', 3), ('internode', 6)]


    plt.figure(figsize=cm2inch(30, 25))
    gs = gridspec.GridSpec(3, 4, hspace=0, wspace=0)

    axes = []
    IC95 = {'density': [], 't': [], 'species': [], 'organ': [], 'metamer': [], 'IC95_up': [], 'IC95_bot': []}

    for density in DENSITIES:
        for stand in stands:
            # CSV files
            PARa_diffus_csv = os.path.join(DATA_DIRPATH, 'Reconstruction_interception\outputs', stand, 'soc_a4z5_1000_plants_Density_{}.csv'.format(density))

            # dfs
            PARa_diffus_df_all_data = pd.read_csv(PARa_diffus_csv)

            # Select data for 5<=t<=19 and listed organs
            PARa_diffus_df = PARa_diffus_df_all_data[(5 <= PARa_diffus_df_all_data['t']) & (PARa_diffus_df_all_data['t'] <= 19)]

            # Computation of theoritical PARa
            ids = []
            PARa_diffus_theoritical_list = []
            IC95_list_bot = []
            IC95_list_up = []
            for group_name, group_data in PARa_diffus_df.groupby(['t', 'species', 'organ', 'metamer']):
                if (group_name[2], group_name[3]) not in organ_select: continue
                ids.append(group_name)
                t = group_name[0]
                PARi_theoritical = PARi_theoritical_df['PARi theorique'][PARi_theoritical_df['hour']==t].iloc[0]
                PARa_diffus_theoritical_list.append(group_data['Eabs'].iloc[0] * PARi_theoritical * 3600 * 1E-6)
                if group_data['IC95'].iloc[0] == '(nan, nan)':
                    IC95_list_bot.append(group_data['Eabs'].iloc[0] * PARi_theoritical * 3600 * 1E-6)
                    IC95_list_up.append(group_data['Eabs'].iloc[0] * PARi_theoritical * 3600 * 1E-6)
                else:
                    IC95_list_bot.append(literal_eval(group_data['IC95'].iloc[0])[0] * PARi_theoritical * 3600 * 1E-6)
                    IC95_list_up.append(literal_eval(group_data['IC95'].iloc[0])[1] * PARi_theoritical * 3600 * 1E-6)

            ids_df = pd.DataFrame(ids, columns=['t', 'species', 'organ', 'metamer'])
            data_df = pd.DataFrame({'PARa theorique': PARa_diffus_theoritical_list, 'IC95_bot': IC95_list_bot, 'IC95_up': IC95_list_up})
            PARa_diffus_df = pd.concat([ids_df, data_df], axis=1)
            PARa_diffus_df.sort_values(['t', 'species', 'organ', 'metamer'], inplace=True)

            # Graphs
            for species in PARa_diffus_df['species'].unique():
                if stand == 'Mixture':
                    cv_map = species + '_Mixed'
                else:
                    cv_map = species

                ax1 = plt.subplot(gs[0, AXES_MAP[density]])
                ax2 = plt.subplot(gs[1, AXES_MAP[density]])
                ax3 = plt.subplot(gs[2, AXES_MAP[density]])
                ## Blade 4
                blade4_diffus = PARa_diffus_df[(PARa_diffus_df['species'] == species) & (PARa_diffus_df['metamer'] == 4) & (PARa_diffus_df['organ'] == 'blade')]
                PARa_blade4_diffus = blade4_diffus['PARa theorique']
                ax1.plot(PARa_diffus_df['t'].unique(), PARa_blade4_diffus, label=label_mapping[cv_map], marker=None, color=colors[cv_map], linestyle=lines[cv_map], linewidth=thickness)
                ax1.fill_between(PARa_diffus_df['t'].unique(), blade4_diffus['IC95_bot'], blade4_diffus['IC95_up'], alpha=alpha, color=colors[cv_map], linewidth=0.0)
                ## Blade 3
                blade3_diffus = PARa_diffus_df[(PARa_diffus_df['species'] == species) & (PARa_diffus_df['metamer'] == 3) & (PARa_diffus_df['organ'] == 'blade')]
                PARa_blade3_diffus  = blade3_diffus['PARa theorique']
                ax2.plot(PARa_diffus_df['t'].unique(), PARa_blade3_diffus, label=label_mapping[cv_map], marker=None, color=colors[cv_map], linestyle=lines[cv_map], linewidth=thickness)
                ax2.fill_between(PARa_diffus_df['t'].unique(), blade3_diffus['IC95_bot'], blade3_diffus['IC95_up'], alpha=alpha, color=colors[cv_map], linewidth=0.0)
                ## Chaff
                chaff_diffus = PARa_diffus_df[(PARa_diffus_df['species'] == species) & (PARa_diffus_df['metamer'] == 6) & (PARa_diffus_df['organ'] == 'internode')]
                PARa_chaff_diffus = chaff_diffus['PARa theorique']
                ax3.plot(PARa_diffus_df['t'].unique(), PARa_chaff_diffus, label=label_mapping[cv_map], marker=None, color=colors[cv_map], linestyle=lines[cv_map], linewidth=thickness)
                ax3.fill_between(PARa_diffus_df['t'].unique(), chaff_diffus['IC95_bot'], chaff_diffus['IC95_up'], alpha=alpha, color=colors[cv_map], linewidth=0.0)

                # Store axes
                for ax in (ax1, ax2, ax3):
                    axes.append(ax)

    IC95_df = pd.DataFrame(IC95)
    mean_IC95 = {}
    for density in DENSITIES:
        mean_IC95[density] = {}
        for organ_id in organ_select:
            mean_IC95[density][organ_id] = {'t': [], 'IC95': []}

    for group_name, group_data in IC95_df.groupby(['density', 't', 'organ', 'metamer']):
        density = group_name[0]
        organ_id = (group_name[2], group_name[3])
        mean_IC95[density][organ_id]['t'].append(group_name[1])
        mean_IC95[density][organ_id]['IC95'].append(group_data['IC95'].mean())

    for ax in axes:
        [i.set_linewidth(thickness) for i in ax.spines.itervalues()]
        ax.tick_params(width=thickness, length=10, labelbottom='off', labelleft='off', top='off', right='off')
        ax.set_yticks(np.arange(0, 4, 1))
        ax.set_xticks([5, 9, 13, 17])

    # Labels
    ax_label_left = (plt.subplot(gs[0, 0]), plt.subplot(gs[1, 0]), plt.subplot(gs[2, 0]))
    for ax in ax_label_left:
        ax.tick_params(labelsize=fontsize+2, labelleft='on')
        if ax !=plt.subplot(gs[0, 0]):
            ax.set_yticks(ax.get_yticks()[:-1])
    ax_label_bottom = (plt.subplot(gs[2, 0]), plt.subplot(gs[2, 1]), plt.subplot(gs[2, 2]), plt.subplot(gs[2, 3]))
    for ax in ax_label_bottom:
        ax.tick_params(labelsize=fontsize+2, labelbottom='on')

    plt.subplot(gs[1, 0]).set_ylabel('PAR absorbed (mol m$^{-2}$ h$^{-1}$)', fontsize=fontsize+5)
    plt.text(30, -3.5, 'Time (hour UTC)', fontsize=fontsize+5)

    plt.text(8, 6.2, 'Density 200', fontsize=fontsize+5)
    plt.text(24, 6.2, 'Density 410', fontsize=fontsize+5)
    plt.text(40, 6.2, 'Density 600', fontsize=fontsize+5)
    plt.text(56, 6.2, 'Density 800', fontsize=fontsize+5)

    plt.subplot(gs[0, 1]).text(17, 2.5, 'Lamina n', fontsize=fontsize+5, backgroundcolor='w', verticalalignment='top')
    plt.subplot(gs[0, 1]).set_zorder(1)
    plt.subplot(gs[1, 1]).text(17, 2.5, 'Lamina n-1', fontsize=fontsize+5, backgroundcolor='w', verticalalignment='top')
    plt.subplot(gs[1, 1]).set_zorder(1)
    plt.subplot(gs[2, 1]).text(18, 2.5, 'Chaff', fontsize=fontsize+5, backgroundcolor='w', verticalalignment='top')
    plt.subplot(gs[2, 1]).set_zorder(1)

    plt.savefig(os.path.join(GRAPHS_DIRPATH, 'Figure_{}_PAR_diffus.TIFF'.format(num)), dpi=300, format='TIFF', bbox_inches='tight')
    plt.close()


def figure_profile_PARa_10_plants(num):

    plt.figure(figsize=cm2inch(30, 15))
    gs = gridspec.GridSpec(3, 3, hspace=0, wspace=0)

    ax1 = plt.subplot(gs[0, 0])
    ax2 = plt.subplot(gs[1, 0])
    ax3 = plt.subplot(gs[2, 0])
    ax4 = plt.subplot(gs[0, 1])
    ax5 = plt.subplot(gs[1, 1])
    ax6 = plt.subplot(gs[2, 1])
    ax7 = plt.subplot(gs[0, 2])
    ax8 = plt.subplot(gs[1, 2])
    ax9 = plt.subplot(gs[2, 2])


    TREATMENTS = {'Plano_diffus': ax1, 'Erect_diffus': ax2, 'Asso_diffus': ax3,
                  'Plano_mixte': ax4, 'Erect_mixte': ax5, 'Asso_mixte': ax6,
                  'Plano_direct': ax7, 'Erect_direct': ax8, 'Asso_direct': ax9}

    for density in DENSITIES:
        for treatement, ax in TREATMENTS.iteritems():
            data_dirpath = os.path.join(DATA_DIRPATH, 'Densite_{}'.format(density), treatement, 'outputs', 'elements_states.csv')
            data_df = pd.read_csv(data_dirpath)
            data_df.dropna(inplace=True)

            if treatement not in ('Asso_diffus', 'Asso_mixte', 'Asso_direct'):
                daily_Eabs_organ_all_plants = {1: [], 2: [], 3: [], 4: []}
                mean_height_organ_all_plants = {1: [], 2: [], 3: [], 4: []}

                data_df_selection = data_df[(5 <= data_df['t']) & (data_df['t'] <= 19) & (data_df['element'] == 'LeafElement1')].groupby(['plant', 'metamer'])
                for element_id, data in data_df_selection:
                    daily_Eabs_organ_all_plants[element_id[1]].append(data['PARa'].sum() * 3600 * 1E-6)
                    mean_height_organ_all_plants[element_id[1]].append(data['height'].mean())

                daily_Eabs_organ_mean = []
                daily_Eabs_organ_IC95 = {'bot': [], 'up': []}
                height_organ_mean = []
                for metamer, data in daily_Eabs_organ_all_plants.iteritems():
                    mean = np.mean(data)
                    daily_Eabs_organ_mean.append(mean)
                    IC95 = stats.t.interval(0.95, len(data)-1, loc=mean, scale=stats.sem(data))
                    daily_Eabs_organ_IC95['bot'].append(IC95[0]), daily_Eabs_organ_IC95['up'].append(IC95[1])

                for metamer, data in mean_height_organ_all_plants.iteritems():
                    height_organ_mean.append(np.mean(data))

                # Plot data
                ax.plot(height_organ_mean, daily_Eabs_organ_mean, color=colors[density], linewidth=thickness, marker=markers[treatement], label=density, markersize=8)
                ax.fill_between(height_organ_mean, daily_Eabs_organ_IC95['bot'], daily_Eabs_organ_IC95['up'], alpha=alpha, color=colors[density], linewidth=0.0)

            else:
                daily_Eabs_organ_all_plants_erectophile = {1: [], 2: [], 3: [], 4: []}
                mean_height_organ_all_plants_erectophile = {1: [], 2: [], 3: [], 4: []}
                daily_Eabs_organ_all_plants_planophile = {1: [], 2: [], 3: [], 4: []}
                mean_height_organ_all_plants_planophile = {1: [], 2: [], 3: [], 4: []}

                data_df_selection_erectophile = data_df[(5 <= data_df['t']) & (data_df['t'] <= 19) & (1 <= data_df['plant']) & (data_df['plant'] <= 5) & (data_df['element'] == 'LeafElement1')].groupby(['plant', 'metamer'])
                data_df_selection_planophile = data_df[(5 <= data_df['t']) & (data_df['t'] <= 19) & (6 <= data_df['plant']) & (data_df['plant'] <= 10) & (data_df['element'] == 'LeafElement1')].groupby(['plant', 'metamer'])

                for element_id, data in data_df_selection_erectophile:
                    daily_Eabs_organ_all_plants_erectophile[element_id[1]].append(data['PARa'].sum() * 3600 * 1E-6)
                    mean_height_organ_all_plants_erectophile[element_id[1]].append(data['height'].mean())

                for element_id, data in data_df_selection_planophile:
                    daily_Eabs_organ_all_plants_planophile[element_id[1]].append(data['PARa'].sum() * 3600 * 1E-6)
                    mean_height_organ_all_plants_planophile[element_id[1]].append(data['height'].mean())

                daily_Eabs_organ_mean_erectophile = []
                daily_Eabs_organ_mean_planophile = []
                daily_Eabs_organ_IC95_erectophile = {'bot': [], 'up': []}
                daily_Eabs_organ_IC95_planophile = {'bot': [], 'up': []}
                height_organ_mean_erectophile = []
                height_organ_mean_planophile = []

                for metamer, data in daily_Eabs_organ_all_plants_erectophile.iteritems():
                    daily_Eabs_organ_mean_erectophile.append(np.mean(data))
                    IC95 = stats.t.interval(0.95, len(data)-1, loc=np.mean(data), scale=stats.sem(data))
                    daily_Eabs_organ_IC95_erectophile['bot'].append(IC95[0]), daily_Eabs_organ_IC95_erectophile['up'].append(IC95[1])

                for metamer, data in mean_height_organ_all_plants_erectophile.iteritems():
                    height_organ_mean_erectophile.append(np.mean(data))

                for metamer, data in daily_Eabs_organ_all_plants_planophile.iteritems():
                    daily_Eabs_organ_mean_planophile.append(np.mean(data))
                    IC95 = stats.t.interval(0.95, len(data)-1, loc=np.mean(data), scale=stats.sem(data))
                    daily_Eabs_organ_IC95_planophile['bot'].append(IC95[0]), daily_Eabs_organ_IC95_planophile['up'].append(IC95[1])

                for metamer, data in mean_height_organ_all_plants_planophile.iteritems():
                    height_organ_mean_planophile.append(np.mean(data))

                ax.plot(height_organ_mean_erectophile, daily_Eabs_organ_mean_erectophile, color=colors[density], linewidth=thickness, marker='^', label=density, markersize=8)
                ax.fill_between(height_organ_mean_erectophile, daily_Eabs_organ_IC95_erectophile['bot'], daily_Eabs_organ_IC95_erectophile['up'], alpha=alpha,
                                color=colors[density], linewidth=0.0)
                ax.plot(height_organ_mean_planophile, daily_Eabs_organ_mean_planophile, color=colors[density], linewidth=thickness, marker='s', label=density, markersize=8)
                ax.fill_between(height_organ_mean_planophile, daily_Eabs_organ_IC95_planophile['bot'], daily_Eabs_organ_IC95_planophile['up'], alpha=alpha,
                                color=colors[density], linewidth=0.0)

    # Formatting
    for ax in TREATMENTS.itervalues():
        ax.set_yticks(np.arange(0, 30, 5))
        ax.set_xticks(np.arange(0, 1, 0.2))
        [i.set_linewidth(thickness) for i in ax.spines.itervalues()]
        ax.tick_params(width=thickness, labelsize=fontsize)

        if ax in (ax1, ax2, ax4, ax5, ax7, ax8):
            ax.tick_params(labelbottom='off')
        else:
            ax.set_xlabel(u'Height (m)', fontsize=fontsize+2)

        if ax in (ax1, ax2, ax3):
            ax.set_ylabel(u'PARa (mol m$^{-2}$ day$^{-1}$)', fontsize=fontsize + 2)
    ax2.set_yticks(ax.get_yticks()[:-1])
    ax3.set_yticks(ax.get_yticks()[:-1])

    ax1.text(0.05, 0.95, 'A - Planophile', transform=ax1.transAxes, fontsize=fontsize+2, verticalalignment='top')
    ax2.text(0.05, 0.95, 'B - Erectophile', transform=ax2.transAxes, fontsize=fontsize+2, verticalalignment='top')
    ax3.text(0.05, 0.95, 'C - Mixture', transform=ax3.transAxes, fontsize=fontsize+2, verticalalignment='top')

    handles, labels = ax1.get_legend_handles_labels()
    plt.legend(handles, labels, prop={'size': 12}, bbox_to_anchor=(0.5, 3.15), loc='upper center', ncol=4, frameon=True, markerscale=0, borderaxespad=0.)
    plt.savefig(os.path.join(GRAPHS_DIRPATH, 'Figure_{}_Profiles_PARa.TIFF'.format(num)), dpi=300, format='TIFF', bbox_inches='tight')
    plt.close()

def figure_profile_PARa_1000_plants(num):

    plt.figure(figsize=cm2inch(15, 30))
    gs = gridspec.GridSpec(3, 1, hspace=0, wspace=0)

    ax1 = plt.subplot(gs[0, 0])
    ax2 = plt.subplot(gs[1, 0])
    ax3 = plt.subplot(gs[2, 0])


    TREATMENTS = {'Plano_direct': ax1, 'Erect_direct': ax2, 'Asso_direct': ax3}

    for density in DENSITIES:
        for treatement, ax in TREATMENTS.iteritems():
            data_dirpath = os.path.join(DATA_DIRPATH, 'Densite_{}'.format(density), treatement, 'outputs', 'elements_states.csv')
            data_df = pd.read_csv(data_dirpath)
            data_df.dropna(inplace=True)

            if treatement not in ('Asso_direct'):
                daily_Eabs_organ_all_plants = {1: [], 2: [], 3: [], 4: []}
                mean_height_organ_all_plants = {1: [], 2: [], 3: [], 4: []}

                data_df_selection = data_df[(5 <= data_df['t']) & (data_df['t'] <= 19) & (data_df['element'] == 'LeafElement1')].groupby(['plant', 'metamer'])
                for element_id, data in data_df_selection:
                    daily_Eabs_organ_all_plants[element_id[1]].append(data['PARa'].sum() * 3600 * 1E-6)
                    mean_height_organ_all_plants[element_id[1]].append(data['height'].mean())

                daily_Eabs_organ_mean = []
                daily_Eabs_organ_IC95 = {'bot': [], 'up': []}
                height_organ_mean = []
                for metamer, data in daily_Eabs_organ_all_plants.iteritems():
                    mean = np.mean(data)
                    daily_Eabs_organ_mean.append(mean)
                    IC95 = stats.t.interval(0.95, len(data)-1, loc=mean, scale=stats.sem(data))
                    daily_Eabs_organ_IC95['bot'].append(IC95[0]), daily_Eabs_organ_IC95['up'].append(IC95[1])

                for metamer, data in mean_height_organ_all_plants.iteritems():
                    height_organ_mean.append(np.mean(data))

                # Plot data
                ax.plot(daily_Eabs_organ_mean, height_organ_mean, color=colors[density], linewidth=thickness, marker=markers[treatement], label=density, markersize=8)

            else:
                daily_Eabs_organ_all_plants_erectophile = {1: [], 2: [], 3: [], 4: []}
                mean_height_organ_all_plants_erectophile = {1: [], 2: [], 3: [], 4: []}
                daily_Eabs_organ_all_plants_planophile = {1: [], 2: [], 3: [], 4: []}
                mean_height_organ_all_plants_planophile = {1: [], 2: [], 3: [], 4: []}

                data_df_selection_erectophile = data_df[(5 <= data_df['t']) & (data_df['t'] <= 19) & (data_df['plant'] == 1) & (data_df['element'] == 'LeafElement1')].groupby(['plant', 'metamer'])
                data_df_selection_planophile = data_df[(5 <= data_df['t']) & (data_df['t'] <= 19) & (data_df['plant'] == 2) & (data_df['element'] == 'LeafElement1')].groupby(['plant', 'metamer'])

                for element_id, data in data_df_selection_erectophile:
                    daily_Eabs_organ_all_plants_erectophile[element_id[1]].append(data['PARa'].sum() * 3600 * 1E-6)
                    mean_height_organ_all_plants_erectophile[element_id[1]].append(data['height'].mean())

                for element_id, data in data_df_selection_planophile:
                    daily_Eabs_organ_all_plants_planophile[element_id[1]].append(data['PARa'].sum() * 3600 * 1E-6)
                    mean_height_organ_all_plants_planophile[element_id[1]].append(data['height'].mean())

                daily_Eabs_organ_mean_erectophile = []
                daily_Eabs_organ_mean_planophile = []
                daily_Eabs_organ_IC95_erectophile = {'bot': [], 'up': []}
                daily_Eabs_organ_IC95_planophile = {'bot': [], 'up': []}
                height_organ_mean_erectophile = []
                height_organ_mean_planophile = []

                for metamer, data in daily_Eabs_organ_all_plants_erectophile.iteritems():
                    daily_Eabs_organ_mean_erectophile.append(np.mean(data))
                    IC95 = stats.t.interval(0.95, len(data)-1, loc=np.mean(data), scale=stats.sem(data))
                    daily_Eabs_organ_IC95_erectophile['bot'].append(IC95[0]), daily_Eabs_organ_IC95_erectophile['up'].append(IC95[1])

                for metamer, data in mean_height_organ_all_plants_erectophile.iteritems():
                    height_organ_mean_erectophile.append(np.mean(data))

                for metamer, data in daily_Eabs_organ_all_plants_planophile.iteritems():
                    daily_Eabs_organ_mean_planophile.append(np.mean(data))
                    IC95 = stats.t.interval(0.95, len(data)-1, loc=np.mean(data), scale=stats.sem(data))
                    daily_Eabs_organ_IC95_planophile['bot'].append(IC95[0]), daily_Eabs_organ_IC95_planophile['up'].append(IC95[1])

                for metamer, data in mean_height_organ_all_plants_planophile.iteritems():
                    height_organ_mean_planophile.append(np.mean(data))

                ax.plot(daily_Eabs_organ_mean_erectophile, height_organ_mean_erectophile, color=colors[density], linewidth=thickness, marker='^', label=density, markersize=8)
                ax.plot(daily_Eabs_organ_mean_planophile, height_organ_mean_planophile, color=colors[density], linewidth=thickness, marker='s', label=density, markersize=8)

    # Formatting
    for ax in TREATMENTS.itervalues():
        ax.set_xticks(np.arange(0, 30, 5))
        ax.set_yticks(np.arange(0, 1, 0.2))
        [i.set_linewidth(thickness) for i in ax.spines.itervalues()]
        ax.tick_params(width=thickness, labelsize=fontsize)

        if ax in (ax1, ax2):
            ax.tick_params(labelbottom='off')

        ax.set_ylabel(u'Height (m)', fontsize=fontsize+2)

        if ax in (ax1, ax2, ax3):
            ax.set_xlabel(u'PARa (mol m$^{-2}$ day$^{-1}$)', fontsize=fontsize + 2)
    ax2.set_yticks(ax.get_yticks()[:-1])
    ax3.set_yticks(ax.get_yticks()[:-1])

    ax1.text(0.05, 0.95, 'A - Planophile', transform=ax1.transAxes, fontsize=fontsize+2, verticalalignment='top')
    ax2.text(0.05, 0.95, 'B - Erectophile', transform=ax2.transAxes, fontsize=fontsize+2, verticalalignment='top')
    ax3.text(0.05, 0.95, 'C - Mixture', transform=ax3.transAxes, fontsize=fontsize+2, verticalalignment='top')

    handles, labels = ax1.get_legend_handles_labels()
    plt.legend(handles, labels, prop={'size': 12}, bbox_to_anchor=(0.5, 3.15), loc='upper center', ncol=4, frameon=True, markerscale=0, borderaxespad=0.)
    plt.savefig(os.path.join(GRAPHS_DIRPATH, 'Figure_{}_Profiles_PARa.TIFF'.format(num)), dpi=300, format='TIFF', bbox_inches='tight')
    plt.close()


def figure_PARa_Cass(num):

    PARa_Cass = {'Density': [], 'Leaf posture': [], 'Stand': [], 'Sky': [], 'plant': [], 'Total_PARa': [], 'Photosynthesis': [], 'Net_Photosynthesis': [], 'Simulation_type': []}

    for density in DENSITIES:
        for treatment in TREATMENTS:
            data_dirpath = os.path.join(DATA_DIRPATH, 'Densite_{}'.format(density), treatment, 'outputs', 'elements_states.csv')
            data_df = pd.read_csv(data_dirpath)
            data_df.dropna(inplace=True)
            data_df['Total_PARa'] = data_df['PARa'] * data_df['green_area'] * 3600 * 1E-6

            if treatment not in ('Asso_diffus', 'Asso_mixte', 'Asso_direct'):
                leaf_posture, sky = treatment.split('_')
                data = data_df.groupby('plant')['Total_PARa', 'Photosynthesis', 'R_Nnit_red', 'R_phloem_loading', 'R_residual']
                for plant, values in data:
                    PARa_Cass['Density'].append(density)
                    PARa_Cass['Leaf posture'].append(leaf_posture)
                    PARa_Cass['Stand'].append('pure')
                    PARa_Cass['Sky'].append(sky)
                    PARa_Cass['plant'].append(plant)
                    PARa_Cass['Total_PARa'].append(values['Total_PARa'].sum())
                    PARa_Cass['Photosynthesis'].append(values['Photosynthesis'].sum() * 1E-3)
                    PARa_Cass['Net_Photosynthesis'].append((values['Photosynthesis'].sum() - values['R_Nnit_red'].sum() -
                                                            values['R_phloem_loading'].sum() - values['R_residual'].sum()) * 1E-3)
                    if density not in ('200', '410', '600', '800'):
                        PARa_Cass['Simulation_type'].append('10_plants')
                    else:
                        PARa_Cass['Simulation_type'].append('mean_plant')

                    del plant, values
                    gc.collect()

            else:
                stand, sky = treatment.split('_')
                if density not in ('200', '410', '600', '800'):
                    data_erect = data_df[(1 <= data_df['plant']) & (data_df['plant'] <= 5)].groupby('plant')['Total_PARa', 'Photosynthesis']
                    data_plano = data_df[(6 <= data_df['plant']) & (data_df['plant'] <= 10)].groupby('plant')['Total_PARa', 'Photosynthesis']
                else:
                    data_erect = data_df[(data_df['plant'] == 1)].groupby('plant')['Total_PARa', 'Photosynthesis']
                    data_plano = data_df[(data_df['plant'] == 2)].groupby('plant')['Total_PARa', 'Photosynthesis']

                for plant, values in data_erect:
                    PARa_Cass['Density'].append(density)
                    PARa_Cass['Stand'].append(stand)
                    PARa_Cass['Leaf posture'].append('Erect')
                    PARa_Cass['Sky'].append(sky)
                    PARa_Cass['plant'].append(plant)
                    PARa_Cass['Total_PARa'].append(values['Total_PARa'].sum())
                    PARa_Cass['Photosynthesis'].append(values['Photosynthesis'].sum() * 1E-3)
                    PARa_Cass['Net_Photosynthesis'].append((values['Photosynthesis'].sum() - values['R_Nnit_red'].sum() -
                                                            values['R_phloem_loading'].sum() - values['R_residual'].sum()) * 1E-3)
                    if density not in ('200', '410', '600', '800'):
                        PARa_Cass['Simulation_type'].append('10_plants')
                    else:
                        PARa_Cass['Simulation_type'].append('mean_plant')

                for plant, values in data_plano:
                    PARa_Cass['Density'].append(density)
                    PARa_Cass['Stand'].append(stand)
                    PARa_Cass['Leaf posture'].append('Erect')
                    PARa_Cass['Sky'].append(sky)
                    PARa_Cass['plant'].append(plant)
                    PARa_Cass['Total_PARa'].append(values['Total_PARa'].sum())
                    PARa_Cass['Photosynthesis'].append(values['Photosynthesis'].sum() * 1E-3)
                    PARa_Cass['Net_Photosynthesis'].append((values['Photosynthesis'].sum() - values['R_Nnit_red'].sum() -
                                                            values['R_phloem_loading'].sum() - values['R_residual'].sum()) * 1E-3)
                    if density not in ('200', '410', '600', '800'):
                        PARa_Cass['Simulation_type'].append('10_plants')
                    else:
                        PARa_Cass['Simulation_type'].append('mean_plant')

                del plant, values
                gc.collect()

    PARa_Cass_df = pd.DataFrame(PARa_Cass)

    data_df_grouped = PARa_Cass_df.groupby(['Density', 'Leaf posture'])

    fig, ax = plt.subplots()

    for group_name, group in data_df_grouped:
        density = group_name[0]
        leaf_inclination = group_name[1]
        if density in ('200', '410', '600', '800'):
            ax.plot(group['Total_PARa'], group['Net_Photosynthesis'], marker=markers[leaf_inclination], markersize=10, color=colors[density], label=density, linewidth=0, zorder=10)
        else:
            ax.plot(group['Total_PARa'], group['Net_Photosynthesis'], marker=markers[leaf_inclination], markersize=5, markerfacecolor='None',
                    markeredgecolor=colors[density], label=density, linewidth=0, zorder=1)

    # linear regression mean_plant
    data = PARa_Cass_df[PARa_Cass_df['Simulation_type'] == 'mean_plant']
    slope, intercept, r_value, p_value, std_err = stats.linregress(data['Total_PARa'], data['Net_Photosynthesis'])
    line = [slope * data['Total_PARa'].min() + intercept, slope * data['Total_PARa'].max() + intercept]
    ax.plot([data['Total_PARa'].min(), data['Total_PARa'].max()], line, 'k-', zorder=2)

    sum_difference = 0
    for i in data.index:
        sum_difference += (data['Net_Photosynthesis'][i] * 1E-3 - (slope * 1E-3 * data['Total_PARa'][i] + intercept * 1E-3))**2

    RMSE = np.sqrt(sum_difference / len(data.index))

    ax.text(0.4, 0.77, '- Average culms \ny = {}x + {} \nR$^2$ = {}, RMSE = {} mol'.format(round(slope, 3), round(intercept, 3), round(r_value**2, 3), round(RMSE, 3)), transform=ax.transAxes, fontsize=fontsize, verticalalignment='top')

    # linear regression 10_plants
    data = PARa_Cass_df[PARa_Cass_df['Simulation_type'] == '10_plants']
    slope, intercept, r_value, p_value, std_err = stats.linregress(data['Total_PARa'], data['Net_Photosynthesis'])
    line = [slope * data['Total_PARa'].min() + intercept, slope * data['Total_PARa'].max() + intercept]
    ax.plot([data['Total_PARa'].min(), data['Total_PARa'].max()], line, 'k--', zorder=2)

    sum_difference = 0
    for i in data.index:
        sum_difference += (data['Net_Photosynthesis'][i] * 1E-3 - (slope * 1E-3 * data['Total_PARa'][i] + intercept * 1E-3))**2

    RMSE = np.sqrt(sum_difference / len(data.index))

    ax.text(0.6, 0.5, '-- 10 subsampled culms \ny = {}x + {} \nR$^2$ = {}, RMSE = {} mol'.format(round(slope, 3), round(intercept, 3), round(r_value**2, 3), round(RMSE, 3)), transform=ax.transAxes, fontsize=fontsize, verticalalignment='top')

    # Format
    [i.set_linewidth(thickness) for i in ax.spines.itervalues()]
    ax.set_xticks(np.arange(0, 12, 2))
    ax.set_yticks(np.arange(0, 350, 50))
    ax.set_xlabel('Total absorbed PAR (mol)', fontsize=fontsize+5)
    ax.set_ylabel('Net C assimilated  (mmol)', fontsize=fontsize+5)

    plt.savefig(os.path.join(GRAPHS_DIRPATH, 'Figure_{}_PARa_Cass.TIFF'.format(num)), dpi=300, format='TIFF', bbox_inches='tight')
    plt.close()

def figure_green_area_proteins(num):
    plt.figure(figsize=cm2inch(30, 15))
    gs = gridspec.GridSpec(2, 4, hspace=0, wspace=0)

    for density in DENSITIES:
        for treatment in TREATMENTS:
            if treatment not in ('Plano_diffus', 'Erect_diffus', 'Asso_diffus'): continue
            # Photosynthetic elements
            ph_elements_path = os.path.join(DATA_DIRPATH, 'Densite_{}'.format(density), treatment, 'outputs\elements_states.csv')
            ph_elements = pd.read_csv(ph_elements_path)

            # 1.1) Lamina area
            ax1 = plt.subplot(gs[0, AXES_MAP[density]])

            if treatment == 'Asso_diffus':
                if density in ('200', '410', '600', '800'):
                    data_erect = (ph_elements[(ph_elements['plant'] == 1)
                                              & (ph_elements['organ'] == 'blade')].groupby(['t', 'plant'])[
                                      'green_area'].aggregate(np.sum) * 10000).groupby(['t'])
                    data_plano = (ph_elements[(ph_elements['plant'] == 2)
                                              & (ph_elements['organ'] == 'blade')].groupby(['t', 'plant'])[
                                      'green_area'].aggregate(np.sum) * 10000).groupby(['t'])
                else:
                    data_erect = (ph_elements[(1 <= ph_elements['plant']) & (ph_elements['plant'] <= 5)
                                                    & (ph_elements['organ'] == 'blade')].groupby(['t', 'plant'])['green_area'].aggregate(np.sum) * 10000).groupby(['t'])
                    data_plano = (ph_elements[(6 <= ph_elements['plant']) & (ph_elements['plant'] <= 10)
                                                    & (ph_elements['organ'] == 'blade')].groupby(['t', 'plant'])['green_area'].aggregate(np.sum) * 10000).groupby(['t'])

                green_area_lamina_mean_erect = data_erect.mean()
                green_area_lamina_mean_plano = data_plano.mean()

                green_area_lamina_IC95_erect = {'bot': [], 'up': []}
                for t, green_area in data_erect:
                    IC95 = sms.DescrStatsW(green_area).tconfint_mean()
                    green_area_lamina_IC95_erect['bot'].append(IC95[0])
                    green_area_lamina_IC95_erect['up'].append(IC95[1])

                green_area_lamina_IC95_plano = {'bot': [], 'up': []}
                for t, green_area in data_plano:
                    IC95 = sms.DescrStatsW(green_area).tconfint_mean()
                    green_area_lamina_IC95_plano['bot'].append(IC95[0])
                    green_area_lamina_IC95_plano['up'].append(IC95[1])

                ax1.plot(green_area_lamina_mean_erect.index, green_area_lamina_mean_erect, colors['Caphorn_Mixed'], label=label_mapping[treatment + '_1'], linestyle=lines['Caphorn_Mixed'], linewidth=thickness)
                ax1.fill_between(green_area_lamina_mean_erect.index, green_area_lamina_IC95_erect['bot'], green_area_lamina_IC95_erect['up'],
                                 alpha=alpha, color=colors['Caphorn_Mixed'], linewidth=0.0)

                ax1.plot(green_area_lamina_mean_plano.index, green_area_lamina_mean_plano, colors['Soissons_Mixed'], label=label_mapping[treatment + '_2'], linestyle=lines['Soissons_Mixed'], linewidth=thickness)
                ax1.fill_between(green_area_lamina_mean_plano.index, green_area_lamina_IC95_plano['bot'], green_area_lamina_IC95_plano['up'],
                                 alpha=alpha, color=colors['Soissons_Mixed'], linewidth=0.0)

            else:
                data = (ph_elements[(ph_elements['organ'] == 'blade')].groupby(['t', 'plant'])['green_area'].aggregate(np.sum)*10000).groupby(['t'])
                green_area_lamina_mean = data.mean()
                green_area_lamina_IC95 = {'bot': [], 'up': []}
                for t, green_area in data:
                    IC95 = sms.DescrStatsW(green_area).tconfint_mean()
                    green_area_lamina_IC95['bot'].append(IC95[0])
                    green_area_lamina_IC95['up'].append(IC95[1])
                ax1.plot(green_area_lamina_mean.index, green_area_lamina_mean, colors[treatment], label=label_mapping[treatment], linewidth=thickness)
                ax1.fill_between(green_area_lamina_mean.index, green_area_lamina_IC95['bot'], green_area_lamina_IC95['up'], alpha=alpha, color=colors[treatment], linewidth=0.0)

            # Formatting
            ax1.set_xticks(np.arange(0, 1300, 200))
            ax1.set_yticks(np.arange(0, 150, 25))
            if AXES_MAP[density] == 0:
                ax1.tick_params(labelbottom='off', width=thickness, labelsize=fontsize)
                ax1.set_ylabel(u'Lamina green area (cm$^{2}$)', fontsize=fontsize+2)
            else:
                ax1.tick_params(labelleft='off', labelbottom='off', width=thickness, labelsize=fontsize)
            ax1.yaxis.set_label_coords(labelx, 0.5)
            [i.set_linewidth(thickness) for i in ax1.spines.itervalues()]
            ax1.axvline(360, color='k', linestyle='--')

            # 1.2) Lamina proteins
            ax2 = plt.subplot(gs[1, AXES_MAP[density]])
            ph_elements['Conc_proteins'] = (ph_elements['proteins'] * 0.145) / ph_elements['mstruct']

            if treatment == 'Asso_diffus':
                if density in ('200', '410', '600', '800'):
                    # La, n
                    data_erect = (ph_elements[(ph_elements['plant'] == 1)
                                              & (ph_elements['organ'] == 'blade') & (
                                                          ph_elements['metamer'] == 4)].groupby(['t'])['Conc_proteins'])
                    data_plano = (ph_elements[(ph_elements['plant'] == 2)
                                              & (ph_elements['organ'] == 'blade') & (
                                                          ph_elements['metamer'] == 4)].groupby(['t'])['Conc_proteins'])
                else:
                    data_erect = (ph_elements[(1 <= ph_elements['plant']) & (ph_elements['plant'] <= 5)
                                              & (ph_elements['organ'] == 'blade') & (
                                                          ph_elements['metamer'] == 4)].groupby(['t', 'plant'])[
                        'Conc_proteins'])
                    data_plano = (ph_elements[(6 <= ph_elements['plant']) & (ph_elements['plant'] <= 10)
                                              & (ph_elements['organ'] == 'blade') & (
                                                          ph_elements['metamer'] == 4)].groupby(['t', 'plant'])[
                        'Conc_proteins'])

                proteins_lamina_mean_erect = data_erect.mean()
                proteins_lamina_mean_plano = data_plano.mean()

                proteins_lamina_IC95_erect = {'bot': [], 'up': []}
                for t, proteins in data_erect:
                    IC95 = sms.DescrStatsW(proteins).tconfint_mean()
                    proteins_lamina_IC95_erect['bot'].append(IC95[0])
                    proteins_lamina_IC95_erect['up'].append(IC95[1])

                proteins_lamina_IC95_plano = {'bot': [], 'up': []}
                for t, proteins in data_plano:
                    IC95 = sms.DescrStatsW(proteins).tconfint_mean()
                    proteins_lamina_IC95_plano['bot'].append(IC95[0])
                    proteins_lamina_IC95_plano['up'].append(IC95[1])

                ax2.plot(proteins_lamina_mean_erect.index, proteins_lamina_mean_erect, colors['Caphorn_Mixed'],
                         label='_nolegend_', linestyle=lines['Caphorn_Mixed'], linewidth=thickness)
                ax2.fill_between(proteins_lamina_mean_erect.index, proteins_lamina_IC95_erect['bot'],
                                 proteins_lamina_IC95_erect['up'],
                                 label='_nolegend_', alpha=alpha, color=colors['Caphorn_Mixed'], linewidth=0.0)

                ax2.plot(proteins_lamina_mean_plano.index, proteins_lamina_mean_plano, colors['Soissons_Mixed'],
                         label='_nolegend_', linestyle=lines['Soissons_Mixed'], linewidth=thickness)
                ax2.fill_between(proteins_lamina_mean_plano.index, proteins_lamina_IC95_plano['bot'],
                                 proteins_lamina_IC95_plano['up'],
                                 label='_nolegend_', alpha=alpha, color=colors['Soissons_Mixed'], linewidth=0.0)

                # La, n-1
                if density in ('200', '410', '600', '800'):
                    data_erect = (ph_elements[(ph_elements['plant'] == 1)
                                              & (ph_elements['organ'] == 'blade') & (
                                                          ph_elements['metamer'] == 3)].groupby(['t'])['Conc_proteins'])
                    data_plano = (ph_elements[(ph_elements['plant'] == 2)
                                              & (ph_elements['organ'] == 'blade') & (
                                                          ph_elements['metamer'] == 3)].groupby(['t'])['Conc_proteins'])
                else:
                    data_erect = (ph_elements[(1 <= ph_elements['plant']) & (ph_elements['plant'] <= 5)
                                              & (ph_elements['organ'] == 'blade') & (
                                                          ph_elements['metamer'] == 3)].groupby(['t', 'plant'])[
                        'Conc_proteins'])
                    data_plano = (ph_elements[(6 <= ph_elements['plant']) & (ph_elements['plant'] <= 10)
                                              & (ph_elements['organ'] == 'blade') & (
                                                          ph_elements['metamer'] == 3)].groupby(['t', 'plant'])[
                        'Conc_proteins'])

                proteins_lamina_mean_erect = data_erect.mean()
                proteins_lamina_mean_plano = data_plano.mean()

                proteins_lamina_IC95_erect = {'bot': [], 'up': []}
                for t, proteins in data_erect:
                    IC95 = sms.DescrStatsW(proteins).tconfint_mean()
                    proteins_lamina_IC95_erect['bot'].append(IC95[0])
                    proteins_lamina_IC95_erect['up'].append(IC95[1])

                proteins_lamina_IC95_plano = {'bot': [], 'up': []}
                for t, proteins in data_plano:
                    IC95 = sms.DescrStatsW(proteins).tconfint_mean()
                    proteins_lamina_IC95_plano['bot'].append(IC95[0])
                    proteins_lamina_IC95_plano['up'].append(IC95[1])

                ax2.plot(proteins_lamina_mean_erect.index, proteins_lamina_mean_erect, colors['Caphorn_Mixed'],
                         label='_nolegend_', linestyle=lines['Caphorn_Mixed'], linewidth=thickness,
                         alpha=alpha * 3)
                ax2.plot(proteins_lamina_mean_plano.index, proteins_lamina_mean_plano, colors['Soissons_Mixed'],
                         label='_nolegend_', linestyle=lines['Soissons_Mixed'], linewidth=thickness,
                         alpha=alpha * 3)

                # # La, n-2
                # if density in ('200', '410', '600', '800'):
                #     data_erect = (ph_elements[(ph_elements['plant'] == 1)
                #                               & (ph_elements['organ'] == 'blade') & (
                #                                           ph_elements['metamer'] == 2)].groupby(['t'])['Conc_proteins'])
                #     data_plano = (ph_elements[(ph_elements['plant'] == 2)
                #                               & (ph_elements['organ'] == 'blade') & (
                #                                           ph_elements['metamer'] == 2)].groupby(['t'])['Conc_proteins'])
                # else:
                #     data_erect = (ph_elements[(1 <= ph_elements['plant']) & (ph_elements['plant'] <= 5)
                #                               & (ph_elements['organ'] == 'blade') & (
                #                                           ph_elements['metamer'] == 2)].groupby(['t', 'plant'])[
                #         'Conc_proteins'])
                #     data_plano = (ph_elements[(6 <= ph_elements['plant']) & (ph_elements['plant'] <= 10)
                #                               & (ph_elements['organ'] == 'blade') & (
                #                                           ph_elements['metamer'] == 2)].groupby(['t', 'plant'])[
                #         'Conc_proteins'])
                #
                # proteins_lamina_mean_erect = data_erect.mean()
                # proteins_lamina_mean_plano = data_plano.mean()
                #
                # proteins_lamina_IC95_erect = {'bot': [], 'up': []}
                # for t, proteins in data_erect:
                #     IC95 = sms.DescrStatsW(proteins).tconfint_mean()
                #     proteins_lamina_IC95_erect['bot'].append(IC95[0])
                #     proteins_lamina_IC95_erect['up'].append(IC95[1])
                #
                # proteins_lamina_IC95_plano = {'bot': [], 'up': []}
                # for t, proteins in data_plano:
                #     IC95 = sms.DescrStatsW(proteins).tconfint_mean()
                #     proteins_lamina_IC95_plano['bot'].append(IC95[0])
                #     proteins_lamina_IC95_plano['up'].append(IC95[1])
                #
                # ax2.plot(proteins_lamina_mean_erect.index, proteins_lamina_mean_erect, colors['Caphorn_Mixed'],
                #          label='_nolegend_', linestyle=lines['Caphorn_Mixed'], linewidth=thickness,
                #          alpha=alpha * 3)
                # ax2.plot(proteins_lamina_mean_plano.index, proteins_lamina_mean_plano, colors['Soissons_Mixed'],
                #          label='_nolegend_', linestyle=lines['Soissons_Mixed'], linewidth=thickness,
                #          alpha=alpha * 3)
                #
                # # La, n-3
                # if density in ('200', '410', '600', '800'):
                #     data_erect = (ph_elements[(ph_elements['plant'] == 1)
                #                               & (ph_elements['organ'] == 'blade') & (
                #                                           ph_elements['metamer'] == 1)].groupby(['t'])['Conc_proteins'])
                #     data_plano = (ph_elements[(ph_elements['plant'] == 2)
                #                               & (ph_elements['organ'] == 'blade') & (
                #                                           ph_elements['metamer'] == 1)].groupby(['t'])['Conc_proteins'])
                # else:
                #     data_erect = (ph_elements[(1 <= ph_elements['plant']) & (ph_elements['plant'] <= 5)
                #                               & (ph_elements['organ'] == 'blade') & (
                #                                           ph_elements['metamer'] == 1)].groupby(['t', 'plant'])[
                #         'Conc_proteins'])
                #     data_plano = (ph_elements[(6 <= ph_elements['plant']) & (ph_elements['plant'] <= 10)
                #                               & (ph_elements['organ'] == 'blade') & (
                #                                           ph_elements['metamer'] == 1)].groupby(['t', 'plant'])[
                #         'Conc_proteins'])
                #
                # proteins_lamina_mean_erect = data_erect.mean()
                # proteins_lamina_mean_plano = data_plano.mean()
                #
                # proteins_lamina_IC95_erect = {'bot': [], 'up': []}
                # for t, proteins in data_erect:
                #     IC95 = sms.DescrStatsW(proteins).tconfint_mean()
                #     proteins_lamina_IC95_erect['bot'].append(IC95[0])
                #     proteins_lamina_IC95_erect['up'].append(IC95[1])
                #
                # proteins_lamina_IC95_plano = {'bot': [], 'up': []}
                # for t, proteins in data_plano:
                #     IC95 = sms.DescrStatsW(proteins).tconfint_mean()
                #     proteins_lamina_IC95_plano['bot'].append(IC95[0])
                #     proteins_lamina_IC95_plano['up'].append(IC95[1])
                #
                # ax2.plot(proteins_lamina_mean_erect.index, proteins_lamina_mean_erect, colors['Caphorn_Mixed'],
                #          label='_nolegend_', linestyle=lines['Caphorn_Mixed'], linewidth=thickness,
                #          alpha=alpha * 3)
                # ax2.plot(proteins_lamina_mean_plano.index, proteins_lamina_mean_plano, colors['Soissons_Mixed'],
                #          label='_nolegend_', linestyle=lines['Soissons_Mixed'], linewidth=thickness,
                #          alpha=alpha * 3)
            else:
                # La, n
                data = (ph_elements[(ph_elements['organ'] == 'blade') & (ph_elements['metamer'] == 4)].groupby(['t'])[
                    'Conc_proteins'])
                proteins_lamina_mean = data.mean()
                proteins_lamina_IC95 = {'bot': [], 'up': []}
                for t, proteins in data:
                    IC95 = sms.DescrStatsW(proteins).tconfint_mean()
                    proteins_lamina_IC95['bot'].append(IC95[0])
                    proteins_lamina_IC95['up'].append(IC95[1])
                ax2.plot(proteins_lamina_mean.index, proteins_lamina_mean, colors[treatment],
                         label='_nolegend_', linewidth=thickness)
                ax2.fill_between(proteins_lamina_mean.index, proteins_lamina_IC95['bot'], proteins_lamina_IC95['up'],
                                 label='_nolegend_', alpha=alpha, color=colors[treatment], linewidth=0.0)

                # La, n-1
                data = (ph_elements[(ph_elements['organ'] == 'blade') & (ph_elements['metamer'] == 3)].groupby(['t'])[
                    'Conc_proteins'])
                proteins_lamina_mean = data.mean()
                proteins_lamina_IC95 = {'bot': [], 'up': []}
                for t, proteins in data:
                    IC95 = sms.DescrStatsW(proteins).tconfint_mean()
                    proteins_lamina_IC95['bot'].append(IC95[0])
                    proteins_lamina_IC95['up'].append(IC95[1])
                ax2.plot(proteins_lamina_mean.index, proteins_lamina_mean, colors[treatment],
                         label='_nolegend_', linewidth=thickness, alpha=alpha * 3)

                # # La, n-2
                # data = (ph_elements[(ph_elements['organ'] == 'blade') & (ph_elements['metamer'] == 2)].groupby(['t'])[
                #     'Conc_proteins'])
                # proteins_lamina_mean = data.mean()
                # proteins_lamina_IC95 = {'bot': [], 'up': []}
                # for t, proteins in data:
                #     IC95 = sms.DescrStatsW(proteins).tconfint_mean()
                #     proteins_lamina_IC95['bot'].append(IC95[0])
                #     proteins_lamina_IC95['up'].append(IC95[1])
                # ax2.plot(proteins_lamina_mean.index, proteins_lamina_mean, colors[treatment],
                #          label='_nolegend_', linewidth=thickness)
                # ax2.fill_between(proteins_lamina_mean.index, proteins_lamina_IC95['bot'], proteins_lamina_IC95['up'],
                #                  label='_nolegend_', alpha=alpha, color=colors[treatment], linewidth=0.0)
                #
                # # La, n-3
                # data = (ph_elements[(ph_elements['organ'] == 'blade') & (ph_elements['metamer'] == 1)].groupby(['t'])[
                #     'Conc_proteins'])
                # proteins_lamina_mean = data.mean()
                # proteins_lamina_IC95 = {'bot': [], 'up': []}
                # for t, proteins in data:
                #     IC95 = sms.DescrStatsW(proteins).tconfint_mean()
                #     proteins_lamina_IC95['bot'].append(IC95[0])
                #     proteins_lamina_IC95['up'].append(IC95[1])
                # ax2.plot(proteins_lamina_mean.index, proteins_lamina_mean, colors[treatment],
                #          label='_nolegend_', linewidth=thickness, alpha=alpha * 3)

            # Formatting
            ax2.set_xticks(np.arange(0, 1300, 200))
            ax2.set_xlabel('Time from flowering (hour)', fontsize=fontsize)
            ax2.set_yticks(np.arange(0, 500, 100))
            ax2.set_yticks(ax2.get_yticks()[:-1])
            if AXES_MAP[density] == 0:
                ax2.tick_params(labelbottom='on', width=thickness, labelsize=fontsize)
                ax2.set_ylabel(u'[proteins] (µmol g$^{-1}$)', fontsize=fontsize + 2)
                ax2.set_xticks(ax2.get_xticks()[:-1])
            else:
                ax2.tick_params(labelleft='off', labelbottom='on', width=thickness, labelsize=fontsize)
                if AXES_MAP[density] != 3:
                    ax2.set_xticks(ax2.get_xticks()[:-1])
            ax2.yaxis.set_label_coords(labelx, 0.5)
            [i.set_linewidth(thickness) for i in ax2.spines.itervalues()]
            ax2.axvline(360, color='k', linestyle='--')

    ax1.legend(prop={'size': 11}, bbox_to_anchor=(0.5, 1.4), ncol=4, frameon=True)
    plt.subplot(gs[0, 0]).text(300, 140, 'Density 200', fontsize=fontsize + 2, backgroundcolor='w',
                               verticalalignment='top')
    plt.subplot(gs[0, 0]).set_zorder(1)
    plt.subplot(gs[1, 0]).text(100, 380, 'Lamina n', fontsize=fontsize, backgroundcolor='w',
                               verticalalignment='top')
    plt.subplot(gs[1, 0]).text(100, 230, 'Lamina n-1', fontsize=fontsize, backgroundcolor='w',
                               verticalalignment='top')
    plt.subplot(gs[0, 1]).text(300, 140, 'Density 410', fontsize=fontsize + 2, backgroundcolor='w',
                               verticalalignment='top')
    plt.subplot(gs[0, 1]).set_zorder(1)
    plt.subplot(gs[0, 2]).text(300, 140, 'Density 600', fontsize=fontsize + 2, backgroundcolor='w',
                               verticalalignment='top')
    plt.subplot(gs[0, 2]).set_zorder(1)
    plt.subplot(gs[0, 3]).text(300, 140, 'Density 800', fontsize=fontsize + 2, backgroundcolor='w',
                               verticalalignment='top')
    plt.subplot(gs[0, 3]).text(1100, 120, 'A', fontsize=fontsize + 2, backgroundcolor='w',
                               verticalalignment='top')
    plt.subplot(gs[0, 3]).set_zorder(1)
    plt.subplot(gs[1, 3]).text(1100, 380, 'B', fontsize=fontsize + 2, backgroundcolor='w',
                               verticalalignment='top')
    plt.subplot(gs[1, 3]).set_zorder(1)

    plt.savefig(os.path.join(GRAPHS_DIRPATH, 'Figure_{}_Culm_scale.TIFF'.format(num)), dpi=300, format='TIFF',
                bbox_inches='tight')
    plt.close()

def figure_Cass_grain_mass(num):

    Cass_grain_mass = {'Density': [], 'Leaf posture': [], 'Stand': [], 'Sky': [], 'plant': [], 'Net_Photosynthesis': [], 'grain_mass': []}

    for density in DENSITIES:
        for treatment in TREATMENTS:
            data_elements_dirpath = os.path.join(DATA_DIRPATH, 'Densite_{}'.format(density), treatment, 'outputs', 'elements_states.csv')
            data_elements_df = pd.read_csv(data_elements_dirpath)
            data_elements_df.dropna(inplace=True)

            data_organs_dirpath = os.path.join(DATA_DIRPATH, 'Densite_{}'.format(density), treatment, 'outputs', 'organs_states.csv')
            data_organs_df = pd.read_csv(data_organs_dirpath)
            data_organs_df['Dry_Mass'] = ((data_organs_df['structure'] + data_organs_df['starch']) * 1E-6 * C_MOLAR_MASS / 0.4) + (
                        (data_organs_df['proteins'] * 1E-6 * N_MOLAR_MASS) / 0.136)

            if treatment not in ('Asso_diffus', 'Asso_mixte', 'Asso_direct'):
                leaf_posture, sky = treatment.split('_')
                data_elements = data_elements_df.groupby('plant')['Photosynthesis', 'R_Nnit_red', 'R_phloem_loading', 'R_residual']
                for plant, data in data_elements:
                    Cass_grain_mass['Density'].append(density)
                    Cass_grain_mass['Leaf posture'].append(leaf_posture)
                    Cass_grain_mass['Stand'].append('pure')
                    Cass_grain_mass['Sky'].append(sky)
                    Cass_grain_mass['plant'].append(plant)
                    Cass_grain_mass['Net_Photosynthesis'].append((data['Photosynthesis'].sum() - data['R_Nnit_red'].sum() - data['R_phloem_loading'].sum() - data['R_residual'].sum()) * C_MOLAR_MASS * 1E-6)
                    Cass_grain_mass['grain_mass'].append(data_organs_df[(data_organs_df['plant'] == plant) & (data_organs_df['organ'] == 'grains')]['Dry_Mass'].iloc[-1])

            else:
                stand, sky = treatment.split('_')
                if density in ('200', '410', '600', '800'):
                    data_erect = data_elements_df[data_elements_df['plant'] == 1].groupby('plant')['Photosynthesis', 'R_Nnit_red', 'R_phloem_loading', 'R_residual']
                    data_plano = data_elements_df[data_elements_df['plant'] == 2].groupby('plant')['Photosynthesis', 'R_Nnit_red', 'R_phloem_loading', 'R_residual']
                else:
                    data_erect = data_elements_df[(1 <= data_elements_df['plant']) & (data_elements_df['plant'] <= 5)].groupby('plant')['Photosynthesis', 'R_Nnit_red', 'R_phloem_loading', 'R_residual']
                    data_plano = data_elements_df[(6 <= data_elements_df['plant']) & (data_elements_df['plant'] <= 10)].groupby('plant')['Photosynthesis', 'R_Nnit_red', 'R_phloem_loading', 'R_residual']

                for plant, data in data_erect:
                    Cass_grain_mass['Density'].append(density)
                    Cass_grain_mass['Stand'].append(stand)
                    Cass_grain_mass['Leaf posture'].append('Erect')
                    Cass_grain_mass['Sky'].append(sky)
                    Cass_grain_mass['plant'].append(plant)
                    Cass_grain_mass['Net_Photosynthesis'].append((data['Photosynthesis'].sum() - data['R_Nnit_red'].sum() - data['R_phloem_loading'].sum() - data['R_residual'].sum()) * C_MOLAR_MASS * 1E-6)
                    Cass_grain_mass['grain_mass'].append(data_organs_df[(data_organs_df['plant'] == plant) & (data_organs_df['organ'] == 'grains')]['Dry_Mass'].iloc[-1])

                for plant, data in data_plano:
                    Cass_grain_mass['Density'].append(density)
                    Cass_grain_mass['Stand'].append(stand)
                    Cass_grain_mass['Leaf posture'].append('Erect')
                    Cass_grain_mass['Sky'].append(sky)
                    Cass_grain_mass['plant'].append(plant)
                    Cass_grain_mass['Net_Photosynthesis'].append((data['Photosynthesis'].sum() - data['R_Nnit_red'].sum() - data['R_phloem_loading'].sum() - data['R_residual'].sum()) * C_MOLAR_MASS * 1E-6)
                    Cass_grain_mass['grain_mass'].append(data_organs_df[(data_organs_df['plant'] == plant) & (data_organs_df['organ'] == 'grains')]['Dry_Mass'].iloc[-1])

            del plant, data
            gc.collect()


    Cass_grain_mass_df = pd.DataFrame(Cass_grain_mass)
    data_df_grouped = Cass_grain_mass_df.groupby(['Density', 'Leaf posture'])

    fig, ax = plt.subplots()

    for group_name, group in data_df_grouped:
        density = group_name[0]
        leaf_inclination = group_name[1]
        if density in ('200', '410', '600', '800'):
            ax.plot(group['Net_Photosynthesis'], group['grain_mass'], marker=markers[leaf_inclination], markersize=10, label=density, linewidth=0, color=colors[density], zorder=10)
        else:
            ax.plot(group['Net_Photosynthesis'], group['grain_mass'], marker=markers[leaf_inclination], markersize=5, markerfacecolor='None',
                        markeredgecolor=colors[density], label=density, linewidth=0, zorder=1)

    z = np.polyfit(Cass_grain_mass_df['Net_Photosynthesis'], Cass_grain_mass_df['grain_mass'], 2)
    print r'y = {}x² + {}x + {}'.format(round(z[0], 3), round(z[1], 3), round(z[2], 3))
    p = np.poly1d(z)
    xp = np.linspace(min(Cass_grain_mass_df['Net_Photosynthesis']), max(Cass_grain_mass_df['Net_Photosynthesis']), 100)
    ax.plot(xp, p(xp), 'k-', zorder=2)
    ax.text(0.4, 0.9, 'y = {}x$^2$ + {}x + {}'.format(round(z[0], 3), round(z[1], 3), round(z[2], 3)), transform=ax.transAxes, fontsize=fontsize, verticalalignment='top')

    [i.set_linewidth(thickness) for i in ax.spines.itervalues()]
    ax.set_xticks(np.arange(0, 5, 1))
    ax.set_yticks(np.arange(0, 5, 1))
    ax.set_xlabel('Net C assimilated (g)', fontsize=fontsize+5)
    ax.set_ylabel('Grain final dry mass (g)', fontsize=fontsize+5)

    plt.savefig(os.path.join(GRAPHS_DIRPATH, 'Figure_{}_Cass_vs_grain_mass.TIFF'.format(num)), dpi=300, format='TIFF', bbox_inches='tight')
    plt.close()


def figure_C_non_struct(num):

    plt.figure(figsize=cm2inch(30, 20))
    gs = gridspec.GridSpec(4, 4, hspace=0, wspace=0)

    for density in DENSITIES:
        for treatment in TREATMENTS:
            if treatment not in ('Plano_diffus', 'Erect_diffus', 'Asso_diffus'): continue
            stand = treatment.split('_')[0]

            axes = []

            # Photosynthetic elements
            ph_elements_path = os.path.join(DATA_DIRPATH, 'Densite_{}'.format(density), treatment, 'outputs\elements_states.csv')
            ph_elements = pd.read_csv(ph_elements_path)
            ph_elements['sum_C'] = (ph_elements['triosesP'] + ph_elements['starch'] + ph_elements['fructan'] + ph_elements['sucrose']) * (1E-3 * C_MOLAR_MASS)
            total_area_shoot = ph_elements.groupby(['t', 'plant'])['green_area'].sum()
            area_stem = ph_elements[ph_elements.organ.isin(['internode', 'sheath', 'peduncle'])].groupby(['t', 'plant'])['green_area'].sum()
            contrib_area_stem = (area_stem / total_area_shoot) * phloem_shoot_root

            # Organs
            organs_path = os.path.join(DATA_DIRPATH, 'Densite_{}'.format(density), treatment, 'outputs\organs_states.csv')
            organs = pd.read_csv(organs_path)
            organs['sum_C'] = (organs['sucrose']) * (1E-3 * C_MOLAR_MASS)
            organs['Dry_Mass'] = ((organs['structure'] + organs['starch']) * 1E-6 * 12 / 0.4) + (
                        (organs['proteins'] * 1E-6 * 14) / 0.136) # for grains

            ax1 = plt.subplot(gs[0, AXES_MAP[density]])
            ax2 = plt.subplot(gs[1, AXES_MAP[density]])
            ax3 = plt.subplot(gs[2, AXES_MAP[density]])
            ax4 = plt.subplot(gs[3, AXES_MAP[density]])
            axes.append(ax1)
            axes.append(ax2)
            axes.append(ax3)
            axes.append(ax4)

            stand_data = {}
            if stand == 'Asso':
                if density in ('200', '410', '600', '800'):
                    # Erectophile
                    grains = (organs[(organs['plant'] == 1) &
                                    (organs['organ'] == 'grains')].groupby(['t', 'plant'])['Dry_Mass'].aggregate(
                        np.sum)).groupby(['t'])
                    stem = (ph_elements[(ph_elements['plant'] == 1) &
                                              (ph_elements.organ.isin(['internode', 'sheath', 'peduncle']))].groupby(['t', 'plant'])['sum_C'].aggregate(np.sum)).groupby(['t'])
                    phloem = (organs[(organs['plant'] == 1) &
                                           (organs['organ'] == 'phloem')].groupby(['t', 'plant'])['sum_C'].aggregate(np.sum)).groupby(['t'])
                    roots = (organs[(organs['plant'] ==1) &
                                           (organs['organ'] == 'roots')].groupby(['t', 'plant'])['sum_C'].aggregate(np.sum)).groupby(['t'])
                    stand_data['Caphorn_Mixed'] = [grains, stem, phloem, roots]

                    # Planophile
                    grains = (organs[(organs['plant'] == 2) &
                                    (organs['organ'] == 'grains')].groupby(['t', 'plant'])['Dry_Mass'].aggregate(
                        np.sum)).groupby(['t'])
                    stem = (ph_elements[(ph_elements['plant'] == 2) &
                                              (ph_elements.organ.isin(['internode', 'sheath', 'peduncle']))].groupby(['t', 'plant'])['sum_C'].aggregate(np.sum)).groupby(['t'])
                    phloem = (organs[(organs['plant'] == 2) &
                                           (organs['organ'] == 'phloem')].groupby(['t', 'plant'])['sum_C'].aggregate(np.sum)).groupby(['t'])
                    roots = (organs[(organs['plant'] == 2) &
                                           (organs['organ'] == 'roots')].groupby(['t', 'plant'])['sum_C'].aggregate(np.sum)).groupby(['t'])
                    stand_data['Soissons_Mixed'] = [grains, stem, phloem, roots]
                else:
                    # Erectophile
                    grains = (organs[(1 <= organs['plant']) & (organs['plant'] <= 5) &
                                    (organs['organ'] == 'grains')].groupby(['t', 'plant'])['Dry_Mass'].aggregate(
                        np.sum)).groupby(['t'])
                    stem = (ph_elements[(1 <= ph_elements['plant']) & (ph_elements['plant'] <= 5) &
                                              (ph_elements.organ.isin(['internode', 'sheath', 'peduncle']))].groupby(['t', 'plant'])['sum_C'].aggregate(np.sum)).groupby(['t'])
                    phloem = (organs[(1 <= organs['plant']) & (organs['plant'] <= 5) &
                                           (organs['organ'] == 'phloem')].groupby(['t', 'plant'])['sum_C'].aggregate(np.sum)).groupby(['t'])
                    roots = (organs[(1 <= organs['plant']) & (organs['plant'] <= 5) &
                                           (organs['organ'] == 'roots')].groupby(['t', 'plant'])['sum_C'].aggregate(np.sum)).groupby(['t'])
                    stand_data['Caphorn_Mixed'] = [grains, stem, phloem, roots]

                    # Planophile
                    grains = (organs[(6 <= organs['plant']) & (organs['plant'] <= 10) &
                                    (organs['organ'] == 'grains')].groupby(['t', 'plant'])['Dry_Mass'].aggregate(
                        np.sum)).groupby(['t'])
                    stem = (ph_elements[(6 <= ph_elements['plant']) & (ph_elements['plant'] <= 10) &
                                              (ph_elements.organ.isin(['internode', 'sheath', 'peduncle']))].groupby(['t', 'plant'])['sum_C'].aggregate(np.sum)).groupby(['t'])
                    phloem = (organs[(6 <= organs['plant']) & (organs['plant'] <= 10) &
                                           (organs['organ'] == 'phloem')].groupby(['t', 'plant'])['sum_C'].aggregate(np.sum)).groupby(['t'])
                    roots = (organs[(6 <= organs['plant']) & (organs['plant'] <= 10) &
                                           (organs['organ'] == 'roots')].groupby(['t', 'plant'])['sum_C'].aggregate(np.sum)).groupby(['t'])
                    stand_data['Soissons_Mixed'] = [grains, stem, phloem, roots]
            else:
                grains = (organs[(organs['organ'] == 'grains')].groupby(['t', 'plant'])['Dry_Mass'].aggregate(np.sum)).groupby(['t'])
                stem = (ph_elements[ph_elements.organ.isin(['internode', 'sheath', 'peduncle'])].groupby(['t', 'plant'])['sum_C'].aggregate(np.sum)).groupby(['t'])
                phloem = (organs[(organs['organ'] == 'phloem')].groupby(['t', 'plant'])['sum_C'].aggregate(np.sum)).groupby(['t'])
                roots = (organs[(organs['organ'] == 'roots')].groupby(['t', 'plant'])['sum_C'].aggregate(np.sum)).groupby(['t'])
                stand_data[treatment] = [grains, stem, phloem, roots]

            for stand, data in stand_data.iteritems():
                grains, stem, phloem, roots = data[0], data[1], data[2],  data[3]
                # Grains
                dry_mass_mean = grains.mean()
                dry_mass_IC95 = {'bot': [], 'up': []}
                for t, dry_mass in grains:
                   IC95 = sms.DescrStatsW(dry_mass).tconfint_mean()
                   dry_mass_IC95['bot'].append(IC95[0])
                   dry_mass_IC95['up'].append(IC95[1])

                ax1.plot(dry_mass_mean.index, dry_mass_mean, colors[stand], label=label_mapping[stand], linestyle=lines[stand], linewidth=thickness)
                ax1.fill_between(dry_mass_mean.index, dry_mass_IC95['bot'],
                                 dry_mass_IC95['up'], alpha=alpha, color=colors[stand], linewidth=0.0)

                ax1.set_yticks(np.arange(0, 6, 2))

                # Stem
                sum_C_corrected_dict = {'t': [], 'sum_C_corrected': []}
                sum_C_IC95 = {'bot': [], 'up': []}
                for t, sum_C in stem:
                    sum_C_corrected = sum_C + (phloem.get_group(t) * contrib_area_stem[t])
                    sum_C_corrected_dict['t'].append(t)
                    sum_C_corrected_dict['sum_C_corrected'].append(sum_C_corrected.mean())
                    IC95 = sms.DescrStatsW(sum_C_corrected).tconfint_mean()
                    sum_C_IC95['bot'].append(IC95[0])
                    sum_C_IC95['up'].append(IC95[1])

                ax2.plot(sum_C_corrected_dict['t'], sum_C_corrected_dict['sum_C_corrected'], colors[stand], label=label_mapping[stand], linestyle=lines[stand], linewidth=thickness)
                ax2.fill_between(sum_C_corrected_dict['t'], sum_C_IC95['bot'],
                                sum_C_IC95['up'], alpha=alpha, color=colors[stand], linewidth=0.0)

                ax2.set_yticks(np.arange(0, 800, 200))

                # Phloem
                sum_C_mean = phloem.mean()
                sum_C_IC95 = {'bot': [], 'up': []}
                for t, sum_C in phloem:
                   IC95 = sms.DescrStatsW(sum_C).tconfint_mean()
                   sum_C_IC95['bot'].append(IC95[0])
                   sum_C_IC95['up'].append(IC95[1])

                ax3.plot(sum_C_mean.index, sum_C_mean, colors[stand], label=label_mapping[stand], linestyle=lines[stand], linewidth=thickness)
                ax3.fill_between(sum_C_mean.index, sum_C_IC95['bot'],
                                sum_C_IC95['up'], alpha=alpha, color=colors[stand], linewidth=0.0)

                ax3.set_yticks(np.arange(0, 800, 200))

                # Roots
                sum_C_corrected_dict = {'t': [], 'sum_C_corrected': []}
                sum_C_IC95 = {'bot': [], 'up': []}
                for t, sum_C in roots:
                    sum_C_corrected = sum_C + (phloem.get_group(t) * (1 - phloem_shoot_root))
                    sum_C_corrected_dict['t'].append(t)
                    sum_C_corrected_dict['sum_C_corrected'].append(sum_C_corrected.mean())
                    IC95 = sms.DescrStatsW(sum_C_corrected).tconfint_mean()
                    sum_C_IC95['bot'].append(IC95[0])
                    sum_C_IC95['up'].append(IC95[1])

                ax4.plot(sum_C_corrected_dict['t'], sum_C_corrected_dict['sum_C_corrected'], colors[stand], label=label_mapping[stand], linestyle=lines[stand], linewidth=thickness)
                ax4.fill_between(sum_C_corrected_dict['t'], sum_C_IC95['bot'],
                                sum_C_IC95['up'], alpha=alpha, color=colors[stand], linewidth=0.0)

                ax4.set_yticks(np.arange(0, 800, 200))

            for ax in axes:
               ax.yaxis.set_label_coords(labelx, 0.5)
               [i.set_linewidth(thickness) for i in ax.spines.itervalues()]
               ax.tick_params(width=thickness, length=10, labelbottom='off', labelleft='off', top='off', right='off', labelsize=fontsize)
               ax.set_xticks(np.arange(0, 1300, 200))
               ax.axvline(360, color='k', linestyle='--')

    ax_label_left = (plt.subplot(gs[0, 0]), plt.subplot(gs[1, 0]), plt.subplot(gs[2, 0]), plt.subplot(gs[3, 0]))
    for ax in ax_label_left:
       ax.tick_params(labelleft='on')
       if ax == plt.subplot(gs[0, 0]):
           ax.set_ylabel('Grain dry mass (g)', fontsize=fontsize+2)
       if ax == plt.subplot(gs[2, 0]):
           ax.set_ylabel('Non structural C mass (mg C)', fontsize=fontsize+2)
       if ax != plt.subplot(gs[0, 0]):
           ax.set_yticks(ax.get_yticks()[:-1])

    ax_label_bottom = (plt.subplot(gs[3, 0]), plt.subplot(gs[3, 1]), plt.subplot(gs[3, 2]), plt.subplot(gs[3, 3]))
    for ax in ax_label_bottom:
       ax.tick_params(labelbottom='on')
       ax.set_xlabel('Time from flowering (hour)', fontsize=fontsize)
       if ax != plt.subplot(gs[2, 3]):
           ax.set_xticks(ax.get_xticks()[:-1])

    xtext = 300
    ytext = 4.75
    plt.subplot(gs[0, 0]).text(xtext, ytext, 'Density 200', fontsize=fontsize+3, backgroundcolor='w', verticalalignment='top')
    plt.subplot(gs[0, 0]).set_zorder(1)
    plt.subplot(gs[0, 1]).text(xtext, ytext, 'Density 410', fontsize=fontsize+3, backgroundcolor='w', verticalalignment='top')
    plt.subplot(gs[0, 1]).set_zorder(1)
    plt.subplot(gs[0, 2]).text(xtext, ytext, 'Density 600', fontsize=fontsize+3, backgroundcolor='w', verticalalignment='top')
    plt.subplot(gs[0, 2]).set_zorder(1)
    plt.subplot(gs[0, 3]).text(xtext, ytext, 'Density 800', fontsize=fontsize+3, backgroundcolor='w', verticalalignment='top')
    plt.subplot(gs[0, 3]).set_zorder(1)

    xtext = 50
    ytext = 575
    plt.subplot(gs[0, 0]).text(xtext, 3.8, 'Grains', fontsize=fontsize+2, verticalalignment='top')
    plt.subplot(gs[0, 0]).set_zorder(1)
    plt.subplot(gs[1, 0]).text(xtext, ytext, 'Stem', fontsize=fontsize+2, backgroundcolor='w', verticalalignment='top')
    plt.subplot(gs[1, 0]).set_zorder(1)
    plt.subplot(gs[2, 0]).text(xtext, ytext, 'Phloem', fontsize=fontsize+2, backgroundcolor='w', verticalalignment='top')
    plt.subplot(gs[2, 0]).set_zorder(1)
    plt.subplot(gs[3, 0]).text(xtext, ytext, 'Roots', fontsize=fontsize+2, backgroundcolor='w', verticalalignment='top')
    plt.subplot(gs[3, 0]).set_zorder(1)

    ax1.legend(prop={'size':11}, bbox_to_anchor=(0.5, 1.5), ncol=4, frameon=True)
    plt.savefig(os.path.join(GRAPHS_DIRPATH, 'Figure_{}_C_dynamics.TIFF'.format(num)), dpi=300, format='TIFF', bbox_inches='tight')
    plt.close()


def table_PARa(num):

    Total_PARa = {'Density': [], 'Leaf posture': [], 'Stand': [], 'Sky': [], 'Total_PARa': [], 'IC95': [], 'Variance': []}

    for density in DENSITIES:
        for treatment in TREATMENTS:
            data_dirpath = os.path.join(DATA_DIRPATH, 'Densite_{}'.format(density), treatment, 'outputs\elements_states.csv')
            data_df = pd.read_csv(data_dirpath)
            data_df.dropna(inplace=True)
            data_df['Total_PARa'] = data_df['PARa'] * data_df['green_area']

            if treatment not in ('Asso_diffus', 'Asso_mixte', 'Asso_direct'):
                Total_PARa['Density'].append(density)
                leaf_posture, sky = treatment.split('_')
                Total_PARa['Leaf posture'].append(leaf_posture)
                Total_PARa['Stand'].append('pure')
                Total_PARa['Sky'].append(sky)
                data = data_df.groupby('plant')['Total_PARa'].sum() * 3600 * 1E-6
                Total_PARa['Total_PARa'].append(data.mean())
                Total_PARa['IC95'].append(sms.DescrStatsW(data).tconfint_mean())
                Total_PARa['Variance'].append(data.var())

            else:
                stand, sky = treatment.split('_')
                Total_PARa['Density'].append(density)
                Total_PARa['Density'].append(density)
                Total_PARa['Stand'].append(stand)
                Total_PARa['Stand'].append(stand)
                Total_PARa['Leaf posture'].append('Erect')
                Total_PARa['Leaf posture'].append('Plano')
                Total_PARa['Sky'].append(sky)
                Total_PARa['Sky'].append(sky)
                if density in ('200', '410', '600', '800'):
                    data_erect = data_df[(data_df['plant'] == 1)].groupby('plant')[
                                     'Total_PARa'].sum() * 3600 * 1E-6
                    data_plano = data_df[(data_df['plant'] == 2)].groupby('plant')[
                                     'Total_PARa'].sum() * 3600 * 1E-6
                else:
                    data_erect = data_df[(1 <= data_df['plant']) & (data_df['plant'] <= 5)].groupby('plant')['Total_PARa'].sum() * 3600 * 1E-6
                    data_plano = data_df[(6 <= data_df['plant']) & (data_df['plant'] <= 10)].groupby('plant')['Total_PARa'].sum() * 3600 * 1E-6
                Total_PARa['Total_PARa'].append(data_erect.mean())
                Total_PARa['Total_PARa'].append(data_plano.mean())
                Total_PARa['IC95'].append(sms.DescrStatsW(data_erect).tconfint_mean())
                Total_PARa['IC95'].append(sms.DescrStatsW(data_plano).tconfint_mean())
                Total_PARa['Variance'].append(data_erect.var())
                Total_PARa['Variance'].append(data_plano.var())

    Total_PARa_df = pd.DataFrame(Total_PARa)
    Total_PARa_df.to_csv(os.path.join(GRAPHS_DIRPATH, 'Table_{}_PARa.csv'.format(num)), index=False)


def table_SumC(num):

    sumC = {'Density': [], 'Leaf posture': [], 'Stand': [], 'Sky': [], 'Photosynthesis': [], 'IC95': [], 'Variance': []}

    for density in DENSITIES:
        for treatment in TREATMENTS:
            data_dirpath = os.path.join(DATA_DIRPATH, 'Densite_{}'.format(density), treatment, 'outputs\elements_states.csv')
            data_df = pd.read_csv(data_dirpath)

            if treatment not in ('Asso_diffus', 'Asso_mixte', 'Asso_direct'):
                sumC['Density'].append(density)
                leaf_posture, sky = treatment.split('_')
                sumC['Leaf posture'].append(leaf_posture)
                sumC['Stand'].append('pure')
                sumC['Sky'].append(sky)
                data = data_df.groupby('plant')['Photosynthesis'].sum() * 1E-3
                sumC['Photosynthesis'].append(data.mean())
                sumC['IC95'].append(sms.DescrStatsW(data).tconfint_mean())
                sumC['Variance'].append(data.var())

            else:
                stand, sky = treatment.split('_')
                sumC['Density'].append(density)
                sumC['Density'].append(density)
                sumC['Stand'].append(stand)
                sumC['Stand'].append(stand)
                sumC['Leaf posture'].append('Erect')
                sumC['Leaf posture'].append('Plano')
                sumC['Sky'].append(sky)
                sumC['Sky'].append(sky)
                if density in ('200', '410', '600', '800'):
                    data_erect = data_df[data_df['plant'] == 1].groupby('plant')[
                                     'Photosynthesis'].sum() * 1E-3
                    data_plano = data_df[data_df['plant'] == 2].groupby('plant')[
                                     'Photosynthesis'].sum() * 1E-3
                else:
                    data_erect = data_df[(1 <= data_df['plant']) & (data_df['plant'] <= 5)].groupby('plant')[
                                    'Photosynthesis'].sum() * 1E-3
                    data_plano = data_df[(6 <= data_df['plant']) & (data_df['plant'] <= 10)].groupby('plant')[
                                    'Photosynthesis'].sum() * 1E-3
                sumC['Photosynthesis'].append(data_erect.mean())
                sumC['Photosynthesis'].append(data_plano.mean())
                sumC['IC95'].append(sms.DescrStatsW(data_erect).tconfint_mean())
                sumC['IC95'].append(sms.DescrStatsW(data_plano).tconfint_mean())
                sumC['Variance'].append(data_erect.var())
                sumC['Variance'].append(data_plano.var())

    sumC_df = pd.DataFrame(sumC)
    sumC_df.to_csv(os.path.join(GRAPHS_DIRPATH, 'Table_{}_Sum_C.csv'.format(num)), index=False)


def table_grains_final_state():
    data_dict_grains_final = {'Density': [], 'Leaf posture': [], 'Stand': [], 'Sky': [], 'Dry_Mass': [], 'Proteins_N_Mass': []}

    for density in DENSITIES:
        for treatment in TREATMENTS:
            if (density, treatment) == (410, 'Plano_mixte'): continue
            data_dirpath = os.path.join(DATA_DIRPATH, 'Densite_{}'.format(density), treatment, 'outputs\organs_states.csv')
            data_df = pd.read_csv(data_dirpath)
            data_df['Dry_Mass'] = ((data_df['structure'] + data_df['starch']) * 1E-6 * C_MOLAR_MASS / 0.4) + ((data_df['proteins'] * 1E-6 * N_MOLAR_MASS) / 0.136)
            data_df['Proteins_N_Mass'] = (data_df['proteins'] * 1E-6 * N_MOLAR_MASS) / 0.136
            data_df_grains_final = data_df[(data_df['organ'] == 'grains') & (data_df['t'] == data_df['t'].unique()[-1])]


            if treatment not in ('Asso_direct', 'Asso_mixte', 'Asso_diffus'):
                data_dict_grains_final['Density'].append(density)
                leaf_posture, sky = treatment.split('_')
                data_dict_grains_final['Leaf posture'].append(leaf_posture)
                data_dict_grains_final['Stand'].append('pure')
                data_dict_grains_final['Sky'].append(sky)
                data_dict_grains_final['Dry_Mass'].append(data_df_grains_final.groupby('plant')['Dry_Mass'].mean().mean())
                data_dict_grains_final['Proteins_N_Mass'].append(data_df_grains_final.groupby('plant')['Proteins_N_Mass'].mean().mean())

            else:
                stand, sky = treatment.split('_')
                data_dict_grains_final['Density'].append(density)
                data_dict_grains_final['Density'].append(density)
                data_dict_grains_final['Stand'].append(stand)
                data_dict_grains_final['Stand'].append(stand)
                data_dict_grains_final['Leaf posture'].append('Erect')
                data_dict_grains_final['Leaf posture'].append('Plano')
                data_dict_grains_final['Sky'].append(sky)
                data_dict_grains_final['Sky'].append(sky)
                if density in ('200', '410', '600', '800'):
                    data_dict_grains_final['Dry_Mass'].append(data_df_grains_final[(data_df_grains_final['plant'] == 1)].groupby('plant')['Dry_Mass'].mean().mean())
                    data_dict_grains_final['Dry_Mass'].append(data_df_grains_final[(data_df_grains_final['plant'] == 2)].groupby('plant')['Dry_Mass'].mean().mean())
                else:
                    data_dict_grains_final['Dry_Mass'].append(data_df_grains_final[(1 <= data_df_grains_final['plant']) & (data_df_grains_final['plant'] <= 5)].groupby('plant')['Dry_Mass'].mean().mean())
                    data_dict_grains_final['Dry_Mass'].append(data_df_grains_final[(6 <= data_df_grains_final['plant']) & (data_df_grains_final['plant'] <= 10)].groupby('plant')['Dry_Mass'].mean().mean())

                data_dict_grains_final['Proteins_N_Mass'].append(data_df_grains_final[(1 <= data_df_grains_final['plant']) & (data_df_grains_final['plant'] <= 5)].groupby('plant')['Proteins_N_Mass'].mean().mean())
                data_dict_grains_final['Proteins_N_Mass'].append(data_df_grains_final[(6 <= data_df_grains_final['plant']) & (data_df_grains_final['plant'] <= 10)].groupby('plant')['Proteins_N_Mass'].mean().mean())

    grains_data_df = pd.DataFrame(data_dict_grains_final)
    grains_data_df.to_csv(os.path.join(GRAPHS_DIRPATH, 'grains_final_state.csv'), index=False)


def table_total_N_uptake():

    sumN = {'Density':[], 'Leaf posture':[], 'Stand':[], 'Sky':[], 'sumN':[]}

    for density in DENSITIES:
        for treatment in TREATMENTS:
            if (density, treatment) == (410, 'Plano_mixte'): continue
            data_dirpath = os.path.join(DATA_DIRPATH, 'Densite_{}'.format(density), treatment, 'outputs\organs_states.csv')
            data_df = pd.read_csv(data_dirpath)
            data_df = data_df[data_df['organ'] == 'roots']


            if treatment not in ('Asso_direct', 'Asso_mixte', 'Asso_diffus'):
                sumN['Density'].append(density)
                leaf_posture, sky = treatment.split('_')
                sumN['Leaf posture'].append(leaf_posture)
                sumN['Stand'].append('pure')
                sumN['Sky'].append(sky)
                sumN['sumN'].append(data_df.groupby('plant')['Uptake_Nitrates'].sum().mean())

            else:
                stand, sky = treatment.split('_')
                sumN['Density'].append(density)
                sumN['Density'].append(density)
                sumN['Stand'].append(stand)
                sumN['Stand'].append(stand)
                sumN['Leaf posture'].append('Erect')
                sumN['Leaf posture'].append('Plano')
                sumN['Sky'].append(sky)
                sumN['Sky'].append(sky)
                sumN['sumN'].append(data_df[(1 <= data_df['plant']) & (data_df['plant'] <= 5)].groupby('plant')['Uptake_Nitrates'].sum().mean())
                sumN['sumN'].append(data_df[(6 <= data_df['plant']) & (data_df['plant'] <= 10)].groupby('plant')['Uptake_Nitrates'].sum().mean())

    sumN_df = pd.DataFrame(sumN)
    sumN_df.to_csv(os.path.join(GRAPHS_DIRPATH, 'Table_total_N_uptake.csv'), index=False)


def RUE_end_grains():

    plt.figure(figsize=cm2inch(30, 25))
    gs = gridspec.GridSpec(3, 4, hspace=0, wspace=0)
    axes = []

    t_end = 900
    phloem_shoot_root = 0.75

    RUE = {'Density': [], 'Leaf posture': [], 'Stand': [], 'Sky': [], 'RG (MJ)': [], 'BM_shoot_cum': [], 'BM_total_cum': [], 'RUE shoot (g MJ-1)': [], 'RUE total (g MJ-1)': []}
    RUE_daily = {'Density': [], 'Leaf posture': [], 'Stand': [], 'Sky': [], 'day': [], 'RUE': [], 'dry_mass_increment': [], 'RG_sum': []}

    for density in DENSITIES:
        for treatment in TREATMENTS:
            data_dirpath_ph = os.path.join(DATA_DIRPATH, 'Densite_{}'.format(density), treatment,
                                        'outputs\elements_states.csv')
            data_dirpath_organs = os.path.join(DATA_DIRPATH, 'Densite_{}'.format(density), treatment,
                                        'outputs\organs_states.csv')
            data_df_ph = pd.read_csv(data_dirpath_ph)
            data_df_ph = data_df_ph[data_df_ph['t'] <= t_end]
            data_df_ph['Total_PARa'] = data_df_ph['PARa'] * data_df_ph['green_area']
            data_df_ph['day'] = data_df_ph['t'] // 24+1
            data_df_ph['sum_dry_mass'] = ((data_df_ph['triosesP'] * 1E-6*C_MOLAR_MASS) / TRIOSESP_MOLAR_MASS_C_RATIO +
                                       (data_df_ph['sucrose'] * 1E-6*C_MOLAR_MASS) / SUCROSE_MOLAR_MASS_C_RATIO +
                                       (data_df_ph['starch'] * 1E-6*C_MOLAR_MASS) / HEXOSE_MOLAR_MASS_C_RATIO +
                                       (data_df_ph['fructan'] * 1E-6*C_MOLAR_MASS) / HEXOSE_MOLAR_MASS_C_RATIO +
                                       (data_df_ph['nitrates'] * 1E-6*N_MOLAR_MASS) / NITRATES_MOLAR_MASS_N_RATIO +
                                       (data_df_ph['amino_acids'] * 1E-6*N_MOLAR_MASS) / AMINO_ACIDS_MOLAR_MASS_N_RATIO +
                                       (data_df_ph['proteins'] * 1E-6*N_MOLAR_MASS) / AMINO_ACIDS_MOLAR_MASS_N_RATIO)

            data_df_organs = pd.read_csv(data_dirpath_organs)
            data_df_organs = data_df_organs[data_df_organs['t'] <= t_end]
            data_df_organs['day'] = data_df_organs['t'] // 24 + 1
            data_df_organs['sum_dry_mass'] = (((data_df_organs.fillna(0)['structure'] + data_df_organs.fillna(0)['starch']) * 1E-6 * 12 / 0.4) +
                                     ((data_df_organs.fillna(0)['proteins'] * 1E-6 * 14) / 0.136) +
                                     (data_df_organs.fillna(0)['sucrose'] * 1E-6 * C_MOLAR_MASS) / SUCROSE_MOLAR_MASS_C_RATIO +
                                     (data_df_organs.fillna(0)['amino_acids'] * 1E-6 * N_MOLAR_MASS) / AMINO_ACIDS_MOLAR_MASS_N_RATIO +
                                     (data_df_organs.fillna(0)['nitrates'] * 1E-6 * N_MOLAR_MASS) / NITRATES_MOLAR_MASS_N_RATIO +
                                              data_df_organs.fillna(0)['mstruct'])

            if treatment not in ('Asso_direct', 'Asso_mixte', 'Asso_diffus'):
                RUE['Density'].append(density)
                leaf_posture, sky = treatment.split('_')
                RUE['Leaf posture'].append(leaf_posture)
                RUE['Stand'].append('pure')
                RUE['Sky'].append(sky)

                # Ph
                sum_dry_mass_org_ph_0 = data_df_ph[data_df_ph['t'] == 0][['plant', 'sum_dry_mass']].groupby('plant').sum().mean().iloc[0]
                sum_dry_mass_org_ph_end = data_df_ph[data_df_ph['t'] == t_end][['plant', 'sum_dry_mass']].groupby('plant').sum().mean().iloc[0]

                # Grains
                sum_dry_mass_grains_0 = data_df_organs[(data_df_organs['organ'] == 'grains') & (data_df_organs['t'] == 0)][['plant', 'sum_dry_mass']].groupby('plant').sum().mean().iloc[0]
                sum_dry_mass_grains_end = data_df_organs[(data_df_organs['organ'] == 'grains') & (data_df_organs['t'] == t_end)][['plant', 'sum_dry_mass']].groupby('plant').sum().mean().iloc[0]

                # Roots
                sum_dry_mass_roots_0 = data_df_organs[(data_df_organs['organ'] == 'roots') & (data_df_organs['t'] == 0)][['plant', 'sum_dry_mass']].groupby('plant').sum().mean().iloc[0]
                sum_dry_mass_roots_end = data_df_organs[(data_df_organs['organ'] == 'roots') & (data_df_organs['t'] == t_end)][['plant', 'sum_dry_mass']].groupby('plant').sum().mean().iloc[0]

                # Phloem
                sum_dry_mass_phloem_0_shoot = data_df_organs[(data_df_organs['organ'] == 'phloem') & (data_df_organs['t'] == 0)][['plant', 'sum_dry_mass']].groupby('plant').sum().mean().iloc[0] * phloem_shoot_root
                sum_dry_mass_phloem_end_shoot = data_df_organs[(data_df_organs['organ'] == 'phloem') & (data_df_organs['t'] == t_end)][['plant', 'sum_dry_mass']].groupby('plant').sum().mean().iloc[0] * phloem_shoot_root

                sum_dry_mass_phloem_0_roots = data_df_organs[(data_df_organs['organ'] == 'phloem') & (data_df_organs['t'] == 0)][['plant', 'sum_dry_mass']].groupby('plant').sum().mean().iloc[0] * (1 - phloem_shoot_root)
                sum_dry_mass_phloem_end_roots = data_df_organs[(data_df_organs['organ'] == 'phloem') & (data_df_organs['t'] == t_end)][['plant', 'sum_dry_mass']].groupby('plant').sum().mean().iloc[0] * (1 - phloem_shoot_root)

                # Total shoot
                total_0 = sum_dry_mass_org_ph_0 + sum_dry_mass_grains_0 + sum_dry_mass_phloem_0_shoot
                total_end = sum_dry_mass_org_ph_end + sum_dry_mass_grains_end + sum_dry_mass_phloem_end_shoot
                dry_mass_increment = total_end - total_0
                RUE['BM_shoot_cum'].append(dry_mass_increment)
                mole_PAR_to_Watt_RG = 1 / 2.02
                RUE['RG (MJ)'].append(data_df_ph[(0 <= data_df_ph['t']) & (data_df_ph['t'] <= t_end)][['plant', 'Total_PARa']].groupby('plant').sum().mean().iloc[0] * 3600 * 1E-6 * mole_PAR_to_Watt_RG)
                RUE['RUE shoot (g MJ-1)'].append(dry_mass_increment / (
                            data_df_ph[(0 <= data_df_ph['t']) & (data_df_ph['t'] <= t_end)][['plant', 'Total_PARa']].groupby('plant').sum().mean().iloc[0] * 3600 * 1E-6 * mole_PAR_to_Watt_RG))

                # Total shoot + roots
                total_0 = sum_dry_mass_org_ph_0 + sum_dry_mass_grains_0 + sum_dry_mass_phloem_0_shoot + sum_dry_mass_phloem_0_roots + sum_dry_mass_roots_0
                total_end = sum_dry_mass_org_ph_end + sum_dry_mass_grains_end + sum_dry_mass_phloem_end_shoot + sum_dry_mass_phloem_end_roots + sum_dry_mass_roots_end
                dry_mass_increment = total_end - total_0
                RUE['BM_total_cum'].append(dry_mass_increment)
                RUE['RUE total (g MJ-1)'].append(dry_mass_increment / (
                            data_df_ph[(0 <= data_df_ph['t']) & (data_df_ph['t'] <= t_end)][['plant', 'Total_PARa']].groupby('plant').sum().mean().iloc[0] * 3600 * 1E-6 * mole_PAR_to_Watt_RG))

            else:
                stand, sky = treatment.split('_')
                RUE['Density'].append(density)
                RUE['Density'].append(density)
                RUE['Leaf posture'].append('Erect')
                RUE['Leaf posture'].append('Plano')
                RUE['Stand'].append(stand)
                RUE['Stand'].append(stand)
                RUE['Sky'].append(sky)
                RUE['Sky'].append(sky)

                # ERECT
                # Ph
                sum_dry_mass_org_ph_0 = data_df_ph[(data_df_ph['t'] == 0) & (1 <= data_df_ph['plant']) & (data_df_ph['plant'] <= 5)][['plant', 'sum_dry_mass']].groupby('plant').sum().mean().iloc[0]
                sum_dry_mass_org_ph_end = data_df_ph[(data_df_ph['t'] == t_end) & (1 <= data_df_ph['plant']) & (data_df_ph['plant'] <= 5)][['plant', 'sum_dry_mass']].groupby('plant').sum().mean().iloc[0]

                # Grains
                sum_dry_mass_grains_0 = data_df_organs[(data_df_organs['organ'] == 'grains') & (data_df_organs['t'] == 0) & (1 <= data_df_organs['plant']) & (data_df_organs['plant'] <= 5)][['plant', 'sum_dry_mass']].groupby('plant').sum().mean().iloc[0]
                sum_dry_mass_grains_end = data_df_organs[(data_df_organs['organ'] == 'grains') & (data_df_organs['t'] == t_end) & (1 <= data_df_organs['plant']) & (data_df_organs['plant'] <= 5)][['plant', 'sum_dry_mass']].groupby('plant').sum().mean().iloc[0]

                # Roots
                sum_dry_mass_roots_0 = data_df_organs[(data_df_organs['organ'] == 'roots') & (data_df_organs['t'] == 0) & (1 <= data_df_organs['plant']) & (data_df_organs['plant'] <= 5)][['plant', 'sum_dry_mass']].groupby('plant').sum().mean().iloc[0]
                sum_dry_mass_roots_end = data_df_organs[(data_df_organs['organ'] == 'roots') & (data_df_organs['t'] == t_end) & (1 <= data_df_organs['plant']) & (data_df_organs['plant'] <= 5)][['plant', 'sum_dry_mass']].groupby('plant').sum().mean().iloc[0]

                # Phloem
                sum_dry_mass_phloem_0_shoot = data_df_organs[(data_df_organs['organ'] == 'phloem') & (data_df_organs['t'] == 0) & (1 <= data_df_organs['plant']) & (data_df_organs['plant'] <= 5)][['plant', 'sum_dry_mass']].groupby('plant').sum().mean().iloc[0] * phloem_shoot_root
                sum_dry_mass_phloem_end_shoot = data_df_organs[(data_df_organs['organ'] == 'phloem') & (data_df_organs['t'] == t_end) & (1 <= data_df_organs['plant']) & (data_df_organs['plant'] <= 5)][['plant', 'sum_dry_mass']].groupby('plant').sum().mean().iloc[0] * phloem_shoot_root

                sum_dry_mass_phloem_0_roots = data_df_organs[(data_df_organs['organ'] == 'phloem') & (data_df_organs['t'] == 0) & (1 <= data_df_organs['plant']) & (data_df_organs['plant'] <= 5)][['plant', 'sum_dry_mass']].groupby('plant').sum().mean().iloc[0] * (1 - phloem_shoot_root)
                sum_dry_mass_phloem_end_roots = data_df_organs[(data_df_organs['organ'] == 'phloem') & (data_df_organs['t'] == t_end) & (1 <= data_df_organs['plant']) & (data_df_organs['plant'] <= 5)][['plant', 'sum_dry_mass']].groupby('plant').sum().mean().iloc[0] * (1 - phloem_shoot_root)

                # Total shoot
                total_0 = sum_dry_mass_org_ph_0 + sum_dry_mass_grains_0 + sum_dry_mass_phloem_0_shoot
                total_end = sum_dry_mass_org_ph_end + sum_dry_mass_grains_end + sum_dry_mass_phloem_end_shoot
                dry_mass_increment = total_end - total_0
                RUE['BM_shoot_cum'].append(dry_mass_increment)
                mole_PAR_to_Watt_RG = 1 / 2.02
                RUE['RG (MJ)'].append(data_df_ph[(0 <= data_df_ph['t']) & (data_df_ph['t'] <= t_end) & (1 <= data_df_ph['plant']) & (data_df_ph['plant'] <= 5)][['plant', 'Total_PARa']].groupby('plant').sum().mean().iloc[0] * 3600 * 1E-6 * mole_PAR_to_Watt_RG)
                RUE['RUE shoot (g MJ-1)'].append(dry_mass_increment / (data_df_ph[(0 <= data_df_ph['t']) & (data_df_ph['t'] <= t_end) & (1 <= data_df_ph['plant']) & (data_df_ph['plant'] <= 5)][['plant', 'Total_PARa']].groupby('plant').sum().mean().iloc[0] * 3600 * 1E-6 * mole_PAR_to_Watt_RG))

                # Total shoot + roots
                total_0 = sum_dry_mass_org_ph_0 + sum_dry_mass_grains_0 + sum_dry_mass_phloem_0_shoot + sum_dry_mass_phloem_0_roots + sum_dry_mass_roots_0
                total_end = sum_dry_mass_org_ph_end + sum_dry_mass_grains_end + sum_dry_mass_phloem_end_shoot + sum_dry_mass_phloem_end_roots + sum_dry_mass_roots_end
                dry_mass_increment = total_end - total_0
                RUE['BM_total_cum'].append(dry_mass_increment)
                RUE['RUE total (g MJ-1)'].append(dry_mass_increment / (data_df_ph[(0 <= data_df_ph['t']) & (data_df_ph['t'] <= t_end) & (1 <= data_df_ph['plant']) & (data_df_ph['plant'] <= 5)][['plant', 'Total_PARa']].groupby('plant').sum().mean().iloc[0] * 3600 * 1E-6 * mole_PAR_to_Watt_RG))

                # PLANO
                # Ph
                sum_dry_mass_org_ph_0 = data_df_ph[(data_df_ph['t'] == 0) & (6 <= data_df_ph['plant']) & (data_df_ph['plant'] <= 10)][['plant', 'sum_dry_mass']].groupby('plant').sum().mean().iloc[0]
                sum_dry_mass_org_ph_end = data_df_ph[(data_df_ph['t'] == t_end) & (6 <= data_df_ph['plant']) & (data_df_ph['plant'] <= 10)][['plant', 'sum_dry_mass']].groupby('plant').sum().mean().iloc[0]

                # Grains
                sum_dry_mass_grains_0 = data_df_organs[(data_df_organs['organ'] == 'grains') & (data_df_organs['t'] == 0) & (6 <= data_df_organs['plant']) & (data_df_organs['plant'] <= 10)][['plant', 'sum_dry_mass']].groupby('plant').sum().mean().iloc[0]
                sum_dry_mass_grains_end = data_df_organs[(data_df_organs['organ'] == 'grains') & (data_df_organs['t'] == t_end) & (6 <= data_df_organs['plant']) & (data_df_organs['plant'] <= 10)][['plant', 'sum_dry_mass']].groupby('plant').sum().mean().iloc[0]

                # Roots
                sum_dry_mass_roots_0 = data_df_organs[(data_df_organs['organ'] == 'roots') & (data_df_organs['t'] == 0) & (6 <= data_df_organs['plant']) & (data_df_organs['plant'] <= 10)][['plant', 'sum_dry_mass']].groupby('plant').sum().mean().iloc[0]
                sum_dry_mass_roots_end = data_df_organs[(data_df_organs['organ'] == 'roots') & (data_df_organs['t'] == t_end) & (6 <= data_df_organs['plant']) & (data_df_organs['plant'] <= 10)][['plant', 'sum_dry_mass']].groupby('plant').sum().mean().iloc[0]

                # Phloem
                sum_dry_mass_phloem_0_shoot = data_df_organs[(data_df_organs['organ'] == 'phloem') & (data_df_organs['t'] == 0) & (6 <= data_df_organs['plant']) & (data_df_organs['plant'] <= 10)][['plant', 'sum_dry_mass']].groupby('plant').sum().mean().iloc[0] * phloem_shoot_root
                sum_dry_mass_phloem_end_shoot = data_df_organs[(data_df_organs['organ'] == 'phloem') & (data_df_organs['t'] == t_end) & (6 <= data_df_organs['plant']) & (data_df_organs['plant'] <= 10)][['plant', 'sum_dry_mass']].groupby('plant').sum().mean().iloc[0] * phloem_shoot_root

                sum_dry_mass_phloem_0_roots = data_df_organs[(data_df_organs['organ'] == 'phloem') & (data_df_organs['t'] == 0) & (6 <= data_df_organs['plant']) & (data_df_organs['plant'] <= 10)][['plant', 'sum_dry_mass']].groupby('plant').sum().mean().iloc[0] * (1 - phloem_shoot_root)
                sum_dry_mass_phloem_end_roots = data_df_organs[(data_df_organs['organ'] == 'phloem') & (data_df_organs['t'] == t_end) & (6 <= data_df_organs['plant']) & (data_df_organs['plant'] <= 10)][['plant', 'sum_dry_mass']].groupby('plant').sum().mean().iloc[0] * (1 - phloem_shoot_root)

                # Total shoot
                total_0 = sum_dry_mass_org_ph_0 + sum_dry_mass_grains_0 + sum_dry_mass_phloem_0_shoot
                total_end = sum_dry_mass_org_ph_end + sum_dry_mass_grains_end + sum_dry_mass_phloem_end_shoot
                dry_mass_increment = total_end - total_0
                RUE['BM_shoot_cum'].append(dry_mass_increment)
                mole_PAR_to_Watt_RG = 1 / 2.02
                RUE['RG (MJ)'].append(data_df_ph[(0 <= data_df_ph['t']) & (data_df_ph['t'] <= t_end) & (6 <= data_df_ph['plant']) & (data_df_ph['plant'] <= 10)][['plant', 'Total_PARa']].groupby('plant').sum().mean().iloc[0] * 3600 * 1E-6 * mole_PAR_to_Watt_RG)
                RUE['RUE shoot (g MJ-1)'].append(dry_mass_increment / (data_df_ph[(0 <= data_df_ph['t']) & (data_df_ph['t'] <= t_end) & (6 <= data_df_ph['plant']) & (data_df_ph['plant'] <= 10)][['plant', 'Total_PARa']].groupby('plant').sum().mean().iloc[0] * 3600 * 1E-6 * mole_PAR_to_Watt_RG))

                # Total shoot + roots
                total_0 = sum_dry_mass_org_ph_0 + sum_dry_mass_grains_0 + sum_dry_mass_phloem_0_shoot + sum_dry_mass_phloem_0_roots + sum_dry_mass_roots_0
                total_end = sum_dry_mass_org_ph_end + sum_dry_mass_grains_end + sum_dry_mass_phloem_end_shoot + sum_dry_mass_phloem_end_roots + sum_dry_mass_roots_end
                dry_mass_increment = total_end - total_0
                RUE['BM_total_cum'].append(dry_mass_increment)
                RUE['RUE total (g MJ-1)'].append(dry_mass_increment / (data_df_ph[(0 <= data_df_ph['t']) & (data_df_ph['t'] <= t_end) & (6 <= data_df_ph['plant']) & (data_df_ph['plant'] <= 10)][['plant', 'Total_PARa']].groupby('plant').sum().mean().iloc[0] * 3600 * 1E-6 * mole_PAR_to_Watt_RG))

            # Graph
            for day, data in data_df_ph.groupby('day'):
                if treatment not in ('Asso_direct', 'Asso_mixte', 'Asso_diffus'):
                    grains_dry_mass = data_df_organs[(data_df_organs['organ'] == 'grains') & (data_df_organs['day'] == day)].groupby('t')
                    grains_dry_mass_begin = grains_dry_mass['sum_dry_mass'].mean().iloc[0]
                    grains_dry_mass_end = grains_dry_mass['sum_dry_mass'].mean().iloc[-1]

                    phloem = data_df_organs[(data_df_organs['organ'] == 'phloem') & (data_df_organs['day'] == day)].groupby('t')
                    phloem_dry_mass_begin = phloem['sum_dry_mass'].mean().iloc[0] * phloem_shoot_root
                    phloem_dry_mass_end = phloem['sum_dry_mass'].mean().iloc[-1] * phloem_shoot_root

                    hourly_ph_dry_mass = data.groupby(['t', 'plant'])
                    ph_dry_mass_begin = hourly_ph_dry_mass['sum_dry_mass'].sum().groupby('t').mean().iloc[0]
                    ph_dry_mass_end = hourly_ph_dry_mass['sum_dry_mass'].sum().groupby('t').mean().iloc[0]

                    dry_mass_increment = (ph_dry_mass_end + grains_dry_mass_end + phloem_dry_mass_end) - (ph_dry_mass_begin + grains_dry_mass_begin + phloem_dry_mass_begin)

                    mole_PAR_to_Watt_RG = 1 / 2.02
                    RG_sum = data[['plant', 'Total_PARa']].groupby('plant').sum().mean().iloc[0] * 3600 * 1E-6 * mole_PAR_to_Watt_RG

                    RUE_daily['Density'].append(density)
                    leaf_posture, sky = treatment.split('_')
                    RUE_daily['Leaf posture'].append(leaf_posture)
                    RUE_daily['Stand'].append('pure')
                    RUE_daily['Sky'].append(sky)
                    RUE_daily['day'].append(day)
                    RUE_daily['RG_sum'].append(RG_sum)
                    RUE_daily['dry_mass_increment'].append(dry_mass_increment)
                    RUE_daily['RUE'].append(dry_mass_increment / RG_sum)

                else:
                    # Erect
                    grains_dry_mass = data_df_organs[(data_df_organs['organ'] == 'grains') & (data_df_organs['day'] == day) & (1 <= data_df_organs['plant']) & (data_df_organs['plant'] <= 5)].groupby('t')
                    grains_dry_mass_begin = grains_dry_mass['sum_dry_mass'].mean().iloc[0]
                    grains_dry_mass_end = grains_dry_mass['sum_dry_mass'].mean().iloc[-1]

                    phloem = data_df_organs[(data_df_organs['organ'] == 'phloem') & (data_df_organs['day'] == day) & (1 <= data_df_organs['plant']) & (data_df_organs['plant'] <= 5)].groupby('t')
                    phloem_dry_mass_begin = phloem['sum_dry_mass'].mean().iloc[0] * phloem_shoot_root
                    phloem_dry_mass_end = phloem['sum_dry_mass'].mean().iloc[-1] * phloem_shoot_root

                    hourly_ph_dry_mass = data[(1 <= data['plant']) & (data['plant'] <= 5)].groupby(['t', 'plant'])
                    ph_dry_mass_begin = hourly_ph_dry_mass['sum_dry_mass'].sum().groupby('t').mean().iloc[0]
                    ph_dry_mass_end = hourly_ph_dry_mass['sum_dry_mass'].sum().groupby('t').mean().iloc[-1]

                    dry_mass_increment = (ph_dry_mass_end + grains_dry_mass_end + phloem_dry_mass_end) - (ph_dry_mass_begin + grains_dry_mass_begin + phloem_dry_mass_begin)

                    mole_PAR_to_Watt_RG = 1 / 2.02
                    RG_sum = data[(1 <= data['plant']) & (data['plant'] <= 5)][['plant', 'Total_PARa']].groupby('plant').sum().mean().iloc[0] * 3600 * 1E-6 * mole_PAR_to_Watt_RG

                    RUE_daily['Density'].append(density)
                    RUE_daily['Leaf posture'].append('erectophile')
                    RUE_daily['Stand'].append('mixture')
                    leaf_posture, sky = treatment.split('_')
                    RUE_daily['Sky'].append(sky)
                    RUE_daily['day'].append(day)
                    RUE_daily['RG_sum'].append(RG_sum)
                    RUE_daily['dry_mass_increment'].append(dry_mass_increment)
                    RUE_daily['RUE'].append(dry_mass_increment / RG_sum)

                    # Plano
                    grains_dry_mass = data_df_organs[(data_df_organs['organ'] == 'grains') & (data_df_organs['day'] == day) & (6 <= data_df_organs['plant']) & (data_df_organs['plant'] <= 10)].groupby('t')
                    grains_dry_mass_begin = grains_dry_mass['sum_dry_mass'].mean().iloc[0]
                    grains_dry_mass_end = grains_dry_mass['sum_dry_mass'].mean().iloc[-1]

                    phloem = data_df_organs[(data_df_organs['organ'] == 'phloem') & (data_df_organs['day'] == day) & (6 <= data_df_organs['plant']) & (data_df_organs['plant'] <= 10)].groupby('t')
                    phloem_dry_mass_begin = phloem['sum_dry_mass'].mean().iloc[0] * phloem_shoot_root
                    phloem_dry_mass_end = phloem['sum_dry_mass'].mean().iloc[-1] * phloem_shoot_root

                    hourly_ph_dry_mass = data[(6 <= data['plant']) & (data['plant'] <= 10)].groupby(['t', 'plant'])
                    ph_dry_mass_begin = hourly_ph_dry_mass['sum_dry_mass'].sum().groupby('t').mean().iloc[0]
                    ph_dry_mass_end = hourly_ph_dry_mass['sum_dry_mass'].sum().groupby('t').mean().iloc[-1]

                    dry_mass_increment = (ph_dry_mass_end + grains_dry_mass_end + phloem_dry_mass_end) - (ph_dry_mass_begin + grains_dry_mass_begin + phloem_dry_mass_begin)

                    mole_PAR_to_Watt_RG = 1 / 2.02
                    RG_sum = data[(6 <= data['plant']) & (data['plant'] <= 10)][['plant', 'Total_PARa']].groupby('plant').sum().mean().iloc[0] * 3600 * 1E-6 * mole_PAR_to_Watt_RG

                    RUE_daily['Density'].append(density)
                    RUE_daily['Leaf posture'].append('planophile')
                    RUE_daily['Stand'].append('mixture')
                    RUE_daily['Sky'].append(sky)
                    RUE_daily['day'].append(day)
                    RUE_daily['RG_sum'].append(RG_sum)
                    RUE_daily['dry_mass_increment'].append(dry_mass_increment)
                    RUE_daily['RUE'].append(dry_mass_increment / RG_sum)

    RUE_df = pd.DataFrame(RUE)
    RUE_daily_df = pd.DataFrame(RUE_daily)
    mean_RUE_daily = RUE_daily_df.groupby('Density')

    for density, data in mean_RUE_daily:
        ax1 = plt.subplot(gs[0, AXES_MAP[density]])
        ax2 = plt.subplot(gs[1, AXES_MAP[density]])
        ax3 = plt.subplot(gs[2, AXES_MAP[density]])
        axes.append(ax1)
        axes.append(ax2)
        axes.append(ax3)

        mean_dry_mass_increment = data.groupby('day')['dry_mass_increment'].mean()
        mean_RG_sum = data.groupby('day')['RG_sum'].mean()
        mean_RUE = data.groupby('day')['RUE'].mean()
        std_dry_mass_increment = data.groupby('day')['dry_mass_increment'].std()
        std_RG_sum = data.groupby('day')['RG_sum'].std()
        std_RUE = data.groupby('day')['RUE'].std()

        ax1.fill_between(data['day'].unique(), mean_dry_mass_increment - std_dry_mass_increment,
                         mean_dry_mass_increment + std_dry_mass_increment, alpha=alpha*3, color='b', linewidth=thickness)
        ax1.set_yticks(np.arange(-0.05, 0.2, 0.05))

        ax2.fill_between(data['day'].unique(), mean_RG_sum - std_RG_sum,
                         mean_RG_sum + std_RG_sum, alpha=alpha*3, color='b', linewidth=thickness)
        ax2.set_yticks(np.arange(0, 0.2, 0.05))

        ax3.fill_between(data['day'].unique(), mean_RUE - std_RUE,
                         mean_RUE + std_RUE, alpha=alpha*3, color='b', linewidth=thickness)

        # Mean RUE
        total_mean_RUE = RUE_df[RUE_df['Density'] == density]['RUE shoot (g MJ-1)'].mean()
        total_std_RUE = RUE_df[RUE_df['Density'] == density]['RUE shoot (g MJ-1)'].std()
        ax3.plot([0, 38], [total_mean_RUE]*2, '--', color='k', linewidth=thickness+1)
        ax3.set_yticks(np.arange(-3, 6, 2))

        xtext = 16
        ytext = 3.5
        ax3.text(xtext, ytext, u'µ = {} \u00B1 {}'.format(round(total_mean_RUE, 2), round(total_std_RUE, 2)), fontsize=fontsize + 2, verticalalignment='top')
        plt.subplot(gs[0, 0]).set_zorder(1)

    for ax in axes:
        ax.set_xticks(np.arange(0, 50, 10))
        [ax.yaxis.set_label_coords(labelx, 0.5) for i in ax.spines.itervalues()]
        ax.axvline(360 / 24., color='k', linestyle='--')
        ax.tick_params(labelbottom='off', labelleft='off')

        if ax in (plt.subplot(gs[2, 0]), plt.subplot(gs[2, 1]), plt.subplot(gs[2, 2]), plt.subplot(gs[2, 3])):
            ax.tick_params(labelbottom='on')
            ax.set_xlabel('Time from flowering (day)', fontsize=fontsize)
            if ax != plt.subplot(gs[2, 3]):
                ax.set_xticks(ax.get_xticks()[:-1])
        if ax in (plt.subplot(gs[0, 0]), plt.subplot(gs[1, 0]), plt.subplot(gs[2, 0])):
            ax.tick_params(labelleft='on')
            if ax == plt.subplot(gs[0, 0]):
                ax.set_ylabel(u'delta dry mass (g)', fontsize=fontsize)
            elif ax == plt.subplot(gs[1, 0]):
                ax.set_ylabel(u'RGa (MJ)', fontsize=fontsize)
                ax.set_yticks(ax.get_yticks()[:-1])
            else:
                ax.set_ylabel(u'RUE (g MJ$^{-1}$)', fontsize=fontsize)
                ax.set_yticks(ax.get_yticks()[:-1])

    xtext = 10
    ytext = 0.175
    plt.subplot(gs[0, 0]).text(xtext, ytext, 'Density 200', fontsize=fontsize + 2, backgroundcolor='w',
                               verticalalignment='top')
    plt.subplot(gs[0, 0]).set_zorder(1)
    plt.subplot(gs[0, 1]).text(xtext, ytext, 'Density 410', fontsize=fontsize + 2, backgroundcolor='w',
                               verticalalignment='top')
    plt.subplot(gs[0, 1]).set_zorder(1)
    plt.subplot(gs[0, 2]).text(xtext, ytext, 'Density 600', fontsize=fontsize + 2, backgroundcolor='w',
                               verticalalignment='top')
    plt.subplot(gs[0, 2]).set_zorder(1)
    plt.subplot(gs[0, 3]).text(xtext, ytext, 'Density 800', fontsize=fontsize + 2, backgroundcolor='w',
                               verticalalignment='top')


    plt.savefig(os.path.join(GRAPHS_DIRPATH, 'RUE_end_grains.TIFF'), dpi=300, format='TIFF',bbox_inches='tight')
    plt.close()

    RUE_df = RUE_df[['Density', 'Leaf posture', 'Stand', 'Sky', 'RG (MJ)', 'BM_shoot_cum', 'BM_total_cum', 'RUE shoot (g MJ-1)', 'RUE total (g MJ-1)']]
    RUE_df.to_csv(os.path.join(GRAPHS_DIRPATH, 'RUE_end_grains.csv'), index=False)

def figure_roots():

    plt.figure(figsize=cm2inch(30, 20))
    gs = gridspec.GridSpec(5, 4, hspace=0, wspace=0)

    for density in DENSITIES:
        for treatment in TREATMENTS:
            if treatment not in ('Plano_diffus', 'Erect_diffus', 'Asso_diffus'): continue
            stand = treatment.split('_')[0]

            axes = []

            # dry mass
            organs_path = os.path.join(DATA_DIRPATH, 'Densite_{}'.format(density), treatment,
                                       'outputs\organs_states.csv')
            organs = pd.read_csv(organs_path)
            organs_roots = organs[organs['organ'] == 'roots']
            organs_roots['Conc_sucrose'] = organs_roots['sucrose'] / organs_roots['mstruct']
            organs_roots['Conc_amino_acids'] = organs_roots['amino_acids'] / organs_roots['mstruct']
            organs_roots['Conc_nitrates'] = organs_roots['nitrates'] / organs_roots['mstruct']

            ax1 = plt.subplot(gs[0, AXES_MAP[density]])
            ax2 = plt.subplot(gs[1, AXES_MAP[density]])
            ax3 = plt.subplot(gs[2, AXES_MAP[density]])
            ax4 = plt.subplot(gs[3, AXES_MAP[density]])
            ax5 = plt.subplot(gs[4, AXES_MAP[density]])
            axes.append(ax1)
            axes.append(ax2)
            axes.append(ax3)
            axes.append(ax4)
            axes.append(ax5)

            stand_data = {}
            if stand == 'Asso':
                if density in ('200', '410', '600', '800'):
                    # Erectophile
                    mstruct = (organs_roots[organs_roots['plant'] == 1].groupby(['t', 'plant'])['mstruct'].aggregate(
                        np.sum)).groupby(['t'])
                    sucrose = (organs_roots[organs_roots['plant'] == 1].groupby(['t', 'plant'])['Conc_sucrose'].aggregate(
                        np.sum)).groupby(['t'])
                    amino_acids = (organs_roots[organs_roots['plant'] == 1].groupby(['t', 'plant'])['Conc_amino_acids'].aggregate(
                        np.sum)).groupby(['t'])
                    nitrates = (organs_roots[organs_roots['plant'] == 1].groupby(['t', 'plant'])['Conc_nitrates'].aggregate(
                        np.sum)).groupby(['t'])
                    Unloading_Sucrose = (organs_roots[organs_roots['plant'] == 1].groupby(['t', 'plant'])['Unloading_Sucrose'].aggregate(
                        np.sum)).groupby(['t'])
                    stand_data['Caphorn_Mixed'] = [mstruct, sucrose, amino_acids, nitrates, Unloading_Sucrose]

                    # Planophile
                    mstruct = (organs_roots[organs_roots['plant'] == 2].groupby(['t', 'plant'])['mstruct'].aggregate(
                        np.sum)).groupby(['t'])
                    sucrose = (organs_roots[organs_roots['plant'] == 2].groupby(['t', 'plant'])['Conc_sucrose'].aggregate(
                        np.sum)).groupby(['t'])
                    amino_acids = (organs_roots[organs_roots['plant'] == 2].groupby(['t', 'plant'])['Conc_amino_acids'].aggregate(
                        np.sum)).groupby(['t'])
                    nitrates = (organs_roots[organs_roots['plant'] == 2].groupby(['t', 'plant'])['Conc_nitrates'].aggregate(
                        np.sum)).groupby(['t'])
                    Unloading_Sucrose = (organs_roots[organs_roots['plant'] == 2].groupby(['t', 'plant'])['Unloading_Sucrose'].aggregate(
                        np.sum)).groupby(['t'])
                    stand_data['Soissons_Mixed'] = [mstruct, sucrose, amino_acids, nitrates, Unloading_Sucrose]

                else:
                    # Erectophile
                    mstruct = (organs_roots[(1 <= organs_roots['plant']) & (organs_roots['plant'] <= 5)].groupby(['t', 'plant'])['mstruct'].aggregate(
                        np.sum)).groupby(['t'])
                    sucrose = (organs_roots[(1 <= organs_roots['plant']) & (organs_roots['plant'] <= 5)].groupby(['t', 'plant'])['Conc_sucrose'].aggregate(
                        np.sum)).groupby(['t'])
                    amino_acids = (organs_roots[(1 <= organs_roots['plant']) & (organs_roots['plant'] <= 5)].groupby(['t', 'plant'])['Conc_amino_acids'].aggregate(
                        np.sum)).groupby(['t'])
                    nitrates = (organs_roots[(1 <= organs_roots['plant']) & (organs_roots['plant'] <= 5)].groupby(['t', 'plant'])['Conc_nitrates'].aggregate(
                        np.sum)).groupby(['t'])
                    Unloading_Sucrose = (organs_roots[(1 <= organs_roots['plant']) & (organs_roots['plant'] <= 5)].groupby(['t', 'plant'])['Unloading_Sucrose'].aggregate(
                        np.sum)).groupby(['t'])
                    stand_data['Caphorn_Mixed'] = [mstruct, sucrose, amino_acids, nitrates, Unloading_Sucrose]

                    # Planophile
                    mstruct = (organs_roots[(6 <= organs_roots['plant']) & (organs_roots['plant'] <= 10)].groupby(['t', 'plant'])['mstruct'].aggregate(
                        np.sum)).groupby(['t'])
                    sucrose = (organs_roots[(6 <= organs_roots['plant']) & (organs_roots['plant'] <= 10)].groupby(['t', 'plant'])['Conc_sucrose'].aggregate(
                        np.sum)).groupby(['t'])
                    amino_acids = (organs_roots[(6 <= organs_roots['plant']) & (organs_roots['plant'] <= 10)].groupby(['t', 'plant'])['Conc_amino_acids'].aggregate(
                        np.sum)).groupby(['t'])
                    nitrates = (organs_roots[(6 <= organs_roots['plant']) & (organs_roots['plant'] <= 10)].groupby(['t', 'plant'])['Conc_nitrates'].aggregate(
                        np.sum)).groupby(['t'])
                    Unloading_Sucrose = (organs_roots[(6 <= organs_roots['plant']) & (organs_roots['plant'] <= 10)].groupby(['t', 'plant'])['Unloading_Sucrose'].aggregate(
                        np.sum)).groupby(['t'])
                    stand_data['Caphorn_Mixed'] = [mstruct, sucrose, amino_acids, nitrates, Unloading_Sucrose]
            else:
                mstruct = (organs_roots.groupby(['t'])['mstruct'].aggregate(np.sum)).groupby(['t'])
                sucrose = (organs_roots.groupby(['t'])['Conc_sucrose'].aggregate(np.sum)).groupby(['t'])
                amino_acids = (organs_roots.groupby(['t'])['Conc_amino_acids'].aggregate(np.sum)).groupby(['t'])
                nitrates = (organs_roots.groupby(['t'])['Conc_nitrates'].aggregate(np.sum)).groupby(['t'])
                Unloading_Sucrose = (organs_roots.groupby(['t'])['Unloading_Sucrose'].aggregate(np.sum)).groupby(['t'])
                stand_data[treatment] = [mstruct, sucrose, amino_acids, nitrates, Unloading_Sucrose]

            for stand, data in stand_data.iteritems():
                mstruct_data, sucrose_data, amino_acids_data, nitrates_data, Unloading_Sucrose_data = data[0], data[1], data[2], data[3], data[4]

                # mstruct
                mstruct_mean = mstruct_data.mean()
                ax1.plot(mstruct_mean.index, mstruct_mean, colors[stand], label=label_mapping[stand],
                         linestyle=lines[stand], linewidth=thickness)

                ax1.set_yticks(np.arange(0, 0.8, 0.2))

                # Sucrose
                sucrose_mean = sucrose_data.mean()
                ax2.plot(sucrose_mean.index, sucrose_mean, colors[stand], label=label_mapping[stand],
                         linestyle=lines[stand], linewidth=thickness)

                ax2.set_yticks(np.arange(0, 8000, 2000))

                # amino_acids
                amino_acids_mean = amino_acids_data.mean()
                ax3.plot(amino_acids_mean.index, amino_acids_mean, colors[stand], label=label_mapping[stand],
                         linestyle=lines[stand], linewidth=thickness)

                ax3.set_yticks(np.arange(0, 200, 50))

                # nitrates
                nitrates_mean = nitrates_data.mean()
                ax4.plot(nitrates_mean.index, nitrates_mean, colors[stand], label=label_mapping[stand],
                         linestyle=lines[stand], linewidth=thickness)

                ax4.set_yticks(np.arange(0, 2000, 500))

                # Unloading_Sucrose
                Unloading_Sucrose_mean = Unloading_Sucrose_data.mean()
                ax5.plot(Unloading_Sucrose_mean.index, Unloading_Sucrose_mean, colors[stand], label=label_mapping[stand],
                         linestyle=lines[stand], linewidth=thickness)

                ax5.set_yticks(np.arange(0, 125, 25))

            for ax in axes:
                ax.yaxis.set_label_coords(labelx, 0.5)
                [i.set_linewidth(thickness) for i in ax.spines.itervalues()]
                ax.tick_params(width=thickness, length=10, labelbottom='off', labelleft='off', top='off',
                               right='off', labelsize=fontsize)
                ax.set_xticks(np.arange(0, 1300, 200))
                ax.axvline(360, color='k', linestyle='--')

    ax_label_left = (plt.subplot(gs[0, 0]), plt.subplot(gs[1, 0]), plt.subplot(gs[2, 0]), plt.subplot(gs[3, 0]), plt.subplot(gs[4, 0]))
    for ax in ax_label_left:
        ax.tick_params(labelleft='on')
        if ax == plt.subplot(gs[0, 0]):
            ax.set_ylabel('Root structural \ndry mass (g)', fontsize=fontsize + 2)
        if ax == plt.subplot(gs[1, 0]):
            ax.set_ylabel(u'[sucrose] \n(µmol g$^{-1}$)', fontsize=fontsize + 2)
        if ax == plt.subplot(gs[2, 0]):
            ax.set_ylabel(u'[amino acids] \n(µmol g$^{-1}$)', fontsize=fontsize + 2)
        if ax == plt.subplot(gs[3, 0]):
            ax.set_ylabel(u'[nitrates] \n(µmol g$^{-1}$)', fontsize=fontsize + 2)
        if ax == plt.subplot(gs[4, 0]):
            ax.set_ylabel(u'Unloading sucrose \n(µmol g$^{-1}$ h$^{-1}$)', fontsize=fontsize + 2)
        if ax != plt.subplot(gs[0, 0]):
            ax.set_yticks(ax.get_yticks()[:-1])

    ax_label_bottom = (plt.subplot(gs[4, 0]), plt.subplot(gs[4, 1]), plt.subplot(gs[4, 2]), plt.subplot(gs[4, 3]))
    for ax in ax_label_bottom:
        ax.tick_params(labelbottom='on')
        ax.set_xlabel('Time from flowering (hour)', fontsize=fontsize)
        if ax != plt.subplot(gs[2, 3]):
            ax.set_xticks(ax.get_xticks()[:-1])

    # xtext = 300
    # ytext = 4.75
    # plt.subplot(gs[0, 0]).text(xtext, ytext, 'Density 200', fontsize=fontsize + 3, backgroundcolor='w',
    #                            verticalalignment='top')
    # plt.subplot(gs[0, 0]).set_zorder(1)
    # plt.subplot(gs[0, 1]).text(xtext, ytext, 'Density 410', fontsize=fontsize + 3, backgroundcolor='w',
    #                            verticalalignment='top')
    # plt.subplot(gs[0, 1]).set_zorder(1)
    # plt.subplot(gs[0, 2]).text(xtext, ytext, 'Density 600', fontsize=fontsize + 3, backgroundcolor='w',
    #                            verticalalignment='top')
    # plt.subplot(gs[0, 2]).set_zorder(1)
    # plt.subplot(gs[0, 3]).text(xtext, ytext, 'Density 800', fontsize=fontsize + 3, backgroundcolor='w',
    #                            verticalalignment='top')
    # plt.subplot(gs[0, 3]).set_zorder(1)
    #
    # xtext = 50
    # ytext = 575
    # plt.subplot(gs[0, 0]).text(xtext, 3.8, 'Grains', fontsize=fontsize + 2, verticalalignment='top')
    # plt.subplot(gs[0, 0]).set_zorder(1)
    # plt.subplot(gs[1, 0]).text(xtext, ytext, 'Stem', fontsize=fontsize + 2, backgroundcolor='w',
    #                            verticalalignment='top')
    # plt.subplot(gs[1, 0]).set_zorder(1)
    # plt.subplot(gs[2, 0]).text(xtext, ytext, 'Phloem', fontsize=fontsize + 2, backgroundcolor='w',
    #                            verticalalignment='top')
    # plt.subplot(gs[2, 0]).set_zorder(1)
    # plt.subplot(gs[3, 0]).text(xtext, ytext, 'Roots', fontsize=fontsize + 2, backgroundcolor='w',
    #                            verticalalignment='top')
    # plt.subplot(gs[3, 0]).set_zorder(1)

    ax1.legend(prop={'size': 11}, bbox_to_anchor=(0.5, 1.5), ncol=4, frameon=True)
    plt.savefig(os.path.join(GRAPHS_DIRPATH, 'Roots.TIFF'), dpi=300, format='TIFF',
                bbox_inches='tight')
    plt.close()

def figure_shoot_root():

    plt.figure(figsize=cm2inch(30, 10))
    gs = gridspec.GridSpec(1, 4, hspace=0, wspace=0)

    for density in DENSITIES:
        for treatment in TREATMENTS:
            if treatment not in ('Plano_diffus', 'Erect_diffus', 'Asso_diffus'): continue

            stand = treatment.split('_')[0]
            # Photosynthetic elements
            ph_elements_path = os.path.join(DATA_DIRPATH, 'Densite_{}'.format(density), treatment, 'outputs\elements_states.csv')
            ph_elements = pd.read_csv(ph_elements_path)
            ph_elements['sum_dry_mass'] = ((ph_elements['triosesP'] * 1E-6*C_MOLAR_MASS) / TRIOSESP_MOLAR_MASS_C_RATIO +
                                       (ph_elements['sucrose'] * 1E-6*C_MOLAR_MASS) / SUCROSE_MOLAR_MASS_C_RATIO +
                                       (ph_elements['starch'] * 1E-6*C_MOLAR_MASS) / HEXOSE_MOLAR_MASS_C_RATIO +
                                       (ph_elements['fructan'] * 1E-6*C_MOLAR_MASS) / HEXOSE_MOLAR_MASS_C_RATIO +
                                       (ph_elements['nitrates'] * 1E-6*N_MOLAR_MASS) / NITRATES_MOLAR_MASS_N_RATIO +
                                       (ph_elements['amino_acids'] * 1E-6*N_MOLAR_MASS) / AMINO_ACIDS_MOLAR_MASS_N_RATIO +
                                       (ph_elements['proteins'] * 1E-6*N_MOLAR_MASS) / AMINO_ACIDS_MOLAR_MASS_N_RATIO +
                                        ph_elements['mstruct'])
            # Organs
            organs_path = os.path.join(DATA_DIRPATH, 'Densite_{}'.format(density), treatment, 'outputs\organs_states.csv')
            organs = pd.read_csv(organs_path)
            organs['sum_dry_mass'] = (((organs.fillna(0)['structure'] + organs.fillna(0)['starch']) * 1E-6 * 12 / 0.4) +
                                     ((organs.fillna(0)['proteins'] * 1E-6 * 14) / 0.136) +
                                     (organs.fillna(0)['sucrose'] * 1E-6 * C_MOLAR_MASS) / SUCROSE_MOLAR_MASS_C_RATIO +
                                     (organs.fillna(0)['amino_acids'] * 1E-6 * N_MOLAR_MASS) / AMINO_ACIDS_MOLAR_MASS_N_RATIO +
                                     (organs.fillna(0)['nitrates'] * 1E-6 * N_MOLAR_MASS) / NITRATES_MOLAR_MASS_N_RATIO +
                                      organs.fillna(0)['mstruct'])

            ax1 = plt.subplot(gs[0, AXES_MAP[density]])

            stand_data = {}
            if stand == 'Asso':
                if density in ('200', '410', '600', '800'):
                    # Erectophile
                    dry_mass_shoot = ((ph_elements[ph_elements['plant'] == 1].groupby(['t'])['sum_dry_mass'].aggregate(np.sum)) +
                                     (organs[(organs['plant'] == 1) & (organs['organ'] == 'grains')].groupby(['t'])['sum_dry_mass'].aggregate(np.sum)))

                    dry_mass_root = (organs[(organs['plant'] == 1) & (organs['organ'] == 'roots')].groupby(['t'])['sum_dry_mass'].aggregate(np.sum))

                    shoot_root_ratio = dry_mass_shoot / dry_mass_root

                    stand_data['Caphorn_Mixed'] = shoot_root_ratio

                    # Planophile
                    dry_mass_shoot = ((ph_elements[ph_elements['plant'] == 2].groupby(['t'])[
                                           'sum_dry_mass'].aggregate(np.sum)) +
                                      (organs[(organs['plant'] == 2) & (organs['organ'] == 'grains')].groupby(
                                          ['t'])['sum_dry_mass'].aggregate(np.sum)))

                    dry_mass_root = (
                        organs[(organs['plant'] == 2) & (organs['organ'] == 'roots')].groupby(['t'])[
                            'sum_dry_mass'].aggregate(np.sum))

                    shoot_root_ratio = dry_mass_shoot / dry_mass_root

                    stand_data['Soissons_Mixed'] = shoot_root_ratio

                else:
                    print 'not implemented'
                    break
                    # # Erectophile
                    # dry_mass_shoot = ((ph_elements[(1 <= organs['plant']) & (organs['plant'] <= 5)].groupby(['t', 'plant'])[
                    #                        'sum_dry_mass'].aggregate(np.sum)) +
                    #                   (organs[(1 <= organs['plant']) & (organs['plant'] <= 5) & (organs['organ'] == 'grains')].groupby(
                    #                       ['t', 'plant'])['sum_dry_mass'].aggregate(np.sum)))
                    #
                    # dry_mass_root = (
                    #     organs[(1 <= organs['plant']) & (organs['plant'] <= 5) & (organs['organ'] == 'roots')].groupby(['t', 'plant'])[
                    #         'sum_dry_mass'].aggregate(np.sum))
                    #
                    # shoot_root_ratio = dry_mass_shoot / dry_mass_root
                    #
                    # stand_data['Caphorn_Mixed'] = [shoot_root_ratio]
                    #
                    # # Planophile
                    # dry_mass_shoot = ((ph_elements[(6 <= organs['plant']) & (organs['plant'] <= 10)].groupby(['t', 'plant'])[
                    #                        'sum_dry_mass'].aggregate(np.sum)) +
                    #                   (organs[(6 <= organs['plant']) & (organs['plant'] <= 10) & (organs['organ'] == 'grains')].groupby(
                    #                       ['t', 'plant'])['sum_dry_mass'].aggregate(np.sum)))
                    #
                    # dry_mass_root = (
                    #     organs[(6 <= organs['plant']) & (organs['plant'] <= 10) & (organs['organ'] == 'roots')].groupby(['t', 'plant'])[
                    #         'sum_dry_mass'].aggregate(np.sum))

                    shoot_root_ratio = dry_mass_shoot / dry_mass_root

                    stand_data['Soissons_Mixed'] = [shoot_root_ratio]
            else:
                dry_mass_shoot = ((ph_elements.groupby(['t'])[
                                       'sum_dry_mass'].aggregate(np.sum)) +
                                  (organs[organs['organ'] == 'grains'].groupby(
                                      ['t'])['sum_dry_mass'].aggregate(np.sum)))

                dry_mass_root = (organs[organs['organ'] == 'roots'].groupby(['t'])[
                                     'sum_dry_mass'].aggregate(np.sum))

                shoot_root_ratio = dry_mass_shoot / dry_mass_root
                stand_data[treatment] = shoot_root_ratio

            for stand, shoot_root_ratio in stand_data.iteritems():

                ax1.plot(shoot_root_ratio.index, shoot_root_ratio, colors[stand], label=label_mapping[stand], linestyle=lines[stand], linewidth=thickness)

            ax1.set_yticks(np.arange(0, 14, 2))
            ax1.set_xticks(np.arange(0, 1300, 200))
            ax1.yaxis.set_label_coords(labelx, 0.5)
            [i.set_linewidth(thickness) for i in ax1.spines.itervalues()]
            ax1.tick_params(width=thickness, length=10, top='off', right='off', labelsize=fontsize)
            ax1.axvline(360, color='k', linestyle='--')


    ax_label = (plt.subplot(gs[0, 0]), plt.subplot(gs[0, 1]), plt.subplot(gs[0, 2]), plt.subplot(gs[0, 3]))
    for ax in ax_label:
        ax.tick_params(labelbottom='on')
        ax.set_xlabel('Time from flowering (hour)', fontsize=fontsize)
        if ax != plt.subplot(gs[0, 3]):
            ax.set_xticks(ax.get_xticks()[:-1])

    for ax in ax_label:
        if ax == plt.subplot(gs[0, 0]):
            ax.tick_params(labelleft='on')
            ax.set_ylabel('Shoot : Root ratio ', fontsize=fontsize)
        else:
            ax.set_yticks(ax.get_yticks()[:-1])
            ax.tick_params(labelleft='off')

    # xtext = 300
    # ytext = 4.75
    # plt.subplot(gs[0, 0]).text(xtext, ytext, 'Density 200', fontsize=fontsize+3, backgroundcolor='w', verticalalignment='top')
    # plt.subplot(gs[0, 0]).set_zorder(1)
    # plt.subplot(gs[0, 1]).text(xtext, ytext, 'Density 410', fontsize=fontsize+3, backgroundcolor='w', verticalalignment='top')
    # plt.subplot(gs[0, 1]).set_zorder(1)
    # plt.subplot(gs[0, 2]).text(xtext, ytext, 'Density 600', fontsize=fontsize+3, backgroundcolor='w', verticalalignment='top')
    # plt.subplot(gs[0, 2]).set_zorder(1)
    # plt.subplot(gs[0, 3]).text(xtext, ytext, 'Density 800', fontsize=fontsize+3, backgroundcolor='w', verticalalignment='top')
    # plt.subplot(gs[0, 3]).set_zorder(1)
    #
    # xtext = 50
    # ytext = 575
    # plt.subplot(gs[0, 0]).text(xtext, 3.8, 'Grains', fontsize=fontsize+2, verticalalignment='top')
    # plt.subplot(gs[0, 0]).set_zorder(1)
    # plt.subplot(gs[1, 0]).text(xtext, ytext, 'Stem', fontsize=fontsize+2, backgroundcolor='w', verticalalignment='top')
    # plt.subplot(gs[1, 0]).set_zorder(1)
    # plt.subplot(gs[2, 0]).text(xtext, ytext, 'Phloem', fontsize=fontsize+2, backgroundcolor='w', verticalalignment='top')
    # plt.subplot(gs[2, 0]).set_zorder(1)
    # plt.subplot(gs[3, 0]).text(xtext, ytext, 'Roots', fontsize=fontsize+2, backgroundcolor='w', verticalalignment='top')
    # plt.subplot(gs[3, 0]).set_zorder(1)

    ax1.legend(prop={'size': 11}, bbox_to_anchor=(0.5, 1.25), ncol=4, frameon=True)
    plt.savefig(os.path.join(GRAPHS_DIRPATH, 'Shoot_root_ratio.TIFF'), dpi=300, format='TIFF', bbox_inches='tight')
    plt.close()


if __name__ == '__main__':

    # figure_inclination('03')
    # figure_PARa_dynamic_direct('04')
    # figure_PARa_dynamic_diffuse('05')
    # figure_profile_PARa_1000_plants('06_direct')
    # figure_PARa_Cass('07')
    # figure_green_area_proteins('08')
    # figure_C_non_struct('09')
    # figure_Cass_grain_mass('10')
    # table_PARa('01_10plants')
    # table_SumC('02')
    # grains_final_state()
    # table_total_N_uptake()
    # RUE_total()
    RUE_end_grains()
    # figure_roots()
    # figure_shoot_root()
