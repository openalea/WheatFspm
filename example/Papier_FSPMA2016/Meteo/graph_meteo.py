# -*- coding: latin-1 -*-

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

##data_csv = r'D:\Documents\PostDoc_Grignon\Modeles\FSPM_Wheat\example\Papier_FSPMA2016\Meteo\data2graph.csv'
##
##data_df = pd.read_csv(data_csv)
##
##
### Daily
##ax1 = plt.subplot(211)
##ax2 = plt.subplot(212)
##
##for year, group in data_df.groupby('aaaa'):
##    cumul_PAR_daily = group.groupby('jjmm')['PAR_corrige'].sum() * (3600/2) / 1E6
##    mean_temp = group.groupby('jjmm')['PAR_corrige'].sum() * (3600/2) / 1E6
##    mean_Rd_PAR_ratio = group.groupby('jjmm')['Rd:PAR'].mean().tolist()
##    t = np.arange(1, len(cumul_PAR_daily)+ 1)
##    ax1.plot(t, cumul_PAR_daily, label=year)
##    ax2.plot(t, mean_Rd_PAR_ratio)
##
##
##ax1.set_ylabel("Daily cumulated PAR  '\n' (mol m-2)")
##ax1.tick_params(labelright='off', labelbottom='off')
##ax1.set_yticks(np.arange(0, 80, 20))
##ax2.set_xlabel('Day (year 2009, Lusignan)')
##ax2.set_ylabel("Mean daily '\n' diffuse:total PAR ratio")
##ax2.set_yticks(np.arange(0, 1, 0.2))
##
##plt.subplots_adjust(hspace = .001)
##plt.savefig('Daily.TIFF', dpi=300, format='TIFF', bbox_inches='tight')
##plt.close()
##
### Monthly
##ax1 = plt.subplot(211)
##ax2 = plt.subplot(212)
##
##for year, group in data_df.groupby('aaaa'):
##    group['PAR_corrige']
##    cumul_PAR_monthly = group.groupby('mm')['PAR_corrige'].sum() * (3600/2) / 1E6
##    mean_Rd_PAR_ratio = group.groupby('mm')['Rd:PAR'].mean().tolist()
##    SD_Rd_PAR_ratio = group.groupby('mm')['Rd:PAR'].std().tolist()
##    t = np.arange(1, len(cumul_PAR_monthly)+ 1)
##    ax1.plot(t, cumul_PAR_monthly.tolist(), label=year, marker='o')
##    ax2.errorbar(t, mean_Rd_PAR_ratio, yerr= SD_Rd_PAR_ratio, marker='o')
##
##
##plt.xticks(t, ['Mai', 'Juin', 'Juillet'])
##ax1.tick_params(labelright='off', labelbottom='off')
##ax1.set_ylabel('Cumulated PAR (mol m-2)')
##ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
##
##ax2.set_ylabel('Mean diffuse:total PAR ratio')
##
##plt.subplots_adjust(hspace = .001)
##plt.savefig('Monthly.TIFF', dpi=300, format='TIFF', bbox_inches='tight')
##plt.close()
##
### Boxplots
##ax = plt.subplot(111)
##
##mean_Rd_PAR_ratio_list = []
##for year, group in data_df.groupby('aaaa'):
##    mean_Rd_PAR_ratio_list.append(group['Rd:PAR'].dropna())
##
##bp = ax.boxplot(mean_Rd_PAR_ratio_list, showmeans=True, showfliers=True)
##
##ax.set_xticklabels(data_df['aaaa'].unique())
##ax.set_ylabel('Diffuse:total PAR ratio')
##ax.set_xlabel('Year')
##
##plt.savefig('Boxplots.TIFF', dpi=300, format='TIFF', bbox_inches='tight')
##plt.close()
##

#2009

def cm2inch(*tupl):
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)

thickness = 1.25
fontsize = 14
data_2009 = r'D:\Documents\PostDoc_Grignon\Modeles\FSPM_Wheat\example\Papier_FSPMA2016\Meteo\Meteo_Lusignan_2009.csv'

data_df = pd.read_csv(data_2009)

fig = plt.figure(figsize=cm2inch(15,30))
gs = gridspec.GridSpec(2, 1, hspace=0, wspace=0)

ax1 = plt.subplot(gs[0,0])
ax1b = ax1.twinx()
ax2 = plt.subplot(gs[1,0])

cumul_PAR_daily = data_df.groupby('jjmm')['PAR_corrige'].sum() * (3600/2) / 1E6
mean_temp = data_df.groupby('jjmm')['Temp air'].mean().tolist()
mean_Rd_PAR_ratio = data_df.groupby('jjmm')['Rd:PAR'].mean().tolist()
t = np.arange(1, len(cumul_PAR_daily)+ 1)
lns1=ax1.plot(t, cumul_PAR_daily, linewidth = thickness, label = 'PAR')
lns2=ax1b.plot(t, mean_temp, color='r', linewidth = thickness, label = 'Temperature')
ax2.plot(t, mean_Rd_PAR_ratio, color='k', linewidth = thickness)


ax1.set_ylabel("Cumulated PAR (mol m$^{-2}$ day$^{-1}$)", multialignment='center', fontsize = fontsize)
ax1.set_yticks(np.arange(0, 80, 20))
ax1.tick_params(labelright='off', labelbottom='off')
ax1b.set_yticks(np.arange(0, 40, 10))
ax1b.set_ylabel(u"Daily mean air temperature (°C)", fontsize = fontsize)
ax2.set_xlabel('Time from flowering (day)', fontsize = fontsize)
ax2.set_ylabel("Mean daily diffuse:total PAR ratio", multialignment='center', fontsize = fontsize)
ax2.set_yticks(np.arange(0, 1, 0.2))
ax1.text(0.05, 0.95, 'A', transform=ax1.transAxes, fontsize=fontsize, verticalalignment='top')
ax2.text(0.05, 0.95, 'B', transform=ax2.transAxes, fontsize=fontsize, verticalalignment='top')

lns = lns1+lns2
labs = [l.get_label() for l in lns]
plt.legend(lns, labs, bbox_to_anchor=(1, 2), frameon=False)


for ax in (ax1, ax1b, ax2):
    [i.set_linewidth(thickness) for i in ax.spines.itervalues()]
    ax.tick_params(width = thickness, labelsize=fontsize)

plt.savefig(r'D:\Documents\PostDoc_Grignon\Modeles\FSPM_Wheat\example\Papier_FSPMA2016\Figures\Figure_S01_Meteo.TIFF', dpi=300, format='TIFF', bbox_inches='tight')
plt.close()