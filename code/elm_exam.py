#!/bin/env python3

import h5py
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as plticker

data_filename = '../data/KSTAR_detach_ctrl_data.h5'
full_filename = '../data/KSTAR_detach_ctrl_data_full_res.h5'

device = 'KSTAR'
shot = 36161
sn = f'{shot}'

fig, axs = plt.subplots(3, 5, sharex='col', sharey='row', gridspec_kw={'width_ratios':[2, 0.5, 0.5, 0.5, 0.5], 'wspace':0.3, 'hspace':0.1, 'right':0.98, 'left':0.1, 'top':0.98}, figsize=(7.5, 6))

tr1 = (7.1, 7.2)
tr2 = (10.1, 10.2)
tr4 = (13.8, 13.9)
tr3 = (12.0, 12.1)
axda = axs[0, :]
axafrac = axs[1, :]
axgas = axs[2, :]
#axnitro = axs[3, :]
axda[0].set_ylim(0, 8)
axs[0, 0].set_xlim(7,16)
colors = ['r', 'm', 'b', 'g']
axda[0].set_ylabel(r'$D_\alpha$ / A.U.')
axafrac[0].set_ylabel('$A_{frac}$')
for ax in axs[-1, 1:]:
    ax.tick_params(axis='x', labelrotation=90)
    loc = plticker.MultipleLocator(base=0.1) # this locator puts ticks at regular intervals
    ax.xaxis.set_major_locator(loc)

afrac_corr_orig = 4.67955
afrac_corr_new = 5.5
afrac_adjustment_factor = afrac_corr_orig / afrac_corr_new

axafrac[0].set_ylim(0, 1.2)
axgas[0].set_ylabel('$PVBLN_2$ \n Cmd / V')
for i, tr in enumerate([tr1, tr2, tr3, tr4]):
    axs[0, 1+i].set_xlim(*tr)
    for ax in axs[:, 0]:
        ax.axvline(tr[0], ls='--', color=colors[i])
        ax.axvline(tr[1], ls='--', color=colors[i])
    for ax in axs[:, 1+i]:
        plt.setp(ax.spines.values(), color=colors[i])
        plt.setp([ax.get_xticklines(), ax.get_yticklines()], color=colors[i])

for ax in axs[1, :]:
    ax.axhline(1, color='k', ls=':')

with h5py.File(data_filename, 'r') as h5, h5py.File(full_filename, 'r') as h5full:
    dat = h5full[sn]['KSTAR']['POL_HA03']['dim0'][:]
    da = h5full[sn]['KSTAR']['POL_HA03']['data'][:]
    for ax in axda:
        ax.plot(dat, da/1e18, color='k')
    afrac = h5[sn]['PCS_KSTAR']['DVSAFRAC']['data'][:] * afrac_adjustment_factor
    afract = h5[sn]['PCS_KSTAR']['DVSAFRAC']['dim0'][:]
    afrac[afract < 6] = np.nan
    for ax in axafrac:
        ax.plot(afract, afrac, color='k')
        
    gas_enable = h5[sn]['PCS_KSTAR']['GVIPVBLN2ON']['data'][:]
    gas_max = h5[sn]['PCS_KSTAR']['GVTPVBLN2MAX']['data'][:]
    gas_min = h5[sn]['PCS_KSTAR']['GVTPVBLN2MIN']['data'][:]
    gas_raw = h5[sn]['PCS_KSTAR']['GVSPVBLN2']['data'][:]
    gast = h5[sn]['PCS_KSTAR']['GVSPVBLN2']['dim0'][:]
    gas = gas_raw * gas_enable
    gas[gas > gas_max] = gas_max[gas > gas_max]
    gas[gas < gas_min] = gas_min[gas < gas_min]
    gas *= gas_enable
    for ax in axgas:
        ax.plot(gast, gas, color='k')

    # nitrogen_vii_spec = h5[sn]['SPECTRO']['VUV_IM_01']['data'][:]
    # nitrogen_vii_spect = h5[sn]['SPECTRO']['VUV_IM_01']['dim0'][:]
    # for ax in axnitro:
    #     ax.plot(nitrogen_vii_spect, nitrogen_vii_spec, color='k')

axs[0, 0].text(0.5, 0.02, 'Time / s', transform=fig.transFigure, ha='center', va='bottom')
axs[0, 0].text(0.99, 0.01, f'{device}#{shot}', transform=fig.transFigure, ha='right', va='bottom', fontsize=8)

fig.savefig(f"../figures/elm_exam_{shot}.pdf")

plt.show()
