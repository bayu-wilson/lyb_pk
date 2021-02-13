#!/usr/bin/env python
import pandas as pd
import numpy as np
import sys
sys.path.append('../pipeline')
import options as opt
import inis
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import os

save_plot = True
show_plot = False

#Load in the measurements
mf_df = pd.read_csv(inis.save_mf_path)
mf_boot = np.loadtxt(inis.save_boot_mf_path)
errorbars = np.reshape(np.cov(mf_boot).diagonal()**0.5,(2,opt.zbinlen))

#Calculate the continuum correction from Faucher-Giguére+2008
C_a_ratio = opt.continuum_correction(opt.zbin_centers)
C_t_ratio = opt.continuum_correction(opt.zbin_centers-0.2)
C_b_ratio = np.abs(np.abs(C_a_ratio)-np.abs(C_t_ratio))

#Continuum correction is done in the pipeline so here we UNCORRECT the mean flux
mfa = mf_df.mf_a #mean lyman alpha flux #continuum corrected
mfa_uncorr = mfa / (1-C_a_ratio)
nans_mft = np.ones_like(opt.zbin_centers)
nans_mft[:2] = np.nan
mft = mf_df.mf_tot*nans_mft #mean total flux in lyman beta forest (includes both lyb and lya)
mft_uncorr = mft / (1-C_t_ratio)*nans_mft #one data point is bad
mfb = mf_df.mf_b # mean lyman beta flux
mfb_uncorr = mfb / (1-C_b_ratio)

#Errorbars
err_mfa = errorbars[0]#mf_df.err_mfa
err_mft = errorbars[1]#mf_df.err_mft
err_mfb = opt.find_err_mf_beta(mfb,mfa,err_mfa,mft,err_mft)

z_mf = mf_df.z

#preparing for the figure
colors = ['red','blue','green','black','gray'] # 3 colors for mfa,mft,mfb
labels = [r"$\overline{F}_{\beta}$",r"$\overline{F}_{T}$",r"$\overline{F}_{\alpha}$",
            r"Model: $\gamma$=1.5", r"Model: $\gamma$=1.0",
            "This work", "Iršič+2017","Becker+2013",
            "Continuum Uncorrected", "Continuum Corrected"]
linestyles = ['-','dashed','-.','dotted']
custom_lines = [Line2D([0], [0], color=colors[2], lw=9),
                Line2D([0], [0], color=colors[1], lw=9),
                Line2D([0], [0], color=colors[0], lw=9),
                Line2D([0], [0], color=colors[3], lw=9),
                Line2D([0], [0], color=colors[4], lw=9),
                [],
                [],
                [],
                Line2D([0], [0], color='k',ls=linestyles[1]),
                Line2D([0], [0], color='k',ls=linestyles[0])]

#plotting routine
fig,ax = plt.subplots(1)
fig.set_size_inches(13,13)
ax.set_ylabel("F(z)",fontsize = 30)
ax.set_xlabel("z",fontsize = 30)
ax.set_ylim(0.15,1.05)#0.31,0.85)
ax.tick_params(axis='both', which='major', labelsize=22.5)
# fontsize = 18


capsize,markersize,lw,capthick,elinewidth = 6,0,1.5,3,4
ax.errorbar(z_mf,mfa,yerr=err_mfa, #F_ALPHA MEAN FLUX, CORRECTED
            color=colors[0],capsize=capsize,markersize=markersize,lw=lw,capthick=capthick,elinewidth=elinewidth)
ax.errorbar(z_mf*1.001,mfa_uncorr,yerr=err_mfa, #F_ALPHA MEAN FLUX, UNCORRECTED
            color=colors[0],capsize=capsize,markersize=markersize,lw=lw,capthick=capthick,ls=linestyles[1],elinewidth=elinewidth)
ax.errorbar(z_mf,mft,yerr=err_mft,color=colors[1], #F_TOTAL MEAN FLUX, CORRECTED
            capsize=capsize,markersize=markersize,lw=lw,capthick=capthick,elinewidth=elinewidth)
ax.errorbar(z_mf*1.000,mft_uncorr,yerr=err_mft,color=colors[1], #F_TOTAL MEAN FLUX, UNCORRECTED
            capsize=capsize,markersize=markersize,lw=lw,capthick=capthick,ls=linestyles[1],elinewidth=elinewidth)
ax.errorbar(z_mf,mfb,yerr=err_mfb, #F_BETA MEAN FLUX, CORRECTED
            color=colors[2],capsize=capsize,markersize=markersize,lw=lw,capthick=capthick,elinewidth=elinewidth)
ax.errorbar(z_mf*0.998,mfb_uncorr,yerr=err_mfb, #F_BETA MEAN FLUX, UNCORRECTED
            color=colors[2],capsize=capsize,markersize=markersize,lw=lw,capthick=capthick,ls=linestyles[1],elinewidth=elinewidth)

#adding error bar to legend
custom_lines[5] = ax.errorbar([],[],yerr=1,color='k',capsize=7,markersize=0,lw=1.5,capthick=3,elinewidth=4)
custom_lines[6] = ax.errorbar([],[],yerr=1,color='k',fmt='^',capsize=4,markersize=7,lw=1,capthick=1)
custom_lines[7] = ax.errorbar([],[],color='k',fmt='o',markersize=7,lw=1,capthick=1)

ax.legend(custom_lines,labels,fontsize=25,loc='lower left',ncol=2,frameon=False)# bbox_to_anchor=(1, 0.5))

figname = "plot_meanflux.pdf"
if save_plot:
    plt.savefig(figname,bbox_inches='tight')
print(figname)
if show_plot:
    plt.show()
plt.clf()
