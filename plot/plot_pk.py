#!/usr/bin/env python
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib
import numpy as np
import pandas as pd
import sys
sys.path.append('../pipeline')
import options as opt
import inis
#Note for Bayu -  This is the same as `paper_pk_v3.py`

##### Control Area ######
save_plot = True
show_plot = False
#########################
colors = ['red', 'green','blue','purple','orange','gold','indigo','black','gray']
marker = ['s','D','^','d','*']
SMALL_SIZE = 12
MEDIUM_SIZE = 15
BIGGER_SIZE = 20
plt.rc('font', size=BIGGER_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=BIGGER_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)
tick_length = 7.5
tick_width = 1.5
matplotlib.rcParams['xtick.major.size'] = tick_length
matplotlib.rcParams['xtick.major.width'] = tick_width
matplotlib.rcParams['xtick.minor.size'] = 5
matplotlib.rcParams['xtick.minor.width'] = 1.5
matplotlib.rcParams['ytick.major.size'] = tick_length
matplotlib.rcParams['ytick.major.width'] = tick_width
matplotlib.rcParams['errorbar.capsize'] = 4
####################################################################################################
labels = [r"$\widehat{P}_{\alpha \alpha}$",r"$\widehat{P}_{TT}$",
        r"$\widehat{P}_{\alpha \beta}$"]#r"$\widehat{\cal{P}}_{\alpha \beta}$"]
                #,, r"$P_{\alpha \beta}$",r"$P_{\beta \beta}$"]
custom_lines = [Line2D([0], [0], color=colors[0], lw=9, marker=None),
                Line2D([0], [0], color=colors[1], lw=9, marker=None),
                Line2D([0], [0], color=colors[2], lw=9, marker=None),
                Line2D([0], [0], lw=0)]
                # Line2D([0], [0], color=colors[3], lw=9, marker=None)]


#pkdata = pd.read_csv(inis.save_pk_with_err_path)
# pkdata = pd.read_csv("../output/continuum_correction/pk_errboot_obs_corrNR_continuum_uncorrected.txt")
# pkdata = pd.read_csv("../final_results/pk_errboot_obs_corrNR.txt")
# pkdata = pd.read_csv("../output/pk_errboot_obs_corrNR_leftcentered_DLATrue_metalTrue_res0_nb1000.csv")
# pkdata = pd.read_csv("../output/pk_errboot_obs_corrNR_DLATrue_metalTrue_res0_nb1000.csv")
pkdata = pd.read_csv(inis.save_pk_with_err_path)
print("Plotting data from {}".format(inis.save_pk_with_err_path))
#"../output/pk_errboot_obs_corrNR_leftcentered_DLATrue_metalTrue_res0_nb1000.csv")#inis.save_pk_with_err_path)
# pkdata_corrcont = pd.read_csv("../output/continuum_correction/pk_errboot_obs_corrNR_continuum_corrected.txt")

# plt.style.use('classic')
fig,ax = plt.subplots(2,3,sharex=True, sharey=True,gridspec_kw={'hspace': 0,'wspace': 0}) #ROW COLUMN

fig.add_subplot(111, frameon=False)
# hide tick and tick label of the big axes
plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
# fig.text(0.5, 0.04, r"$log_{10}(k/[km^{-1}s])$", ha='center', va='center',fontsize=20)
fig.text(0.02, 0.5, r"$log_{10}(k P(k) / \pi)$", ha='center', va='center', rotation='vertical',fontsize=25)
plt.grid(False)
# plt.xlabel(r"log$_{10}(k/[km^{-1}s])$")
# plt.ylabel(r"$log_{10}(k P(k) / \pi)$")

# fig.delaxes(ax[-1][2])
fig.set_size_inches(15,9)
ax[-1,1].set_xlabel(r"$log_{10}(k/$[km$^{-1}$s])",fontsize=25)
#[i.set_xlabel(r"log$_{10}(k/[km^{-1}s])$") for i in ax[-1]] #xlabels on bottom row
# ax[1,1].set_ylabel(r"$k P(k) / \pi$")

zidx = 2
for i in range(2):
    for j in range(3):
        if zidx<7:
            pkmask = (pkdata.z == opt.zbin_centers[zidx])
            t = pkdata[pkmask]
            k,paa,ptt,pab,pbb = t.k.values,t.paa.values,t.ptt.values,t.pab.values,t.pbb.values
            err_paa,err_pab,err_ptt,err_pbb= t.err_paa.values, t.err_pab.values, t.err_ptt.values,t.err_pbb.values

            ### INCLUDING RESOLUTION UNCERTAINTY
            sigma_res_aa = opt.find_res_uncertainty(k,t.z.values,paa)
            err_paa = np.sqrt(err_paa**2+sigma_res_aa**2)
            sigma_res_tt = opt.find_res_uncertainty(k,t.z.values,ptt)
            err_ptt = np.sqrt(err_ptt**2+sigma_res_tt**2)
            sigma_res_ab = opt.find_res_uncertainty(k,t.z.values,pab)
            err_pab = np.sqrt(err_pab**2+sigma_res_ab**2)
            sigma_res_bb = opt.find_res_uncertainty(k,t.z.values,pbb)
            err_pbb = np.sqrt(err_pbb**2+sigma_res_bb**2)


            log_err_paa = [np.abs(np.log10((paa-err_paa))-np.log10(paa)),
                           np.abs(np.log10((paa+err_paa))-np.log10(paa))] #0.434*err_paa/paa
            log_err_ptt = [np.abs(np.log10((ptt-err_ptt))-np.log10(ptt)),
                           np.abs(np.log10((ptt+err_ptt))-np.log10(ptt))] #0.434*err_ptt/ptt
            log_err_pab = [np.abs(np.log10(np.abs(pab-err_pab))-np.log10(pab)),
                           np.abs(np.log10((pab+err_pab))-np.log10(pab))] #0.434*err_pab/pab
            k_x = np.log10(k)
            # for test_idx in range(len(log_err_pab)):
            #     #if np.log10(k[test_idx]*pab[test_idx]/np.pi)-(log_err_pab[0])[test_idx]>0:
            #     if pab[test_idx]-err_pab[test_idx]<0:
            #         (log_err_pab[0])[test_idx] = 0
            #         print("!")

            if (i == 0)&(j==2):
                pass
            else:
                #ax[i,j].set_yscale('log')
                ax[i,j].errorbar(k_x,np.log10(k*paa/np.pi),yerr=log_err_paa,color=colors[0], fmt='.')
                ax[i,j].errorbar(k_x*0.99,np.log10(k*ptt/np.pi),yerr=log_err_ptt,color=colors[1], fmt='.')
                ax[i,j].errorbar(k_x,np.log10(k*pab/np.pi),yerr=log_err_pab,color=colors[2], fmt='.')
                #ax[i,j].errorbar(k_x,np.log10(k*pab/np.pi))
                #ax[i,j].errorbar(k_x,k*pbb/np.pi,yerr=err_pbb*k/np.pi,color=colors[3], fmt='.')
                ax[i,j].text(0.82, 0.10,"z={0}".format(opt.zbin_centers[zidx]) #0.80, 0.95
                                                          , ha='center', va='center',
                                                          transform=ax[i,j].transAxes,fontsize=25)
                ax[i,j].xaxis.set_ticks_position('both')
                ax[i,j].yaxis.set_ticks_position('both')
                ax[i,j].xaxis.set_tick_params(direction='in')#, which='top')
                ax[i,j].yaxis.set_tick_params(direction='in')#, which='top')
                #ax[i,j].set_yticks(ax[i,j].get_yticks()[::1])
                #print(zidx,i,j)
                zidx += 1

                ################################################################
                ################################################################
                ################################################################
                # #CONTINUUM CORRECTION
                # t2 = pkdata_corrcont[pkmask]
                # k,paa,ptt,pab,pbb = t2.k.values,t2.paa.values,t2.ptt.values,t2.pab.values,t2.pbb.values
                # err_paa,err_pab,err_ptt,err_pbb= t2.err_paa.values, t2.err_pab.values, t2.err_ptt.values,t2.err_pbb.values

                # ### INCLUDING RESOLUTION UNCERTAINTY
                # sigma_res_aa = opt.find_res_uncertainty(k,t2.z.values,paa)
                # err_paa = np.sqrt(err_paa**2+sigma_res_aa**2)
                # sigma_res_tt = opt.find_res_uncertainty(k,t2.z.values,ptt)
                # err_ptt = np.sqrt(err_ptt**2+sigma_res_tt**2)
                # sigma_res_ab = opt.find_res_uncertainty(k,t2.z.values,pab)
                # err_pab = np.sqrt(err_pab**2+sigma_res_ab**2)
                # sigma_res_bb = opt.find_res_uncertainty(k,t2.z.values,pbb)
                # err_pbb = np.sqrt(err_pbb**2+sigma_res_bb**2)
                #
                #
                # log_err_paa = [np.abs(np.log10((paa-err_paa))-np.log10(paa)),
                #                np.abs(np.log10((paa+err_paa))-np.log10(paa))] #0.434*err_paa/paa
                # log_err_ptt = [np.abs(np.log10((ptt-err_ptt))-np.log10(ptt)),
                #                np.abs(np.log10((ptt+err_ptt))-np.log10(ptt))] #0.434*err_ptt/ptt
                # log_err_pab = [np.abs(np.log10(np.abs(pab-err_pab))-np.log10(pab)),
                #                np.abs(np.log10((pab+err_pab))-np.log10(pab))] #0.434*err_pab/pab
                # k_x = np.log10(k)*1.01
                #
                # ax[i,j].errorbar(k_x,np.log10(k*paa/np.pi),yerr=log_err_paa,color=colors[0], fmt='.',alpha=0.2)
                # ax[i,j].errorbar(k_x*0.99,np.log10(k*ptt/np.pi),yerr=log_err_ptt,color=colors[1], fmt='.',alpha=0.2)
                # ax[i,j].errorbar(k_x,np.log10(k*pab/np.pi),yerr=log_err_pab,color=colors[2], fmt='.',alpha=0.2)
                ################################################################
                ################################################################
                ################################################################

# ax[1,2].set_ylim(np.log10(8e-3),np.log10(4e-1))
ax[0,0].set_xlim(-2.7,-0.9)
ax[0,0].set_ylim(-2.24,-0.35)
fig.delaxes(ax[0][2])

# if plot_inis:
#     #labels.append("Wilson+19")
#     custom_lines.append(Line2D([0], [0], color='k',lw=0,marker='.'))
# box = ax[-1,1].get_position()
# ax[-1,1].set_position([box.x0, box.y0, box.width, box.height])
# ax[-1,1].legend(custom_lines,labels,fontsize=25,loc='center left', bbox_to_anchor=(1.1, 0.5),frameon=False)
ax[1,0].legend(custom_lines[:2],labels[:2],loc='lower left',ncol=1,frameon=False)
ax[1,1].legend(custom_lines[2:4],labels[2:4],loc='lower left',ncol=1,frameon=False)


# for i in ax[0,2].get_xticklabels():
#     i.set_visible(True)
ax[0,2].xaxis.set_tick_params(labelbottom=True)

# ax[0,1].set_xticklabels(['-2.5','-2.0','-1.5','-1.0'])

plt.tight_layout()
if save_plot:
    plt.savefig("figures/plot_pk.pdf")
if show_plot:
    plt.show()
else:
    plt.clf()
