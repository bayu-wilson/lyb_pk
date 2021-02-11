#!/usr/bin/env python

import os
########################################################################################
######################## Control Center for entire pipeline ############################
########################################################################################
cat_name = "obs/XQ-100_catalogue" #"mocks/XQ-100_catalogue_n100"
tag = "lyb_nocorr_leftcentered" #"lyb_wNR_leftcentered" #badseeing_addp2"
mock_or_obs = cat_name.split("/")[0]

subtract_metal_power = True
continuum_correction = False
M = 1000 # Bootstrap samples, line 19 in boot_indo.py
remove_dla = True
wR2 = False # If True use original R2, if False, use wR column (11 km/s)
log_kbinning = True
carswell_res = 0 #1 == use the resolution from Carswell+18, 0== use my resolution
no_overlap = True  #currently only uses VIS arm in z=3.6 Lyman-alpha
SEEING_MAXIMUM = 1.5 #limit to reasonably good seeing

individual_qso_kmax = 0   #This is a vestigal option that looks like it is irrelevant (kmax is enforced eslewhere)

#flags for mocks
add_beta = False # only works on "noB" mocks, see line 188 in QuasarSpectrum.py
add_ovi = False # also only works on mocks
add_sithree = False # also only works on mocks

#bootstraps
BayusBootstrap = 0 #0=mine is much faster but otherwise identical

########################################################################################
#See line 86 in QuasarSpectrum.py for when `cat_name` gets used
#Common `cat_name`: "obs/XQ-100_catalogue", "mocks/XQ-100_catalogue_n100","obs/XQ-100_catalogue", "mocks/XQ-100_catalogue_n100", "mocks/XQ-100_catalogue_n5000"
#Most important use of `tag` in load_qso_data method in QuasarSpectrum.py. Line 113
#Common `tag`: "lyb_nocorr","corrNR","lyb_nocorr_n5000","gaussian_n600","noB_onlyA", "noB_n5000", "uncorr", "gaussian"


ntag = "n100" # default
nqso = 100
if "n5000" in tag:
    ntag = "n5000"
    nqso = 5000

#Choosing reshift and wavenumber binning
zmin = 3.0 #3.4
zmax = 4.2
zbinlen = 7 #5
if log_kbinning:
    kmin = 10**-2.5 #3e-3
    kmax = 10**-1.2 #6e-2
    kbinlen = 13
else: #VID BINS
    kmin = 3e-3 #0.0
    kmax = 0.06-3e-3#0.06 #0.06
    kbinlen = 19

# Rescaling
rescale_flux = 1.0

#####################
###  Saving Data  ###
#####################
save_mf_path = "../output/mf_{0}_{1}.csv".format(mock_or_obs,tag)
save_mf_with_err = True #get_errorbars.py
save_mf_with_err_path = "../output/mf_errboot_{0}_{1}.txt".format(mock_or_obs,tag)  #bootstrapped value

save_pk = True #main.py
save_pk_with_err = True #plot_pk.py
if mock_or_obs == 'obs':
    save_pk_path = "../output/pk_{0}_{1}_DLA{2}_metal{3}_res{4}.csv".format(mock_or_obs,tag, remove_dla, subtract_metal_power, carswell_res)
    save_pk_with_err_path = "../output/pk_errboot_{0}_{1}_DLA{2}_metal{3}_res{4}_nb{5}".format(mock_or_obs,tag, remove_dla, subtract_metal_power, carswell_res, M)
else:
    save_pk_path = "../output/pk_{0}_{1}.csv".format(mock_or_obs,tag)
    save_pk_with_err_path = "../output/pk_errboot_{0}_{1}.txt".format(mock_or_obs,tag)

save_boot_mf = True #bootstrap_pk.py
save_boot_mf_path = "../output/mf_boot_{0}_{1}_nb{2}_DLA{3}.csv".format(mock_or_obs,tag, M, remove_dla)

save_boot_pk = True #bootstrap_pk.py
save_boot_pk_path = "../output/pk_boot_{0}_{1}_DLA{2}_metal{3}_res{4}_nb{5}.csv".format(mock_or_obs,tag, remove_dla, subtract_metal_power, carswell_res, M)


#save in k,z,qidx bins July 18  (These are for the bootstrap samples)
save_kzq_mf_path = "../output/qsos_mf_{0}_{1}.txt".format(mock_or_obs, tag)
if mock_or_obs == 'obs':
    save_kzq_pk_path = "../output/qsos_pk_{0}_{1}_DLA{2}_metal{3}_res{4}.txt".format(mock_or_obs,tag, remove_dla, subtract_metal_power, carswell_res)
else:
    save_kzq_pk_path = "../output/qsos_pk_{0}_{1}.txt".format(mock_or_obs,tag)
#save masks bins July 31
save_masks_path = "../output/qsos_masks_{0}.txt".format(tag)


save_flux_pdf = False
save_flux_pdf_path = "../output/flux_pdf_lya.txt"

########################
###  Saving Figures  ###
########################
#save_mf_fig = True #plot_mf.py
#save_mf_fig_path = "figures/mf_{0}_{1}.pdf".format(mock_or_obs,tag)
#save_pk_fig = True  #plot_pk.py
#
#save_pk_fig_path = "figures/pk_{0}_{1}.pdf".format(mock_or_obs,tag)
#save_pk_fig_path = "figures/pk_{0}_{1}.pdf".format(mock_or_obs,tag)
#
#save_pk_ratio_fig = True #plot_ratio_pk.py
#save_pk_ratio_fig_path = "figures/pk_{0}_ratio_{1}.pdf".format(mock_or_obs,ntag)
#
#save_wR_ratio_fig = True #plot_ratio_wR.py
#save_wR_ratio_fig_path = "figures/pk_wR_ratio_{0}.pdf".format(ntag)
#
#save_nocorr_ratio_fig = True #plot_nocorr_ratio.py
#save_nocorr_ratio_fig_path = "figures/pk_metals_ratio_{0}.pdf".format(ntag)
#
#save_boot_fig = True #plot_bootstrap.py
#save_boot_fig_path = "figures/boot_{0}_{1}.pdf".format(mock_or_obs,tag)
#
#include_metal_model = False #plot_nocorr_ratio.py

#save_res_fig = True #plot_res.py
#save_res_fig_path = "figures/res_comparison.pdf"
#
#save_ratio_dla_fig = True #plot_ratio_dla.py
#save_dla_fig = True #plot_specific.py
#save_dla_fig_path = "figures/pk_REMOVING_DLAs_all.pdf"


# PAPER PLOTS AND STUFF
#save_paper_mf = True
#save_paper_mf_path = "figures/paper_mf.pdf"
#
#save_paper_mf_v2 = True
#save_paper_mf_v2_path = "figures/paper_mf_v2.pdf"
#
#save_paper_mf_v3 = True
#save_paper_mf_v3_path = "figures/paper_mf_v3.pdf"
#
#save_paper_err_pk = True
#save_paper_err_pk_path = "figures/paper_err_pk.pdf"
#
#save_paper_realization_err_pk = "True"
#save_paper_realization_err_pk_path = "figures/paper_realization_err_pk.pdf"

#save_paper_pk_v0 = True
#save_paper_pk_path_v0 = "figures/paper_pk_v0.pdf"

#save_paper_pk_v1 = True
#save_paper_pk_v1_path = "figures/paper_pk_v1.pdf"

#save_paper_pk_v2 = True
#save_paper_pk_v2_path = "figures/paper_pk_v2.pdf"

#save_paper_pk_v3 = True
#save_paper_pk_v3_path = "figures/paper_pk_v3.pdf"

#save_paper_pk_v4 = True
#save_paper_pk_v4_path = "figures/paper_pk_v4.pdf"
#
#save_paper_pk_v4_metal_power_ratio = True
#save_paper_pk_v4_metal_power_ratio_path = "figures/paper_metal_power_ratio.pdf"
#
#save_paper_pk_v5 = True
#save_paper_pk_v5_path = "figures/paper_pk_v5.pdf"
#
save_paper_covmatrix = True
save_paper_covmatrix_path = "figures/paper_covmatrix_{0}_{1}_DLA{2}_metal{3}_res{4}_nb{5}.pdf".format(mock_or_obs,tag, remove_dla, subtract_metal_power, carswell_res, M)
#
#save_paper_ratio = True
#save_paper_ratio_path = "figures/paper_ratio.pdf"
#
#save_paper_metal_power = True
#save_paper_metal_power_path = "figures/paper_metal_power.pdf"
#
#save_paper_tau = True
#save_paper_tau_path = "figures/paper_tau.pdf"
#
#save_paper_pixel_dist = True
#save_paper_pixel_dist_path = "figures/paper_pixel_distribution_{0}.pdf".format(tag)
#
#save_paper_OVI = True
#save_paper_OVI_path = "figures/paper_OVI_ratio.pdf"
#
#save_dla_ratio = True
#save_dla_ratio_path = "figures/paper_dla_ratio.pdf"
#
#save_paper_compare_vid = True
#save_paper_compare_vid_path = "figures/paper_compare_vid.pdf"
#
#save_paper_compare_pk_bv = True
#save_paper_compare_pk_bv_path = "figures/paper_compare_pk_bv.pdf"
#
#save_paper_compare_err_vid = True
#save_paper_compare_err_vid_path = "figures/paper_compare_err_vid.pdf"
#
#save_compare_mocks_bv = True
#save_compare_mocks_bv_path = "figures/compare_mocks_bv.pdf"
#
#save_composite_spectra = True
#save_composite_spectra_path = "figures/composite_spectra.pdf"
