#!/usr/bin/env python

from QuasarSpectrum import QuasarSpectrum
import inis
import numpy as np
import pandas as pd
import options as opt
from scipy.interpolate import griddata

cat_name = inis.cat_name
tag = inis.tag
rescale_flux = inis.rescale_flux #typically I am not rescaling the flux so this is one. If I am, it would be on mocks.
print("Catalogue:", cat_name)
print("Tag:", tag, "\n")

# In line ~109 and ~335 there is 'RuntimeWarning: invalid value encountered in true_divide'
np.seterr(divide='ignore', invalid='ignore') # I am ignoring this error because it is not important

##############################
######## Loading data ########
##############################
QuasarSpectrum.load_cat(cat_name) # loading catalog using QuasarSpectrum class method
nqso = QuasarSpectrum.nqso
print("Loading Data")
qso_arr = []
for i in range(nqso): #looping through each quasar filepath to load it into an array of nqso objects
    q = QuasarSpectrum.load_qso_data(i,tag=tag,rescale_flux=rescale_flux)
    #####################################
    # Adding more lines to mock spectra #
    #####################################
    #Lyb
    if ("noB" in tag)&(inis.add_beta): # adding Lyb to mocks
        q.get_new_forest(rescale_flux=rescale_flux,wrange = (opt.lyb_min,opt.lyb_max,opt.lyb_rest, opt.xs_beta))
    #OVI
    if inis.cat_name.startswith('mocks')&(inis.add_ovi):
        q.get_new_forest(rescale_flux=rescale_flux,wrange = (opt.ovi_min,opt.ovi_max,opt.ovi_rest_d1, opt.xs_ovi)) #ovi_rest_d1 = 1032
        q.get_new_forest(rescale_flux=rescale_flux,wrange = (opt.ovi_min,opt.ovi_max,opt.ovi_rest_d2, opt.xs_ovi*opt.ovi_factor)) #ovi_rest_d2 = 1038
    #SiIII
    if inis.cat_name.startswith('mocks')&(inis.add_sithree):
        q.get_new_forest(rescale_flux=rescale_flux,wrange = (opt.sithree_min,opt.sithree_max,opt.sithree_rest_d1, opt.xs_sithree))
    qso_arr.append(q)
print("Done!\n")

###################################
###### Mean flux calculation ######
###################################
mf_msrmnts = opt.mf_msrmnts #list of columns names for the mean flux data table
n_mf_msrmnts = len(mf_msrmnts)
zbin_msr_matrix = np.zeros((opt.zbinlen,n_mf_msrmnts)) # 7 zbins by 12 measurements
lya_flux_pdf = [[] for tmp in range(opt.zbinlen)]
pdf_bins = opt.pdf_bins
print("Mean Flux Calculation")
f = open(inis.save_kzq_mf_path,'w') #july10 #opening a file that I'll write to... relating to mean flux measurements. (This makes is easier to do bootstraps.)
for zidx in range(opt.zbinlen):
    #looping through redshift bins
    msrmnt_in_zbin =  np.zeros(n_mf_msrmnts) #measurements in the z'th bin
    count_in_zbin = np.zeros(n_mf_msrmnts)   #counts in the zidx'th bin
    opt.updt(opt.zbinlen, zidx)
    for i in range(nqso):
        #looping through each quasar
        name = qso_arr[i].name
        zpix_a = qso_arr[i].get_zpix(opt.lya_rest) #redshifts of LYA forest for the i'th quasar
        mask = qso_arr[i].get_zmask(forest=(opt.lya_min,opt.lya_max,opt.lya_rest),zpix=zpix_a,zidx=zidx,zedges=opt.zbin_edges,name=name) #masks forest and zbin
        zpix_b = qso_arr[i].get_zpix(opt.lyb_rest) #redshifts of LYB forest for the i'th quasar
        mask_b = qso_arr[i].get_zmask(forest=(opt.lyb_min,opt.lyb_max,opt.lyb_rest),
                                     zpix=zpix_b,zidx=zidx,zedges=opt.zbin_edges,name=name)
        za = zpix_a[mask] #redshifts of LYA forest for the i'th quasar AND in zidx'th bin
        ztot = zpix_b[mask_b] #redshifts of LYB forest for the i'th quasar AND in zidx'th bin
        try:
            new_af_mask = (za>np.min(ztot))&(za<np.max(ztot)) #where the LYA forest z-pixels are within the upper and lower redshift bounds of the LYB forest
            new_bf_mask = (ztot>np.min(za))&(ztot<np.max(za)) #where the LYB forest z-pixels are within the upper and lower redshift bounds of the LYA forest
            ferra = qso_arr[i].err_flux[mask][new_af_mask]
            ferrtot = qso_arr[i].err_flux[mask_b][new_bf_mask]

            #Interpolating to the forest with the least number of pixels
            if len(ferrtot)<=len(ferra):
                ferra = griddata(za[new_af_mask],ferra,ztot[new_bf_mask],method='linear')
                ferra = ferra[np.isfinite(ferra)] #if you don't do this, something goes wrong... I think
            else:
                ferrtot = griddata(ztot[new_bf_mask],ferrtot,za[new_af_mask],method='linear')
                ferrtot = ferrtot[np.isfinite(ferrtot)]
            msrmnt_in_zbin[10]+= np.sum(ferra*ferrtot)*0 #CHANGED 5/3/19 after meeting with Matt. Flux covariance between alpha and total should be zero -- That's correct
            count_in_zbin[10] += len(ferra)                         # len lya_tot 10
        except:
            #sometimes the two forests don't overlap
            pass

        msrmnt_in_zbin[0]+= np.sum(qso_arr[i].flux[mask])          # mf lya 0
        count_in_zbin[0] += np.sum(mask)                           # len lya 0
        msrmnt_in_zbin[1]+= np.sum(qso_arr[i].err_flux[mask]**2)   # var lya 1
        count_in_zbin[1] += np.sum(mask)                           # len lya 1
        msrmnt_in_zbin[2]+= np.sum(qso_arr[i].dloglambda[mask])    # dloglam lya 2
        count_in_zbin[2] += np.sum(mask)                           # len lya 2
        msrmnt_in_zbin[4]+= np.sum(qso_arr[i].flux[mask_b])        # mf tot 4
        count_in_zbin[4] += np.sum(mask_b)                         # len tot 4
        msrmnt_in_zbin[5]+= np.sum(qso_arr[i].err_flux[mask_b]**2) # var tot 5
        count_in_zbin[5] += np.sum(mask_b)                         # len tot 5
        msrmnt_in_zbin[6]+= np.sum(qso_arr[i].dloglambda[mask_b])  # dloglam tot 6
        count_in_zbin[6] += np.sum(mask_b)                         # len tot 6

        #writing results per quasar to an external file so we can calculate bootstraps faster
        s = "{0:g} {1:g} {2:g} {3:g} {4:g} {5:g}".format(i,opt.zbin_centers[zidx],np.sum(qso_arr[i].flux[mask]), np.sum(mask),np.sum(qso_arr[i].flux[mask_b]),np.sum(mask_b))
        f.write(s+'\n') #july10

        lya_flux_pdf[zidx].append(np.histogram(qso_arr[i].flux[mask],bins=pdf_bins)[0]/np.nansum(mask))#FLUX PDF
    zbin_msr_matrix[zidx] = msrmnt_in_zbin/count_in_zbin #since we are want the MEAN of these measurements per zbin, I SUMMED them up and then DIVIDED them by the counts
opt.updt(opt.zbinlen, opt.zbinlen)
f.close()
print("Done!\n")

zbin_msr_matrix.T[3] = list(QuasarSpectrum.get_npow(mf=zbin_msr_matrix.T[0], nvar=zbin_msr_matrix.T[1], dloglambda=zbin_msr_matrix.T[2]))
zbin_msr_matrix.T[7] = list(QuasarSpectrum.get_npow(mf=zbin_msr_matrix.T[4], nvar=zbin_msr_matrix.T[5],dloglambda=zbin_msr_matrix.T[6]))
zbin_msr_matrix.T[8] = opt.zbin_centers # zbins 8
zbin_msr_matrix.T[11] = ((zbin_msr_matrix.T[4]*zbin_msr_matrix.T[0])**(-1)*zbin_msr_matrix.T[10] * np.pi / (opt.kmax-opt.kmin)) #npow lya-tot 11... I think this is zero.

mf_output_df = pd.DataFrame(zbin_msr_matrix) #I prefer working with pandas dataframes
mf_output_df.columns = mf_msrmnts

zab_centers = opt.find_za(opt.zbin_centers) #transforming the lyb zbins to the equivalent, lower, lya zbins.
len_zab = len(zab_centers)

#Gives corresponding lya bin for each lyb bin. organized by increasing z.
bin_zab=np.ones(len_zab)*np.nan
for i in range(len_zab):
    for j in range(len_zab):
        if (zab_centers[i]>opt.zbin_edges[j])&(zab_centers[i]<opt.zbin_edges[j+1]):
            bin_zab[i] = (opt.zbin_centers[j])

#Calulating mean flux for LYB forest by dividing TOT mean flux with corresponding, lower z, LYA mean flux
mf_lyb = np.ones(len_zab)*np.nan #nan until proven otherwise
for i in range(len_zab):
    if bin_zab[i] in mf_output_df.z.values:
        za_idx = mf_output_df.z == bin_zab[i]
        ztot_idx = i
        mf_lyb[i] = mf_output_df.mf_tot[ztot_idx]/mf_output_df.mf_a[za_idx]
mf_output_df['mf_b'] = mf_lyb

# CONTINUUM CORRECTION. F_true = F_est * (1-deltaC/C_true)
if inis.continuum_correction:
    C_a_ratio = opt.continuum_correction(opt.zbin_centers)
    C_t_ratio = opt.continuum_correction(opt.zbin_centers-0.2)  #Matt M:  the 0.2 is from eyeballing where the mean flux is about the same in total relative to alpha
    C_b_ratio = np.abs(np.abs(C_a_ratio)-np.abs(C_t_ratio))
    mf_output_df['mf_a'] = mf_output_df['mf_a']*(1-C_a_ratio)
    mf_output_df['mf_tot'] = mf_output_df['mf_tot']*(1-C_t_ratio)
    mf_output_df['mf_b'] = mf_output_df['mf_b']*(1-C_b_ratio)

#Flux pdf for lya (for fun)... possibly a bug
if inis.save_flux_pdf:
    lya_fpdf_arr = np.reshape(np.concatenate((lya_flux_pdf)),(nqso,opt.zbinlen,len(pdf_bins)-1))
    fpdf_x = np.concatenate([np.nanmean(lya_fpdf_arr,0)[i] for i in range(opt.zbinlen)])
    np.savetxt(inis.save_flux_pdf_path,fpdf_x)

########################
#### Power Spectrum ####
########################

pk_msrmnts = opt.pk_msrmnts
n_pk_msrmnts = len(pk_msrmnts)
chunklengthlist_a = [np.array([0])]*opt.zbinlen; chunklengthlist_t = [np.array([0])]*opt.zbinlen; chunklengthlist_c = [np.array([0])]*opt.zbinlen   #to output list of pixel lengths in each redshift bin
znk_matrix = np.zeros((opt.zbinlen,n_pk_msrmnts,opt.kbinlen)) #  7 zbins, 10 measurements, 20 kbins
print("Power Spectra Calculation")
f = open(inis.save_kzq_pk_path,'w') #july1
for zidx in range(opt.zbinlen): #aug20
    #looping through the zidx'th zbin
    opt.updt(opt.zbinlen, zidx)
    msrmnt_in_kbin = np.zeros((n_pk_msrmnts,opt.kbinlen)) #measurements in the k'th bin
    count_in_kbin = np.zeros((n_pk_msrmnts,opt.kbinlen)) #counts in the k'th bin
    msrmnt_in_kbin[0] = opt.kbin_centers
    count_in_kbin[0] = np.ones_like(opt.kbin_centers) # these have to be ones. I just divide by one later.
    count_in_kbin[6] = np.ones_like(opt.kbin_centers) #counts for npix_aa
    count_in_kbin[7] = np.ones_like(opt.kbin_centers) #counts for npix_tt
    count_in_kbin[8] = np.ones_like(opt.kbin_centers) #counts for npix_ab
    msrmnt_in_kbin[-1] = np.ones_like(opt.kbin_centers) * opt.zbin_centers[zidx]
    count_in_kbin[-1] = np.ones_like(opt.kbin_centers) #counts for redshift bins
    for qidx in range(nqso): #aug20
        #looping through the qidx'th quasar
        ###################################################################
        #################### LYA FOREST: P ALPHA ALPHA ####################
        ###################################################################
        name = qso_arr[qidx].name
        #[kmaxa, kmaxb, kmaxab] = qso_arr[qidx].kmax_arr  #Matt M:added this

        zpix_a = qso_arr[qidx].get_zpix(opt.lya_rest)
        zmask_a = qso_arr[qidx].get_zmask(forest=(opt.lya_min,opt.lya_max,opt.lya_rest),
                                        zpix=zpix_a,zidx=zidx,zedges=opt.zbin_edges,name=name)
        zpix_tot = qso_arr[qidx].get_zpix(opt.lyb_rest)
        zmask_tot = qso_arr[qidx].get_zmask(forest=(opt.lyb_min,opt.lyb_max,opt.lyb_rest),
                                        zpix=zpix_tot,zidx=zidx,zedges=opt.zbin_edges,name=name)
        nchunks_a = opt.how_many_chunks(zmask_a) #usually 1
        nchunks_tot = opt.how_many_chunks(zmask_tot)

        Lya_str = 'VIS' if opt.overlap_maxwav < opt.lya_rest*(1+opt.zbin_centers[zidx]) else 'UV'
        Lyb_str = 'VIS' if opt.overlap_maxwav < opt.lyb_rest*(1+opt.zbin_centers[zidx]) else 'UV'
        kmax = [qso_arr[qidx].kmax[Lya_str], qso_arr[qidx].kmax[Lyb_str], 2**.5*(1./qso_arr[qidx].kmax[Lya_str]**2 +1./qso_arr[qidx].kmax[Lyb_str]**2)**-0.5]
        #print("kmax = ", qso_arr[qidx].kmax,  kmax[2])
        #exit()

        for chunk_idx_a in range(nchunks_a):
            #looping through each chunk. Probably 1. maybe 2. not more than 3.
            zmask_a_sub = opt.get_chunks(zmask_a)[chunk_idx_a]
            chunk_length = len(zmask_a_sub) #opt.get_chunk_length(zmask_a)

            #if chunk_length >1200:
            #    print("very long segment: ", np.min(zmask_a), np.max(zmask_a))
            #print("chunck = ", len(zmask_a_sub), chunk_length)

            if chunk_length >opt.min_pix:
                chunklengthlist_a[zidx] = np.append(chunklengthlist_a[zidx],  chunk_length) #make a list of number of pixels used
                #auto-power if passes minimum pized constraint
                kpix,pk = qso_arr[qidx].get_autopower(mf_output_df.mf_a[zidx],zmask_a_sub)
                for kidx in range(opt.kbinlen):
                    npow = mf_output_df.npow_a.values[zidx]
                    kmask = qso_arr[qidx].get_kmask(kpix=kpix,kidx=kidx,kedges=opt.kbin_edges, kmax=kmax[0])
                    if inis.individual_qso_kmax == 1 and qso_arr[qidx].get_kmax() < opt.kbin_centers[kidx]:  #Matt M: added this and next line to not add highest k
                            break
                    pk_sub = qso_arr[qidx].get_pk_subsets(kpix=kpix,pk=pk,zmask=zmask_a_sub,kmask=kmask,corr_tag=tag,npow=npow) #pk in the kidx'th kbin                                                                      #Matt M: This is where resolution correction occurs.
                    msrmnt_in_kbin[1,kidx] += np.sum(pk_sub) #adding paa in k,z bin for the qidx'th quasar
                    count_in_kbin[1,kidx] += len(pk_sub) #adding number of paa pixels
                    msrmnt_in_kbin[6,kidx] += len(pk_sub) ##adding number of npix_aa pixels

                    s = "{4:g} 0 {0:g} {1:g} {2:g} {3:g}".format(opt.zbin_centers[zidx],opt.kbin_centers[kidx], np.sum(pk_sub),len(pk_sub),qidx) #july1 #writing to external file
                    #0 means paa
                    f.write(s+'\n') #july1

        ###################################################################
        #################### LYB FOREST: P TOTAL TOTAL ####################
        ###################################################################
        for chunk_idx_tot in range(nchunks_tot):
            zmask_tot_sub = opt.get_chunks(zmask_tot)[chunk_idx_tot]
            chunk_length = len(zmask_tot_sub) #opt.get_chunk_length(zmask_tot)
            #print("chunck ttot", chunk_length)
            if (chunk_length >opt.min_pix):
                chunklengthlist_t[zidx] = np.append(chunklengthlist_t[zidx],  chunk_length)
                kpix,pk = qso_arr[qidx].get_autopower(mf_output_df.mf_tot[zidx],zmask_tot_sub)
                for kidx in range(opt.kbinlen):
                    npow = mf_output_df.npow_tot.values[zidx]
                    kmask = qso_arr[qidx].get_kmask(kpix=kpix,kidx=kidx,kedges=opt.kbin_edges, kmax=kmax[1])
                    if inis.individual_qso_kmax == 1 and kmaxb < opt.kbin_centers[kidx]:  #Matt M: added this and next line to not add highest k
                            break
                    pk_sub = qso_arr[qidx].get_pk_subsets(kpix=kpix,pk=pk,zmask=zmask_tot_sub,kmask=kmask,corr_tag=tag,npow=npow)
                                       #Matt M: This is where resolution correction occurs.
                    msrmnt_in_kbin[2,kidx] += np.sum(pk_sub) #adding ptot in k,z bin for the qidx'th quasar
                    count_in_kbin[2,kidx] += len(pk_sub) #adding number of ptot pixels
                    msrmnt_in_kbin[7,kidx] += len(pk_sub) #adding number of npix_tt pixels

                    s = "{4:g} 1 {0:g} {1:g} {2:g} {3:g}".format(opt.zbin_centers[zidx],opt.kbin_centers[kidx],np.sum(pk_sub),len(pk_sub),qidx)
                    #1 means ptt
                    f.write(s+'\n') #july1

        # Loop through more chunks here. Check zmask_tot and zmask_a CROSS POWER
        ####################################################################
        #################### CROSS-POWER: P ALPHA-TOTAL ####################
        ####################################################################
        idx_a = 0
        idx_tot = 0
        while (idx_tot < nchunks_tot)&(idx_a < nchunks_a):
            #Keep looping until there are no more chunks!
            if (nchunks_a>0)&(nchunks_tot>0):
                #only proceed if both forests have non-zero chunks!
                za = zpix_a[opt.get_chunks(zmask_a)[idx_a]]
                #chunk_length_a = opt.get_chunk_length(zmask_a)
                ztot = zpix_tot[opt.get_chunks(zmask_tot)[idx_tot]]
                #chunk_length_tot = opt.get_chunk_length(zmask_tot)
                za_min = np.min(za)
                za_max = np.max(za)
                ztot_min = np.min(ztot)
                ztot_max = np.max(ztot)
                mask_chk_a = (zpix_a>np.max([ztot_min,za_min]))&(zpix_a<np.min([ztot_max,za_max])) #mask for alpha chunk
                mask_chk_tot = (zpix_tot>np.max([ztot_min,za_min]))&(zpix_tot<np.min([ztot_max,za_max])) #mask for beta chunk
                #print("chunks", chunk_length_a, chunk_length_tot, np.sum(mask_chk_a), np.sum(mask_chk_tot), mf_output_df.mf_tot[zidx])
                if (np.sum(mask_chk_a) >opt.min_pix)&(np.sum(mask_chk_tot)>opt.min_pix):
                    chunklengthlist_c[zidx] = np.append(chunklengthlist_c[zidx],  np.min([np.sum(mask_chk_a), np.sum(mask_chk_tot)]))
                    kpix,pab,qab,dlama, dlamb,resa,resb = qso_arr[qidx].cross_pk_fft(mask_lya=mask_chk_a,mask_lyb=mask_chk_tot,
                                          mf_lya=mf_output_df.mf_a[zidx],
                                          mf_lyb=mf_output_df.mf_tot[zidx])
                    npow = mf_output_df.npow_atot.values[zidx]
                    for kidx in range(opt.kbinlen): #aug20
                        kmask = qso_arr[qidx].get_kmask(kpix=kpix,kidx=kidx,kedges=opt.kbin_edges, kmax=kmax[2])
                        if inis.individual_qso_kmax == 1 and kmaxab < opt.kbin_centers[kidx]:  #Matt M: added this and next line to not add highest k
                            break
                        pab_sub,qab_sub = qso_arr[qidx].get_xpk_subsets(kpix,pab,qab,dlama,dlamb, resa,resb,tag,npow,kmask)

                        s = "{4:g} 2 {0:g} {1:g} {2:g} {3:g}".format(opt.zbin_centers[zidx],opt.kbin_centers[kidx],
                                                        np.sum(pab_sub),len(pab_sub), qidx) #july1 #average of pab_sub??
                        f.write(s+'\n')
                        msrmnt_in_kbin[3,kidx] += np.sum(pab_sub)
                        count_in_kbin[3,kidx] += len(pab_sub)
                        msrmnt_in_kbin[4,kidx] += np.sum(qab_sub)
                        count_in_kbin[4,kidx] += len(qab_sub)
                        msrmnt_in_kbin[8,kidx] += len(pab_sub) #npix_ab #sept25
                if za_max<ztot_min: # no overlap, increase alpha index to try to match redshift of ztot
                    idx_a +=1
                    idx_tot +=0
                elif ztot_max<za_min: # opposite way of no overlap, increase beta index to match redshift of z-alpha
                    idx_a +=0
                    idx_tot +=1
                else: # lya dla will always be in lyb forest. but lyb dla could possibly not be in lya forest (absorbs at low-z)
                    idx_tot +=1
            else:
                break
    znk_matrix[zidx] = msrmnt_in_kbin/count_in_kbin
f.close() #july1
opt.updt(opt.zbinlen, opt.zbinlen)
print("Done!\n")

#with open('Lyaskewerlengths.txt') as f:
#    for i in range(opt.zbinlen):
#        #f.write
#        f.write(str(chunklengthlist[i]))
#Matt M: for my own use to know lengths of skewers: should delete
# np.savez("Lyaskewerlengths", chunklengthlist_a[0], chunklengthlist_a[1], chunklengthlist_a[2], chunklengthlist_a[3], chunklengthlist_a[4], chunklengthlist_a[5],chunklengthlist_a[6], kwds=['3.0','3.2', '3.4', '3.6', '3.8', '4.0', '4.2']) # for i in range(opt.zbinlen)], header=str())
# np.savez("LyTskewerlengths", chunklengthlist_t[0], chunklengthlist_t[1], chunklengthlist_t[2], chunklengthlist_t[3], chunklengthlist_t[4], chunklengthlist_t[5],chunklengthlist_t[6], kwds=['3.0','3.2', '3.4', '3.6', '3.8', '4.0', '4.2'])
# np.savez("LyCskewerlengths", chunklengthlist_c[0], chunklengthlist_c[1], chunklengthlist_c[2], chunklengthlist_c[3], chunklengthlist_c[4], chunklengthlist_c[5],chunklengthlist_c[6], kwds=['3.0','3.2', '3.4', '3.6', '3.8', '4.0', '4.2'])

#Finding Lyman beta power
for i in range(len_zab):
    if bin_zab[i] in opt.zbin_centers:
        za_idx = np.where(opt.zbin_centers == bin_zab[i])[0][0]
        znk_matrix[i][5] = znk_matrix[i][2]-znk_matrix[za_idx][1]

#Making 3d pk matrix into 2d pk data frame
x = pd.DataFrame(znk_matrix[0].T,columns=pk_msrmnts)
for i in range(1,opt.zbinlen):
    x = x.append(pd.DataFrame(znk_matrix[i].T,columns=pk_msrmnts))


########## SUBTRACTING REDSIDE METAL POWER ##########
if inis.subtract_metal_power and inis.use_obs:
    metal_power = np.concatenate([np.loadtxt('../data/obs/pk_xs_avg.txt')]*opt.zbinlen) #subtracting metals again sep25
    x.paa = x.paa.values-metal_power
    x.ptt = x.ptt.values-metal_power
    #x.Pbb = x.Pbb.values-metal_power #redside metals won't affect beta power. contribution cancels out

###########################################################
#Saving figures
##########################################################
mf_output_df.to_csv(inis.save_mf_path,index=False)
x.to_csv(inis.save_pk_path,index=False)


print("OVI added?                           : ", inis.add_ovi)
print("SiIII added?                         : ", inis.add_sithree)
print("LYB added?                           : ", inis.add_beta)
print("Carswell Resolution?                 : ", bool(inis.carswell_res))
print("Don't use the overlap region at all? : ", inis.no_overlap)
if (inis.remove_dla)&inis.cat_name.startswith('obs'):
    print("DLA's removed?                       : ", inis.remove_dla)
else:
    print("DLA's removed?                       : False")
print("Metals Subtracted?                   : ", bool(inis.subtract_metal_power and inis.use_obs)) #see line 310
print("Continuum Corrected?                 : ", inis.continuum_correction)
print("Log-binning for wavenumber (k)?      : ", inis.log_kbinning)
print("saved pre-boot mf file here          : ", inis.save_kzq_mf_path)
print("saved pre-boot pk file here          : ", inis.save_kzq_pk_path)
print("saved mf here                        : ", inis.save_mf_path)
print("saved pf here                        : ", inis.save_pk_path)
print("---- Script complete ----")
