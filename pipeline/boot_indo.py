#!/usr/bin/env python

import numpy as np
import inis
import options as opt
import matplotlib.pyplot as plt
import pandas as pd

np.seterr(divide='ignore', invalid='ignore') # I am ignoring this error because it is not important

print("reading in ", inis.save_kzq_mf_path, inis.save_kzq_pk_path, " for bootstraps")
table_mf = np.loadtxt(inis.save_kzq_mf_path)
qtable_pk = np.loadtxt(inis.save_kzq_pk_path)
#quasar index, z-index, sum of lya flux in z-bin, sum/num of lya pixels in z-bin, sum of lyb (total) flux in z-bin, sum/num of lyb (total) pixels in z-bin
#quasar index, type of power measurement (numbered from 0-2 for paa,ptt, and pab), redshift, wavenumber (k), sum of power in k,z bin, sum of pixels in k,z bin.

M = inis.M #number of bootstrap samples
np.random.seed(1)
nqso = inis.nqso
rand_matrix = [np.floor(np.random.rand(nqso)*nqso).astype(int) for i in range(M)]
boot_mf_arr = np.zeros((M,opt.zbinlen*2)) #empty matrix to be filled with mf values for each bootstrap sample
boot_pk_arr = np.zeros((M,opt.zbinlen*3,opt.kbinlen))#empty matrix to be filled with pk values for each bootstrap sample

print("Doing mean flux bootraps")
for m in range(M):
    rand_qidxs = rand_matrix[m] #array of 100 random integers with replacement
    mf_arr = np.zeros((opt.zbinlen*2)) #empty matrix to be filled with lya and lya-lyb mf values in z bins
    mf_N_arr = np.zeros((opt.zbinlen*2)) #empty matrix to be filled with total number of pixels in each z bin
    for qidx in rand_qidxs: #looping through these 100 integers (quasar indices)
        qtable_mf = table_mf[table_mf[:,0] == qidx] #choosing portion of table matching the quasar index
        for zidx in range(opt.zbinlen):
            ztable_mf = qtable_mf[qtable_mf[:,1] == opt.zbin_centers[zidx]] #choosing portion of table matching zbin
            mf_arr[zidx] += np.sum(ztable_mf[:,2]) #summing alpha mf into zbin
            mf_N_arr[zidx] += np.sum(ztable_mf[:,3])  #summing alpha pixels into zbin
            mf_arr[opt.zbinlen+zidx] += np.sum(ztable_mf[:,4]) #summing alpha-beta mf into zbin
            mf_N_arr[opt.zbinlen+zidx] += np.sum(ztable_mf[:,5]) #summing alpha-beta pixels into zbin
    boot_mf_arr[m] = mf_arr/mf_N_arr # sum(mf)/num(mf) = <mf> in zbin for the mth boot sample

print("Doing power spectrum bootstraps")


if inis.BayusBootstrap:
    for m in range(M):
        opt.updt(M,m)
        rand_qidxs = rand_matrix[m] #array of 100 random integers with replacement
        pk_arr = np.zeros((opt.zbinlen*3,opt.kbinlen)) #empty matrix to be filled with pk values in k,z bins
        N_arr = np.zeros((opt.zbinlen*3,opt.kbinlen)) #empty matrix to be filled with total number of pixels for in each k,z bins
        for pidx in range(3):
            for qidx in rand_qidxs: #looping through these 100 integers (quasar indices)
                qmask = (qtable_pk[:,0] == qidx)&(qtable_pk[:,1] == pidx)
                qsub_table = qtable_pk[qmask].T #choosing portion of table matching the quasar index
                for zidx in range(opt.zbinlen):
                    zmask = qsub_table[2] == opt.zbin_centers[zidx] #choosing portion of table matching zbin
                    zsub_table = qsub_table.T[zmask].T
                    for kidx in range(opt.kbinlen):
                        kmask = np.round(zsub_table[3],5) == np.round(opt.kbin_centers[kidx],5) #this is fine
                        pk_arr[pidx*opt.zbinlen+zidx][kidx] += np.sum(zsub_table[4][kmask])
                        N_arr[pidx*opt.zbinlen+zidx][kidx] += np.sum(zsub_table[5][kmask])
        boot_pk_arr[m] = pk_arr/(N_arr+1e-10)

else: #mine
    pk_qso = np.zeros((nqso, opt.zbinlen*3,opt.kbinlen))
    N_qso = np.zeros((nqso, opt.zbinlen*3,opt.kbinlen))
    for qidx in range(nqso): #looping through these 100 integers (quasar indices)
        for pidx in range(3):
            qmask = (qtable_pk[:,0] == qidx)&(qtable_pk[:,1] == pidx)
            qsub_table = qtable_pk[qmask].T #choosing portion of table matching the quasar index
            for zidx in range(opt.zbinlen):
                zmask = qsub_table[2] == opt.zbin_centers[zidx] #choosing portion of table matching zbin
                zsub_table = qsub_table.T[zmask].T
                for kidx in range(opt.kbinlen):
                    kmask = np.round(zsub_table[3],5) == np.round(opt.kbin_centers[kidx],5) #this is fine
                    pk_qso[qidx][pidx*opt.zbinlen+zidx][kidx] = np.sum(zsub_table[4][kmask])
                    N_qso[qidx][pidx*opt.zbinlen+zidx][kidx] = np.sum(zsub_table[5][kmask])
                    #print(qidx, pidx, zidx, kidx, " num pix ", np.sum(zsub_table[5][kmask]), zsub_table[5][kmask])
                #Need to average Pk in some way
    for m in range(M):
        opt.updt(M,m)
        rand_qidxs = rand_matrix[m] #array of 100 random integers with replacement
        pk_arr = np.zeros((opt.zbinlen*3,opt.kbinlen)) #empty matrix to be filled with pk values in k,z bins
        N_arr = np.zeros((opt.zbinlen*3,opt.kbinlen)) #empty matrix to be filled with total number of pixels for in each k,z bins
        for pidx in range(3):
            for qidx in rand_qidxs: #looping through these 100 integers (quasar indices)
                for zidx in range(opt.zbinlen):
                    for kidx in range(opt.kbinlen):
                        pk_arr[pidx*opt.zbinlen+zidx][kidx] += pk_qso[qidx][pidx*opt.zbinlen+zidx][kidx]
                        N_arr[pidx*opt.zbinlen+zidx][kidx] += N_qso[qidx][pidx*opt.zbinlen+zidx][kidx]
        boot_pk_arr[m] = pk_arr/(N_arr+1e-10)
                    
        
mf_summed = pd.read_csv(inis.save_mf_path)
mf_boot = np.reshape(np.mean(boot_mf_arr,axis=0),(2,opt.zbinlen))

lya_factor = (np.nansum(mf_summed.mf_a)*np.nansum(mf_summed.mf_a)) / (np.nansum(mf_boot[0])*np.nansum(mf_boot[0]))
boot_pk_corr = lya_factor * boot_pk_arr
#print("lyafactor = ", lya_factor)
#print("boot = ", boot_pk_arr)
#exit(5)



#P_XY = [(<F_X> * <F_y>)_bootstrapped / (<F_X> * <F_Y>)_summed up]^2 * P_XY_bootstrapped
###########################################################################
################# NEW ERRORBARS AND SAVING THE BOOTSTRAPS #################
###########################################################################

print("saving")
mfdata = pd.read_csv(inis.save_mf_path)
pkdata = pd.read_csv(inis.save_pk_path)
if inis.save_boot_mf:
    np.savetxt(inis.save_boot_mf_path,boot_mf_arr.T)
    mf = np.loadtxt(inis.save_boot_mf_path)
    err_mfa = np.cov(mf).diagonal()[0:opt.zbinlen]**0.5
    err_mft = np.cov(mf).diagonal()[opt.zbinlen:opt.zbinlen*2]**0.5
    err_mfb = opt.find_err_mf_beta(mfdata.mf_b.values,mfdata.mf_a.values,err_mfa,
                                                      mfdata.mf_tot.values,err_mft)
    if inis.save_mf_with_err:
        columns = ["z","mfa","err_mfa","mft","err_mft","mfb","err_mfb"]
        mf_everything = np.column_stack((mfdata.z,
                        mfdata.mf_a, err_mfa,
                        mfdata.mf_tot, err_mft,
                        mfdata.mf_b, err_mfb))
        df_meanflux = pd.DataFrame(mf_everything,columns=columns)
        #df_meanflux.mft[1] = np.nan
        df_meanflux.to_csv(inis.save_mf_with_err_path, index=False)
        #pd.savetxt(inis.save_mf_with_err_path,mf_everything)

########################################################################
#add uncertainty to mean flux in bootstap errors
########################################################################
#if add_meanfluxerror_bootstrap:
####################################################################
#add 10% uncertainty for resolution
####################################################################
#if add_error_resolution:
#    sigma_dict = {'UV': 15., 'VIS': 8.} #effective resolutions



        
if inis.save_boot_pk:
    np.savetxt(inis.save_boot_pk_path,np.reshape(boot_pk_corr,(M,opt.zbinlen*3*opt.kbinlen)).T)
    p = np.loadtxt(inis.save_boot_pk_path)
    N_kz = opt.zbinlen*opt.kbinlen
    covPk = np.cov(p)
    #covPk =  pd.DataFrame(p).cov
    pk_err_diag = np.diagonal(covPk)**0.5
    err_paa = pk_err_diag[0:N_kz]
    err_ptt = pk_err_diag[N_kz:N_kz*2]
    err_pab = pk_err_diag[N_kz*2:N_kz*3]

    print("err_paa = ", np.cov(p).diagonal()[0:N_kz])
    print("err_ptt = ", np.cov(p).diagonal()[N_kz:2*N_kz])
    print("err_pab = ", np.cov(p).diagonal()[2*N_kz:3*N_kz])
      
    err_paa_sub = err_paa[0*opt.kbinlen:3*opt.kbinlen]
    err_ptt_sub = err_ptt[4*opt.kbinlen:7*opt.kbinlen]
    err_pbb_sub = np.sqrt(err_paa_sub**2+err_ptt_sub**2)
    err_pbb = np.ones(N_kz)*np.nan
    err_pbb[4*opt.kbinlen:7*opt.kbinlen] = err_pbb_sub

    if inis.save_pk_with_err:
        columns = ['k','z','paa','err_paa','N_aa','ptt', 'err_ptt','N_tt','pab','err_pab','N_ab','pbb','err_pbb']
        pk_everything = np.column_stack((pkdata.k,pkdata.z,
                        pkdata.paa, err_paa, pkdata.npix_aa,
                        pkdata.ptt, err_ptt, pkdata.npix_tt,
                        pkdata.pab,err_pab, pkdata.npix_ab,
                        pkdata.pbb, err_pbb,))
        df_pk = pd.DataFrame(pk_everything,columns=columns)
        df_pk.to_csv(inis.save_pk_with_err_path + '.csv', index=False)

opt.updt(M, M)
print("Saving bootstraps here:\n{0}\n{1}".format(inis.save_boot_mf_path,inis.save_boot_pk_path))
print("Saving new datatables here:\n{0}\n{1}".format(inis.save_mf_with_err_path,inis.save_pk_with_err_path+ '.csv'))

# ### INCLUDING RESOLUTION UNCERTAINTY
# sigma_res_aa = opt.find_res_uncertainty(pkdata.k,pkdata.z,pkdata.Paa)
# err_paa = np.sqrt(err_paa**2+sigma_res_aa**2)
# sigma_res_tt = opt.find_res_uncertainty(pkdata.k,pkdata.z,pkdata.Ptot)
# err_ptt = np.sqrt(err_ptt**2+sigma_res_tt**2)
# sigma_res_ab = opt.find_res_uncertainty(pkdata.k,pkdata.z,pkdata.Pab)
# err_pab = np.sqrt(err_pab**2+sigma_res_ab**2)
# sigma_res_bb = opt.find_res_uncertainty(pkdata.k,pkdata.z,pkdata.Pbb)
# err_pbb = np.sqrt(err_pbb**2+sigma_res_bb**2)

# N = 91
# boot_mat = np.reshape(boot_arr,(M,N))
# # [np.mean(boot_mat[:,i]) for i in range(91)]
# # np.reshape(boot_arr,(10,91))[:,90]
# # np.array([paa[:,i] - np.mean(paa,axis=1) for i in range(3000)]).T
# # np.array([boot_mat.T[:,i] - np.mean(boot_mat.T,axis=1) for i in range(10)]).T
#
# paa_minus_paamean = np.array([boot_mat[i] - np.mean(boot_mat,axis=0) for i in range(M)]).T
# Cij = np.zeros((N,N))
# for i in range(N):
#     for j in range(N):
#         Cij[i][j] = np.mean(paa_minus_paamean[i]*paa_minus_paamean[j])
# rij = np.zeros((N,N))
# for i in range(N):
#     for j in range(N):
#         rij[i][j] = Cij[i][j]/np.sqrt(np.abs(Cij[i][i]*Cij[j][j]))
# plt.imshow(rij)
# plt.show()
