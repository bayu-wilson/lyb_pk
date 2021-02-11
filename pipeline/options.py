# least unique -> most unique

import numpy as np
import sys
import inis

#Cosmological Parameters
c_kms = 2.99792458e5 #km/s
h = 0.678
H_0 = h*100 #km/s/Mpc
omega_b = 0.0484
omega_cdm = 0.2596
omega_m0 = omega_b + omega_cdm

#LYMAN ALPHA
lya_min = 1045 # 1159.11 # 1045 #aug2
lya_max = 1185 # 1201.78 # 1185 #aug2
lya_rest =  1215.67
f_osc_lya = 0.4162
xs_alpha = f_osc_lya * lya_rest
v_alpha = 2271 #km/s Matt M: appears in a function that is never used

#LYMAN BETA
lyb_min = 978 # 881.717 # 978 #aug2
lyb_max = 1014 # 1014 # 999.8 # 1014 #aug2
lyb_rest = 1025.72
f_osc_lyb = 0.0791
xs_beta = f_osc_lyb * lyb_rest
v_beta = -1801 #km/s Matt M: appears in a function that is never used


#O VI
ovi_min = 978


ovi_max = 1014
ovi_rest_d1 = 1031.91 # 1032  #Matt M changed to add sig fig
fosc_ovi_d1 = 0.13290
ovi_rest_d2 = 1037.61 #Matt M changed to add sig fig
fosc_ovi_d2 = 0.06609
xs_ovi = 10**(-2.5) * xs_alpha
ovi_factor = fosc_ovi_d2/fosc_ovi_d1


#Si III tau_SIII = 10^-3.5 to 10^-3  tau_HI"  SiIII, lambda= 1206.5 \AA
sithree_min = 1045 # 1159.11 # 1045 #aug2
sithree_max = 1185 # 1201.78 # 1185 #aug2
sithree_rest_d1 = 1206.5
sithree_rest_d2 = 1206.6
xs_sithree = 10**(-3.5) * xs_alpha
# https://physics.nist.gov/cgi-bin/ASD/lines1.pl?spectra=Si%20III&limits_type=0&low_w=1200&upp_w=1250&unit=0&submit=Retrieve%20Data&de=0&format=0&line_out=0&en_unit=0&output=0&bibrefs=1&page_size=15&show_obs_wl=1&show_calc_wl=1&unc_out=1&order_out=0&show_av=2&tsb_value=0&A_out=0&intens_out=on&allowed_out=1&forbid_out=1&conf_out=on&term_out=on&enrg_out=on&J_out=on&level_id=014003.000007


#REDSHIFT BINS
zmin = inis.zmin
zmax = inis.zmax
zbinlen = inis.zbinlen
zbin_centers = np.round(np.linspace(zmin,zmax,zbinlen),2)
dz = zbin_centers[1]-zbin_centers[0]
zbin_edges = np.linspace(zmin-dz/2,zmax+dz/2,zbinlen+1)

#flux pdf
pdf_bins = np.arange(-0.025,1.05,.05)

mf_msrmnts = ['mf_a', 'nvar_a','dloglambda_a', 'npow_a',
           'mf_tot', 'nvar_tot','dloglambda_tot', 'npow_tot','z','mf_b',
           'var_atot','npow_atot']
pk_msrmnts = ['k','paa','ptt','pab','qab','pbb','npix_aa','npix_tt','npix_ab','z'] #sept25

#Wavenumbers
if inis.log_kbinning:
    kmin = inis.kmin #3e-3
    kmax = inis.kmax #6e-2
    kbinlen = inis.kbinlen #13?
    kbin_centers = np.logspace(np.log10(kmin),np.log10(kmax),kbinlen)
    dk = (np.log10(kmax) - np.log10(kmin))/(kbinlen-1) #feb4
    kbin_edges = np.logspace(np.log10(kmin)-dk/2,np.log10(kmax)+dk/2,kbinlen+1)

else:
    kmin = inis.kmin #3e-3 #0.0
    kmax = inis.kmax #- 3e-3#0.06 #0.06
    kbinlen = inis.kbinlen
    kbin_centers = np.linspace(kmin,kmax,kbinlen)
    dk = (kmax - kmin)/(kbinlen-1) #feb4
    kbin_edges = np.linspace(kmin-dk/2,kmax+dk/2,kbinlen+1)

#DLAs
DLAcat_file = "Data/XQ-100_DLA_catalogue.txt"
DLA_cut_factor = 1. #measured in equivalent widths (goes out half this distance in each direction from center

#Bad Pixels
min_pix = 75  #changed from 100 to allow for z=3.6 bin
min_trans = -100

#SIGMA_RMS (resolution) in XQ-100 files
#R_UV = 20 #RMS RESOLUTION (sigma_R) = FWHM/(2*sqrt(2ln(2))).
#R_VIS = 11 #converting FWHM to sigma_R.

#CARSWELL+18 values
R_UV_carswell = 41.52075368/(2*np.sqrt(2*np.log(2)))
R_VIS_carswell = 23.60134842/(2*np.sqrt(2*np.log(2)))
overlap_maxwav = 5599.14
overlap_minwav = 5499.85

ADC_off_qsos  = ['J1013+0650', 'J1416+1811', 'J1524+2123', 'J1542+0955', 'J1552+1005', 'J1621-0042', 'J1723+2243'] #these are the qsos for which the wavelength calibration is likely to be worst

def updt(total, progress):
    """
    Displays or updates a console progress bar.

    Original source: https://stackoverflow.com/a/15860757/1391441
    """
    barLength, status = 20, ""
    progress = float(progress) / float(total)
    if progress >= 1.:
        progress, status = 1, "\r\n"
    block = int(round(barLength * progress))
    text = "\r[{}] {:.0f}% {}".format(
        "#" * block + "-" * (barLength - block), round(progress * 100, 0),
        status)
    sys.stdout.write(text)
    sys.stdout.flush()

def find_za(lyb_z):
    """
    Gives redshift of pixel in lyb forest from lya transition
    INPUT:
    redshift of pixel from lyb transition
    OUTPUT:
    redshift of the same pixel but at lya transition (lower redshift)
    i.e.)
    z_beta:          3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2
    z_alpha(z_beta): 2.4, 2.5, 2.7, 2.9, 3.0, 3.2, 3.4
    """
    return (lyb_z + 1)*(lyb_rest/lya_rest)-1

def find_zb(lya_z):
    """
    Gives redshift of pixel in lya forest from lyb transition
    INPUT:
    redshift of pixel from lyb transition
    OUTPUT:
    redshift of the same pixel but at lya transition (lower redshift)
    i.e.)
    z_beta:          3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2
    z_alpha(z_beta): 2.4, 2.5, 2.7, 2.9, 3.0, 3.2, 3.4
    """
    return (lya_z + 1)*(lya_rest/lyb_rest)-1

def window(k,p,R):
    """
    Deconvolution kernel.
    k: velocity wavenumber
    p: pixel width (in velocity, km/s)
    R: resolution
    """
    gaussian = np.exp(-0.5*k**2*R**2)
    tophat =  np.sinc(0.5*k*p/np.pi) # THIS PI NEEDS TO BE HERE#np.sin(k*p/2)/(k*p/2)
    deconvolution_kernel = gaussian*tophat
    return deconvolution_kernel

def metal_model_pab(a,b,k,pab):
    n = 1 + a*np.cos(k*v_alpha) + b*np.cos(k*v_beta) * a*b*np.cos(k*(v_beta-v_alpha))
    return n*pab

def comoving_distance(v):
    z = v/c_kms
    return v/(H_0*np.sqrt(omega_m0*(1+z)))

def find_EW(NHI,line,z_dla):
    # Galaxy Formation & Evolution, Frank van den Boschm Pg. 712, Eq. 16.113
    W_alpha = 7.3 * (10**NHI/10**20)**(0.5)*(1+z_dla)#angstroms
    W_beta = xs_beta/xs_alpha * W_alpha
    if line == 'alpha':
        return  W_alpha
    if line == 'beta':
        return W_beta #angstroms

def how_many_chunks(mask):
    #See line ~185 in main.py
    increases,decreases = 0,0
    idx_inc = []
    idx_dec = []
    for i in range(len(mask)-1):
        if mask[i+1]>mask[i]: #DLA on left edge only
            increases += 1
            idx_inc.append(i+1)
        if mask[i]>mask[i+1]: #DLA on right edge only
            decreases += 1
            idx_dec.append(i)
    if np.sum(mask) == 1:
        return 1
    if increases != decreases: #DLA on one edge
        return np.max([increases,decreases])
    if np.all(mask==0): #No data
        return 0
    if not idx_inc: #No DLA
        return 1
    if (increases == decreases)&(idx_inc[0]<idx_dec[0]): #DLA on both edges
        return increases
    if (increases == decreases)&(idx_inc[0]>=idx_dec[0]): #DLA in middle area
        return increases+1
    #return np.max([increases,decreases])

def get_chunk_length(mask):
    chunk_edges = []
    for i in range(len(mask)-1):
        if mask[i+1]>mask[i]:
            chunk_edges.append(i+1)
        if mask[i]>mask[i+1]:
            chunk_edges.append(i+1)
    chunk_edges.sort()
    
    if len(chunk_edges) == 2:
        if chunk_edges[1] - chunk_edges[0] >1200:
            print("ce1", chunk_edges[1], chunk_edges[0])
        return chunk_edges[1] - chunk_edges[0]
    else:
        if chunk_edges[0] >1200:
            print("ce2", chunk_edges)
        return chunk_edges[0]
        
def get_chunks(mask):
    #Used in line ~190 of main.py
    chunk_edges = []
    for i in range(len(mask)-1):
        if mask[i+1]>mask[i]:
            chunk_edges.append(i+1)
        if mask[i]>mask[i+1]:
            chunk_edges.append(i+1)
    chunk_edges.sort()

    splitted = np.split(mask, chunk_edges)
    splitted_idx = np.split(np.arange(len(mask)), chunk_edges)
    subset_chunk_edges = []
    for s in range(len(splitted)):
        if np.all(splitted[s]==True):
            subset_chunk_edges.append(splitted_idx[s])
    #if len( chunk_edges)>2:
    #    print("subset chunk edges", subset_chunk_edges)
    return subset_chunk_edges

def continuum_correction(z):
    # returns (delta C)/(C_true) where C stands for continuum
    return 1.58e-5*(1+z)**5.63

def find_err_mf_beta(mf_beta,mf_alpha,err_mf_alpha,mf_total,err_mf_total):
    return mf_beta*np.sqrt((err_mf_alpha/mf_alpha)**2+(err_mf_total/mf_total)**2)

def find_res_uncertainty(k,z,pk):
    res_uncertainty_array = np.zeros_like(z)
    for i in range(len(z)):
        if z[i]<=3.6: #UVB ARM
            rel_err_sigmaR = 0.05369698
            sigmaR = 17.63224063
        elif z[i]>3.6: #VIS ARM
            rel_err_sigmaR = 0.04838776
            sigmaR = 10.02256986
        res_uncertainty_array[i]=(2*(k[i]*sigmaR)**2*rel_err_sigmaR*pk[i])
    return res_uncertainty_array

a=2.88314e-26#cm^2
b=2.02161e-8
# c=3e+10 #cm/s
def get_dla_profile(v_over_c,N_HI):
  tau = a*N_HI/((v_over_c)**2+b**2)
  return np.exp(-tau)

ab=6.2164025654e-27#cm^2
bb=5.16495523076e-9
def get_dlb_profile(v_over_c,N_HI):
  tau = ab*N_HI/((v_over_c)**2+bb**2)
  return np.exp(-tau)
