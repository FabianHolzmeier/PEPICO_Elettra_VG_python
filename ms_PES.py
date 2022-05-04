# -*- coding: utf-8 -*-
"""
Created on Fri Jan 21 15:02:02 2022

@author: holzme33
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.optimize import curve_fit

PATH = "B:\\analysis\\tBMA\PEPICO\\TBMA_automatic_factor_53_88_93\\"
ide = "TBMA_combined"
hv = 92     # photon energy in eV

tof = np.loadtxt(PATH+ide+"_tof_full.txt")
en = np.loadtxt(PATH+ide+"_en_full.txt")
mat = np.loadtxt(PATH+ide+"_PEPICO_matrix_full.txt")


# mass spectrum
ms = np.sum(mat,axis=1)

fig,ax = plt.subplots()
ax.plot(tof,ms,'k')
ax.set_xlabel("TOF ($\mu$s)")
ax.set_ylabel("Integrated TOF counts")
ax.set_title(ide+" TOF spectrum")

#%% mass calibration

def calibrate(m1,m2,t1,t2):
    # mass calibration of TOF spectra
    # m = alpha*t**2+beta (mass-to-charge is proportional to the square of the TOF)
    alpha = (m1-m2)/(t1**2-t2**2)
    beta = m1-alpha*t1**2
    return alpha*tof**2+beta

mass = calibrate(29,142,3.154,6.748)        # reference calibrations in OneNote

fig,ax = plt.subplots()
ax.plot(mass,ms,'r')
ax.set_xlabel("m/z (amu)")
ax.set_ylabel("Integrated TOF counts")
ax.set_title(ide+" mass spectrum")

#%% find peaks and integration limits

# threshold counts (only masses with higher peak maxima will be evaluated)
threshold = 5
# find peaks
peaks, _ = find_peaks(ms, height = threshold, distance=10)

# Gauss fitting function
def fit_Gauss(x,a,b,c):
    # a = amplitude, b = center, c = standard deviation sigma, d = offset
    return a*np.exp(-((x-b)**2)/(2*c**2))


# control figure showing the selected peaks and their Gaussian fits
fig,ax = plt.subplots()
ax.plot(mass,ms,'k')
ax.plot(mass[peaks],ms[peaks],"rx")
ax.set_xlabel("m/z (amu)")
ax.set_ylabel("Integrated TOF counts")
ax.set_title(ide+" mass spectrum")

# initialize peak centers, lower and upper integration limits
center,ll,ul = np.zeros((len(peaks),)),np.zeros((len(peaks),)),np.zeros((len(peaks),)) 
for i in range(len(peaks)):
    # if the fit is not good, play with the initial guesses (p0) or the boundary conditions (bounds)
    # poptx is an array containing the optimized parameters, pcovx the covariance
    poptx, pcovx = curve_fit(fit_Gauss, mass, ms, p0=[ms[peaks[i]],mass[peaks[i]],0.2], bounds = ([ms[peaks[i]]*0.75,mass[peaks[i]]-0.05,-np.inf],[ms[peaks[i]]*1.25,mass[peaks[i]]+0.05,1]))
    
    center[i] = poptx[1] # peak centers in m/z
    # plot fit result
    ax.plot(mass,fit_Gauss(mass,*poptx),'r:')
    # define integration limits
    ll[i] = poptx[1]-3*poptx[2]     # lower limit: center-3*std
    ul[i] = poptx[1]+3*poptx[2]     # upper limit: center+3*std

ax.legend(['mass spectrum','selected peaks','fit'])

#%% get mass-selected PES
msPES = np.zeros((len(en),len(peaks)))

for i in range(len(peaks)):
    # get indeces of the integration limits   
    ill = np.min(np.where(mass>ll[i]))     #lower limit
    iul = np.max(np.where(mass<ul[i]))     #upper limit
    
    # integrate PEPICO matrix over the selected mass range
    msPES[:,i] = np.sum(mat[ill:iul,:],axis=0)
    
masses = np.round(center)

# plot results
fig,ax = plt.subplots(figsize=(12, 6))
for i in range(len(masses)):
    ax.plot(hv-en,msPES[:,i],label=str(masses[i]))
    ax.fill_between(hv-en,msPES[:,i],alpha=0.3)
ax.plot(hv-en,np.sum(mat,axis=0),'k',linewidth=2,label='all masses')
ax.plot(hv-en,np.sum(msPES,axis=1),'k--',label='sum of sel. masses')
ax.legend(loc='right')
ax.set_xlabel("Binding Energy (eV)")
ax.set_ylabel("ms-PES Signal (arb. u.)")
ax.set_title(ide+" mass-selected PES")

#%% save results

np.savetxt(PATH+ide+'_ms_PES.out',msPES,delimiter='\t')
np.savetxt(PATH+ide+'_masses.out',masses,delimiter='\t')
np.savetxt(PATH+ide+'_binding_en.out',hv-en,delimiter='\t')