# -*- coding: utf-8 -*-
"""
Created on Fri Jan  7 09:32:27 2022

@author: holzme33
"""
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import time

# PATH = "B:\\analysis\\BSn\\PEPICO\\BSn_168\\"
PATH = "B:\\analysis\\BSn\\PEPICO\\BSn_168_fac0.65\\"
# PATH = "B:\\analysis\\tBMA\\PEPICO\\TBMA_088_tot\\"

ide = 'BSn_168'
MaxTof = '9'
tof_bin = '2'
x_step = '2'

# work function in eV
wf = 4.1

# reference PES (no coincidences)
path_ref = "B:\\VG\\rawdata\\"
file_ref = "BSn_165Sum"
DT_on = 1           # was the reference PES recorded with drift tube on?

px2eV = 0.85        # calibration factor px to eV
DT_shift = 0.6      # energy shift with drift tube on (eV)


# kinetic energies at which the measurement was done (eV) are retrieved automatically
path_en = "B:\\PEPICO\\rawdata\\"
file_en = "Energy_"+ide+".txt"
energy = np.loadtxt(path_en+file_en, skiprows=1, usecols=(2))
central_en = np.flip(np.unique(energy))


#******************************************************************************
# no need to change anything below this line
#******************************************************************************

files = []
# get the filenames of the raw data saved as .mat
for file in os.listdir(PATH):
    if file.endswith(".dat"):
        files.append(file)

for i in range(len(files)):
    tmp = np.loadtxt(PATH+ide+'_KE'+str(i)+'_PEPICO_matrix_'+MaxTof+'_'+tof_bin+'_'+x_step+'.dat')
    if i==0:
        comb_m = np.zeros((tmp.shape[0],tmp.shape[1],len(files)))
                        
    comb_m[:,:,i] = tmp
print("loaded all PEPICO matrices           "+str(time.ctime()))

# load reference
ref = np.loadtxt(path_ref+file_ref,skiprows=25)

#%% plot spectrum integrated over all masses for each KE and compare to PES


xx = np.zeros((comb_m.shape[0],len(central_en)))

for i in range(len(files)):
    xx[:,i] = np.arange(comb_m.shape[0]/2,-comb_m.shape[0]/2,-1)*0.085+central_en[i]-DT_shift-wf

fig,ax = plt.subplots()
ax.plot(ref[:,0]-(DT_shift*DT_on)-wf,ref[:,1],'k')
ax.legend(['reference PES'],loc = 'upper left')
ax.set_xlabel('Kinetic Energy (eV)')
ax.set_ylabel('PES counts')
ax2 = ax.twinx()
for i in range(len(files)):
    ax2.plot(xx[:,i],np.sum(comb_m [:,:,i],axis=1),'.-')
ax2.legend(['KE0','KE1','etc.'],loc='upper center')
ax2.set_ylabel('PEPICO counts')

    
#%% add all PEPICO matrices together to obtain final result

# energy points to cut on the high energy side for each KE (MCP prob)
cut = 2

# maximum intensity in TOF matrix
z_max = 5

# make new energy axis with spacing en_step in eV
en_step = 0.1
en_new = np.arange(np.floor(np.min(xx)),np.ceil(np.max(xx)),en_step)

# TOF axis
tof = np.arange(0,float(MaxTof),float(tof_bin)/1000)

interp_m = np.zeros((len(en_new),comb_m.shape[1],len(files)))


for i in range(len(files)):
    # original energies
    en_old = np.flip(xx[:,i])
    for j in range(comb_m.shape[1]):
        msPES = np.flip(comb_m[:,j,i])
        # make interpolation
        interp_m[:,j,i] = np.interp(en_new,en_old,msPES)
    # set values outside of original energy range to Nan
    for k in range(len(en_new)):
        if (en_new[k]<min(en_old) or en_new[k]>en_old[-1-cut]):
            interp_m[k,:,i] = np.nan
            
# take mean of all PEPICO matrices to obtain complete matrix            
res = np.nanmean(interp_m,axis=2)
res[np.isnan(res)] = 0
# plot complete matrix 
x,y = np.meshgrid(en_new,tof)
cm = 1/2.54

fig = plt.figure(figsize=(15*cm, 15*cm))
gs = GridSpec(4,4,figure=fig)

ax = fig.add_subplot(gs[:3,:])
ax2 = fig.add_subplot(gs[3,:],sharex = ax)

ax.pcolor(x,y,res.T,shading='auto',cmap='jet',vmin=0,vmax=z_max)
ax.set_xlim([min(en_new)-2,max(en_new)+2])
ax.set_ylim([0,float(MaxTof)])
plt.setp(ax.get_xticklabels(), visible=False)
ax.set_ylabel('TOF ($\mu$s)')
ax.set_title(ide+' full PEPICO matrix')


ax2.plot(ref[:,0]-(DT_on*DT_shift)-wf,ref[:,1]*0.5,'k--')
ax2.set_xlabel('Kinetic Energy (eV)')
ax2.set_ylabel('PES counts')
ax2.legend(['PES'],loc='upper left')
ax3 = ax2.twinx()

ax3.plot(en_new,np.sum(res,axis=1),'r')
ax3.set_ylabel('PEPICO counts')
ax3.legend(['PEPICO'],loc='upper right')
plt.savefig(PATH+ide+'_PEPICO_matrix_full.png',format='png')

#%% save final result

np.savetxt(PATH+ide+'_PEPICO_matrix_full.txt', res.T, delimiter='\t')
np.savetxt(PATH+ide+'_en_full.txt',en_new,delimiter='\t')
np.savetxt(PATH+ide+'_tof_full.txt',tof,delimiter='\t')

