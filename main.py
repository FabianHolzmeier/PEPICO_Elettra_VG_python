# -*- coding: utf-8 -*-
"""
Created on Sun Dec 12 08:56:51 2021

@author: holzme33
"""
import glob, os
import os.path
from os import path
import time
import numpy as np
import inspect
import scipy.io as sio
from VG_tof import VG_tof
from VG_Tof_only import VG_Tof_only
from VG_figures import VG_figures
from VG_Movie_Tof import VG_Movie_Tof

start_t = time.time()

VG_tof = VG_tof()
VG_Tof_only = VG_Tof_only()
VG_figures = VG_figures()
VG_Movie_Tof = VG_Movie_Tof()

# directory with the .mat files to be analyzed
PATH = "B:\\PEPICO\\rawdata\\BSn_157\\KE9\\"
group = "KE9"

MaxTof = 9                                      # TOF range in microseconds
show_singles = 1                                # show histograms of singles                
show_doubles = 0                                # show histograms of doubles (for future release)

# results_folder
path_res = "B:\\analysis\\BSn\\"

#******************************************************************************
# Do not change anything below this line (unless you know what you're doing)

files = []
# get the filenames of the raw data saved as .mat
for file in os.listdir(PATH):
    if file.endswith(".mat"):
        files.append(file)

# all events and identifiers on the delay lines
All_X1, All_Y1, All_E_id = [],[],[]
# all ion signals
All_E_id_2, All_TS_2, All_TX_2, All_TX1_2, All_TY1_2 = [],[],[],[],[]
# SRS trigger = real start for the TOF
All_E_id_4, All_TS_4, All_TX_4, All_TX1_4, All_TY1_4 = [],[],[],[],[]
# 'true' triggers, similar to 4 except for pulses coming during the TOF (ringing or real)
All_E_id_7, All_TS_7, All_TX_7, All_TX1_7, All_TY1_7 = [],[],[],[],[]
# 'true' with random delayed; coincidence indicates true&random = source unknown
All_E_id_8, All_TS_8, All_TX_8, All_TX1_8, All_TY1_8 = [],[],[],[],[]
# ion times on 2 with trigger on 4
All_TX_2_4 = []

# ion identifiers
valid_start = 1
All_Start_ID_2, All_Tof2, All_Start_ID_8, All_Tof8 = [],[],[],[]


# loop over all .mat files in the directory
for i in range(len(files)):
# for i in range(1):
    
    # load .mat file
    dat = sio.loadmat(PATH+files[i])
    
    print('loaded file: '+files[i]+'     '+str(time.ctime()))
    
    channel = dat['channel'].flatten()              # list of channels in which the event was detected
    start = dat['start'].flatten()                  # start of clock (1 Hz)
    stat_times1 = dat['stat_times'][:,1]            # start of trigger
    stoptime = dat['stoptime'].flatten()            # stop times
    
    # read data from file and sort it by channel 
    VG_tof.prepare(channel,stat_times1,start,stoptime)
    
    
    #**************************************************************************
    #                        ELECTRONS
    #**************************************************************************
    
    # find the four electron signals on the delay line anode that belong together
    # TX, TY: arrival time on X and Y delay line average of X1,X2 and Y1,Y2, respectively
    # E_id: running number of event
    TX,X1,X2,TY,Y1,Y2,E_id = VG_tof.find_t_xy()     
    
    # find ions that are detected in the various channels that could belong to the electron
    # TS is the start time of the electron with the identifier E_id
    # TX_2_4 is the difference of ions detected in ch2 and ch4 with the same identifier (electron start), i.e. real TOF
    TX_8,TX1_8,TY1_8,TS_8,E_id_8,TX_7,TX1_7,TY1_7,TS_7,E_id_7,TX_4,TX1_4,TY1_4,TS_4,E_id_4,TX_2,TX1_2,TY1_2,TS_2,E_id_2,TX_2_4 = VG_tof.find_channels(MaxTof)
    
    E_id+=len(All_E_id)                      # increase event id by number of events from previous files
    
    # concatenate different files at the same energy
    All_X1 = np.concatenate((All_X1,X1))                        # all signals on X delay line
    All_Y1 = np.concatenate((All_Y1,Y1))                        # all signals on Y delay line
    All_E_id = np.concatenate((All_E_id,E_id))                  # all event identifiers
    
    All_TX_8 = np.concatenate((All_TX_8, TX_8))                 # all ion stop times
    All_TX_7 = np.concatenate((All_TX_7,TX_7))
    All_TX_4 = np.concatenate((All_TX_4,TX_4))
    All_TX_2 = np.concatenate((All_TX_2,TX_2))
    All_TX_2_4 = np.concatenate((All_TX_2_4,TX_2_4))            # all real ion time of flights
    
    All_TX1_8 = np.concatenate((All_TX1_8,TX1_8))               # all coordinates x and y of the electrons
    All_TY1_8 = np.concatenate((All_TY1_8,TY1_8))
    All_TX1_7 = np.concatenate((All_TX1_7,TX1_7))
    All_TY1_7 = np.concatenate((All_TY1_7,TY1_7))
    All_TX1_4 = np.concatenate((All_TX1_4,TX1_4))
    All_TY1_4 = np.concatenate((All_TY1_4,TY1_4))
    All_TX1_2 = np.concatenate((All_TX1_2,TX1_2))
    All_TY1_2 = np.concatenate((All_TY1_2,TY1_2))
    
    All_TS_8 = np.concatenate((All_TS_8,TS_8))                  # all starts triggered by electron
    All_TS_7 = np.concatenate((All_TS_7,TS_7))
    All_TS_4 = np.concatenate((All_TS_4,TS_4))
    All_TS_2 = np.concatenate((All_TS_2,TS_2))
    
    E_id_8+=len(All_E_id_8)
    E_id_7+=len(All_E_id_7)
    E_id_4+=len(All_E_id_4)
    E_id_2+=len(All_E_id_2)
    
    All_E_id_8 = np.concatenate((All_E_id_8,E_id_8))            # all event identifiers
    All_E_id_7 = np.concatenate((All_E_id_7,E_id_7))
    All_E_id_4 = np.concatenate((All_E_id_4,E_id_4))
    All_E_id_2 = np.concatenate((All_E_id_2,E_id_2))
    
    print('processed electrons from file: '+files[i]+'     '+str(time.ctime()))
    
    
    # #**************************************************************************
    # #                           IONS
    # #**************************************************************************
    
    VG_Tof_only.prepare(channel,stat_times1,start,stoptime)     # (maybe) to be optimized to not call twice the Gstops
    # get ion events of single files
    Start_ID_2,Tof2,Start_ID_8,Tof8 = VG_Tof_only.Tof_only(MaxTof,valid_start)
    
    # concatenate files
    All_Start_ID_2 = np.concatenate((All_Start_ID_2,Start_ID_2))
    All_Tof2 = np.concatenate((All_Tof2,Tof2))
    All_Start_ID_8 = np.concatenate((All_Start_ID_8,Start_ID_8))
    All_Tof8 = np.concatenate((All_Tof8,Tof8))
    
    valid_start = max(max(All_Start_ID_2),max(All_Start_ID_8))
    
    print('processed ions from file: '+files[i]+'          '+str(time.ctime()))
    

# find orders of electron events (how many stops)
for i in [2,4,7,8]:
    exec('tmp = All_TS_'+str(i))
    if len(tmp)>0:
        T_Order = VG_tof.find_order(i,tmp)
        exec('All_TS_'+str(i)+'_Order = T_Order')
    
# convert x and y times to pix
bins = 100
for i in [2,4,7,8]:    
    exec('All_X_'+str(i)+', All_Y_'+str(i)+' = VG_tof.GetXY(bins,All_TX1_'+str(i)+', All_TY1_'+str(i)+')')
    
#%% clean up electrons

# remove all events with 2 or more stops on ch 8 (HV pulse source unwknown)
# remove all events with more than 2 ions detected on ch 2 (improbable)
# remove all events with more than 2 stops on ch 4 (source unknown)
# remove all events with an ion (2) without a corresponding trigger on 4 (still have an infinite time)

# list of all events to remove
Ids = VG_tof.clean_up1(All_E_id,All_E_id_2,All_E_id_4,All_E_id_7,All_E_id_8,All_TS_2_Order,All_TS_4_Order,All_TS_8_Order,All_TX_2_4)
# remove the events above
All_TX_2_4 = VG_tof.clean_up3(Ids,All_E_id_2,All_TX_2_4)
for i in [2,4,7,8]:
    exec('All_E_id_'+str(i)+', All_TS_'+str(i)+', All_TS_'+str(i)+'_Order, All_TX_'+str(i)+', All_TX1_'+str(i)+', All_TY1_'+str(i)+
         ' = VG_tof.clean_up2(Ids,All_E_id_'+str(i)+',All_TS_'+str(i)+', All_TS_'+str(i)+'_Order, All_TX_'+str(i)+', All_TX1_'+str(i)+', All_TY1_'+str(i)+')')


# channel 9 = start-corrected 2_4
All_E_id_9,All_TS_9,All_TS_9_Order,All_TX_9,All_TX1_9,All_TY1_9,All_TX_9 = All_E_id_2,All_TS_2,All_TS_2_Order,All_TX_2_4,All_TX1_2,All_TY1_2, All_TX_2_4

# Get XYs also for channel 9
bins = 100
for i in [2,4,7,8,9]:    
    exec('All_X_'+str(i)+', All_Y_'+str(i)+' = VG_tof.GetXY(bins,All_TX1_'+str(i)+', All_TY1_'+str(i)+')')

# histograms and electron MCP image

if show_singles:
    # histograms
    for i in [2,4,7,8,9]:
        hist = VG_figures.fig_singles(i,MaxTof*1000,1,eval('All_TX_'+str(i)),eval('All_TS_'+str(i)+'_Order'))
    # electron image    
    VG_figures.fig_ele_image(All_X_9,All_Y_9)
        
# if show_doubles:
    

#%% clean up tof

All_Start_ID_8, All_Tof8, All_Tof2_Start, All_Start_ID_2, All_Tof2 = VG_Tof_only.clean_up_tof(All_Start_ID_2,All_Start_ID_8,All_Tof2,All_Tof8)

# histograms of all ions (Tof 2) and of triggers by anode or pulser (Tof 8) 

int_lim,pulser_starts = VG_figures.ion_hist(All_Tof2,All_Tof8,MaxTof)

# Correlate ions to start from anode or pulser
First_Tof2_Start,First_Tof2,First_Start_ID_2, Second_Tof2_Start,Second_Tof2,Second_Start_ID_2 = VG_Tof_only.first_second(int_lim[0],int_lim[1],int_lim[2],int_lim[3])

# TOF spectra for ions from anode and from pulser

bin_size = 2            # bin size in ns

First_STof_Hist = VG_Tof_only.Fst_single_tof(MaxTof*1000,bin_size)
Second_STof_Hist =  VG_Tof_only.Snd_single_tof(MaxTof*1000,bin_size)
VG_figures.STof(First_STof_Hist,Second_STof_Hist)

#%% coincidences

# select ROI on electron image in px
xmin,xmax = 5,95        # x limits
x_step = 2              # step size in x (electron energy)
ymin,ymax = 20,65       # y limits
VG_Movie_Tof.ROI(xmin,xmax,x_step,ymin,ymax,MaxTof*1000,bin_size)

# factor for subtraction of randoms
fac = 0.65

# test MXY
MXYSingles_ch4,MXYch4_hist = VG_Movie_Tof.show_singlesMXY(All_TX_4,All_TS_4_Order,All_X_4,All_Y_4,xmin,xmax)
MXYSingles_ch9,MXYch9_hist = VG_Movie_Tof.show_singlesMXY(All_TX_9,All_TS_9_Order,All_X_9,All_Y_9,xmin,xmax)

# make the loop, returns the spectrum of Starts and the PEPICO matrix (res)
Starts, res = VG_Movie_Tof.Tof_loop(All_TX_4,All_TS_4_Order,All_X_4,All_Y_4,All_TX_9,All_TS_9_Order,All_X_9,All_Y_9,Second_STof_Hist,pulser_starts,fac)
# normalize the result matrix on the number of files
res = res/len(files)

# save result
path_pepico = path_res+'PEPICO\\'

if path.isdir(path_pepico) == False:
    os.mkdir(path_pepico)

res_file = files[0][:-15]+"_"+group+"_PEPICO_matrix_"+str(MaxTof)+"_"+str(bin_size)+"_"+str(x_step)

np.savetxt(path_pepico+res_file+".dat", res, delimiter='\t')

# figure of PEPICO matrix
VG_figures.PEPICO_matrix(res,MaxTof,xmin,xmax,x_step,path_pepico,res_file)

VG_figures.tof_PEPICO(res,MaxTof,bin_size)

run_time = time.time()-start_t
print()
print("analysis finished after {:.1f} min".format(run_time/60))
print(str(time.ctime()))













