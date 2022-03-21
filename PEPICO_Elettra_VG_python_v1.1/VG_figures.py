# -*- coding: utf-8 -*-
"""
Created on Sat Dec 18 19:51:04 2021

@author: Fabian
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.signal import find_peaks

class VG_figures:
    # produces figures for PEPICO experiment with VG electron analyzer
    
    def fig_singles(self,ch,tof_ns,Bin,All_TX_n,All_TS_n_Order):
        points = tof_ns/(Bin+1)
        Ts = All_TX_n[All_TS_n_Order==1]    # select only events with Order = 1
        hist = np.histogram(Ts,bins=round(points),range=(0,tof_ns*1e3)) #returns a tuple with the histogram (0) and the x limits of the bins (1)
        
        fig,ax = plt.subplots()
        ax.plot(hist[1][:-1]/1e6,hist[0],'r-')
        ax.set_xlim([0,tof_ns/1e3])
        ax.set_title('Channel '+str(ch)+' - Singles') 
        ax.set_xlabel('TOF ($\mu$s)')
        ax.set_ylabel('counts')
        
    def fig_ele_image(self,x9,y9):
        
        hist = np.histogram2d(x9,y9,bins=(100,100))
        
        fig = plt.figure()
        gs = GridSpec(3,3,figure=fig)
        
        ax = fig.add_subplot(gs[:2,:])
        ax2 = fig.add_subplot(gs[2,:],sharex=ax)
        
        # histogram
        ax.pcolor(hist[2],hist[1],hist[0].T,shading='auto')
        ax.set_ylabel('y (px)')
        ax.set_title('electron image')
        
        # electron spectrum as a function of electron energy in px
        ax2.plot(hist[2][:-1],np.sum(hist[0],axis=1),'r')
        ax2.set_xlabel('x (px)')
        ax2.set_ylabel('PES sginal')
        
    def ion_hist(self,All_Tof2,All_Tof8,MaxTof):
        # histograms for ions only
        All_Tof2_hist = np.histogram(All_Tof2,bins=round(MaxTof*1000/2),range=(0,MaxTof*1e6))
        All_Tof8_hist = np.histogram(All_Tof8,bins=round(MaxTof*1000/2),range=(0,MaxTof*1e6))
        
        # play with height if peaks are not found (or too many)
        peaks, _ = find_peaks(All_Tof8_hist[0], height = 0.02*max(All_Tof8_hist[0]))

        anode_starts = sum(All_Tof8_hist[0][peaks[0]-5:peaks[0]+5])
        pulser_starts = sum(All_Tof8_hist[0][peaks[1]-5:peaks[1]+5])
        all_starts = sum(All_Tof8_hist[0])
        print("anode starts: ", anode_starts) 
        print("pulser starts: ", pulser_starts)
        print("all starts: ", all_starts)
        if abs(all_starts-(anode_starts+pulser_starts)) > 2:
            print("check peak picking!")

        fig,(ax1,ax2) = plt.subplots(2,1)
        # channel 2 (all ions)
        ax1.plot(All_Tof2_hist[1][:-1]/1e6,All_Tof2_hist[0],'k-')
        ax1.legend(['All_Tof2'])
        ax1.set_title("Eventually play with 'height' in VG_figures.ion_hist line 53")
        ax2.plot(All_Tof8_hist[1][:-1]/1e6,All_Tof8_hist[0],'k-')
        ax2.plot(All_Tof8_hist[1][peaks]/1e6,All_Tof8_hist[0][peaks],"rx")
        ax2.legend(['All_Tof8'])
        ax2.set_xlabel('TOF ($\mu$s)')
        fig.suptitle("Check that peaks in lower subplot are picked correctly!", fontweight = 'bold')
        
        # return integration limits for future use
        int_limits = [All_Tof8_hist[1][peaks[0]-5],All_Tof8_hist[1][peaks[0]+5],All_Tof8_hist[1][peaks[1]-5],All_Tof8_hist[1][peaks[1]+5]]
        
        return int_limits, pulser_starts, anode_starts
    
    def STof(self,first_tof,second_tof,bg):
        
        #% background range in channels
        bg_ch = [np.min(np.where(first_tof[1]>bg[0]*1e6)), np.min(np.where(first_tof[1]>bg[1]*1e6))]
        
        fig,(ax1,ax2) = plt.subplots(2,1)
        ax1.plot(first_tof[1][:-1]/1e6,first_tof[0],'k-')
        ax1.plot(first_tof[1][bg_ch[0]:bg_ch[1]]/1e6,first_tof[0][bg_ch[0]:bg_ch[1]],'r-')
        ax1.legend(['anode'])
        ax1.set_title("TOF from anode and pulser")
        ax2.plot(second_tof[1][:-1]/1e6,second_tof[0],'k-')
        ax2.plot(second_tof[1][bg_ch[0]:bg_ch[1]]/1e6,second_tof[0][bg_ch[0]:bg_ch[1]],'r-')
        ax2.legend(['pulser'])
        ax2.set_xlabel('TOF ($\mu$s)')
        
        return bg_ch
        
    def PEPICO_matrix(self,res,MaxTof,xmin,xmax,x_step,path,file):
        
        TOF = np.arange(0,MaxTof,2e-3)
        px_y = np.arange(xmin,xmax,x_step)
        
        x,y = np.meshgrid(px_y,TOF)
        
        cm = 1/2.54
        
        fig = plt.figure(figsize=(15*cm, 15*cm))
        gs = GridSpec(4,4,figure=fig)
        
        ax = fig.add_subplot(gs[:3,:])
        ax2 = fig.add_subplot(gs[3,:],sharex = ax)

        p = ax.pcolor(x,y,res.T,shading='auto',cmap='jet')
        ax.set_title(file)
        ax.set_ylabel('TOF ($\mu$s)')
        # plt.colorbar(p)
        ax2.plot(px_y,np.sum(res,axis=1))
        ax2.set_xlabel('x (px)')
        
        plt.setp(ax.get_xticklabels(), visible=False)
        plt.savefig(path+file+'.png',format='png')
        
    def tof_PEPICO(self,res,MaxTof,bin_size):
        
        tof = np.arange(0,MaxTof,bin_size/1000)
        
        fig,ax = plt.subplots()
        ax.plot(tof,np.sum(res,axis=0),'k')
        ax.set_xlabel('TOF ($\mu$s)')
        ax.set_ylabel('PEPICO ion counts')
        ax.set_title("Decrease factor in 'main.py' line 231 if negative peaks!")

        