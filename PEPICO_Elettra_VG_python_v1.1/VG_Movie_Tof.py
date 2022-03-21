# -*- coding: utf-8 -*-
"""
Created on Wed Jan  5 18:05:01 2022

@author: holzme33
"""
import numpy as np


class VG_Movie_Tof:
    
    def __init__(self):
        return
        
    def ROI(self,xmin,xmax,x_step,ymin,ymax,tof_ns,tof_bin):
        self.xmin = xmin
        self.xmax = xmax
        self.x_step = x_step
        self.ymin = ymin
        self.ymax = ymax
        self.tof_ns = tof_ns
        self.tof_bin = tof_bin
    
    def show_singlesMXY(self,channel_stops,channel_order,channel_X,channel_Y,x_min,x_max):
        # x limits change in loop, y limits are always constant
        Ts = np.zeros(len(channel_stops))
        for i in range(len(channel_stops)):
            if (channel_order[i]==1 and channel_X[i]<x_max and channel_X[i]>x_min and channel_Y[i]<self.ymax and channel_Y[i]>self.ymin):
                Ts[i] = channel_stops[i]
        Ts = Ts[Ts!=0]
        hist = np.histogram(Ts,bins=round(self.tof_ns/self.tof_bin),range=(0,self.tof_ns*1e3))
        
        return Ts,hist
        
    def Tof_loop(self,All_TX_4,All_TS_4_Order,All_X_4,All_Y_4,All_TX_9,All_TS_9_Order,All_X_9,All_Y_9,Second_STof_Hist,scaling,fac):
        
        LTOF3 = Second_STof_Hist[0]     # TOF for all random ions
        
        x_nst = round((self.xmax-self.xmin)/self.x_step)
        Starts = np.zeros((x_nst,))
        LNN = np.zeros((x_nst,len(LTOF3)))
        
        for i in range(x_nst):
            
            x_min = self.xmin+i*self.x_step
            x_max = x_min+self.x_step
            
            Ts4,hist4 = self.show_singlesMXY(All_TX_4,All_TS_4_Order,All_X_4,All_Y_4,x_min,x_max)
            Ts9,hist9 = self.show_singlesMXY(All_TX_9,All_TS_9_Order,All_X_9,All_Y_9,x_min,x_max)
            
            Starts[i] = sum(hist4[0])   # all starts of the pulser for the region between x_min and x_max (by anode or random trigger)
            
            #******************************************************************
            # arbitrary factor (around 0.9) that nobody knows where it is coming from!!!
            vv = fac*Starts[i]/scaling
            #******************************************************************
            
            LNN[i,:] = hist9[0]-vv*LTOF3        # TOF9 - (fac * Starts * randomTOF/total_randoms)
        
        return Starts, LNN
        
        
        