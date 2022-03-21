# -*- coding: utf-8 -*-
"""
Created on Wed Dec 15 10:18:49 2021

@author: holzme33
"""
import numpy as np
from VG_tof import VG_tof
import time


class VG_Tof_only:
    
    def __init__(self):
        # vgt = VG_tof()
        return
    
    def prepare(self,channel,stat_times1,start,stoptime):

        # assigns to all stop signals the channel they occured in
        trigger = 1e6*max(stat_times1) 

        self.Gstop = stoptime+trigger*(start-1)             # absolute time for (starting at 0 for each file)
        self.Gstop1 = self.Gstop[np.where(channel==1)]      # x1 electron
        self.Gstop2 = self.Gstop[np.where(channel==2)]      # all ions
        self.Gstop3 = self.Gstop[np.where(channel==3)]      # x2 electron
        self.Gstop4 = self.Gstop[np.where(channel==4)]      # true start of ions (signal on e-anode)
        self.Gstop5 = self.Gstop[np.where(channel==5)]      # y1 electron
        self.Gstop6 = self.Gstop[np.where(channel==6)]      # y2 electron
        self.Gstop7 = self.Gstop[np.where(channel==7)]      # SRS trigger (true or random)
        self.Gstop8 = self.Gstop[np.where(channel==8)]      # true and random starts with delay D1 and D2, respectively
        self.Gstop9 = self.Gstop[np.where(channel==9)]
        
    def find_stops(self,mtof,start_id):         # mtof in in picoseconds
        
        max_ch = 2*max(len(self.Gstop2),len(self.Gstop8))
        self.Tof2,self.Start_ID_2,self.Tof8,self.Start_ID_8 = np.zeros(max_ch),np.zeros(max_ch),np.zeros(max_ch),np.zeros(max_ch)
        self.Tof2[:], self.Tof8[:] = np.inf, np.inf                 # all ions (2) and delayed ion signal (8)
        self.Start_ID_2[:], self.Start_ID_8[:] = -np.inf, -np.inf        
        
        
        kk2,kk8,firstK2,firstK8,jj2,jj8 = 0,0,0,0,0,0
        for i in range(len(self.Gstop4)):       # we look for ions corresponding to each detected electron
            start_x = self.Gstop4[i]
            kk2 = firstK2
            kk8 = firstK8
            stop_x = self.Gstop8[kk8]
            while (kk8<len(self.Gstop8) and (start_x+mtof)>stop_x):
                stop_x = self.Gstop8[kk8]
                tof_dif = stop_x-start_x                #TOF difference corresponds to the delay between the true start and delayed start
                if (tof_dif<mtof and tof_dif>0):
                    self.Tof8[jj8] = tof_dif                #Tof8: helps identifying true from random signal
                    self.Start_ID_8[jj8] = i+start_id
                    
                    jj8+=1
                    firstK8=kk8
                kk8+=1
                
            stop_x = self.Gstop2[kk2]    
            while (kk2<len(self.Gstop2) and (start_x+mtof)>stop_x):
                stop_x = self.Gstop2[kk2]
                tof_dif = stop_x-start_x            #TOF difference corresponds to the delay between the true start and the ion signal
                if (tof_dif<mtof and tof_dif>0):
                    self.Tof2[jj2] = tof_dif                #Tof2: real ion time-of-flight
                    self.Start_ID_2[jj2] = i+start_id
                    jj2+=1
                    firstK2=kk2
                kk2+=1
        
        # get rid of empty events        
        self.Tof2 = self.Tof2[self.Start_ID_2 != -np.inf]; self.Start_ID_2 = self.Start_ID_2[self.Start_ID_2 != -np.inf]
        self.Tof8 = self.Tof8[self.Start_ID_8 != -np.inf]; self.Start_ID_8 = self.Start_ID_8[self.Start_ID_8 != -np.inf]

    
    def Tof_only(self,mtof_us,last_start):
        
        # look for stop signals if channels 2,4,8 have counts (i.e. there are ions which don't come from trigger)
        if len(self.Gstop4)>0 and (len(self.Gstop2)>2 and len(self.Gstop8)>0):
                self.find_stops(mtof_us*1e6,last_start)       
        
        return self.Start_ID_2,self.Tof2,self.Start_ID_8,self.Tof8
    
    def clean_up_tof(self,All_Start_ID_2,All_Start_ID_8,All_Tof2,All_Tof8):
        
        self.All_Start_ID_2 = All_Start_ID_2
        self.All_Tof2 = All_Tof2
        
        # sort all_Start_ID, All_Tof8
        mywave = np.full((len(All_Start_ID_8),1),1)
        
        for i in range(len(All_Start_ID_8)-1):
            if All_Start_ID_8[i+1] == All_Start_ID_8[i]:
                mywave[i+1],mywave[i] = 2,2
        
        All_Start_ID_8 = All_Start_ID_8[mywave[:,0]==1]
        All_Tof8 = All_Tof8[mywave[:,0]==1]
        
        
        mywave2 = np.full((int(max(All_Start_ID_8))+1,1),np.inf)
        for i in range (len(All_Start_ID_8)):
            mywave2[round(All_Start_ID_8[i])] = All_Tof8[i]
            
        Tof2_Start = np.full((len(self.All_Start_ID_2),1),np.inf)
        for i in range(len(self.All_Start_ID_2)):
            if round(self.All_Start_ID_2[i])<len(mywave2):                  # dirty bug fix 13/01/22
                Tof2_Start[i] = mywave2[round(self.All_Start_ID_2[i])]
            
        self.All_Tof2_Start = np.copy(Tof2_Start)
        
        self.All_Tof2_Start = self.All_Tof2_Start[Tof2_Start[:,0] != np.inf]
        self.All_Start_ID_2 = self.All_Start_ID_2[Tof2_Start[:,0] != np.inf]
        self.All_Tof2 = self.All_Tof2[Tof2_Start[:,0] != np.inf]
        
        
        return All_Start_ID_8, All_Tof8, self.All_Tof2_Start, self.All_Start_ID_2, self.All_Tof2
    
    def first_second(self,x1,x2,x3,x4):
        # finds which ion belongs to first trigger signal (anode) and which one to second trigger signal (pulser)
        
        self.First_Tof2_Start, self.First_Tof2, self.First_Start_ID_2 = np.zeros((len(self.All_Tof2_Start),1)),np.zeros((len(self.All_Tof2_Start),1)),np.zeros((len(self.All_Tof2_Start),1))
        self.Second_Tof2_Start, self.Second_Tof2, self.Second_Start_ID_2 = np.zeros((len(self.All_Tof2_Start),1)),np.zeros((len(self.All_Tof2_Start),1)),np.zeros((len(self.All_Tof2_Start),1))
        
        for i in range(len(self.All_Tof2_Start)):
            if (self.All_Tof2_Start[i]>x1 and self.All_Tof2_Start[i]<x2):
                self.First_Tof2_Start[i] = self.All_Tof2_Start[i]
                self.First_Tof2[i] = self.All_Tof2[i]
                self.First_Start_ID_2[i] = self.All_Start_ID_2[i]
            elif (self.All_Tof2_Start[i]>x3 and self.All_Tof2_Start[i]<x4):
                self.Second_Tof2_Start[i] = self.All_Tof2_Start[i]
                self.Second_Tof2[i] = self.All_Tof2[i]
                self.Second_Start_ID_2[i] = self.All_Start_ID_2[i]
                
        self.First_Tof2_Start = self.First_Tof2_Start[self.First_Tof2_Start != 0]
        self.First_Tof2 = self.First_Tof2[self.First_Tof2 != 0]
        self.First_Start_ID_2 = self.First_Start_ID_2[self.First_Start_ID_2 != 0]
        self.Second_Tof2_Start = self.Second_Tof2_Start[self.Second_Tof2_Start != 0]
        self.Second_Tof2 = self.Second_Tof2[self.Second_Tof2 != 0]
        self.Second_Start_ID_2 = self.Second_Start_ID_2[self.Second_Start_ID_2 != 0]
                
        return self.First_Tof2_Start,self.First_Tof2,self.First_Start_ID_2, self.Second_Tof2_Start,self.Second_Tof2,self.Second_Start_ID_2
        
    
    def Fst_single_tof(self,max_tof_ns,bin_size_ns):
        # TOF with one ion only
        
        First_order_single = np.full((len(self.First_Start_ID_2),),1)
        for i in range(len(self.First_Start_ID_2)-1):
            if self.First_Start_ID_2[i+1] == self.First_Start_ID_2[i]:
                First_order_single[i+1],First_order_single[i] = 0,0
        
        
        points = max_tof_ns/bin_size_ns
        F_S_H = np.histogram(self.First_Tof2,bins=round(points),range=(0,max_tof_ns*1e3), weights = First_order_single)
        
        return F_S_H
    
    def Snd_single_tof(self,max_tof_ns,bin_size_ns):
        # TOF with one ion only
        
        Second_order_single = np.full((len(self.Second_Start_ID_2),),1)
        for i in range(len(self.Second_Start_ID_2)-1):
            if self.Second_Start_ID_2[i+1] == self.Second_Start_ID_2[i]:
                Second_order_single[i+1],Second_order_single[i] = 0,0
        
        
        points = max_tof_ns/bin_size_ns
        S_S_H = np.histogram(self.Second_Tof2,bins=round(points),range=(0,max_tof_ns*1e3), weights = Second_order_single)
        
        return S_S_H
        
        