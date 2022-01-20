# -*- coding: utf-8 -*-
"""
Created on Sun Dec 12 08:57:31 2021

@author: holzme33
"""

import numpy as np
import time

class VG_tof:
    # electron signal
    
    def __init__(self):
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
        
    def find_t_xy(self):
        # goes through the stop signals and checks which four ones belong together, i.e. their times add up correctly
        how_many = max(len(self.Gstop1),len(self.Gstop3),len(self.Gstop5),len(self.Gstop6))

        self.TX = np.zeros((how_many,))      # hit time X
        self.X1 = np.zeros((how_many,))      # t(X1)
        self.X2 = np.zeros((how_many,))      # t(X2)
        self.TY = np.zeros((how_many,))      # hit time Y
        self.Y1 = np.zeros((how_many,))      # t(Y1)
        self.Y2 = np.zeros((how_many,))      # t(Y2)
        self.E_id = np.zeros((how_many,))    # number of event

        i1,i3,i5,i6,ii = 0,0,0,0,0

        # find the four delay line signals that belong together
        while (i1<len(self.Gstop1) and i3<len(self.Gstop3) and i5<len(self.Gstop5) and i6<len(self.Gstop6)):
            maybe, found = 0,0                      # boolean variables if electrons belong together
            time1 = self.Gstop1[i1]+150             # corrections for different lengths of delay lines
            time3 = self.Gstop3[i3]-150
            time5 = self.Gstop5[i5]-5817
            time6 = self.Gstop6[i6]-5117
            DtimeX = abs(time3-time1)
            DtimeY = abs(time5-time6)
            
            
            if (DtimeX<=28600 and DtimeY<=31200):   # 28600 and 31200 are presumably the cable lengths
                maybe = 1
            
            if maybe:
                TimeX = (time1+time3-28600)/2           # averaged arrival time for x and y
                TimeY = (time5+time6-31200)/2
                if abs(TimeX-TimeY)<3000:               # check if x and y could come from the same electron
                    self.TX[ii] = TimeX
                    self.X1[ii] = -self.TX[ii]+time1
                    self.X2[ii] = -self.TX[ii]+time3
                    self.TY[ii] = TimeY
                    self.Y1[ii] = -self.TY[ii]+time5
                    self.Y2[ii] = -self.TY[ii]+time6
                    self.E_id[ii] = ii+1
                    ii+=1
                    found=1
            
            mintime = min(time1,time3,time5,time6)
            
            if found:
                i1+=1; i3+=1; i5+=1; i6+=1
            else:
                if mintime == time1: i1+=1 
                if mintime == time3: i3+=1 
                if mintime == time5: i5+=1 
                if mintime == time6: i6+=1 


        # delete all zeros at the end
        self.TX = self.TX[self.TX != 0]; self.X1 = self.X1[self.X1 != 0]; self.X2 = self.X2[self.X2 != 0]; self.TY = self.TY[self.TY != 0]; self.Y1 = self.Y1[self.Y1 != 0]; self.Y2 = self.Y2[self.Y2 != 0]; self.E_id = self.E_id[self.E_id != 0]
        
        return self.TX, self.X1, self.X2, self.TY, self.Y1, self.Y2, self.E_id
    
    def find_Ch(self,which,tof):
        LGS = eval('self.Gstop'+str(which))
        max_CH = 5*len(LGS)
        tof=1e6*abs(tof)    # tof range in picoseconds (tof is in us)
        LTXW, LTYW, LTX1W, LTY1W, LTSW, LEidW = np.zeros((max_CH)), np.zeros((max_CH)), np.zeros((max_CH)), np.zeros((max_CH)), np.zeros((max_CH)), np.zeros((max_CH)),
        kk, firstK, jj = 0,0,0
        for i in range(len(self.TX)):
            start_x = self.TX[i]        # electron i is the start
            start_id = self.E_id[i]
            kk = firstK                 # kk is the event to be evaluated
            stop_x = LGS[kk]
            
            while (kk<len(LGS) and (start_x+tof)>stop_x):   # the event must occur within the defined max. tof
                stop_x = LGS[kk]
                tof_dif = stop_x-start_x
                if (tof_dif<tof and tof_dif>0):     # the ion must be detected after the electron
                    LTXW[jj]=tof_dif                # ion time of flight
                    LTX1W[jj] = self.X1[i]          # x of the electron
                    LTY1W[jj] = self.Y1[i]          # y of the electron
                    LTSW[jj] = start_x              # start
                    LEidW[jj] = start_id            # start identifier
                    jj+=1
                    firstK = kk
                kk+=1
            
        LTXW = LTXW[LTXW != 0]; LTYW = LTYW[LTYW != 0]; LTX1W = LTX1W[LTX1W != 0]; LTY1W = LTY1W[LTY1W != 0]; LTSW = LTSW[LTSW != 0]; LEidW = LEidW[LEidW != 0]
        return LTXW, LTYW, LTX1W, LTY1W, LTSW, LEidW

    # true coincidences are detected both on ch 4 and ch 2
    def Match_42(self):
        ii, jj = 0,0
        LTX_2_4 = np.zeros((len(self.TX_2)))
        while (ii<len(self.E_id_2) and jj<len(self.E_id_4)):
            # check if the identifiers of canal 2 and 4 match
            if self.E_id_4[jj]<self.E_id_2[ii]:
                jj+=1
            else:
                if self.E_id_4[jj]>self.E_id_2[ii]:
                    ii+=1
                else:
                    LTX_2_4[ii] = self.TX_2[ii]-self.TX_4[jj]
                    ii+=1
        return LTX_2_4 # tof difference between canal 2 and 4
    
    def find_channels(self,MaxTof):
        # find the events of the "ion" channels
        if len(self.Gstop8)>0:
            TX_8, TY_8, TX1_8, TY1_8, TS_8, E_id_8 = self.find_Ch(8,MaxTof)
        if len(self.Gstop7)>0:
            TX_7, TY_7, TX1_7, TY1_7, TS_7, E_id_7 = self.find_Ch(7,MaxTof)
        if len(self.Gstop4)>0:
            self.TX_4, TY_4, TX1_4, TY1_4, TS_4, self.E_id_4 = self.find_Ch(4,MaxTof)
        if len(self.Gstop2)>0:
            self.TX_2, TY_2, TX1_2, TY1_2, TS_2, self.E_id_2 = self.find_Ch(2,MaxTof)

        if (len(self.Gstop4)>0 and len(self.Gstop2)>0):
            TX_2_4 = self.Match_42()
        
        return TX_8, TX1_8, TY1_8, TS_8, E_id_8, TX_7, TX1_7, TY1_7, TS_7, E_id_7, self.TX_4, TX1_4, TY1_4, TS_4, self.E_id_4, self.TX_2, TX1_2, TY1_2, TS_2, self.E_id_2, TX_2_4
    

    
    def find_order(self,ch,T0):
        # takes the channel number and the start times in this channel as arguments
        # returns how many ions were detected in coincidence in this channel
        ii = 0
        T_Order = np.full((len(T0)),1)
        while (ii<len(T0)-2): ## FH: was: (ii<len(T0)-1)
            start_1 = T0[ii]
            start_2 = T0[ii+1]
            while (start_1>start_2):    # if start_2 is before start_1, one can move on
                start_1 = T0[ii]
                start_2 = T0[ii+1]
                if start_1>start_2:     
                    ii+=1
            
            if start_2>start_1:         # if start_2 is after start_1, the order is at least 1
                T_Order[ii]=1
                ii+=1
            else:
                start_3 = T0[ii+2]
                if start_3 > start_2:                   # 2 events at once
                    T_Order[ii], T_Order[ii+1] = 2,2
                    ii+=2
                else:
                    if ii+3<len(T0):        # added by FH to make sure that the loop does not run out of range
                        start_4 = T0[ii+3]
                    if start_4>start_3:                 # 3 events at once
                        T_Order[ii], T_Order[ii+1], T_Order[ii+2] = 3,3,3
                        ii+=3
                    else:
                        kk = ii+4                                       # everything with order >3 is set to 999
                        while (start_4 == start_1 and kk<len(T0)):
                            start_4 = T0[kk]
                            kk+=1
                        j = ii
                        while (j<kk and j<len(T0)):
                            T_Order[j] = 999
                            j+=1
                        ii=kk-1
        return T_Order
    
    def GetXY(self,nbins,x1,y1):
        XXY, YXY = np.zeros(len(x1)),np.zeros(len(y1))
        for i in range(len(x1)):
            XXY[i] = (x1[i]/28600)*nbins
            YXY[i] = (y1[i]/31200)*nbins   
        return XXY, YXY
    
    def clean_up1(self,All_E_id,All_E_id_2,All_E_id_4,All_E_id_7,All_E_id_8,All_TS_2_Order,All_TS_4_Order,All_TS_8_Order,All_TX_2_4):
        To_rem_8 = All_E_id_8[All_TS_8_Order !=1]           # more than one stop on ch8 (unknown wheter event triggered by electron or pulser)
        To_rem_2 = All_E_id_2[All_TS_2_Order>2]             # more than two ions (improbable event for inner-valence ionization)
        To_rem_4 = All_E_id_4[All_TS_4_Order !=1]           # extraction field started more than once (unknown whether extracted started by anode or pulser)
        To_rem_2_4 = All_E_id_2[All_TX_2_4 == 0]            # event with an ion, but no trigger
        
        To_rem = np.concatenate((To_rem_8,To_rem_2,To_rem_4,To_rem_2_4))
        To_rem = np.sort(To_rem)
        
        TR = np.empty((len(To_rem),))
        TR[:] = np.nan
        
        vv,jj = 0,0
        
        for i in range(len(To_rem)):
            if To_rem[i] > vv:
                TR[jj] = To_rem[i]
                vv = TR[jj]
                jj+=1
        TR = TR[~np.isnan(TR)]
        TR = TR.astype(int)
        
        hi_e_id = max(max(All_E_id),max(All_E_id_8),max(All_E_id_2),max(All_E_id_4),max(All_E_id_7))
        
        # make list with events to be removed (nan) and to be kept (1)
        Ids = np.full((int(hi_e_id)+1,),1.0)
        for i in range(len(TR)):
            Ids[TR[i]] = np.nan
        Ids[0] = np.nan
        
        return Ids
    
    def clean_up2(self,Ids,All_E_Id_n,All_TS_n,All_TS_n_Order,All_TX_n,All_TX1_n,All_TY1_n):
        # use list of Ids to throw out 'bad' events
        All_E_id_n = All_E_Id_n.astype(int)                     # All_E_Id is needed as integer and float array (note capital I in variable name)
        All_E_Id_n = All_E_Id_n.astype(float)
        All_TS_n_Order = All_TS_n_Order.astype(float)
        for i in range(len(All_E_id_n)):
            All_TS_n[i] = All_TS_n[i]*Ids[All_E_id_n[i]]
            All_TS_n_Order[i] = All_TS_n_Order[i]*Ids[All_E_id_n[i]]
            All_TX_n[i] = All_TX_n[i]*Ids[All_E_id_n[i]]
            All_TX1_n[i] = All_TX1_n[i]*Ids[All_E_id_n[i]]
            All_TY1_n[i] = All_TY1_n[i]*Ids[All_E_id_n[i]]
            All_E_Id_n[i] = All_E_Id_n[i]*Ids[All_E_id_n[i]]
            
        All_E_Id_n = All_E_Id_n[~np.isnan(All_E_Id_n)]; All_TS_n = All_TS_n[~np.isnan(All_TS_n)]; All_TS_n_Order = All_TS_n_Order[~np.isnan(All_TS_n_Order)]
        All_TX_n = All_TX_n[~np.isnan(All_TX_n)]; All_TX1_n = All_TX1_n[~np.isnan(All_TX1_n)]; All_TY1_n = All_TY1_n[~np.isnan(All_TY1_n)]
        
        return All_E_Id_n,All_TS_n,All_TS_n_Order,All_TX_n,All_TX1_n,All_TY1_n
            
    def clean_up3(self,Ids,All_E_Id_2,All_TX_2_4):
        All_E_id_2 = All_E_Id_2.astype(int)
        All_E_Id_2 = All_E_Id_2.astype(float)
        for i in range(len(All_E_id_2)):
            All_TX_2_4[i] = All_TX_2_4[i]*Ids[All_E_id_2[i]]
        
        All_TX_2_4 = All_TX_2_4[~np.isnan(All_TX_2_4)]
        
        return All_TX_2_4
    
                        