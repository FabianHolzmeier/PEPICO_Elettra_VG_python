# -*- coding: utf-8 -*-
"""
Created on Sat Mar 19 17:06:02 2022

@author: holzme33
"""
	
from PyQt5.QtWidgets import*
from PyQt5.uic import loadUi
from PyQt5.QtCore import Qt
from PyQt5 import QtCore as qtc

from matplotlib.backends.backend_qt5agg import (NavigationToolbar2QT as NavigationToolbar)

import numpy as np
import random
import sys
     
class MainWindow(QMainWindow):
    
    def __init__(self):
        
        QMainWindow.__init__(self)
        loadUi("msPES_GUI.ui",self)
        self.setWindowTitle("mass-selected PES GUI")


        # buttons
        self.file_btn.clicked.connect(self.dialog)
        self.tof_ms_btn.clicked.connect(self.tof_ms)
        self.save_msPES_btn.clicked.connect(self.save_msPES)
        self.save_esTOF_btn.clicked.connect(self.save_esTOF)
        
        # status variables
        self.matrix_loaded = 0
        self.show_tof = 1       # if true TOF is plotted, else MASS
        self.mass = None        # initialize mass array
        
        # figure toolbars
        self.addToolBar(NavigationToolbar(self.matrix_fig.canvas,self))
        
        # poller
        self.poller = qtc.QTimer()
        self.poller.timeout.connect(self.updater)
        self.poller.start(200)
    
        
    def dialog(self):
        # opens file dialog to load PEPICO matrix and calls the plot functions
        self.file, check = QFileDialog.getOpenFileName(None, 
                                                  "select PEPICO matrix","","Text Files (*.txt);;All Files (*)"
                                                  )
        self.file = self.file.replace("/","\\\\")
        if check:
            self.file_lbl.setText(self.file)
            self.matrix = np.loadtxt(self.file)
            en_file = self.file.replace("PEPICO_matrix","en")
            tof_file = self.file.replace("PEPICO_matrix","tof")
            self.en = np.loadtxt(en_file)
            self.tof = np.loadtxt(tof_file)
            self.plot_matrix()
            
            self.PES = np.sum(self.matrix,axis=0)
            self.msPES = self.PES
            self.TOF = np.sum(self.matrix,axis=1)
            self.esTOF = self.TOF
            # self.plot_PES()
            # self.plot_TOF()
            
    def tof_ms(self):
        # do mass calibration with input parameters
        m1 = float(self.m1_ent.text())
        m2 = float(self.m2_ent.text())
        t1 = float(self.t1_ent.text())
        t2 = float(self.t2_ent.text())
        self.mass = self.calibrate(m1,m2,t1,t2)
        
        # change button text
        if self.show_tof:
            self.show_tof=0
            self.tof_ms_btn.setText("Plot TOF")
        else:
            self.show_tof=1
            self.tof_ms_btn.setText("Plot MS")
        # replot matrix and TOF_MS    
        self.plot_matrix()
        self.plot_TOF()
            
            
    def calibrate(self,m1,m2,t1,t2):
        # mass calibration of TOF spectra
        # m = alpha*t**2+beta (mass-to-charge is proportional to the square of the TOF)
        alpha = (m1-m2)/(t1**2-t2**2)
        beta = m1-alpha*t1**2
        return alpha*self.tof**2+beta
    
    
    def plot_matrix(self):
        # plot the PEPICO matrix
        
        if self.show_tof:
            x,y = np.meshgrid(self.en,self.tof)
        else:
            x,y = np.meshgrid(self.en,self.mass)
        
        self.matrix_fig.canvas.axes.clear()
        self.matrix_fig.canvas.axes.pcolormesh(x,y,self.matrix,shading='gouraud',cmap='gist_earth_r')
        self.matrix_fig.canvas.axes.set_xlabel('Kinetic Energy (eV)',fontsize=14)
        
        # plot on TOF or mass scale
        if self.show_tof:
            self.matrix_fig.canvas.axes.set_ylabel('TOF (us)',fontsize=14)
            self.matrix_fig.canvas.axes.set_ylim([0,9])           
        else:
            self.matrix_fig.canvas.axes.set_ylabel('m/z (amu)',fontsize=14)
            self.matrix_fig.canvas.axes.set_ylim([0,230])
        
        self.matrix_fig.canvas.axes.set_title("PEPICO matrix")
        self.matrix_fig.canvas.draw()
        
        self.matrix_loaded = 1
        
    def plot_PES(self):
        # plot the (mass-selected msPES)
        self.msPES_fig.canvas.axes.clear()
        self.msPES_fig.canvas.axes.plot(self.en,self.msPES,'r',label="ms-PES")
        self.msPES_fig.canvas.axes.fill_between(self.en,self.PES,0,alpha=0.3,color='r',label="PES")
        self.msPES_fig.canvas.axes.set_xlim(self.matrix_fig.canvas.axes.get_xlim()) #use  same x limits as for matrix
        self.msPES_fig.canvas.axes.set_title("Mass-selected Photoelectron Spectrum")
        self.msPES_fig.canvas.axes.legend()
        # plot integration limits for msPES
        ylim = self.msPES_fig.canvas.axes.get_ylim()
        self.msPES_fig.canvas.axes.plot([self.en_lim[0],self.en_lim[0]],[ylim[0],ylim[1]],'b')
        self.msPES_fig.canvas.axes.plot([self.en_lim[1],self.en_lim[1]],[ylim[0],ylim[1]],'b')
        
        self.msPES_fig.canvas.draw()
        
    def plot_TOF(self):
        self.TOF_fig.canvas.axes.clear()
        if self.show_tof:
            self.TOF_fig.canvas.axes.plot(self.TOF,self.tof,'b',label='TOF',alpha=0.3)
            self.TOF_fig.canvas.axes.plot(self.esTOF,self.tof,'b',label='es-TOF')
            self.TOF_fig.canvas.axes.set_title("Energy-selected TOF")
            self.int_lbl.setText("integration range (us):")
        else:
            self.TOF_fig.canvas.axes.plot(self.TOF,self.mass,'b',label='MS',alpha=0.3)
            self.TOF_fig.canvas.axes.plot(self.esTOF,self.mass,'b',label='es-MS')
            self.TOF_fig.canvas.axes.set_title("Energy-selected Mass Spectrum")
            self.int_lbl.setText("integration range (amu):")
            
        self.TOF_fig.canvas.axes.set_ylim(self.matrix_fig.canvas.axes.get_ylim()) #use  same y limits as for matrix
        # autoscale x axis
        ylim = self.matrix_fig.canvas.axes.get_ylim()
        if self.show_tof:
            idx = [np.min(np.where(self.tof>ylim[0])),np.max(np.where(self.tof<ylim[1]))]
        else:
            idx = [np.min(np.where(self.mass>ylim[0])),np.max(np.where(self.mass<ylim[1]))]
        imax = np.max(self.esTOF[idx[0]:idx[1]])
        self.TOF_fig.canvas.axes.set_xlim(-imax*0.02,imax*1.02)
        
        
        # plot integration limits for msPES
        xlim = self.TOF_fig.canvas.axes.get_xlim()
        self.TOF_fig.canvas.axes.plot([xlim[0],xlim[1]],[self.tof_lim[0],self.tof_lim[0]],'r')
        self.TOF_fig.canvas.axes.plot([xlim[0],xlim[1]],[self.tof_lim[1],self.tof_lim[1]],'r')
        
        self.TOF_fig.canvas.axes.legend()
        
        self.TOF_fig.canvas.draw()
        
    def get_msPES(self):
        #returns mass-selected photoelectron spectrum according to set integration limits
        #get integration limits
        self.tof_lim = np.array([float(self.tof_lim1_ent.text()),float(self.tof_lim2_ent.text())])
        if self.show_tof:
            idx = [np.min(np.where(self.tof>self.tof_lim[0])),np.max(np.where(self.tof<self.tof_lim[1]))]
        else:
            idx = [np.min(np.where(self.mass>self.tof_lim[0])),np.max(np.where(self.mass<self.tof_lim[1]))]

        # get mass-selected PES        
        self.msPES = np.sum(self.matrix[idx[0]:idx[1],:],axis=0)
        
    def get_esTOF(self):
        # returns energy-selected TOF spectrum according to set integration limits
        
        # get integration limits
        self.en_lim = np.array([float(self.en_lim1_ent.text()),float(self.en_lim2_ent.text())])
        idx = [np.min(np.where(self.en>self.en_lim[0])),np.max(np.where(self.en<self.en_lim[1]))]
        
        self.esTOF = np.sum(self.matrix[:,idx[0]:idx[1]],axis=1)
        
    def save_msPES(self):
        tmp = np.zeros((len(self.en),2))
        tmp[:,0] = self.en
        tmp[:,1] = self.msPES
        if self.show_tof:
            txt = "msPES_{:.1f}-{:.1f}_us".format(self.tof_lim[0],self.tof_lim[1])
        else:
            txt = "msPES_{:.1f}-{:.1f}_amu".format(self.tof_lim[0],self.tof_lim[1])
        file_new = self.file.replace("PEPICO_matrix_full",txt)
        np.savetxt(file_new,tmp,delimiter='\t')
        
    def save_esTOF(self):
        tmp = np.zeros((len(self.tof),3))
        tmp[:,0] = self.tof
        if self.mass is not None:
            tmp[:,1] = self.mass
        
        tmp[:,2] = self.esTOF
        txt = "esTOF_{:.1f}-{:.1f}_eV".format(self.en_lim[0],self.en_lim[1])
        file_new = self.file.replace("PEPICO_matrix_full",txt)
        np.savetxt(file_new,tmp,delimiter='\t')
        
    
    def updater(self):
        if self.matrix_loaded:
            self.get_msPES()
            self.get_esTOF()
            self.plot_PES()
            self.plot_TOF()
            
        
        
if __name__ == '__main__':
    app = QApplication(sys.argv)
    win = MainWindow()
    win.show()
    app.quitOnLastWindowClosed()
    sys.exit(app.exec_())
    app.exec_()