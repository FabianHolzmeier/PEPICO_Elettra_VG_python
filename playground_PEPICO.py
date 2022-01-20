import numpy as np


PATH = "B:\\PEPICO\\rawdata\\"
    
file = "Energy_TBMA_088.txt"

energy = np.loadtxt(PATH+file, skiprows=1, usecols=(2))

en = np.unique(energy)