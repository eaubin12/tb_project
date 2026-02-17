#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os 
from thunderboltz import ThunderBoltz # imports the main simulation object 
from thunderboltz.input import He_TB # b uilt in He model preset 

import matplotlib.pyplot as plt 

#read single calculations and return thunderboltz objects
from thunderboltz import read 
from thunderboltz import ThunderBoltz # imports the main simulation object 
from thunderboltz import CrossSections 
from thunderboltz import Process 
import numpy as np
import matplotlib.pyplot as plt 
from thunderboltz import read 
import pandas as pd
import thunderboltz as tb
from thunderboltz.parameters import WrapParameters


xs = CrossSections() #initalize an empty cross sections object
#xs.from_LXCat("cs_LXCat_all_mtcs.txt")
xs.from_LXCat("cs_LXCat_all_N.txt")
xs.set_fixed_background(fixed=True)

mask = xs.table.rtype == "Ionization"
xs.table.loc[mask, "rtype"] = "IonizationNoEGen"

### viewing cross sections ###
print(xs.table)
print(xs.data) # view cs associated with each process 

xs.plot_cs() # plot of cs data 
plt.show()

#Fixed-Particle2
#CrossSections.set_fixed_background(fixed=True)

#and add wrap parameters 
#choose relatively reasonable electric field strength 50 - 100 and run the jobs for that one field strength. 
#number of particles should be the same throughout
#if decide to go away from default option, then alot of degrees of freedom can play with. best to keep NP the same, same stastics used throughout
# E_red = [E V/m] / [n m-3] 


# In[3]:


# In[5]:


### run multiple calculations in sequence ###

os.makedirs("Nitrogen-Aubin-1000")

calc = tb.ThunderBoltz(
    cs=xs,
    #DT = 1e-11, 
    NS=500000, # Number of time steps s
    L=1e-6, #The cell size (m).
    EP_0 = 23.9795,
    Nmin=300,
    #DE = 0.05,
    NN=3000,
    pct_ion=0.3, #from helium example and needed for autostep, may need to be changed
   # DT=4.6982187312778576e-11, # Time step
    #NP= [30000, 3000], #Number of electrons, macroparticles
    eesd="uniform", # Use the uniform electron energy sharing ionization model.
    eadf="N_Aubin", # Use isotropic elastic scattering.
   # autostep=True, # (bool) Flag to calculate DT / NP / E from Ered / L /pct_ion / DE / EP_0, default is False.
    egen=False, # Do not generate secondary electrons in ionization events.
   # MEM = 10, # in GB

)
#fields = [0.001, 0.01, 0.1] #, 0.3, 1, 3, 10, 30, 100, 300, 1000]
fields = [1000]

#loop through the field values
for field in fields: 

    # Create a new directory for this reduced-field calculation
    subdir = os.path.join("Nitrogen-Aubin-1000", f"{field}Td")
    os.makedirs(subdir)

    # Update simulation parameters
    calc.set_(Ered=field, autostep=True, directory=subdir)
   # calc.set_(Ered=field, directory=subdir)

    # Run the calculation 
    calc.run()




