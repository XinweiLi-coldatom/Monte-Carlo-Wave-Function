# -*- coding: utf-8 -*-
"""
Created on Sun Apr 08 10:35:03 2018

@author: Xinwei
"""

import numpy as np
import mcwf

global ntot, gamma, c
ntot = 10
m = 0
hdim = (ntot-abs(m))/2+1 "dimension of initial Hilbert space"
c = -1
ntraj = 3000

for nidx in range(1,ntraj+1):
    mcwf(ntot)
    
n0ave = n0ave/ntraj
n1ave = n1ave/ntraj
nm1ave = nm1ave/ntraj
jjp1ave = jjp1ave/ntraj
jzsqave = jzsqave/ntraj
jjzave = jjzave/ntraj
vjzave = vjzave/ntraj
nmeanave = nmeanave/ntraj
entg1ave = entg1ave/ntraj
numeratorave = numeratorave/ntraj