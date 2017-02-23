# -*- coding: utf-8 -*-
"""
Created on Thu Mar 13 16:05:59 2014

@author: rttn
"""

from math import *
import numpy as np
from pylab import *
import matplotlib.pyplot as plt
import random
from os import sys

mp = 938.28
me = 0.511

qe = 1.602176487e-19        #C
na = 6.02214179e23          #mol^-1
eps0 = 8.854187817e-12      #F/m
c = 299792458               #m/s

def beta(v):
    return v/c
    
def gamma(v):
    return 1./sqrt(1-v*v/c/c)
    
def v_of_Ekin_m0(Ekin, m0):
    b2 = 1. - 1./(1+Ekin/m0)**2
    return sqrt(b2)*c

def dEdx(ZA, m0, Ekin, I): # Bethe Bloch
    v = v_of_Ekin_m0(Ekin, m0)
    b2 = beta(v)**2
    g2 = gamma(v)**2
    mem = me / m0
    Tmax = 2*me * b2 * g2 / (1 + 2*mem * sqrt(1 + b2*g2) + mem**2)
    K = 0.307
    ln_term = log(2. * me * b2*g2 * Tmax / (I**2))
    paran = 0.5 * ln_term - b2
    dedx = K * ZA / b2 * paran
    return dedx    

E0 = 188

if sys.argv:
    if 0.01 < int(sys.argv[1]) < 1000:
        E = int(sys.argv[1])
    else:
        E = E0

m0 = mp                   #mass of projectile
I_water = 75e-6              #75 eV in J, i.e., average ionization energy of water
I_air = 85.7e-6               # ionization energy of air
I_W = 727e-6
I_Si = 173e-6
I_Al = 166e-6
#I_eff = 703e-6 # W
#I_eff = 165.51e-6 # Al
I_eff = 87.4e-6 # PMMA
#ZA_eff = 0.392 # W
#ZA_eff = 0.481 # Al
ZA_eff = 0.497 # PMMA
ZA_W = 0.403
ZA_Si = 0.5
ZA_water = 0.555
rho_Al = 2.6989
rho_W = 19.25 # g/cm3
rho_Si = 2.329
rho_air = 0.0012
rho_water = 1
rho_pmma = 1.19
#weights = [0.981, 0.0015, 0.0164, 0.00096] # W
#weights = [0.8835, 0.0095, 0.1012, 0.0059] # Al
weights = [0.777, 0.0182, 0.1936, 0.0113] # PMMA
#rhos = [rho_w, rho_air, rho_Si, rho_water] # W
#rhos = [rho_Al, rho_air, rho_Si, rho_water] # Al
rhos = [rho_pmma, rho_air, rho_Si, rho_water] # PMMA
avg_rho = sum([weights[i]*rhos[i] for i in range(4)])

dx = 5e-4 # 1 mm
thisE = E
x = 0 # in mm
dE = 0.
while thisE > 0:
    dE_norho = dEdx(ZA_eff, m0, thisE, I_eff)*(dx/10.)
    if dE_norho < 0: break
    dE = sum([dE_norho * weights[i] * rhos[i] for i in range(4)])
    x += dx
    thisE -= dE

print "Final range ({} MeV) with rho weighting (PMMA): \t".format(E), x, " mm."

"""
x = 0.
thisE = E
while thisE > 0:
    dE = rho_W * dEdx(ZA_W, m0, thisE, I_W)/10. *(dx)
    if dE<0: break
    x += dx
    thisE -= dE

print "Final range ({} MeV) with tungsten: \t\t".format(E), x, " mm."


x = 0.
thisE = E
while thisE > 0:
    dE = rho_water * dEdx(ZA_water, m0, thisE, I_water)*(dx/10.)
    if dE<0: break
    x += dx
    thisE -= dE

print "Final range ({} MeV) with water: \t\t".format(E), x, " mm."
"""
