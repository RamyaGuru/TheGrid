#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 11:44:05 2019

@author: ramyagurunathan

STGB scattering

Steps:
    1) Set up a step function in the misorientation angle theta
    2) Set up the strain field 
    3) 
    
    Si-Si twist boundary
"""


import ArrayScattering as AS
from AngularVs import ElasticProps, v_long, v_shear, v_sound, v_xy
import Callaway as Cal
import math
import numpy as np
import matplotlib.pyplot as plt
from math import asin, acos

np.seterr(divide='raise', invalid="raise")


'''
Crystal Properties
'''
hbar = 6.626e-34/(2*math.pi)
# Crystal properties
vs = 6084.       # Speed of sound [m/s]
V = 2E-29       # Volume per atom [m^3]
N = 2           # Number of atoms per primitive unit cell
def gamma(k_vector):
    return 1           # Gruneissen parameter
nu = 0.27       # Poisson's ratio
k_max = (6 * math.pi**2 / (V * N))**(1 / 3)  # maximum k-vector (spherical BZ, only acoustic brances)
omega_D = vs * k_max  # Debye frequency [Hz]
def omega_k(k_vector):
    if np.size(k_vector)>1:
        k = AS.k_mag(k_vector)
    else:
        k = k_vector
    return vs * k
def vg_k(k_vector):
    return vs


'''
Microstructural Features
'''

## Microstructural parameters
b = (V * N) ** (1 / 3)  # Burger's vector [m]
d_GS = 350E-9           # Average grain size [m]
n_1D = 3 / d_GS         # Number density of GBs [m^-1]
D = 1E-9                # Linear defects spacing [m]



'''
Angular dependence of the speed of sound 
'''

    

'''
Rotation component:
'''

'''
theta here is the misorientation angle
inc: angle of incidence determined from the k-vector
direction normal to the plane is the x-axis
'''
def V_twiddle_sq_R(k_vector):
    k = AS.k_mag(k_vector)
    #q_vector = np.asarray(kprime_vector)- np.asarray(k_vector)
    #calculate the angle of incidence
    inc = acos(k_vector[0]/(k))
    vs_new = v_sound(inc)
    # k will cancel out
    return abs(hbar*(vs_new - vs)/(2*math.pi))**2


'''
Strain components of the "n" array: Dislocaiton line in z, spacing in y
'''

def Vn_twiddle_sq_E13(k_vector, kprime_vector):
    k = AS.k_mag(k_vector)
    q_vector = np.asarray(kprime_vector) - np.asarray(k_vector)
    return abs(hbar*omega_k(k_vector)*gamma(k)*(b*q_vector[1])/(2* (q_vector[0]**2 + q_vector[1]**2)))**2
    
def Vn_twiddle_sq_E23(k_vector, kprime_vector):
    k = AS.k_mag(k_vector)
    q_vector = np.asarray(kprime_vector) - np.asarray(k_vector)
    return abs(hbar*omega_k(k_vector)*gamma(k)*(b*q_vector[0]/(2*(q_vector[0]**2 + q_vector[1]**2))))**2

def V_twiddle_sq_n(k_vector, kprime_vector):
    return Vn_twiddle_sq_E13(k_vector, kprime_vector) + Vn_twiddle_sq_E23(k_vector, kprime_vector)
    

'''
Strain components of the "m" array: Dislocation line in y, spacing in z
'''
def Vm_twiddle_sq_E12(k_vector, kprime_vector):
    k = AS.k_mag(k_vector)
    q_vector = np.asarray(kprime_vector) - np.asarray(k_vector)
    return abs(hbar*omega_k(k_vector)*gamma(k)*(b*q_vector[2])/(2* (q_vector[0]**2 + q_vector[2]**2)))**2
    
def Vm_twiddle_sq_E23(k_vector, kprime_vector):
    k = AS.k_mag(k_vector)
    q_vector = np.asarray(kprime_vector) - np.asarray(k_vector)
    return abs(hbar*omega_k(k_vector)*gamma(k)*(b*q_vector[0]/(2*(q_vector[0]**2 + q_vector[2]**2))))**2

def V_twiddle_sq_m(k_vector, kprime_vector):
    return Vm_twiddle_sq_E12(k_vector, kprime_vector) + Vm_twiddle_sq_E23(k_vector, kprime_vector)

'''
STGB scattering potential
'''
#Still need to include the rotation term
def Gamma_GBS(k_vector, kprime_yvectors, kprime_zvectors, vg, n_1D, D):
   return AS.GammaArray(k_vector, kprime_yvectors, V_twiddle_sq_n, vg, n_1D, D, 1) \
          + AS.GammaArray(k_vector, kprime_zvectors, V_twiddle_sq_m, vg, n_1D, D, 2) 


def Gamma_GBS_rot(k_vector, kprime_yvectors, kprime_zvectors, vg, n_1D, D):
   return AS.GammaArray(k_vector, kprime_yvectors, V_twiddle_sq_n, vg, n_1D, D, 1) \
          + AS.GammaArray(k_vector, kprime_zvectors, V_twiddle_sq_m, vg, n_1D, D, 2) \
          + V_twiddle_sq_R(k_vector)
          
def Gamma(k_vector, vg):
    return Gamma_GBS(k_vector, AS.kprimes_y(k_vector, D), AS.kprimes_z(k_vector, D), vg, n_1D, D) * 1E-9 #what's this function for?

def Gamma_rot(k_vector, vg):
    return Gamma_GBS_rot(k_vector, AS.kprimes_y(k_vector, D), AS.kprimes_z(k_vector, D), vg, n_1D, D) * 1E-9
#NEED TO MODIFY
k_vector = [0.2 * k_max, 0, 0]


# Plot Gamma_GBS(k_vector) for normal incidence
n_k = 1.E2
dk = k_max / n_k
k_mags = np.arange(dk, k_max, dk)
k_norm = k_mags / k_max
k_vectors = []
for k in k_mags:
    k_vectors.append([k, 0, 0]) #the non-zero term in this could be a problem? perpendicular to axis?

Gamma_GBS_list = []
for k_vector in k_vectors:
    #Gamma_GBS_list.append(Gamma(k_vector, vg_k(k_vector)))
    Gamma_GBS_list.append(Gamma_rot(k_vector, vg_k(k_vector)))
        
plt.figure()
plt.xlim((0, 1))
plt.ylim((0, 20))
plt.xlabel(r'$k/k_{\mathrm{max}}$', fontsize=16)
plt.ylabel(r'$\Gamma \; \mathrm{(ns^{-1})}$', fontsize=16)
plt.plot(k_norm, Gamma_GBS_list)
plt.savefig('twist_diffD1e-09.pdf', dpi=400, bbox_inches='tight')
plt.show(block=False)
# Convergence of tau_spectral, n_angle=100 is sufficient.
# n_angle_list = np.arange(4, 100, 2)
# tau_nlist = []
# for n_angle in n_angle_list:
#     tau_nlist.append(AS.tau_spectral(Gamma, k_max / 5., vg_k, n_angle))

# plt.figure()
# plt.xlabel('n', fontsize=16)
# plt.ylabel(r'$\tau(k)^{-1} \; \mathrm{(ns^{-1})}$', fontsize=16)
# plt.plot(n_angle_list, tau_nlist)
# plt.show(block=False)

#%%
# Calculation of spectral tau and kappa
omega_list = []
vg_list = []
tau_list = []
kappa_list = []
trans_list = []
T = 300
for k in k_mags:
    omega_list.append(omega_k([k,0,0])) # omega and vg are supposed to be only a function of k, not k_vector. This is tacky and needs to be fixed!
    vg_list.append(vg_k([k,0,0]))
    tau_list.append(AS.tau_spectral(Gamma_rot, k, vg_k, 50))
    trans_list.append(AS.transmissivity(k, vg_k, n_1D, Gamma_rot, 50))
    #kappa_list.append(AS.kL_spectral(Gamma_rot, k, vg_k, omega_k, T, 50))


plt.figure()
plt.xlabel(r'$k \; \mathrm{(m^{-1})}$', fontsize=16)
plt.ylabel(r'$\tau \; \mathrm{(ns)}$', fontsize=16)
plt.plot(k_mags, tau_list)
plt.savefig('twistBoundary_D1e-9.pdf', dpi=400, bbox_inches = 'tight')
plt.show(block=False)

#plt.figure
#plt.xlabel(r'$k \; \mathrm{(m^{-1})}$', fontsize=16)
#plt.ylabel(r'$\kappa_\mathrm{L} \; \mathrm{(W/m/K)}$', fontsize=16)
#plt.plot(k_mags, kappa_list)
#plt.savefig('twistKappa_D1e-9.pdf', dpi=400, bbox_inches = 'tight')
#plt.show(block=False)

plt.show()


#kappaT = []
#temps = np.linspace(100, 500, 5)
#for T in temps:
#    kappaT.append(AS.kL_T(Gamma_rot, k_max, dk, vg_k, omega_k, T, 50))
#plt.figure()
#plt.loglog(temps, kappaT)
#plt.xlabel(r'$\mathbf{x}$ in $\mathbf{(Mg_2Sn_{1-x})Si_x}$')
#plt.ylabel(r'$\kappa_L$ $\mathrm{(W/m/K)}$')
    

    