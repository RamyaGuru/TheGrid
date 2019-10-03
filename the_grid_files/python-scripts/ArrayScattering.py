import math
import numpy as np
import matplotlib.pyplot as plt
from math import pi as pi

hbar = 6.6260696 * 10**-34 / (2 * math.pi)
kB = 1.38064852e-23

def k_mag(k_vector):
    return (k_vector[0]**2 + k_vector[1]**2 + k_vector[2]**2)**(1 / 2)


def qm(m, D):
    return 2 * math.pi * m / D


def kxprime_msigma(kx, kD, m, sign, D):
    '''
    kxprime derived by asserting the conservation laws:
    k'=k, kz'=kz, and ky'=ky-qm.
    '''
#    print(qm(m,D)**2)
#    print(kx**2 + (2 * ky * qm(m, D)))
    return (sign * (kx**2 + (2 * kD * qm(m, D)) - qm(m, D)**2)**(1 / 2))


def qx_msigma(kx, ky, m, sign, D):
    '''
    qx = kxprime - kx
    '''
    return kxprime_msigma(kx, ky, m, sign, D) - kx


def kprimes_y(k_vector, D):
    '''
    Calculates all possible kprime_vectors for a given k_vector.
    This applies the conservation laws, k'=k, kz'=kz, and ky'=ky-qm.
    '''
    kx = k_vector[0]
    ky = k_vector[1]
    kz = k_vector[2]

    # largest poisive m value set by requiring that
    # (kx**2 + (2 * ky * qm(m, D)) - qm(m, D)**2) >= 0.
    m_maxPos = math.modf((D * (ky + (kx**2 + ky**2)**(1 / 2))) / (2 * math.pi))[1]
    # largest negative m value set by the same condition
    m_maxNeg = math.modf((D * (ky - (kx**2 + ky**2)**(1 / 2))) / (2 * math.pi))[1]

    # make list of possible m values which are non-zero... so this is ignoring m=0 term?
    m_values = []
    for m in range(int(m_maxNeg), 0):
        m_values.append(m)
    for m in range(1, int(m_maxPos) + 1):
        m_values.append(m)

    # make list of possible k' vectors
    kprime_list = []
    # Add back-scattering, sign=-1 and m=0. Omit forward scattering
    # sign = 1 and m=0 or (k_vector=kprime_vector) as this causes a
    # 0/0 indeterminate in the V1_twiddle_sqs.
    kprime_list.append([-kx, ky, kz])
    for m in m_values:
        for sign in [-1, 1]:
            kprime_list.append([kxprime_msigma(kx, ky, m, sign, D), ky - qm(m, D), kz])

    return kprime_list

def kprimes_z(k_vector, D):
    '''
    Calculates all possible kprime_vectors for a given k_vector.
    This applies the conservation laws, k'=k, kz'=kz, and ky'=ky-qm.
    '''
    kx = k_vector[0]
    ky = k_vector[1]
    kz = k_vector[2]

    # largest poisive m value set by requiring that
    # (kx**2 + (2 * ky * qm(m, D)) - qm(m, D)**2) >= 0.
    m_maxPos = math.modf((D * (kz + (kx**2 + kz**2)**(1 / 2))) / (2 * math.pi))[1]
    # largest negative m value set by the same condition
    m_maxNeg = math.modf((D * (kz - (kx**2 + kz**2)**(1 / 2))) / (2 * math.pi))[1]

    # make list of possible m values which are non-zero
    m_values = []
    for m in range(int(m_maxNeg), 0):
        m_values.append(m)
    for m in range(1, int(m_maxPos) + 1):
        m_values.append(m)

    # make list of possible k' vectors
    kprime_list = []
    # Add back-scattering, sign=-1 and m=0. Omit forward scattering
    # sign = 1 and m=0 or (k_vector=kprime_vector) as this causes a
    # 0/0 indeterminate in the V1_twiddle_sqs.
    kprime_list.append([-kx, ky, kz])
    for m in m_values:
        for sign in [-1, 1]:
            kprime_list.append([kxprime_msigma(kx, kz, m, sign, D), ky, kz-qm(m,D)])

    return kprime_list

def kprimes_plot(k_vector, D):
    '''
    Plots the scattering diagram. (WORKS BUT NEEDS SOME LOVE TO MAKE IT PRETTY)
    '''
    k = k_mag(k_vector)

    dtheta = math.pi / 100.
    kx_cirlist = []
    ky_cirlist = []
    for theta in np.arange(0, 2 * math.pi + dtheta, dtheta):
        kx_cirlist.append(k * math.cos(theta))
        ky_cirlist.append(k * math.sin(theta))
    plt.figure()
    plt.plot(kx_cirlist, ky_cirlist)
    plt.plot(k_vector[0], k_vector[1], 'ko', markersize=8)
    plt.plot(np.asarray(kprimes_y(k_vector, D))[:, 0], \
             np.asarray(kprimes_y(k_vector, D))[:, 1], 'ro', markersize=4)
    plt.axes().set_aspect(aspect=1)
    plt.show(block=False)


def GammaArray(k_vector, kprime_vectors, V1_twiddle_sq, vg, n_1D, D, ax):
    '''
    Performs sum over all possible k' states.
    Requires the magnitude squared of the Fourier transform
    of the scattering potential, V1Twidle2(k_vector, kprime_vector).
    ax: axis of the dislcolation spacing
    '''
    k = k_mag(k_vector)
    kx = k_vector[0]
    kD = k_vector[ax]
    # sum over all possible kprime_vectors
    i=0
    running_sum = 0
    for kprime_vector in kprime_vectors:
        kxprime = kprime_vector[0]
        kDprime = kprime_vector[ax]
        qx = kxprime - kx
        qD = kDprime - kD
#        if i < 10:
#            print(k_vector)
#            print(kprime_vector)
        running_sum = running_sum + \
        V1_twiddle_sq(k_vector, kprime_vector) * (-qx * kx - qD * kD) * abs(kxprime) ** -1  
        i+=1
    return (n_1D / (hbar ** 2 * D ** 2 * vg * k)) * running_sum 


def tau_spectral(Gamma, k, vg_func, n):
    '''
    Calculates a spectral tau which can be applied to the
    single mode, isotropic Callaway model. This relaxation 
    time is for kappa_xx.
    
    Don't think I need to adjust for the second grid. 
     
    '''
    d_angle = (math.pi / 2) / n
    running_integrand = 0
    theta_list = np.arange(d_angle, math.pi / 2 + d_angle, d_angle)
    phi_list = np.arange(d_angle, math.pi / 2, d_angle)
    for theta in theta_list:
        for phi in phi_list:
            #this is giving you the circle of k points # conversion to sphereical coordinates for the integral
            k_vector_int = [k * np.sin(theta - (d_angle / 2)) * np.cos(phi - (d_angle / 2)), \
                            k * np.sin(theta - (d_angle / 2)) * np.sin(phi - (d_angle / 2)), \
                            k * np.cos(theta - (d_angle / 2))]
            running_integrand = running_integrand + \
            ((np.sin(theta - (d_angle / 2))**3 * np.cos(phi - (d_angle / 2))**2) * Gamma(k_vector_int, vg_func(k_vector_int))**(-1))
    
    return 8 * running_integrand * d_angle**2 * (3 / (4 * math.pi))


def transmissivity(k, vg_k, n_1D, Gamma, n):
    
    '''
    calculate the phonon spectral transmissivity
    '''
    vg = vg_k(k)
    tau = tau_spectral(Gamma, k, vg_k, 50)
    a = vg*tau*n_1D
    return a/((3/4)+a)

def Cv(k, T, omega_k):
    '''
    spectral heat capacity. there are no factors of group or phase velocity because it is in terms of k rather than omega
    '''
    k_vec = [k,0,0]
    BEexp = hbar*omega_k(k_vec)/(kB*T)
    return (3*kB/(2*pi**2))*BEexp**2*k**2*math.exp(BEexp)/(math.exp(BEexp) - 1)**2
    
def kL_spectral(Gamma, k, vg_k, omega_k, T, n):
    '''
    calculate the spectral thermal conductivity
    '''
    vg = vg_k(k)
    tau = tau_spectral(Gamma, k, vg_k, n)
# "spectral heat capacity. There is no factors of group and phase velocity because it is in terms of k rather than omega.
    #print(Cv(k, T, omega_k))
    kL = (1/3)*Cv(k,T, omega_k)*vg**2*tau
    return kL
    

def kL_T(Gamma, kmax, dk, vg_k, omega_k, T, n):
    '''
    Calculate the thermal conductivity versus temperature
    Note: this performs a sum; do a numerical integration too?
    '''
    kLint = 0
    for k in np.arange(dk, kmax, dk):
        vg = vg_k(k)
        tau = tau_spectral(Gamma, k, vg_k, n)
        kLint = kLint + Cv(k,T, omega_k)*vg**2*tau*dk
    return (1/3)*kLint

def tbc_spectral(k, vg_k, omega_k, T, Gamma, n_1D, n):
    '''
    Calculate the spectral thermal boundary conductance
    '''
    vg = vg_k(k)
    a = transmissivity(k, vg_k, n_1D, Gamma, n)
    tbc = a*vg*Cv(k, T, omega_k)
    return (1/4)*tbc

def tbc_T(kmax, dk, vg_k, omega_k, T, n_1D, Gamma, n):
    '''
    Calculate the thermal boundary conductance versus T
    Note: 
    '''
    tbc_int = 0
    for k in np.arange(dk, kmax, dk):
        vg = vg_k(k)
        a = transmissivity(k, vg_k, n_1D, Gamma, n)
        tbc_int = tbc_int + a*vg*Cv(k,T, omega_k)*dk
    return (1/4)*tbc_int