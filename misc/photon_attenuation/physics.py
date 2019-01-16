import math as m
import numpy as np

import constants
import conversions

def comptonRatio(energy, theta):

    return 1.0 / (1.0 + (energy/(constants.me*m.pow(constants.c, 2)))*(1.0 - m.cos(theta)))

def comptonDiffCrossSection(energy, theta):

    P = comptonRatio(energy, theta)

    # # energy of the scattered photon
    # E_p = energy * P
    # # energy of the scattered photon in in keV
    # E_p_keV = 0.001*energy_J_to_eV(E_p)
    # # energy of the recoil electron
    # E_e = energy * (1 - P)
    # # energy of the scattered photon in in keV
    # E_e_keV = 0.001*energy_J_to_eV(E_e)

    return 0.5 * m.pow(P, 2) * m.pow(constants.r_e, 2) * (P + 1.0/P - m.pow(m.sin(theta), 2))

def comptonCrossSection(energy):

    total_cross_section = 0

    step = conversions.deg2rad(0.1) # [rad]

    for theta in np.arange(-m.pi, m.pi+step, step):

        # Klein-Nishina formula
        diff_cross_section = comptonDiffCrossSection(energy, theta)

        # target solid angle
        solid_angle = m.pi * m.sin(m.fabs(theta)) * step

        total_cross_section += diff_cross_section * solid_angle

    return total_cross_section

def pe_cs_gavrila_pratt_simplified(material, energy):

    e_rest_mass_energy = 511000.0
    k = energy/e_rest_mass_energy

    return (16.0/3.0)*m.sqrt(2.0)*m.pi*m.pow(constants.r_e, 2)*m.pow(constants.alpha, 4)*(m.pow(material.atomic_number, 5)/m.pow(k, 3.5))
    # return conversions.barn2m2(3.0e12*(m.pow(material.atomic_number, 4)/m.pow(energy, 3.5)))

def pe_cs_gavrila_pratt(material, energy):

    electron_rest_mass = conversions.energy_J_to_eV(constants.me*m.pow(constants.c*100, 2))

    constant = ((4.0*m.pi*m.pow(m.e, 4))/(m.pow(electron_rest_mass, 2)))*(1.0/constants.alpha)
    constant = 1.367e-22

    return constant*m.pow(constants.alpha*material.atomic_number, 5)*(electron_rest_mass/energy)*0.0001

# photoelectric effect coefficients
pe_a = [1.6268e-9, 1.5274e-9, 1.1330e-9, -9.12e-11]
pe_b = [-2.683e-12, -5.110e-13, -2.177e-12, 0.0]
pe_c = [4.173e-2, 1.027e-2, 2.013e-2, 0.0]
pe_p = [1.0, 2.0, 3.5, 4.0]

def pe_cs_hubell_k_shell(material, energy):

    summ = 0
    for i in range(0, 4):
        summ += ((pe_a[i] + pe_b[i]*material.atomic_number)/(1 + pe_c[i]*material.atomic_number))*m.pow(energy/1000000.0, -pe_p[i])

    return conversions.barn2m2(m.pow(material.atomic_number, 5)*summ)

def pe_cs_hubell(material, energy):

    boundary = 1500000.0
    
    if energy < boundary:
        return peeDavisson(material, energy)
    else:
        davisson_boundary = peeDavisson(material, boundary)
        hubell_boundary = pe_cs_hubell_k_shell(material, boundary)

        return (davisson_boundary/hubell_boundary)*pe_cs_hubell_k_shell(material, energy)
