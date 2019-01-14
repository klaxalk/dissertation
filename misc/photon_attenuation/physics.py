import math as m
import numpy as np

import constants
import conversions

def comptonRatio(energy, theta):

    return 1 / (1 + (energy/(constants.me*m.pow(constants.c, 2)))*(1 - m.cos(theta)))

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

# photoelectric effect coefficients
pe_a = [1.6268e-9, 1.5274e-9, 1.1330e-9, -9.12e-11]
pe_b = [-2.683e-12, -5.110e-13, -2.177e-12, 0]
pe_c = [4.173e-2, 1.027e-2, 2.013e-2, 0]
pe_p = [1, 2, 3.5, 4]

def photoelectricEffectCrossSection(material, energy):

    e_rest_mass_energy = 5.11e5

    k = energy/e_rest_mass_energy
    
    if True or k < 0.9:
        return (16/3)*m.sqrt(2)*m.pi*m.pow(constants.r_e, 2)*m.pow(constants.alpha, 4)*(m.pow(material.atomic_number, 5)/m.pow(k, 3.5))
    else:
        summ = 0
        for i in range(0, 4):
            summ += ((pe_a[i] + pe_b[i]*material.atomic_number)/(1 + pe_c[i]*material.atomic_number))*m.pow(k, -pe_p[i])
            print("summ: {}".format(summ))

        return m.pow(material.atomic_number, 5) * summ
