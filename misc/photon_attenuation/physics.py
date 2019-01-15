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

# photoelectric effect coefficients
pe_a = [1.6268e-9, 1.5274e-9, 1.1330e-9, -9.12e-11]
pe_b = [-2.683e-12, -5.110e-13, -2.177e-12, 0.0]
pe_c = [4.173e-2, 1.027e-2, 2.013e-2, 0.0]
pe_p = [1.0, 2.0, 3.5, 4.0]

def peeScofield(material, energy):

    e_rest_mass_energy = 511000.0

    k = energy/e_rest_mass_energy

    return (16.0/3.0)*m.sqrt(2.0)*m.pi*m.pow(constants.r_e, 2)*m.pow(constants.alpha, 4)*(m.pow(material.atomic_number, 5)/m.pow(k, 3.5))

def peeSauter(material, energy):

    E_e = 511000.0 # ev
    Thomson = 0.66526 # barn
    Eb = 4.26 # ev
    Z = material.atomic_number
    gamma = (energy - 4.26 + 511000.0)/511000.0

    square_bracket = 4/3

    return (3.0/2.0)*Thomson*m.pow(constants.alpha, 4)*m.pow((Z*E_e)/energy, 5)*m.pow(m.pow(gamma, 2) - 1, 3/2.0)*square_bracket;

def peeHubell(material, energy):

    summ = 0
    for i in range(0, 4):
        summ += ((pe_a[i] + pe_b[i]*material.atomic_number)/(1 + pe_c[i]*material.atomic_number))*m.pow(energy/1000000.0, -pe_p[i])

    return conversions.barn2m2(m.pow(material.atomic_number, 5)*summ)

def peeCrossSection(material, energy):
    
    boundary = 1000000.0
    
    if energy < boundary:
        return peeScofield(material, energy)
    else:
        # pee_scofield = peeScofield(material, boundary)
        pee_scofield = conversions.barn2m2(0.5)
        pee_hubell = peeHubell(material, boundary)

        return (pee_scofield/pee_hubell)*peeHubell(material, energy)
