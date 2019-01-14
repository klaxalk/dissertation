import math as m
import numpy as np

import constants
import conversions

def comptonRatio(energy, theta):

    return 1 / (1 + (energy/(constants.me*m.pow(constants.c, 2)))*(1 - m.cos(theta)))

def comptonDiffCrossSection(energy, theta):

    P = comptonRatio(energy, theta)

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
pe_a = [1.6268e-9, -2.683e-12, 4.173e-2, 1]
pe_b = [1.5274e-9, -5.110e-13, 1.027e-2, 2]
pe_c = [1.1330e-9, -2.177e-12, 2.013e-2, 3.5]
pe_p = [-9.1e10-1, 0, 0, 4]

def photoelectricEffectCrossSection(energy):
    
   pass 
