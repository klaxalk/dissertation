#!/bin/python

import math as m

global h, hr, c, me

h = 6.6207015e-34 # Planc constant [J * s]
hr = 6.6207015e-34 / (2*m.pi) # reduced Planc constant [J * s]
c = 299792458 # speed of light in vacuum [m/s]
# me = 510.9989461 # mass of an electron [keV/c^2]
me = 9.10938356e-31 # mass of an electron [kg]
alpha = 1 / 137.04 # fine structure constant [-]
r_e = 2.8179e-15 # classical electron radius [m]

cw = h / (me * c) # Compton wavelength[m]
rcw = hr / (me * c) # reduced Compton wavelength [m]

# #{ conversions

# return energy in Jouls
def wavelength_to_J(wavelength):

    return (h * c) / wavelength

# return energy in Jouls
def wavelength_um_to_ev(wavelength):

    return 1.2398 / wavelength

# return energy in eV
def energy_J_to_eV(energy):

    return 6.242e18 * energy

# return energy in eV
def energy_ev_to_J(energy):

    return energy / 6.242e18

# #} end of conversions

# energy of the incoming photon
E_0 = energy_ev_to_J(662000) # [J]

# scattering angle
theta = 30 * (m.pi / 180.0) # [rad]

print("alpha*r_e: {}".format(alpha*r_e))

# ratio of incoming/scattered photon energy
P = 1 / (1 + (E_0/(me*m.pow(c, 2)))*(1 - m.cos(theta)))

# Klein-Nishina formula
KN = 0.5 * m.pow(r_e, 2) * m.pow(P, 2) * (P + 1.0/P - m.pow(m.sin(theta), 2))

# energy of the scatter photon
E_p = E_0 * P

# energy of the scattered photon in in keV
E_p_keV = 0.001*energy_J_to_eV(E_p)

# energy of the recoil electron
E_e = E_0 * (1 - P)

# energy of the scattered photon in in keV
E_e_keV = 0.001*energy_J_to_eV(E_e)

print("E_0: {0:1.3f} keV".format(0.001*energy_J_to_eV(E_0)))
print("P: {0:1.3f}".format(P))
print("E_p_keV: {0:1.3f} keV".format(E_p_keV))
print("E_e_keV: {0:1.3f} keV".format(E_e_keV))
print("KN: {} cm^2 / sr".format(KN*10000))
