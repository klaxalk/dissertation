#!/bin/python

import math as m
import matplotlib
import numpy as np
import matplotlib.pyplot as plt

global h, hr, c, me

h = 6.6207015e-34 # Planc constant [J * s]
hr = 6.6207015e-34 / (2*m.pi) # reduced Planc constant [J * s]
c = 299792458 # speed of light in vacuum [m/s]
# me = 510.9989461 # mass of an electron [keV/c^2]
me = 9.10938356e-31 # mass of an electron [kg]
alpha = 1 / 137.04 # fine structure constant [-]
r_e = 2.8179e-15 # classical electron radius [m]
si_atom = 0.222e-9 # diameter of an silicone atom [m]
N_a = 6.02214086e23 # Avogadro's number [mol^-1]

Si_ro = 2.32 # Density of Si [kg/m^3]
Si_atomar = 14 # Atom number of Si [-]
Si_molar_mass = 0.02808550 # Molar mass of Si [kg/mol]
Si_n_kg = N_a / Si_molar_mass # number of Si atoms in 1 kg [1/kg]
Si_e_density = Si_n_kg * Si_atomar * Si_ro # electron denSity of Si [1/m^3]

Cd_ro = 8.65 # Density of Cd [kg/m^3]
Cd_atomar = 48 # Atom number of Cd [-]
Cd_molar_mass = 0.112411 # Molar mass of Cd [kg/mol]
Cd_n_kg = N_a / Cd_molar_mass # number of Cd atoms in 1 kg [1/kg]
Cd_e_density = Cd_n_kg * Cd_atomar * Cd_ro # electron density of Cd [1/m^3]

Te_ro = 6.24 # Density of Te [kg/m^3]
Te_atomar = 52 # Atom number of Te [-]
Te_molar_mass = 0.12760 # Molar mass of Te [kg/mol]
Te_n_kg = N_a / Te_molar_mass # number of Te atoms in 1 kg [1/kg]
Te_e_density = Te_n_kg * Te_atomar * Te_ro # electron density of Te [1/m^3]

CdTe_ro = 5.85 # Density of CdTe [kg/m^3]
CdTe_atomar = Cd_atomar + Te_atomar # Atom number of CdTe [-]
CdTe_molar_mass = 0.2400110 # Molar mass of CdTe [kg/mol]
CdTe_n_kg = N_a / CdTe_molar_mass # number of CdTe atoms in 1 kg [1/kg]
CdTe_e_density = CdTe_n_kg * CdTe_atomar * CdTe_ro # electron density of CdTe [1/m^3]

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

# scatterer thickness
thickness = 0.001 # [m]

angles = []
diff_cross_section = []
total_cross_section = 0

step = 0.1

for i in np.arange(-180, 180, step):

  # scattering angle
  theta = i * (m.pi / 180.0) # [rad]
  
  # ratio of incoming/scattered photon energy
  P = 1 / (1 + (E_0/(me*m.pow(c, 2)))*(1 - m.cos(theta)))
  
  # Klein-Nishina formula
  KN = 0.5 * m.pow(r_e, 2) * m.pow(P, 2) * (P + 1.0/P - m.pow(m.sin(theta), 2))
  
  # energy of the scattered photon
  E_p = E_0 * P
  
  # energy of the scattered photon in in keV
  E_p_keV = 0.001*energy_J_to_eV(E_p)
  
  # energy of the recoil electron
  E_e = E_0 * (1 - P)
  
  # energy of the scattered photon in in keV
  E_e_keV = 0.001*energy_J_to_eV(E_e)
  
  # target solid angle
  omega_1 = m.pow(step * (m.pi / 180), 2) # TODO: should not be here ^2?
  # omega_1 = 4*m.pow(m.pi, 3)/3
  
  d_sigma_Si = Si_e_density * KN * omega_1 * thickness
  d_sigma_CdTe = CdTe_e_density * KN * omega_1 * thickness

  angles.append(theta)
  diff_cross_section.append(KN)
  total_cross_section += d_sigma_CdTe

# create the radial data
# xs = [diff_cross_section[i]*m.cos(angles[i]) for i in range(0, len(angles))]
# ys = [diff_cross_section[i]*m.sin(angles[i]) for i in range(0, len(angles))]
# fig, ax = plt.subplots()

ax = plt.subplot(111, projection='polar')
ax.plot(angles, diff_cross_section)
# ax.set_rmax(2)
# ax.set_rticks([0.5, 1, 1.5, 2])  # less radial ticks
# ax.set_rlabel_position(-22.5)  # get radial labels away from plotted line
ax.grid(True)

plt.show()

print("E_0: {0:1.3f} keV".format(0.001*energy_J_to_eV(E_0)))
print("P: {0:1.3f}".format(P))
print("E_p_keV: {0:1.3f} keV".format(E_p_keV))
print("E_e_keV: {0:1.3f} keV".format(E_e_keV))
print("KN: {} cm^2 / sr".format(KN*10000)) # convert to cm^2 / sr
print("d_sigma_Si: {}".format(d_sigma_Si))
print("d_sigma_CdTe: {}".format(d_sigma_CdTe))
print("total_cross_section: {}".format(total_cross_section))
