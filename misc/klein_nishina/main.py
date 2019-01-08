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
N_a = 6.02214086e23 # Avogadro's number [1/mol]

Si_ro = 2320 # Density of Si [kg/m^3]
Si_atomic = 14 # Atom number of Si [-]
Si_molar_mass = 0.02808550 # Molar mass of Si [kg/mol]
Si_n_kg = N_a / Si_molar_mass # number of Si atoms in 1 kg [1/kg]
Si_e_density = Si_n_kg * Si_atomic * Si_ro # electron denSity of Si [1/m^3]

Cd_ro = 8650 # Density of Cd [kg/m^3]
Cd_atomic = 48 # Atom number of Cd [-]
Cd_molar_mass = 0.112411 # Molar mass of Cd [kg/mol]
Cd_n_kg = N_a / Cd_molar_mass # number of Cd atoms in 1 kg [1/kg]
Cd_e_density = Cd_n_kg * Cd_atomic * Cd_ro # electron density of Cd [1/m^3]

Na_ro = 968 # Density of Na [kg/m^3]
Na_atomic = 11 # Atom number of Na [-]
Na_molar_mass = 0.022989770 # Molar mass of Na [kg/mol]
Na_n_kg = N_a / Na_molar_mass # number of Na atoms in 1 kg [1/kg]
Na_e_density = Na_n_kg * Na_atomic * Cd_ro # electron density of Na [1/m^3]

I_ro = 4933 # Density of I [kg/m^3]
I_atomic = 53 # Atom number of I [-]
I_molar_mass = 0.253808940 # Molar mass of I [kg/mol]
I_n_kg = N_a / I_molar_mass # number of I atoms in 1 kg [1/kg]
I_e_density = I_n_kg * I_atomic * Cd_ro # electron density of I [1/m^3]

NaI_ro = 3670 # Density of NaI [kg/m^3]
NaI_atomic = Na_atomic + I_atomic # Atom number of NaI [-]
NaI_molar_mass = 0.149894239 # Molar mass of NaI [kg/mol]
NaI_n_kg = N_a / NaI_molar_mass # number of NaI atoms in 1 kg [1/kg]
NaI_e_density = NaI_n_kg * NaI_atomic * NaI_ro # electron density of NaI [1/m^3]
print("NaI_e_density: {}".format(NaI_e_density/(100*100*100)))

Te_ro = 6240 # Density of Te [kg/m^3]
Te_atomic = 52 # Atom number of Te [-]
Te_molar_mass = 0.12760 # Molar mass of Te [kg/mol]
Te_n_kg = N_a / Te_molar_mass # number of Te atoms in 1 kg [1/kg]
Te_e_density = Te_n_kg * Te_atomic * Te_ro # electron density of Te [1/m^3]

CdTe_ro = 5850 # Density of CdTe [kg/m^3]
CdTe_atomic = Cd_atomic + Te_atomic # Atom number of CdTe [-]
CdTe_molar_mass = 0.2400110 # Molar mass of CdTe [kg/mol]
CdTe_n_kg = N_a / CdTe_molar_mass # number of CdTe atoms in 1 kg [1/kg]
CdTe_e_density = CdTe_n_kg * CdTe_atomic * CdTe_ro # electron density of CdTe [1/m^3]
print("CdTe_e_density: {}".format(CdTe_e_density/(100*100*100)))

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
E_0 = energy_ev_to_J(1000000) # [J]

# scatterer thickness
thickness = 0.001 # [m]

angles = []
prob_const_area = []
prob_angle = []
total_cross_section = 0
total_area = 0

step = 0.01 # [deg]
step_theta = step * (m.pi / 180.0) # [rad]

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
  omega_1 = m.pi*m.sin(m.fabs(theta))*step_theta
  
  d_sigma_CdTe = CdTe_e_density * KN * omega_1 * thickness

  try:
      prob_theta = d_sigma_CdTe/step_theta
  except:
      prob_theta = prob_theta

  try:
      prob_area = (d_sigma_CdTe/omega_1)*(m.pi / 180)
  except:
      prob_area = prob_area

  angles.append(theta)
  prob_angle.append(prob_theta)
  prob_const_area.append(prob_area)
  total_cross_section += d_sigma_CdTe
  total_area += omega_1

# # normalize the distribution
# prob_const_area = [x/total_cross_section for x in prob_const_area]

# create the radial data
# xs = [diff_cross_section[i]*m.cos(angles[i]) for i in range(0, len(angles))]
# ys = [diff_cross_section[i]*m.sin(angles[i]) for i in range(0, len(angles))]
# fig, ax = plt.subplots()

ax = plt.subplot(121, projection='polar')
ax.plot(angles, prob_const_area)
ax.grid(True)
ax = plt.subplot(122, projection='polar')
ax.plot(angles, prob_angle)
ax.grid(True)
plt.show()

print("E_0: {0:1.3f} keV".format(0.001*energy_J_to_eV(E_0)))
print("P: {0:1.3f}".format(P))
print("E_p_keV: {0:1.3f} keV".format(E_p_keV))
print("E_e_keV: {0:1.3f} keV".format(E_e_keV))
print("KN: {} cm^2 / sr".format(KN*10000)) # convert to cm^2 / sr
print("d_sigma_CdTe: {}".format(d_sigma_CdTe))
print("total_area: {}".format(total_area))
print("total_cross_section: {}".format(total_cross_section))
