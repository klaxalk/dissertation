#!/bin/python

import math as m
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import os

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

global h, hr, c, me

# #{ physical constants

h = 6.6207015e-34 # Planc constant [J * s]
hr = 6.6207015e-34 / (2*m.pi) # reduced Planc constant [J * s]
c = 299792458 # speed of light in vacuum [m/s]
me = 9.10938356e-31 # mass of an electron [kg]
alpha = 1 / 137.04 # fine structure constant [-]
r_e = 2.8179e-15 # classical electron radius [m]
si_atom = 0.222e-9 # diameter of an silicone atom [m]
N_a = 6.02214086e23 # Avogadro's number [1/mol]

# #} end of physical constants

# #{ materials

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

# #} end of materials

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

cw = h / (me * c) # Compton wavelength[m]
rcw = hr / (me * c) # reduced Compton wavelength [m]

# energy of the incoming photon
E_0_keV = [0.0027, 60, 511, 1460, 10000]
E_0 = [energy_ev_to_J(E * 1000) for E in E_0_keV] # [J]

# scatterer thickness
thickness = 0.003 # [m]

a_list = []
angles = []
klein_nishina = []
prob_angle_normalized = []
prob_angle = []
total_cross_section = []
total_area = 0

step = 0.1 # [deg]
step_theta = step * (m.pi / 180.0) # [rad]
print("step_theta: {}".format(step_theta))

for energy_idx,energy in enumerate(E_0):

    sublist1 = []
    klein_nishina.append(sublist1)
    sublist2 = []
    prob_angle_normalized.append(sublist2)
    sublist3 = []
    prob_angle.append(sublist3)
    total_cross_section.append(0)

    for i in np.arange(-180, 180, step):

      # scattering angle
      theta = i * (m.pi / 180.0) # [rad]
      
      # ratio of scattered/incoming photon energy
      P = 1 / (1 + (energy/(me*m.pow(c, 2)))*(1 - m.cos(theta)))
      
      # Klein-Nishina formula
      KN = 0.5 * m.pow(P, 2) * m.pow(r_e, 2) * (P + 1.0/P - m.pow(m.sin(theta), 2))
      klein_nishina[energy_idx].append(KN)
      
      # energy of the scattered photon
      E_p = energy * P
      
      # energy of the scattered photon in in keV
      E_p_keV = 0.001*energy_J_to_eV(E_p)
      
      # energy of the recoil electron
      E_e = energy * (1 - P)
      
      # energy of the scattered photon in in keV
      E_e_keV = 0.001*energy_J_to_eV(E_e)
      
      # target solid angle
      omega_1 = m.pi*m.sin(m.fabs(theta))*step_theta
      
      d_sigma = CdTe_e_density * KN * omega_1 * thickness

      try:
          # prob_theta = d_sigma/step_theta
          prob_theta = d_sigma
      except:
          prob_theta = prob_theta

      prob_angle[energy_idx].append(prob_theta)
      # prob_const_area[energy_idx].append(prob_area)
      total_cross_section[energy_idx] += d_sigma

      if energy_idx == 0:
          angles.append(theta)
          total_area += omega_1

# # normalize the distribution
print("total_cross_section[0]: {}".format(total_cross_section[0]))
prob_angle_normalized = [[x/total_cross_section[sublist_idx] for x in sublist] for sublist_idx,sublist in enumerate(prob_angle)]
# prob_angle_normalized = [[x for x in sublist] for sublist_idx,sublist in enumerate(prob_angle)]

summ = 0
for idx,value in enumerate(prob_angle_normalized[0]):
    summ += value

print("summ: {}".format(summ))

# create the radial data
# xs = [diff_cross_section[i]*m.cos(angles[i]) for i in range(0, len(angles))]
# ys = [diff_cross_section[i]*m.sin(angles[i]) for i in range(0, len(angles))]
# fig, ax = plt.subplots()

def plot_everything(*args):

    # plt.figure(1)
    # ax = plt.subplot(111, projection='polar')
    # for energy_idx,energy in enumerate(E_0):
    #   ax.plot(angles, klein_nishina[energy_idx], label="{} keV".format(E_0_keV[energy_idx]))
    # ax.legend()
    # ax.grid(True)
    # plt.title("Compton scattering diff. cross section for $\Theta \in [0, \pi]$".format())
    # plt.savefig("klein_nishina_1.png", bbox_inches="tight")

    plt.figure(2)
    ax = plt.subplot(111, projection='polar')
    for energy_idx,energy in enumerate(E_0):
      ax.plot(angles, prob_angle[energy_idx], label="{} keV".format(E_0_keV[energy_idx]))
    ax.grid(True)
    ax.legend()
    plt.title('Abs. prob. of scattering for $\Theta \in [0, \pi]$, {} mm CdTe'.format(thickness*1000))
    plt.savefig("klein_nishina_2.png", bbox_inches="tight")

    # plt.figure(3)
    # ax = plt.subplot(111, projection='polar')
    # for energy_idx,energy in enumerate(E_0):
    #   ax.plot(angles, prob_angle_normalized[energy_idx], label="{} keV".format(E_0_keV[energy_idx]))
    # ax.grid(True)
    # ax.legend()
    # plt.title('Marginalized prob. for radial angle $\Theta \in [0, \pi]$, {} mm CdTe'.format(thickness*1000))
    # plt.savefig("klein_nishina_3.png", bbox_inches="tight")
    plt.show()

pid = os.fork()
if pid == 0:
    plot_everything()

for energy_idx,energy in enumerate(E_0_keV):
    print("{0:2.1f}% of photons are scattered for {1:6.0f} keV".format(total_cross_section[energy_idx]*100, energy))
print("total area: {0:2.2f} sr".format(total_area))
