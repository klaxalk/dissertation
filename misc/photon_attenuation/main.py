#!/bin/python

import math as m
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import os

import materials
import conversions
import constants
import physics

# latex in pyplot
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

# energy of the incoming photon
E_0_keV = [0.0027, 60, 511, 1460, 10000]
E_0 = [conversions.energy_ev_to_J(E * 1000) for E in E_0_keV] # [J]

# scatterer thickness
thickness = 0.001 # [m]

a_list = []
angles = []
klein_nishina = []
sigma_normalized_normalized = []
sigma_normalized = []
abs_prob = []
total_cross_section = []
total_area = 0

step = 0.1 # [deg]
step_theta = conversions.deg2rad(step) # [rad]
print("step_theta: {}".format(step_theta))

for energy_idx,energy in enumerate(E_0):

    sublist1 = []
    klein_nishina.append(sublist1)
    sublist2 = []
    sigma_normalized_normalized.append(sublist2)
    sublist3 = []
    sigma_normalized.append(sublist3)
    sublist4 = []
    abs_prob.append(sublist4)
    total_cross_section.append(physics.comptonCrossSection(energy))

    for i in np.arange(-180, 180+step, step):

      # scattering angle
      theta = conversions.deg2rad(i) # [rad]
      
      # Klein-Nishina formula
      compton_diff_cross_section = physics.comptonDiffCrossSection(energy, theta)
      klein_nishina[energy_idx].append(compton_diff_cross_section)
      
      # # energy of the scattered photon
      # E_p = energy * P
      # # energy of the scattered photon in in keV
      # E_p_keV = 0.001*energy_J_to_eV(E_p)
      # # energy of the recoil electron
      # E_e = energy * (1 - P)
      # # energy of the scattered photon in in keV
      # E_e_keV = 0.001*energy_J_to_eV(E_e)
      
      # target solid angle
      omega_1 = m.pi * m.sin(m.fabs(theta)) * step_theta

      # did not work well for thick materials
      # prob = CdTe_e_density * compton_diff_cross_section * omega_1 * thickness
      # prob = 1 - np.exp(-CdTe_e_density * compton_diff_cross_section * omega_1 * thickness)

      prob = compton_diff_cross_section * omega_1

      sigma_normalized[energy_idx].append(prob)

      if energy_idx == 0:
          angles.append(theta)
          total_area += omega_1

# normalize sigma_normalized
for energy_idx,energy in enumerate(E_0):

    total = sum(sigma_normalized[energy_idx])
    sigma_normalized[energy_idx] = [x/total for x in sigma_normalized[energy_idx]]

# convert to abs. prob
for energy_idx,energy in enumerate(E_0):

    total_prob = 1 - np.exp(-materials.Si.electron_density * total_cross_section[energy_idx] * thickness)
    abs_prob[energy_idx] = [x*total_prob for x in sigma_normalized[energy_idx]]

sigma_marginalized = [[x/total_cross_section[sublist_idx] for x in sublist] for sublist_idx,sublist in enumerate(sigma_normalized)]
# sigma_normalized_normalized = [[x/total_cross_section[sublist_idx] for x in sublist] for sublist_idx,sublist in enumerate(klein_nishina)]
# sigma_normalized_normalized = [[x for x in sublist] for sublist_idx,sublist in enumerate(sigma_normalized)]

for energy_idx,energy in enumerate(E_0):
    summ = 0
    for idx,value in enumerate(sigma_normalized[energy_idx]):
        summ += value
    print("prob_sum for {} keV: {}".format(energy, summ))

def plot_everything(*args):

    multiplot = False

    fig = plt.figure(1)
    fig.canvas.set_window_title("1")
    if multiplot:
        ax = plt.subplot(131, projection='polar')
    else:
        ax = plt.subplot(111, projection='polar')
    for energy_idx,energy in enumerate(E_0):
      ax.plot(angles, klein_nishina[energy_idx], label="{} keV".format(E_0_keV[energy_idx]))
    ax.legend()
    ax.grid(True)
    plt.title("Compton scattering diff. cross section for $\heta \in [0, \pi]$".format())
    plt.savefig("klein_nishina_1.png", bbox_inches="tight")

    if multiplot:
        ax = plt.subplot(132, projection='polar')
    else:
        fig = plt.figure(2)
        ax = plt.subplot(111, projection='polar')
    for energy_idx,energy in enumerate(E_0):
      ax.plot(angles, abs_prob[energy_idx], label="{} keV".format(E_0_keV[energy_idx]))
    ax.grid(True)
    ax.legend()
    plt.title('Posterior prob. of scattering by a radial angle $\Theta \in [0, \pi]$'.format())
    plt.savefig("klein_nishina_2.png", bbox_inches="tight")

    if multiplot:
        ax = plt.subplot(133, projection='polar')
    else:
        fig = plt.figure(3)
        ax = plt.subplot(111, projection='polar')
    for energy_idx,energy in enumerate(E_0):
      ax.plot(angles, sigma_normalized[energy_idx], label="{} keV".format(E_0_keV[energy_idx]))
    ax.grid(True)
    ax.legend()
    plt.title('The likelihood of scattering by a radial angle $\Theta$, {} mm Si'.format(thickness*1000))
    plt.savefig("klein_nishina_3.png", bbox_inches="tight")
    plt.show()

for energy_idx,energy in enumerate(E_0_keV):
    print("")
    print("total_cross_section[energy_idx]: {0:1.2f} re^2".format(total_cross_section[energy_idx]/m.pow(constants.r_e, 2)))

    prob = 1 - np.exp(-materials.Si.electron_density * total_cross_section[energy_idx] * thickness)
    # prob = CdTe_e_density * total_cross_section[energy_idx] * thickness

    print("{0:2.5f}% of photons are scattered for {1:6.1f} keV".format(prob, energy))

print("")
print("total area: {0:2.3f} sr".format(total_area))

pid = os.fork()
if pid == 0:
    plot_everything()
