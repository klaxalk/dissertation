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
E_0_J = [conversions.energy_ev_to_J(E * 1000) for E in E_0_keV] # [J]

# scatterer thickness
thickness = 1.0

a_list = []
angles = []
klein_nishina = []
sigma_normalized_normalized = []
sigma_normalized = []
abs_prob = []
total_cross_section = []
total_area = 0

step = conversions.deg2rad(0.1) # [rad]

material = materials.Si

print("{} electron density: {} 1/cm^3 \n".format(material.name, material.electron_density/(100*100*100)))

for energy_idx,energy in enumerate(E_0_J):

    sublist1 = []
    klein_nishina.append(sublist1)
    sublist2 = []
    sigma_normalized_normalized.append(sublist2)
    sublist3 = []
    sigma_normalized.append(sublist3)
    sublist4 = []
    abs_prob.append(sublist4)
    total_cross_section.append(physics.comptonCrossSection(energy))

    for theta in np.arange(-m.pi, m.pi+step, step):

      # Klein-Nishina formula
      compton_diff_cross_section = physics.comptonDiffCrossSection(energy, theta)
      klein_nishina[energy_idx].append(compton_diff_cross_section)
      
      # target solid angle
      omega_1 = m.pi * m.sin(m.fabs(theta)) * step

      # did not work well for thick materials
      # prob = CdTe_e_density * compton_diff_cross_section * omega_1 * thickness
      # prob = 1 - np.exp(-CdTe_e_density * compton_diff_cross_section * omega_1 * thickness)

      prob = compton_diff_cross_section * omega_1

      sigma_normalized[energy_idx].append(prob)

      if energy_idx == 0:
          angles.append(theta)
          total_area += omega_1

# normalize sigma_normalized
for energy_idx,energy in enumerate(E_0_J):

    total = sum(sigma_normalized[energy_idx])
    sigma_normalized[energy_idx] = [x/total for x in sigma_normalized[energy_idx]]

# convert to abs. prob
for energy_idx,energy in enumerate(E_0_J):

    total_prob = 1 - np.exp(-material.electron_density * total_cross_section[energy_idx] * thickness)
    abs_prob[energy_idx] = [x*total_prob for x in sigma_normalized[energy_idx]]

sigma_marginalized = [[x/total_cross_section[sublist_idx] for x in sublist] for sublist_idx,sublist in enumerate(sigma_normalized)]
# sigma_normalized_normalized = [[x/total_cross_section[sublist_idx] for x in sublist] for sublist_idx,sublist in enumerate(klein_nishina)]
# sigma_normalized_normalized = [[x for x in sublist] for sublist_idx,sublist in enumerate(sigma_normalized)]

for energy_idx,energy in enumerate(E_0_J):
    summ = 0
    for idx,value in enumerate(sigma_normalized[energy_idx]):
        summ += value
    print("prob_sum for {} keV: {}".format(conversions.energy_J_to_eV(energy)/1000, summ))

# #{ plot_everything()
    
    def plot_everything(*args):
    
        multiplot = False
        plot_compton = False
    
        if plot_compton:

          fig = plt.figure(1)
          fig.canvas.set_window_title("1")
          if multiplot:
              ax = plt.subplot(141, projection='polar')
          else:
              ax = plt.subplot(111, projection='polar')
          for energy_idx,energy in enumerate(E_0_J):
            ax.plot(angles, klein_nishina[energy_idx], label="{} keV".format(E_0_keV[energy_idx]))
          ax.legend()
          ax.grid(True)
          plt.title("Compton scattering diff. cross section for $\heta \in [0, \pi]$".format())
          plt.savefig("klein_nishina_1.png", bbox_inches="tight")
          
          if multiplot:
              ax = plt.subplot(142, projection='polar')
          else:
              fig = plt.figure(2)
              ax = plt.subplot(111, projection='polar')
          for energy_idx,energy in enumerate(E_0_J):
            ax.plot(angles, abs_prob[energy_idx], label="{} keV".format(E_0_keV[energy_idx]))
          ax.grid(True)
          ax.legend()
          plt.title('Posterior prob. of scattering by a radial angle $\Theta \in [0, \pi]$'.format())
          plt.savefig("klein_nishina_2.png", bbox_inches="tight")
          
          if multiplot:
              ax = plt.subplot(143, projection='polar')
          else:
              fig = plt.figure(3)
              ax = plt.subplot(111, projection='polar')
          for energy_idx,energy in enumerate(E_0_J):
            ax.plot(angles, sigma_normalized[energy_idx], label="{} keV".format(E_0_keV[energy_idx]))
          ax.grid(True)
          ax.legend()
          plt.title('The likelihood of scattering by a radial angle $\Theta$, {} mm {}'.format(thickness*1000, material.name))
          plt.savefig("klein_nishina_3.png", bbox_inches="tight")
          plt.show()

        fig = plt.figure(10)
        ax = plt.subplot(111)
        ax.plot(pe_energies, prob_1, label="Gavrila-Pratt simplified Si, 0.0003 mm".format())
        ax.plot(pe_energies, prob_2, label="Gavrila-Pratt Si, 0.0003 mm".format())
        ax.grid(True)
        ax.legend()
        plt.show()
    
# #} end of plot_everything()

for energy_idx,energy in enumerate(E_0_keV):
    print("")
    print("total_cross_section[energy_idx]: {0:1.2f} re^2".format(total_cross_section[energy_idx]/m.pow(constants.r_e, 2)))

    prob = 1 - np.exp(-material.electron_density * total_cross_section[energy_idx] * thickness)
    # prob = material.electron_density * total_cross_section[energy_idx] * thickness

    print("{0:2.1f}% of photons are scattered for {1:6.1f} keV".format(prob*100, energy))

print("")
print("total area: {0:2.3f} sr".format(total_area))

prob_1 = []
prob_2 = []
prob_3 = []
pe_energies = []

for e in range(1, 50, 1):
    pe_energies.append(e)

    pe_cs_1 = physics.pe_cs_gavrila_pratt_simplified(materials.Si, e*1000.0)
    print("pe_cs_1: {}".format(pe_cs_1))
    pe_cs_2 = physics.pe_cs_gavrila_pratt(materials.Si, e*1000.0)
    print("pe_cs_2: {}".format(pe_cs_2))

    prob = 1 - np.exp(-materials.Si.atomic_density * pe_cs_1 * 0.0003)
    prob_1.append(prob)

    prob = 1 - np.exp(-materials.Si.atomic_density * pe_cs_2 * 0.0003)
    prob_2.append(prob)

pid = os.fork()
if pid == 0:
    plot_everything()
