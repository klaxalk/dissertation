#!/bin/python

import math as m
import matplotlib
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

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
scatterer_z = 0.001
scatterer_material = materials.Si

# absorber thickness
absorber_z = 0.001
absorber_material = materials.CdTe

a_list = []
angles = []
klein_nishina = []
sigma_normalized_normalized = []
sigma_normalized = []
abs_prob = []
total_cross_section = []
total_area = 0

step = conversions.deg2rad(0.1) # [rad]

print("scatterer {} electron density: {} 1/cm^3 \n".format(scatterer_material.name, scatterer_material.electron_density/(100*100*100)))
print("absorber {} electron density: {} 1/cm^3 \n".format(absorber_material.name, absorber_material.electron_density/(100*100*100)))

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
      # prob = CdTe_e_density * compton_diff_cross_section * omega_1 * scatterer_z
      # prob = 1 - np.exp(-CdTe_e_density * compton_diff_cross_section * omega_1 * scatterer_z)

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

    total_prob = 1 - np.exp(-scatterer_material.electron_density * total_cross_section[energy_idx] * scatterer_z)
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
      plt.title('The likelihood of scattering by a radial angle $\Theta$, {} mm {}'.format(scatterer_z*1000, scatterer_material.name))
      plt.savefig("klein_nishina_3.png", bbox_inches="tight")
      plt.show()

    if multiplot:
        fig = plt.figure(10)
        ax = plt.subplot(211)
    else:
        fig = plt.figure(10)
        ax = plt.subplot(111)
    plt.yscale('log')
    for axis in [ax.xaxis, ax.yaxis]:
        axis.set_major_formatter(ScalarFormatter())
    ax.plot(pe_energies, prob_pe_scatterer, label="Photoelectric effect prob., {}, {} mm".format(scatterer_material.name, scatterer_z/1000.0))
    ax.plot(pe_energies, prob_cs_scatterer, label="Compton scattering prob., {}, {} mm".format(scatterer_material.name, scatterer_z/1000.0))
    ax.plot(pe_energies, prob_scatterer_attenuation, label="Total attenutaion by PE and CS, {}, {} mm".format(scatterer_material.name, scatterer_z/1000.0), linestyle="dashed")
    ax.set_xlabel("Photon energy [keV]")
    ax.set_ylabel("Probability [-]")
    ax.grid(True)
    ax.legend()
    plt.savefig("scatterer_attenuation.png", bbox_inches="tight")

    if multiplot:
        ax = plt.subplot(212)
    else:
        fig = plt.figure(11)
        ax = plt.subplot(111)
    ax.plot(pe_energies, prob_pe_absorber, label="Photoelectric effect prob., {}, {} mm".format(absorber_material.name, absorber_z/1000.0))
    ax.plot(pe_energies, prob_cs_absorber, label="Compton scattering prob., {}, {} mm".format(absorber_material.name, absorber_z/1000.0))
    ax.plot(pe_energies, prob_absorber_attenuation, label="Total attenutaion by PE and CS, {}, {} mm".format(absorber_material.name, absorber_z/1000.0), linestyle="dashed")
    ax.set_xlabel("Photon energy [keV]")
    ax.set_ylabel("Probability [-]")
    ax.legend()
    ax.grid(True)
    plt.savefig("absorber_attenuation.png", bbox_inches="tight")
    plt.show()
    
# #} end of plot_everything()

for energy_idx,energy in enumerate(E_0_keV):
    print("")
    print("total_cross_section[energy_idx]: {0:1.2f} re^2".format(total_cross_section[energy_idx]/m.pow(constants.r_e, 2)))

    prob = 1 - np.exp(-scatterer_material.electron_density * total_cross_section[energy_idx] * scatterer_z)
    # prob = scatterer_material.electron_density * total_cross_section[energy_idx] * scatterer_z

    print("{0:2.1f}% of photons are scattered for {1:6.1f} keV".format(prob*100, energy))

print("")
print("total area: {0:2.3f} sr".format(total_area))

prob_pe_scatterer = []
prob_cs_scatterer = []

prob_pe_absorber = []
prob_cs_absorber = []

prob_absorber_attenuation = []
prob_scatterer_attenuation = []

pe_energies = []

for e in range(1, 1000, 5): # over keV

    pe_energies.append(e)

    pe_cross_section_scatterer = physics.pe_cs_gavrila_pratt_simplified(scatterer_material, e*1000.0)
    cs_cross_section_scatterer = physics.comptonCrossSection(conversions.energy_ev_to_J(e*1000.0))

    pe_cross_section_absorber = [physics.pe_cs_gavrila_pratt_simplified(mat, e*1000.0) for mat in absorber_material.elements]
    cs_cross_section_absorber = physics.comptonCrossSection(conversions.energy_ev_to_J(e*1000.0))

    prob = 1 - np.exp(-scatterer_material.atomic_density * pe_cross_section_scatterer * scatterer_z)
    prob_pe_scatterer.append(prob)

    prob = 1 - np.exp(-scatterer_material.electron_density * cs_cross_section_scatterer * scatterer_z)
    prob_cs_scatterer.append(prob)

    prod = 1
    for index,element in enumerate(pe_cross_section_absorber):
        prod *= np.exp(-absorber_material.atomic_density * pe_cross_section_absorber[index] * absorber_z)
    prob = 1 - prod
    prob_pe_absorber.append(prob)

    prob = 1 - np.exp(-absorber_material.electron_density * cs_cross_section_absorber * absorber_z)
    prob_cs_absorber.append(prob)

    prob_absorber_attenuation = [(1 - (1-x[0])*(1-x[1])) for x in zip(prob_cs_absorber, prob_pe_absorber)]
    prob_scatterer_attenuation = [(1 - (1-x[0])*(1-x[1])) for x in zip(prob_cs_scatterer, prob_pe_scatterer)]

pid = os.fork()
if pid == 0:
    plot_everything()
