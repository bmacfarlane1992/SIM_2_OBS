'''

outburst_thesis.py

Program to generate a variety of SED plots, primarily for Ch. 5 of thesis.
Sequently, plots generated are:

    - 10 L_sol, unscaled model at i = 30 deg., for models with mass scaled discs
    - As above, for varying envelope mass scaling
    -
All model SEDs computed from emission from inner 10 000 AU.

Author: Benjamin MacFarlane
Date: 03/09/2018
Contact: bmacfarlane@uclan.ac.uk

'''
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - MODULE IMPORTS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
    # Standard Python modules

import numpy as np
import os
import math
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

    # Append location of constants and functions modules into pythonpath

import sys
sys.path.insert(0,"./../")

    # Import local modules

import constants as cs
import functions as fs
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - VARIABLE DEFINITIONS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
dist = 140.0
ymin = 1.e-15 ; ymax = 1.e-7

cwd = os.getcwd()
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
		# # # - - - MAIN PROGRAM - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
    # Generate SEDs, comparing unscaled 10 L_sol, i = 30 deg. model with
    # models where disc mass has been scaled to x0.1, or x2 initial mass

runs_dir = cwd+"/../../runs/OUTBURST/1686/1e4AU/POST_10L_cavity10"
pstar_dir = cwd+"/../../protostar"
plt_dir = cwd+"/../../runs/OUTBURST/plots_analysis"
plt_form = "png"

file = pstar_dir+"/spectrum5934.out"

f = open(file, "r")

wav_star, flam_star, lam_flam_star, \
 nu_star, fnu_star, nu_fnu_star = fs.SED_read(file)

lam_flam_star /= dist**(2.0)

labels = ["_disc01","","_disc2"]
plt_labels = ["Disc x0.1","Unscaled","Disc x2"]
linestyle = ["--","-",":"]

wav = [[] for i in range(3)]
lam_flam = [[] for i in range(3)]

for i in range(3):

    file = runs_dir+labels[i]+"/dat/spectrum5934_ISRF_30i.out"

    wav[i], flam_tmp, lam_flam[i], \
     nu_tmp, fnu_tmp, nu_fnu_tmp = fs.SED_read(file)

    lam_flam[i] /= dist**(2.0)

fig = plt.figure(1)
ax1 = plt.subplot(111)

plt.plot(wav_star, lam_flam_star, color = "g", linewidth = 1, \
 label = "Protostar")

for i in range(3):
    plt.plot(wav[i], lam_flam[i], color = "k", linewidth = 1, \
     linestyle = linestyle[i], label = plt_labels[i])

plt.legend(loc = "upper left", fontsize=10)
plt.xlabel("Wavelength ("+(r"$\mu$m")+")", fontsize = 18, labelpad=0.5)
plt.xticks(fontsize = 15) ;  ax1.set_xscale("log")
ax1.set_xlim(5.e-3,10000.)
ax1.set_yscale("log") ; ax1.set_ylim( ymin , ymax)
plt.ylabel(r"$\lambda$ F$_{\lambda}$ (erg cm$^{-2}$ s$^{-1}$)", \
  fontsize = cs.fontsize, labelpad=0.5)
plt.tight_layout()
plt.savefig(plt_dir+"/sed_discrescale."+plt_form, \
  format = plt_form)
plt.clf()
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
    # Generate SEDs, comparing unscaled 10 L_sol, i = 30 deg. model with
    # models where disc mass has been scaled to x0.1, or x2 initial mass

labels = ["_env01","","_env2"]
plt_labels = ["Envelope x0.1","Unscaled","Envelope x2"]
linestyle = ["--","-",":"]

wav = [[] for i in range(3)]
lam_flam = [[] for i in range(3)]

for i in range(3):

    file = runs_dir+labels[i]+"/dat/spectrum5934_ISRF_30i.out"

    wav[i], flam_tmp, lam_flam[i], \
     nu_tmp, fnu_tmp, nu_fnu_tmp = fs.SED_read(file)

    lam_flam[i] /= dist**(2.0)

fig = plt.figure(1)
ax1 = plt.subplot(111)

plt.plot(wav_star, lam_flam_star, color = "g", linewidth = 1, \
 label = "Protostar")

for i in range(3):
    plt.plot(wav[i], lam_flam[i], color = "k", linewidth = 1, \
     linestyle = linestyle[i], label = plt_labels[i])

plt.legend(loc = "upper left", fontsize=10)
plt.xlabel("Wavelength ("+(r"$\mu$m")+")", fontsize = 18, labelpad=0.5)
plt.xticks(fontsize = 15) ;   ax1.set_xscale("log")
ax1.set_xlim(5.e-3,10000.)
ax1.set_yscale("log") ; ax1.set_ylim( ymin , ymax)
plt.ylabel(r"$\lambda$ F$_{\lambda}$ (erg cm$^{-2}$ s$^{-1}$)", \
  fontsize = cs.fontsize, labelpad=0.5)
plt.tight_layout()
plt.savefig(plt_dir+"/sed_envrescale."+plt_form, \
  format = plt_form)
plt.clf()
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
    # Generate SEDs, comparing unscaled RT models with various values of
    # protostellar luminosity (2, 10, 20, 50, and 80 L_sol). SEDs
    # calculated at i = 30 deg.

wav = [[] for i in range(5)]
lam_flam = [[] for i in range(5)]
files = [ runs_dir+"/../../../1611/1e4AU/POST_cavity10/dat/spectrum3974_ISRF_30i.out", \
 runs_dir+"/dat/spectrum5934_ISRF_30i.out", \
 runs_dir+"/../POST_20L_cavity10/dat/spectrum7056_ISRF_30i.out", \
 runs_dir+"/../POST_50L_cavity10/dat/spectrum8871_ISRF_30i.out", \
 runs_dir+"/../POST_80L_cavity10/dat/spectrum9978_ISRF_30i.out" ]
plt_labels = [r"$L_{*,q} = 2$  L$_\odot$", r"$L_{*,o} = 10$  L$_\odot$", \
 r"$L_{*,o} = 20$  L$_\odot$", r"$L_{*,o} = 50$  L$_\odot$", \
 r"$L_{*,o} = 80$  L$_\odot$"]
colors = ["k","r","g","b","k"]
linestyles = ["--","-","-","-","-"]

for i in range(len(files)):

    file = files[i]

    wav[i], flam_tmp, lam_flam[i], \
     nu_tmp, fnu_tmp, nu_fnu_tmp = fs.SED_read(file)

    lam_flam[i] /= dist**(2.0)

fig = plt.figure(1)
ax1 = plt.subplot(111)

for i in range(len(files)):
    plt.plot(wav[i], lam_flam[i], color = colors[i], linewidth = 1, \
     linestyle = linestyles[i], label = plt_labels[i])

plt.legend(loc = "upper right", fontsize=10)
plt.xlabel("Wavelength ("+(r"$\mu$m")+")", fontsize = 18, labelpad=0.5)
plt.xticks(fontsize = 15) ;   ax1.set_xscale("log")
ax1.set_xlim(1.e-1,10000.)
ax1.set_yscale("log") ; ax1.set_ylim(1.e-12,2.e-7)
plt.ylabel(r"$\lambda$ F$_{\lambda}$ (erg cm$^{-2}$ s$^{-1}$)", \
  fontsize = cs.fontsize, labelpad=0.5)
plt.tight_layout()
plt.savefig(plt_dir+"/sed_lumscale."+plt_form, format = plt_form)
plt.clf()

plt_labels = [r"$L_{*} = 10$  L$_\odot$", r"$L_{*} = 20$  L$_\odot$",\
 r"$L_{*} = 50$  L$_\odot$", r"$L_{*} = 80$  L$_\odot$"]
temps = ["5934","7056","8871","9978"]
colors = ["r","g","b","k"]
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
    # Generate SEDs, comparing models with envelope mass scaled to x0.1
    # intitial value, computed at i = 0 deg. Compare to equivalently mass
    # scaled models, with protostellar luminosities of 10, 20, 50, and 80 L_sol.
    # Compute, and plot the location of the wavelength corresponding to
    # flux weighted mean frequency of each SED. 

wav = [[] for i in range(4)]
lam_flam = [[] for i in range(4)]
mean_lam = []

files = [ runs_dir+"/../POST_10L_cavity10_env01/dat/spectrum5934_ISRF_0i.out", \
 runs_dir+"/../POST_20L_cavity10_env01/dat/spectrum7056_ISRF_0i.out", \
 runs_dir+"/../POST_50L_cavity10_env01/dat/spectrum8871_ISRF_0i.out", \
 runs_dir+"/../POST_80L_cavity10_env01/dat/spectrum9978_ISRF_0i.out" ]

for i in range(len(files)):

    file = files[i]

    wav[i], flam_tmp, lam_flam[i], \
     nu_tmp, fnu_tmp, nu_fnu_tmp = fs.SED_read(file)

    lam_flam[i] /= dist**(2.0)

    T_bol = fs.Tbol_calc(nu_tmp, fnu_tmp, nu_fnu_tmp)
    mean_lam.append(2900. / T_bol)

fig = plt.figure(1)
ax1 = plt.subplot(111)

for i in range(len(files)):
    plt.plot(wav[i], lam_flam[i], color = colors[i], linewidth = 1, \
     linestyle = "-", label = plt_labels[i])
    plt.axvline(x = mean_lam[i], linestyle=":", color = colors[i])

plt.legend(loc = "upper right", fontsize=10)
plt.xlabel("Wavelength ("+(r"$\mu$m")+")", fontsize = 18, labelpad=0.5)
plt.xticks(fontsize = 15) ;   ax1.set_xscale("log")
ax1.set_xlim(1.e-1,10000.)
ax1.set_yscale("log") ; ax1.set_ylim(1.e-12,2.e-7)
plt.ylabel(r"$\lambda$ F$_{\lambda}$ (erg cm$^{-2}$ s$^{-1}$)", \
  fontsize = cs.fontsize, labelpad=0.5)
plt.tight_layout()
plt.savefig(plt_dir+"/tbol_01env."+plt_form, format = plt_form)
plt.clf()
