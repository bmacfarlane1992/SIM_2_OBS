'''

ec53_thesis.py

Script to generate thesis specific plots on the modeling of EC53. In sequence,
the ploted data is:
    - Quiescent models with varying r_core
    - Quiescent models for best fitting r_core, over varying inclinations
    - For best fitting r_core, and inclinations, the varying outburst
        luminosities.

All models have L_bol, T_bol and ratio of sub-mm L_bol vs. L_bol computed. For
all bar outbursting case comparison, SED is overplotted with c2d Spitzer,
Herschel and JCMT data of Dunham+ (2015) and Yoo+ (2017)

Author: Benjamin MacFarlane
Date: 20/08/2018
Contact: bmacfarlane@uclan.ac.uk

'''
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - MODULE IMPORTS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Standard Python modules
import os
import shutil
import numpy as np
import math
import random
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

    # Append location of constants and functions modules into pythonpath

import sys
sys.path.insert(0,"./../")

    # Import local modules

import constants as cs
import functions as fs
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - VARIABLE DEFINITIONS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
dist = 436.0        # Float: Distance (pc) to source
best_cavity = 20
best_inclin = 30         # Int. Designated inclination ("best" fit)
best_rcore = "1e4"       # Str. Designated core radius ("best" fit)
#
bolometrics = False     # Bool.: Are bolometric propertes of SEDs computed?
plt_form = "png"    # Str.: Format of output plots ["png","eps"]
#
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
		# # # - - - MAIN PROGRAM - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

cwd = os.getcwd()
run_dir = cwd + "/../../runs/EC53"
plt_dir = run_dir + "/plots_analysis"
#
    # Before data read, store EC 53 photometry data
#
f = run_dir + "/EC53_D15_2876_flux.txt"
ec53_wav, ec53_fnu, ec53_err = np.loadtxt(f, skiprows=1, unpack=True)
ec53_fnu *= 0.001 * cs.Jy_cgs * (cs.c_cgs / (ec53_wav*cs.cm_per_micron))
ec53_err *= 0.001 * cs.Jy_cgs * (cs.c_cgs / (ec53_wav*cs.cm_per_micron))


### ------------------------------------------------------------------------ ###


print("\nComparing SEDs for R_CORE = "+best_rcore+", with varying inclinations")
print("and cavity opening angles\n")
#
wav = [[] for i in range(4)]
flam = [[] for i in range(4)]
lam_flam = [[] for i in range(4)]
#
cwd = os.getcwd()
dat_dir = run_dir + "/dat_" + best_rcore + "_cavity" + str(best_cavity)
#
inclin = [0,30,60,90]
#
leg_names = ["i = 0 deg.", \
  "i = 30 deg.", \
  "i = 60 deg.", \
  "i = 90 deg."]
color = ["r","g","b","k"]
linestyle = ["-","-","-","-"]
linewidth = [1,1,1,1]

for i in range(4):

    # Read SED data

    fsed = dat_dir+"/spectrum6L_ISRF_"+str(inclin[i])+"i.out"

    wav[i], flam[i], lam_flam[i], nu, fnu, nu_fnu = fs.SED_read(fsed)

    for j in range(len(lam_flam[i])):
        lam_flam[i][j] /= dist**(2.0)

    if bolometrics:

        # Model bolometric luminosity

        L_bol = fs.Lbol_calc(nu, fnu)

        # Model bolometric temperature

        T_bol = fs.Tbol_calc(nu,fnu, nu_fnu)

        # Model sub-mm to bolometric luminosity ratio

        L_ratio = fs.L_ratio_calc(nu, fnu)

        # Print results to terminal

        print("\nFor theta = {0} deg., {1} model: \n".format(best_cavity, \
         leg_names[i]))

        print("Bolometric luminosity is: {0} L_sol\n".format(L_bol))
        print("Ratio of Bolometric and sub-mm "+ \
          "luminosity is {0} %\n".format(L_ratio) )
        print("Bolometric temperature is: {0} K\n".format(T_bol))

fig = plt.figure(1)
ax1 = plt.subplot(111)
for i in range(len(inclin)):
    plt.plot( wav[i], lam_flam[i], label = leg_names[i], color=color[i], \
       linestyle = linestyle[i], linewidth = linewidth[i])
plt.errorbar(ec53_wav, ec53_fnu, yerr=ec53_err, \
 fmt='o', mfc='k', ecolor='k', ms=5)
ymax = max(lam_flam[0])*10. ; ymin = 8.0e-14
plt.xlabel("Wavelength ("+(r"$\mu$m")+")", fontsize = 18, labelpad=0.5)
plt.xticks(fontsize = 15) ;   ax1.set_xscale("log")
ax1.set_xlim(8.e-2,2500.)
plt.ylabel(r"$\lambda$ F$_{\lambda}$", fontsize = 18, labelpad=0.5)
ax1.set_yscale("log") ; ax1.set_ylim( ymin, ymax )
plt.legend(loc = "upper left", fontsize=cs.leg_fontsize)
plt.tight_layout()
plt.savefig(plt_dir + "/incl_ec53_"+str(best_rcore)+ \
 "cavity"+str(best_cavity)+"."+plt_form, format = plt_form)
plt.clf()

### ------------------------------------------------------------------------ ###

print("\nComparing SEDs for radially varying models (0 deg. inclination)\n")

wav = [[] for i in range(3)]
flam = [[] for i in range(3)]
lam_flam = [[] for i in range(3)]
#
r_core = ["1e4","3e4","5e4"]
inclin = 0
#
leg_names = [r"$R_{CORE} = 10 000$ AU", \
  r"$R_{CORE} = 30 000$ AU", \
  r"$R_{CORE} = 50 000$ AU"]
color = ["r","g","b"]
linestyle = ["-","-","-"]
linewidth = [1,1,1]

for i in range(3):

    # Read SED data

    fsed = run_dir + "/dat_" + r_core[i] + "_cavity"+str(best_cavity)+ \
     "/spectrum6L_ISRF_"+str(inclin)+"i.out"

    wav[i], flam[i], lam_flam[i], nu, fnu, nu_fnu = fs.SED_read(fsed)

    for j in range(len(lam_flam[i])):
        lam_flam[i][j] /= dist**(2.0)

    if bolometrics:

        # Model bolometric luminosity

        L_bol = fs.Lbol_calc(nu, fnu)

        # Model bolometric temperature

        T_bol = fs.Tbol_calc(nu,fnu, nu_fnu)

        # Model sub-mm to bolometric luminosity ratio

        L_ratio = fs.L_ratio_calc(nu, fnu)

        # Print results to terminal

        print("\nFor theta = {0} deg., {1} model: \n".format(best_cavity, \
         leg_names[i]))

        print("Bolometric luminosity is: {0} L_sol\n".format(L_bol))
        print("Ratio of Bolometric and sub-mm "+ \
          "luminosity is {0} %\n".format(L_ratio) )
        print("Bolometric temperature is: {0} K\n".format(T_bol))

fig = plt.figure(1)
ax1 = plt.subplot(111)
for i in range(3):
    plt.plot( wav[i], lam_flam[i], label = leg_names[i], color=color[i], \
       linestyle = linestyle[i], linewidth = linewidth[i])
plt.errorbar(ec53_wav, ec53_fnu, yerr=ec53_err, \
 fmt='o', mfc='k', ecolor='k', ms=5)
ymax = max(lam_flam[0])*10. ; ymin = 8.0e-14
plt.xlabel("Wavelength ("+(r"$\mu$m")+")", fontsize = 18, labelpad=0.5)
plt.xticks(fontsize = 15) ;   ax1.set_xscale("log")
ax1.set_xlim(8.e-2,2500.)
plt.ylabel(r"$\lambda$ F$_{\lambda}$", fontsize = 18, labelpad=0.5)
ax1.set_yscale("log") ; ax1.set_ylim( ymin, ymax )
plt.legend(loc = "upper left", fontsize=cs.leg_fontsize)
plt.tight_layout()
plt.savefig(plt_dir + "/rad_ec53_cavity"+str(best_cavity)+"."+plt_form, \
 format = plt_form)
plt.clf()


### ------------------------------------------------------------------------ ###

print("\nComparing SEDs with varying outburst luminosity, in models where ")
print("R_CORE  = "+best_rcore+", theta = "+str(best_cavity)+ \
 "deg., and i = "+str(best_inclin)+"deg.")#
#
cwd = os.getcwd()
dat_dir = run_dir + "/dat_"+best_rcore+"_cavity"+str(best_cavity)
#
lums = [10,20,25,30]
leg_names = [r"$L_{*,q} = 6 $ L$_\odot$", \
  r"$L_{*,o} = 10 $ L$_\odot$", \
  r"$L_{*,o} = 20 $ L$_\odot$", \
  r"$L_{*,o} = 25 $ L$_\odot$", \
  r"$L_{*,o} = 30 $ L$_\odot$"]
color = ["k","r","g","b","cyan"]
linestyle = ["--","-","-","-","-"]
linewidth = [1,1,1,1,1]

wav = [[] for i in range(len(lums)+1)]
flam = [[] for i in range(len(lums)+1)]
lam_flam = [[] for i in range(len(lums)+1)]

for i in range(len(lums)+1):

    # Read SED data

    if (i == 0):
        fsed = dat_dir + "/spectrum6L_ISRF_"+str(best_inclin)+"i.out"
    else:
        fsed = dat_dir + "/spectrum"+str(lums[i-1])+"L_ISRF_"+ \
         str(best_inclin)+"i.out"

    wav[i], flam[i], lam_flam[i], nu, fnu, nu_fnu = fs.SED_read(fsed)

    for j in range(len(lam_flam[i])):
        lam_flam[i][j] /= dist**(2.0)

    if bolometrics:

        # Model bolometric luminosity

        L_bol = fs.Lbol_calc(nu, fnu)

        # Model bolometric temperature

        T_bol = fs.Tbol_calc(nu,fnu, nu_fnu)

        # Model sub-mm to bolometric luminosity ratio

        L_ratio = fs.L_ratio_calc(nu, fnu)

        # Print results to terminal

        print("\nFor theta = {0} deg., {1} model: \n".format(best_cavity, \
         leg_names[i]))

        print("Bolometric luminosity is: {0} L_sol\n".format(L_bol))
        print("Ratio of Bolometric and sub-mm "+ \
          "luminosity is {0} %\n".format(L_ratio) )
        print("Bolometric temperature is: {0} K\n".format(T_bol))

flam_int_q = interp1d(wav[0], flam[0])
for i in range(len(lums)):
    flam_int_o= interp1d(wav[i+1], flam[i+1])
    ratio850 = str(round( ( flam_int_o(850) / flam_int_q(850) ), 2))
    leg_names[i+1] += r" ; F$_{\lambda, o}$ / F$_{\lambda, q}$ = "+ratio850

fig = plt.figure(1)
ax1 = plt.subplot(111)
for i in range(len(lums)+1):
    plt.plot( wav[i], lam_flam[i], label = leg_names[i], color=color[i], \
       linestyle = linestyle[i], linewidth = linewidth[i])
ymax = max(lam_flam[3])*10. ; ymin = 8.0e-14
plt.xlabel("Wavelength ("+(r"$\mu$m")+")", fontsize = 18, labelpad=0.5)
plt.xticks(fontsize = 15) ;   ax1.set_xscale("log")
ax1.set_xlim(8.e-2,2500.)
plt.ylabel(r"$\lambda$ F$_{\lambda}$", fontsize = 18, labelpad=0.5)
ax1.set_yscale("log") ; ax1.set_ylim( ymin, ymax )
plt.axvline(x=850.0, linestyle="dashed", color="k")
plt.legend(loc = "upper left", fontsize=cs.leg_fontsize-2)
plt.tight_layout()
plt.savefig(plt_dir + "/burst_ec53."+plt_form, format = plt_form)
plt.clf()
