'''

ec53.py

Routine to check both SED and temperature profiles of two models with
different protostellar luminosity, with/without ISRF in RT calculation. Also
compares to RT run with ISRF only. Uses data to then evaluate flux ratio
at 850 micron.

SEREN snapshot DE05...00788... used, which is scaled to EC 53 envelope mass of
1.84 M_sol (Enoch et al., 2009). Simulation with RT model parameterised by
protostellar temperature defined (in K) as XXXX in spectrumXXXX.out files.
Data with ISRF of Andre et al. (2003) denoted by _ISRF tag.

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
inclin = 30          # Int. Inclination being probed by SED analyses
lums = [6, 20]
rcore = "4e4"
#
plt_form = "png"    # Str.: Format of output plots ["png","eps"]
#
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
		# # # - - - MAIN PROGRAM - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =


print("\nRUNNING ANALYSES FOR EC 53 MODELLING\n")
#
names = [str(lums[0])+"L", str(lums[0])+"L_ISRF", \
  str(lums[1])+"L", str(lums[1])+"L_ISRF"]
#
leg_names = [str(lums[0])+r" L$_\odot$ (Quiescent)", \
  str(lums[0])+r" L$_\odot$ + ISRF (Quiescent )", \
  str(lums[1])+r" L$_\odot$ (Outbursting)", \
  str(lums[1])+r" L$_\odot$ + ISRF (Outbursting)"]
#
cwd = os.getcwd()
run_dir = cwd + "/../../runs/EC53"
dat_dir = run_dir + "/dat_" + rcore
plt_dir = run_dir + "/plots_analysis/" + rcore

if not os.path.isdir(plt_dir):
    os.makedirs(plt_dir)

linewidth = [1,1,2,2]
s = [0.25,0.25,0.5,0.5]
color = ["b","b","r","r"]
linestyle = ["-",":","-",":"]
#
    # Before data read, store EC 53 photometry data
#
f = run_dir + "/EC53_D15_2876_flux.txt"
ec53_wav, ec53_fnu, ec53_err = np.loadtxt(f, skiprows=1, unpack=True)
ec53_fnu *= 0.001 * cs.Jy_cgs * (cs.c_cgs / (ec53_wav*cs.cm_per_micron))
ec53_err *= 0.001 * cs.Jy_cgs * (cs.c_cgs / (ec53_wav*cs.cm_per_micron))

wav = [[] for i in range(4)]
flam = [[] for i in range(4)]
lam_flam = [[] for i in range(4)]
#
    # First, read and store SED data for PRE-, then POST- cut model
#
for i in range(len(names)):

    fsed = dat_dir + "/spectrum"+names[i]+"_"+str(inclin)+"i.out"

    wav[i], flam[i], lam_flam[i], nu, fnu, nu_fnu = fs.SED_read(fsed)
    for j in range(len(lam_flam[i])):
        lam_flam[i][j] /= dist**(2.0)


### ------------------------------------------------------------------------ ###

    # Model bolometric luminosity

    L_bol = fs.Lbol_calc(nu, fnu)

    # Model bolometric temperature

    T_bol = fs.Tbol_calc(nu, fnu, nu_fnu)

    # Model sub-mm to bolometric luminosity ratio

    L_ratio = fs.L_ratio_calc(nu, fnu)

    # Print results to terminal

    print("\nFor theta = {0} deg., {1} model: \n".format(best_cavity, \
     leg_names[i]))

    print("Bolometric luminosity is: {0} L_sol\n".format(L_bol))
    print("Ratio of Bolometric and sub-mm "+ \
      "luminosity is {0} %\n".format(L_ratio) )
    print("Bolometric temperature is: {0} K\n".format(T_bol))


### ------------------------------------------------------------------------ ###


#   Multiple cut method analysis:
#
    # Now plot SED. Zoom into SED in region of interest, using second panel
#
# SED
#
fig = plt.figure(1)
ax1 = plt.subplot(111)
for i in range(4):
    plt.plot( wav[i], lam_flam[i], label = leg_names[i], color=color[i], \
       linestyle = linestyle[i], linewidth = linewidth[i])
plt.errorbar(ec53_wav, ec53_fnu, yerr=ec53_err, \
 fmt='o', mfc='k', ecolor='k', ms=5)
ymax = max(lam_flam[3])*10. ; ymin = 8.0e-14
plt.xlabel("Wavelength ("+(r"$\mu$m")+")", fontsize = 18, labelpad=0.5)
plt.xticks(fontsize = 15) ;   ax1.set_xscale("log")
ax1.set_xlim(8.e-1,2500.)
plt.ylabel(r"$\lambda$ F$_{\lambda}$", fontsize = 18, labelpad=0.5)
ax1.set_yscale("log") ; ax1.set_ylim( ymin, ymax )
#plt.axvline(x=850.0, linestyle="dashed", color="k")
plt.legend(loc = "upper left", fontsize=8)
plt.tight_layout()
plt.savefig(plt_dir + "/SED_"+str(inclin)+"i."+plt_form, format = plt_form)
plt.clf()
#
fig = plt.figure(1)
ax1 = plt.subplot(111)
for i in range(4):
    plt.plot( wav[i], lam_flam[i], label = leg_names[i], color=color[i], \
       linestyle = linestyle[i], linewidth = linewidth[i] )
ymax = 2.00e-9 ; ymin = 2.00e-12
plt.xlabel("Wavelength ("+(r"$\mu$m")+")", fontsize = 18, labelpad=0.5)
plt.xticks(fontsize = 15) ;   ax1.set_xscale("log")
ax1.set_xlim(300.,1000.)
plt.ylabel(r"$\lambda$ F$_{\lambda}$", fontsize = 18, labelpad=0.5)
ax1.set_yscale("log") ; ax1.set_ylim( ymin, ymax )
plt.axvline(x=450.0, linestyle="dashed", color="k")
plt.axvline(x=850.0, linestyle="dashed", color="k")
plt.legend(loc = "upper right", fontsize=10)
plt.tight_layout()
plt.savefig(plt_dir + "/SED_ZOOM._"+str(inclin)+"i."+plt_form, \
  format = plt_form)
plt.clf()

### ------------------------------------------------------------------------ ###


#   Singular cut method analysis:
#
#   450 and 850 micron flux ratio calculation (lambda F lambda)
#
L_ratio = float(lums[1]) / float(lums[0])
print("Ratio of protostellar luminosities, assuming constant "+ \
  "radius is: {0}\n".format(L_ratio) )

print("\nComputing (lambda F_lambda) ratio for "+ \
  "{0} and {1}\n".format(leg_names[0], leg_names[2]) )
#
flam_int1 = interp1d(wav[0], flam[0])
flam_int2 = interp1d(wav[2], flam[2])
ratio450 = round( ( flam_int2(450) / flam_int1(450) ), 3)
print("Ratio 450 micron flux is: {0}\n".format(ratio450) )
ratio850 = round( ( flam_int2(850) / flam_int1(850) ), 3)
print("Ratio 850 micron flux is: {0}\n".format(ratio850) )

print("\nComputing (lambda F_lambda) ratio for "+ \
  "{0} and {1}\n".format(leg_names[1], leg_names[3]) )
#
flam_int1 = interp1d(wav[1], flam[1])
flam_int2 = interp1d(wav[3], flam[3])
ratio450 = round( ( flam_int2(450) / flam_int1(450) ), 3)
print("Ratio 450 micron flux is: {0}\n".format(ratio450) )
ratio850 = round( ( flam_int2(850) / flam_int1(850) ), 3)
print("Ratio 850 micron flux is: {0}\n".format(ratio850) )

### ------------------------------------------------------------------------ ###

#
    # Compare the beam-dependent flux for E1 snapshots with ISRF heating
#

f_beam = dat_dir+"/MULTI/lam_flam_multi6L.dat"
f_beam_isrf = dat_dir+"/MULTI/lam_flam_multi6L_ISRF.dat"

if os.path.isfile(f_beam) and os.path.isfile(f_beam_isrf):

    r_beam1, lam_flam_a1, lam_flam_b1 = np.loadtxt(f_beam, unpack=True)
    r_beam2, lam_flam_a2, lam_flam_b2 = np.loadtxt(f_beam_isrf, unpack=True)

    f_beam = dat_dir+"/MULTI/lam_flam_multi20L.dat"
    f_beam_isrf = dat_dir+"/MULTI/lam_flam_multi20L_ISRF.dat"

    r_beam3, lam_flam_a3, lam_flam_b3  = np.loadtxt(f_beam, unpack=True)
    r_beam4, lam_flam_a4, lam_flam_b4  = np.loadtxt(f_beam_isrf, unpack=True)

    ratio_a = lam_flam_a3 / lam_flam_a1
    ratio_isrf_a = lam_flam_a4 / lam_flam_a2

    ratio_b = lam_flam_b3 / lam_flam_b1
    ratio_isrf_b = lam_flam_b4 / lam_flam_b2

    cols = ["g","k"]
    lines = ["-",":"]
    phase = ["No ISRF","ISRF"]

    fig = plt.figure(1)
    ax1 = plt.subplot(111)
    plt.plot(r_beam1, ratio_a, color = cols[0], linestyle = lines[0], \
     label = r"450 $\mu$m ("+phase[0]+")")
    plt.plot(r_beam1, ratio_isrf_a, color = cols[0], linestyle = lines[1], \
     label = r"450 $\mu$m ("+phase[1]+")")
    plt.plot(r_beam1, ratio_b, color = cols[1], linestyle = lines[0], \
     label = r"850 $\mu$m ("+phase[0]+")")
    plt.plot(r_beam1, ratio_isrf_b, color = cols[1], linestyle = lines[1], \
     label = r"850 $\mu$m ("+phase[1]+")")

    plt.xlabel("Radius (AU)", fontsize = cs.fontsize, labelpad=0.5)
    ax1.set_xlim(0, max(r_beam1))
    ax1.set_ylim(1.4, 1.8)
    ax1.axvline(1580, linestyle="dashed", color="k")
    ax1.axvline(6060, linestyle="dashed", color="k")
    plt.ylabel(r"F$_{\lambda, o}$ / F$_{\lambda, q}$", \
     fontsize = cs.fontsize, labelpad=0.5)
    plt.legend(loc = "upper left", fontsize=cs.leg_fontsize)
    plt.tight_layout()
    plt.savefig(plt_dir+"/lamflam_beam"+"."+plt_form, format = plt_form)
    plt.clf()

    global lum_rcore, wav_probe

    max_r = lum_rcore

### ------------------------------------------------------------------------ ###


    # Estimate luminosity from emissivity of tau = 1 -> domain edge

os.chdir(dat_dir)

    # First, read in temeperature data

f = open(dat_dir+"/dust_temperature6L.dat","r")
trash = f.readline()
ngrid = int(f.readline())
nrspecies = int(f.readline())
temp = [[] for s in range(nrspecies)]
for s in range(nrspecies):
    for j in range(ngrid):
        temp[s].append(float(f.readline()))
f.close()
temp = np.array(temp)

    # Now, read loc_grid.dat, finding maximum radial limits

xloc, yloc, zloc, rloc, xwid = \
 np.loadtxt(dat_dir+"/loc_grid.dat", unpack = True)
rloc *= cs.au_per_cm

    # Copy data file under inspection to radmc3d input, ready for tau test

shutil.copy2(dat_dir+"/dust_temperature6L.dat", dat_dir+"/dust_temperature.dat")

    # Loop over inclinations

sdat_ext = "6L_"+str(inclin)+"i.out"

    # Now run tau test

file = open(dat_dir+"/spectrum"+sdat_ext)
wav, flam, lam_flam, nu, fnu, nu_fnu = fs.SED_read(file)

wav_probe = np.average(wav, weights = lam_flam)

exec_loc = ""
max_r = max(rloc)
npix = 1000
os.system('''{0}radmc3d tausurf 1 lambda {1} incl {2} \\
    zoomau {3} {4} {3} {4} npix {5} secondorder'''.format(exec_loc, \
    wav_probe, inclin, - max_r, max_r, npix) )

    # Find average radial location from resultant tausurface_3d

r_im  = []
f = open(dat_dir+"/image.out","r")
for l in range(6):
    trash = f.readline()
for y in range(npix):
    yloc = -max_r + (max_r / npix) + y*(2*max_r / npix)

    for x in range(npix):

        xloc = -max_r + (max_r / npix) + x*(2*max_r / npix)

        z_tau = float(f.readline())

        if (z_tau == -1.e91):
            continue
        else:
            r_im.append( np.sqrt(xloc**2.0 + yloc**2.0 + \
             (z_tau * cs.au_per_cm)**2.0) )
f.close()
r_tau = np.mean(r_im)

    # Now compute average temperature of regions from
    # tau = 1 surface to domain limits

count = 0.0 ; temp_ave = 0.0
for s in range(nrspecies):
    for j in range(ngrid):
        if (rloc[j] > r_tau) and (rloc[j] < max_r):
            count += 1.0
            temp_ave += temp[s][j]
temp_ave /= count

    # With radial limits, and average temperature, estimate
    # luminosity of emissive envelope

L = ( (max_r - r_tau) * cs.m_per_au )**2.0 * (temp_ave**4.0)
L /= ( (cs.m_per_rsol)**2.0 * (5780.0)**4.0 )


print("\n Luminosity of i={0} deg. model is {1} L_sol\n".format(inclin, L))
