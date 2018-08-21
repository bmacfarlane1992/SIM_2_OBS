'''

evolution.py

Program to read in SED data for EVOLUTION models, and compare between
    - Single snapshots for multiple inclinations
    - Multiple snapshots, for inclination limits

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
import numpy as np
import math
import matplotlib.pyplot as plt

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
snaps = ['00101','00601','01101','01601','02101']   # Str. Arr.: DE05 tags
snaps_leg = ["79 kyr","84 kyr","89 kyr","94 kyr","99 kyr"] # Str. Arr.: Snap t's
inclins = [0,45,60,90]      # Int. Arr.: Inclinations used for snapshots
#
sed_norm = True         # Bool: Normalise to stellar flux, or to distance?
dist = 140.0            # Float: Distance to source (pc)
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
		# # # - - - MAIN PROGRAM - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
cwd = os.getcwd()

run_dir = cwd + "/../../runs/EVOLUTION"
plt_dir = run_dir + "/plots_analysis"

if not os.path.isdir(plt_dir):
    os.system('mkdir '+plt_dir)
#
t_star = [] ; r_star = []
wavs = [ [] for i in range(len(snaps))]
flam = [ [ [] for i in range(len(inclins))] for i in range(len(snaps))]
wavs_s = [ [] for i in range(len(snaps))]
flam_s = [ [] for i in range(len(snaps))]
#
leg_cols = ['r','g','b','c']
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
for i in range(len(snaps)):
#
    # Read snapshot specific protostellar properties/SED
#
    fname = run_dir + "/" + snaps[i] + "/dat/stars.inp"
    f = open(fname, 'r')
    trash = f.readline() ; trash = f.readline()
    r_star.append(float(f.readline().split(' ')[0]) * cs.rsol_per_cm)
    t_star.append(-float(f.readlines()[-1]))
    f.close()
#
    fname = run_dir + "/" + snaps[i] + "/dat/spectrum"+str(int(t_star[i]))+".out"

    wavs_s[i], flam_s[i], lam_flam, nu, fnu, nu_fnu = fs.SED_read(fname)

    if sed_norm:
        fnorm = ( (r_star[i] * cs.cm_per_rsol)**2. * cs.sigma_sb_cgs * \
           (t_star[i])**4. ) / (cs.cm_per_pc)**2.
        for j in range(len(flam_s[i])):
            flam_s[i][j] /= fnorm
    else:
        fnorm = dist**2.0
        for j in range(len(flam_s[i])):
            flam_s[i][j] /= fnorm
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
    # Loop over inclinations, and read in snapshot specific SEDs
#
    for j in range(len(inclins)):
#
        fname = run_dir+"/"+snaps[i]+"/dat/spectrum_"+str(inclins[j])+"i.out"

        wavs[i], flam[i][j], lam_flam, nu, fnu, nu_fnu = fs.SED_read(fname)

        if sed_norm:
            fnorm = ( (r_star[i] * cs.cm_per_rsol)**2. * cs.sigma_sb_cgs * \
               (t_star[i])**4. ) / (cs.cm_per_pc)**2.
            for k in range(len(flam[i][j])):
                flam[i][j][k] /= fnorm
        else:
            fnorm = dist**2.0
            for k in range(len(flam[i][j])):
                flam[i][j][k] /= fnorm
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
    # Plot sinlge snaphot SED, for multiple inclinations
#
for i in range(len(snaps)):
    fig = plt.figure(1)
    ax1 = plt.subplot(111)
    plt.plot( wavs_s[i], flam_s[i], \
       label = "Protostar: "+str(int(t_star[i]))+" K", \
       color = 'k', linestyle = ":")
    for j in range(len(inclins)):
        plt.plot( wavs[i], flam[i][j], label = str(inclins[j])+" deg.", \
           color=leg_cols[j] )
    if sed_norm:
        ymax = 25.0 ; ymin = 1.e-12
        plt.ylabel(r'($\lambda$ F$_{\lambda}$) / F$_{*}$', \
           fontsize = 18, labelpad=0.5)
    else:
        ymax = max(flam[-1][0])*10. ; ymin = ymax * 1.e-9
        plt.ylabel(r'$\lambda$ F$_{\lambda}$', fontsize = 18, labelpad=0.5)
    plt.xlabel('Wavelength ('+(r'$\mu$m')+')', fontsize = 18, labelpad=0.5)
    plt.xticks(fontsize = 15) ;   ax1.set_xscale('log')
    ax1.set_xlim(1.e-1,10000.)
    ax1.set_yscale('log') ; ax1.set_ylim( ymin, ymax )
    plt.legend(loc = 'upper right', fontsize=8, scatterpoints = 20)
    plt.savefig(plt_dir + "/SED_"+snaps[i]+".png") ; plt.clf()
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
    # Plot limit snapshots, for limit inclinations
#
fig = plt.figure(1)
ax1 = plt.subplot(111)
#
plt.plot( wavs_s[0], flam_s[0], \
   label = "Protostar: "+str(int(t_star[0]))+" K", \
   color = leg_cols[0], linestyle = ":")
plt.plot( wavs[0], flam[0][0], \
   label = snaps_leg[0]+" "+str(inclins[0])+" deg.", color=leg_cols[0], \
   linestyle = "-" )
plt.plot( wavs[0], flam[0][-1], \
   label = snaps_leg[0]+" "+str(inclins[-1])+" deg.", color=leg_cols[0], \
   linestyle = "--" )
#
plt.plot( wavs_s[-1], flam_s[-1], \
   label = "Protostar: "+str(int(t_star[-1]))+" K", \
   color = leg_cols[-1], linestyle = ":")
plt.plot( wavs[-1], flam[-1][0], \
   label = snaps_leg[-1]+" "+str(inclins[0])+" deg.", color=leg_cols[-1], \
   linestyle = "-" )
plt.plot( wavs[-1], flam[-1][-1], \
   label = snaps_leg[-1]+" "+str(inclins[-1])+" deg.", color=leg_cols[-1], \
   linestyle = "--" )
if sed_norm:
    ymax = 25.0 ; ymin = 1.e-12
    plt.ylabel(r'($\lambda$ F$_{\lambda}$) / F$_{*}$', \
       fontsize = 18, labelpad=0.5)
else:
    ymax = max(flam[-1][0])*10. ; ymin = ymax * 1.e-9
    plt.ylabel(r'$\lambda$ F$_{\lambda}$', fontsize = 18, labelpad=0.5)
plt.xlabel('Wavelength ('+(r'$\mu$m')+')', fontsize = 18, labelpad=0.5)
plt.xticks(fontsize = 15) ;   ax1.set_xscale('log')
ax1.set_xlim(1.e-1,10000.)
ax1.set_yscale('log') ; ax1.set_ylim( ymin, ymax )
plt.legend(loc = 'upper right', fontsize=8, scatterpoints = 20)
plt.savefig(plt_dir + "/SED_compare.png") ; plt.clf()
