'''

alma_proposal.py

Routine to check both SED and temperature profiles of two models with
different protostellar luminosity, with/without ISRF in RT calculation. Also
compares to RT run with ISRF only. Uses data to then evaluate flux ratio
at 850 micron.

SEREN snapshot DE05...00788... used, which is scaled to EC 53 envelope mass of
1.84 M_sol (Enoch et al., 2009). Simulation with RT model parameterised by
protostellar temperature defined (in K) as XXXX in spectrumXXXX.out files.
Data with ISRF of Andre et al. (2003) denoted by _ISRF tag.

Author: Benjamin MacFarlane
Date: 14/02/2018
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
import random
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import constants as cs
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - VARIABLE DEFINITIONS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
restricted = "POST"
#
dist = 436.0        # Float: Distance (pc) to source
#
plt_form = "png"    # Str.: Format of output plots ["png","eps"]
#
isrf = False
#
phase = ["Quiescent","Outbursting"]
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
		# # # - - - MAIN PROGRAM - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =


l_tags = ["6L","20L"]
restricted = ["POST"]
isrf_ext = ["_0ISRF","_ISRF"]

wav = [ [ [ [] for i in range(len(isrf_ext)) ] for i in range(len(restricted)) ] \
  for i in range(2) ]
flam = [ [ [ [] for i in range(len(isrf_ext)) ] for i in range(len(restricted)) ] \
  for i in range(2) ]

for i in range(len(restricted)):
    for j in range(len(isrf_ext)):
#
        l_star = []
        for t in range(len(temps)):
            l = 4.0 * math.pi * (cs.m_per_rsol)**2.0 * cs.sigma_sb * \
              float(temps[t])**4.0
            lsol = l / cs.Lsol_SI
            l_star.append( np.round( lsol, 2) )

        names = [str(temps[0])+isrf_ext[j], str(temps[1])+isrf_ext[j]]
#
        leg_names = [r"$L_*$ = "+str(l_star[0])+r" $L_\odot$ (Quiescent)", \
          r"$L_*$ = "+str(l_star[1])+r" $L_\odot$ (Outbursting)" ]
#
        cwd = os.getcwd()
        run_dir = cwd+"/../.."
        if (restricted[i] == "POST"):
            dat_dir = run_dir + "/POST/dat"
            plt_dir = cwd + "/../plots/POST"
        elif (restricted[i] == "POST_1000"):
            dat_dir = run_dir + "/POST_1000/dat"
            plt_dir = cwd + "/../plots/POST_1000"
        if not os.path.isdir(plt_dir):
            os.system("mkdir {0}".format(plt_dir))
#
### ------------------------------------------------------------------------ ###

#
    # Read in SED data
#
        for t in range(len(temps)):
            fsed = dat_dir + "/spectrum"+names[t]+"_0i.out"
#
            f = open(fsed,"r")
            for lines in range(3):
                header = f.readline()
            for lines in f:
                lines = lines.strip() ; columns = lines.split()
                wav[t][i][j].append(float(columns[0]))
                flam[t][i][j].append( float(columns[1]) * \
                  ( cs.c_cgs / (float(columns[0])*cs.cm_per_micron) ) / \
                  dist**(2.0))
            f.close()

### ------------------------------------------------------------------------ ###


#   Singular cut method analysis:
#
    # Now plot SED. Zoom into SED in region of interest, using second panel
#
# SED
#
        color = ["b","r","g","k","c"]
#
        fig = plt.figure(1)
        ax1 = plt.subplot(111)
        for t in range(len(temps)):
            plt.plot( wav[t][i][j], flam[t][i][j], label = leg_names[t], \
              color=color[t] )
        ymax = max(flam[1][0][0])*10. ; ymin = 2.0e-19
        plt.xlabel("Wavelength ("+(r"$\mu$m")+")", fontsize = 18, labelpad=0.5)
        plt.xticks(fontsize = 15) ;   ax1.set_xscale("log")
        ax1.set_xlim(1.e-0,10000.)
        plt.ylabel(r"$\lambda$ F$_{\lambda}$ (erg/s/cm$^2$)", \
          fontsize = 18, labelpad=0.5)
        ax1.set_yscale("log") ; ax1.set_ylim( ymin, ymax )
        plt.axvline(x=850.0, linestyle="dashed", color="k")
        plt.axvline(x=450.0, linestyle="dashed", color ="k")
        plt.legend(loc = "lower center", fontsize=12, scatterpoints = 20)
        plt.tight_layout()
        plt.savefig(plt_dir + "/SED"+str(isrf_ext[j])+"."+plt_form, \
          format = plt_form)
        plt.clf()
#
        fig = plt.figure(1)
        ax1 = plt.subplot(111)
        for t in range(len(temps)):
            plt.plot( wav[t][i][j], flam[t][i][j], label = leg_names[t], \
              color=color[t] )
        ymax = 1.00e-10 ; ymin = 1.00e-13
        plt.xlabel("Wavelength ("+(r"$\mu$m")+")", fontsize = 18, labelpad=0.5)
        plt.xticks(fontsize = 15) ;   ax1.set_xscale("log")
        ax1.set_xlim(300.,1200.)
        plt.ylabel(r"$\lambda$ F$_{\lambda}$ (erg/s/cm$^2$)", \
          fontsize = 18, labelpad=0.5)
        ax1.set_yscale("log") ; ax1.set_ylim( ymin, ymax )
        plt.axvline(x=850.0, linestyle="dashed", color="k")
        plt.axvline(x=450.0, linestyle="dashed", color ="k")
        plt.tight_layout()
        plt.savefig(plt_dir + "/SED_ZOOM"+str(isrf_ext[j])+"."+plt_form, \
          format = plt_form)
        plt.clf()

    # Plot all ISRF models for single restriction in both full, and zoom SEDs

    lines = ["-","--",":"]
    leg_names = ["No ISRF","Andre+ (2003) ISRF", "10x Andre+ (2003) ISRF"]

    fig = plt.figure(1)
    ax1 = plt.subplot(111)

    for j in range(len(isrf_ext)):
        for t in range(len(temps)):
            plt.plot( wav[t][i][j], flam[t][i][j], \
              label = phase[t] + ": " + leg_names[j], color=color[t], \
              linestyle = lines[j])

    ymax = max(flam[1][0][0])*10. ; ymin = 2.0e-19
    plt.xlabel("Wavelength ("+(r"$\mu$m")+")", fontsize = 18, labelpad=0.5)
    plt.xticks(fontsize = 15) ;   ax1.set_xscale("log")
    ax1.set_xlim(1.e-0,10000.)
    plt.ylabel(r"$\lambda$ F$_{\lambda}$ (erg/s/cm$^2$)", \
      fontsize = 18, labelpad=0.5)
    ax1.set_yscale("log") ; ax1.set_ylim( ymin, ymax )
    plt.axvline(x=850.0, linestyle="dashed", color="k")
    plt.axvline(x=450.0, linestyle="dashed", color ="k")
    plt.legend(loc = "lower center", fontsize=10, scatterpoints = 20)
    plt.tight_layout()
    plt.savefig(plt_dir + "/SED_comp."+plt_form, format = plt_form)
    plt.clf()

    fig = plt.figure(1)
    ax1 = plt.subplot(111)

    for j in range(len(isrf_ext)):
        for t in range(len(temps)):
            plt.plot( wav[t][i][j], flam[t][i][j], \
              label = phase[t] + ": " + leg_names[j], color=color[t], \
              linestyle = lines[j])
    ymax = 1.00e-10 ; ymin = 1.00e-13
    plt.xlabel("Wavelength ("+(r"$\mu$m")+")", fontsize = 18, labelpad=0.5)
    plt.xticks(fontsize = 15) ;   ax1.set_xscale("log")
    ax1.set_xlim(300.,1200.)
    plt.ylabel(r"$\lambda$ F$_{\lambda}$ (erg/s/cm$^2$)", \
      fontsize = 18, labelpad=0.5)
    ax1.set_yscale("log") ; ax1.set_ylim( ymin, ymax )
    plt.axvline(x=850.0, linestyle="dashed", color="k")
    plt.axvline(x=450.0, linestyle="dashed", color ="k")
    plt.tight_layout()
    plt.savefig( plt_dir + "/SED_ZOOM_comp."+plt_form, \
      format = plt_form)
    plt.clf()



### ------------------------------------------------------------------------ ###


#   Singular cut method analysis:
#
#   450 and 850 micron flux ratio calculation (lambda F lambda)
#
L_ratio = round((float(temps[1])/float(temps[0]))**(4.0),3)
print("Ratio of protostellar luminosities, assuming constant "+ \
  "radius is: {0}\n".format(L_ratio) )
#
    # Single beam restriction
for i in range(len(restricted)):
    for j in range(len(isrf_ext)):
#
        print("\nLAMBDA F_LAMBDA RATIOS FOR "+ \
          "{0} CUT, \"{1}\" ISRF RUN OF EC 53 MODELLING\n".format(restricted[i], \
            isrf_ext[j]) )

        flam_int1 = interp1d(wav[0][i][j], flam[0][i][j])
        flam_int2 = interp1d(wav[1][i][j], flam[1][i][j])

        ratio450 = round( ( flam_int2(450) / flam_int1(450) ), 3)
        print("Ratio 450 micron flux is: {0}\n".format(ratio450) )
        ratio850 = round( ( flam_int2(850) / flam_int1(850) ), 3)
        print("Ratio 850 micron flux is: {0}\n".format(ratio850) )
#
    # Cross comparison on flux ratios between beam sizes
#
print("\n\n\n")
for j in range(len(isrf_ext)):
    for t in range(len(temps)):
#
        print("LAMBDA F_LAMBDA RATIOS FOR "+ \
          "{0} PHASE, \"{1}\" ISRF RUN OF EC 53 MODELLING\n".format(phase[t], \
          isrf_ext[j]))

        flam_int1 = interp1d(wav[t][1][j], flam[t][1][j])
        flam_int2 = interp1d(wav[t][0][j], flam[t][0][j])

        ratio450 = round( ( flam_int2(450) / flam_int1(450) ), 3)
        print("Ratio 450 micron flux is: {0}\n".format(ratio450) )
        ratio850 = round( ( flam_int2(850) / flam_int1(850) ), 3)
        print("Ratio 850 micron flux is: {0}\n".format(ratio850) )
