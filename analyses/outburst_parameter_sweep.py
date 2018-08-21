'''

outburst_parameter_sweep.py

Program to evaluate flux levels for varying model component/outburst
characteristics. All data stored in 1686/ directory (i.e. for E2 outburst), all
POST cut, with ISRF heating, with parameters as:

10 L_sol
20 L_sol
50 L_sol
80 L_sol

all with varied cavity opening angles of 0, 5 and 10 deg. at 10 000 AU.
20 L_sol and 80 L_sol runs also have rescaled disc and envelope masses, of
x0.1 and x2 respectively.

Fluxes then compared to quiescent phase snapshot (1611/), with self-similar
cut, ISRF heating and component rescaling.

Using SEDs of all models, bolometric luminosity vs. bolometric temperature
computed (numerical integration) and plotted, as comparison to observational
equivalents.

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
from matplotlib import rc
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
lums = [10, 20, 50, 80]   # Int. Arr.: Central luminsoity of runs
lums_label = [5, 10, 25, 40]
temps = [5934, 7056, 8871, 9978] # Int. Arr.: Temperature of central source
#
cavity_theta = [10] # Int. Arr.: Cavity opening angles
#
r_cores = ["1e3AU" ,"1e4AU"]  # Radius of cores RRT carried out over
#
mult_tags = ["01", "2"] # Str. Arr: Tags for unique component rescales
#
inclins = [0,30,60,90] # Int. Arr.: Inclinations used in RRT
#
wavs_compare = [70, 250, 850, 1300]   # Int. Arr.: Wavelengths
#                                              at which flux ratios are compared
plt_form = "png"        # Str.: Output plot format ["png","eps"]
#
#   Indices:
#       s - for snaps
#       c - for cut
#       r - for isrf
#       i - for inclins
#       l - for llambda
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
		# # # - - - MAIN PROGRAM - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
cwd = os.getcwd()
print cwd
arch_dir = cwd + "/../.."
dat_dir_o = arch_dir + "/runs/OUTBURST/1686"
dat_dir_q = arch_dir + "/runs/OUTBURST/1611"
plt_dir = arch_dir + "/runs/OUTBURST/plots_analysis"

if not os.path.isdir(plt_dir):
    os.system("mkdir {0}/".format(plt_dir) )
#
# flam_ratio[LUMS][CAVITY_THETA][MULT_TAGS*2][INCLINS][WAVS_COMPARE]
#
flam_ratio = [ [ [ [ [ [ [0.0] for l in range(len(wavs_compare)) ]
  for i in range(len(inclins)) ]
  for r in range(len(mult_tags)*2 + 1) ]
  for c in range(len(cavity_theta)) ]
  for s in range(len(lums)) ]
  for cut in range(len(r_cores)) ]
L_bol = [ [ [ [ [ [] for i in range(len(inclins)) ]
  for r in range(len(mult_tags)*2 + 1) ]
  for c in range(len(cavity_theta)) ]
  for s in range(len(lums) + 1) ]
  for cut in range(len(r_cores)) ]
T_bol = [ [ [ [ [ [] for i in range(len(inclins)) ]
  for r in range(len(mult_tags)*2 + 1) ]
  for c in range(len(cavity_theta)) ]
  for s in range(len(lums) + 1) ]
  for cut in range(len(r_cores)) ]
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
    # First, read cavity-only data (even if theta = 0)
#
for cut in range(len(r_cores)):

    print("\nCarrying out analyses for {0} cut models".format(r_cores[cut]))

    print "\nNow processing cavity-only models\n"

    for s in range(len(lums)):
#
        for i in range(len(inclins)):
#
            for c in range(len(cavity_theta)):
#
                f_sed_o = dat_dir_o+"/"+r_cores[cut]+"/POST_"+ \
                 str(lums[s])+"L_cavity"+ str(cavity_theta[c])+ \
                 "/dat/spectrum"+str(int(temps[s]))+ \
                 "_ISRF_"+str(inclins[i])+"i.out"

                wavs_o, flam_o, lam_flam_o, \
                 nu_o, fnu_o, nu_fnu_o = fs.SED_read(f_sed_o)
#
                f_sed_q = dat_dir_q+"/"+r_cores[cut]+ \
                 "/POST_cavity"+str(cavity_theta[c])+ \
                 "/dat/spectrum3974_ISRF_"+str(inclins[i])+"i.out"
#
                wavs_q, flam_q, lam_flam_q, \
                 nu_q, fnu_q, nu_fnu_q = fs.SED_read(f_sed_q)
#
                flam_int_o = interp1d(wavs_o,flam_o, kind="linear")
                flam_int_q = interp1d(wavs_q,flam_q, kind="linear")

                for w in range(len(wavs_compare)):
                    flam_ratio[cut][s][c][0][i][w] = \
                     flam_int_o(wavs_compare[w]) / flam_int_q(wavs_compare[w])

    # Model bolometric luminosity/temperature - quiescent

                if (s == 0):
#
                    wavs_q *= cs.cm_per_micron
#
                    L_tmp = fs.Lbol_calc(nu_q, fnu_q)
                    L_bol[cut][0][c][0][i].append( L_tmp )
#
                    T_tmp = fs.Tbol_calc(nu_q, fnu_q, nu_fnu_q)
                    T_bol[cut][0][c][0][i].append( T_tmp )

    # Model bolometric luminosity/temperature - outbursting

                wavs_o *= cs.cm_per_micron
#
                L_tmp = fs.Lbol_calc(nu_o, fnu_o)
                L_bol[cut][1 + s][c][0][i].append( L_tmp )
#
                T_tmp = fs.Tbol_calc(nu_o, fnu_o, nu_fnu_o)
                T_bol[cut][1 + s][c][0][i].append( T_tmp )
#
    # Now models with rescaled component mass - disc
#
    print "\nNow processing component mass rescaled models\n"

    for s in range(len(lums)):
#
        for i in range(len(inclins)):
#
            for c in range(len(cavity_theta)):
#
                for r in range(len(mult_tags)):

#
                    f_sed_o = dat_dir_o+"/"+r_cores[cut]+ \
                     "/POST_"+str(lums[s])+"L_cavity"+ \
                     str(cavity_theta[c])+"_disc"+str(mult_tags[r])+ \
                     "/dat/spectrum"+str(int(temps[s]))+"_ISRF_"+ \
                     str(inclins[i])+"i.out"
#
                    wavs_o, flam_o, lam_flam_o, \
                     nu_o, fnu_o, nu_fnu_o = fs.SED_read(f_sed_o)
#
                    f_sed_q = dat_dir_q+"/"+r_cores[cut]+ \
                     "/POST_cavity"+str(cavity_theta[c])+ \
                     "_disc"+str(mult_tags[r])+"/dat/spectrum3974_ISRF_"+ \
                     str(inclins[i])+"i.out"
#
                    wavs_q, flam_q, lam_flam_q, \
                     nu_q, fnu_q, nu_fnu_q = fs.SED_read(f_sed_q)
#
                    flam_int_o = interp1d(wavs_o,flam_o, kind="linear")
                    flam_int_q = interp1d(wavs_q,flam_q, kind="linear")

                    for w in range(len(wavs_compare)):
                        flam_ratio[cut][s][c][1 + r][i][w] = \
                         flam_int_o(wavs_compare[w]) / flam_int_q(wavs_compare[w])

    # Model bolometric luminosity/temperature - quiescent

                    if (s == 0):
#
                        wavs_q *= cs.cm_per_micron
#
                        L_tmp = fs.Lbol_calc(nu_q, fnu_q)
                        L_bol[cut][0][c][1 + r][i].append( L_tmp )
#
                        T_tmp = fs.Tbol_calc(nu_q, fnu_q, nu_fnu_q)
                        T_bol[cut][0][c][1 + r][i].append( T_tmp )

    # Model bolometric luminosity/temperature - outbursting

                    wavs_o *= cs.cm_per_micron
#
                    L_tmp = fs.Lbol_calc(nu_o, fnu_o)
                    L_bol[cut][1 + s][c][1 + r][i].append( L_tmp )
#
                    T_tmp = fs.Tbol_calc(nu_o, fnu_o, nu_fnu_o)
                    T_bol[cut][1 + s][c][1 + r][i].append( T_tmp )
#
    # Now models with rescaled component mass - envelope

    for s in range(len(lums)):
#
        for i in range(len(inclins)):
#
            for c in range(len(cavity_theta)):
#
                for r in range(len(mult_tags)):

#
                    f_sed_o = dat_dir_o+"/"+r_cores[cut]+ \
                     "/POST_"+str(lums[s])+"L_cavity"+ \
                     str(cavity_theta[c])+"_env"+str(mult_tags[r])+ \
                     "/dat/spectrum"+str(int(temps[s]))+"_ISRF_"+ \
                     str(inclins[i])+"i.out"
#
                    wavs_o, flam_o, lam_flam_o, \
                     nu_o, fnu_o, nu_fnu_o = fs.SED_read(f_sed_o)
#
                    f_sed_q = dat_dir_q+"/"+r_cores[cut]+ \
                     "/POST_cavity"+str(cavity_theta[c])+ \
                     "_env"+str(mult_tags[r])+"/dat/spectrum3974_ISRF_"+ \
                     str(inclins[i])+"i.out"
#
                    wavs_q, flam_q, lam_flam_q, \
                     nu_q, fnu_q, nu_fnu_q = fs.SED_read(f_sed_q)
#
                    flam_int_o = interp1d(wavs_o,flam_o, kind="linear")
                    flam_int_q = interp1d(wavs_q,flam_q, kind="linear")

                    for w in range(len(wavs_compare)):
                        flam_ratio[cut][s][c][1 + len(mult_tags) + r][i][w] = \
                         flam_int_o(wavs_compare[w]) / flam_int_q(wavs_compare[w])

    # Model bolometric luminosity/temperature - quiescent

                    if (s == 0):
#
                        wavs_q *= cs.cm_per_micron
#
                        L_tmp = fs.Lbol_calc(nu_q, fnu_q)
                        L_bol[cut][0][c][1 + len(mult_tags) + r][i].append( L_tmp )
#
                        T_tmp = fs.Tbol_calc(nu_q, fnu_q, nu_fnu_q)
                        T_bol[cut][0][c][1 + len(mult_tags) + r][i].append( T_tmp )
#
    # Model bolometric luminosity/temperature - outbursting
#
                    wavs_o *= cs.cm_per_micron
#
                    L_tmp = fs.Lbol_calc(nu_o, fnu_o)
                    L_bol[cut][1 + s][c][1 + len(mult_tags) + r][i].append( L_tmp )
#
                    T_tmp = fs.Tbol_calc(nu_o, fnu_o, nu_fnu_o)
                    T_bol[cut][1 + s][c][1 + len(mult_tags) + r][i].append( T_tmp )
#
    # Finally, read in models with simulation protostellar luminosities
#
    print "\nNow processing hydro output models\n"

    sim_tags = [["0901","0951"],["1611","1686"]]
    sim_temps = [["3421","16678"],["3974","19248"]]
    sim_L_bol = [ [ [ [] for i in range(len(inclins)) ]
     for i in range(len(sim_tags[0])) ] \
     for i in range(len(sim_tags)) ]
    sim_T_bol = [ [ [ [] for i in range(len(inclins)) ]
     for i in range(len(sim_tags[0])) ] \
     for i in range(len(sim_tags)) ]

    for e in range(len(sim_tags)):
        for t in range(len(sim_tags[0])):
            for i in range(len(inclins)):
#
                f_sed = dat_dir_o+"/../"+sim_tags[e][t]+"/POST/dat/spectrum"+ \
                 sim_temps[e][t]+"_ISRF_"+str(inclins[i])+"i.out"
#
                wavs, flam, lam_flam, nu, fnu, nu_fnu = fs.SED_read(f_sed)

    # Model bolometric luminosity/temperature
#
                wavs *= cs.cm_per_micron
#
                L_tmp = fs.Lbol_calc(nu, fnu)
                sim_L_bol[e][t][i].append( L_tmp )

                T_tmp = fs.Tbol_calc(nu, fnu, nu_fnu)
                sim_T_bol[e][t][i].append( T_tmp )
#
    print "\nPlotting\n"

    colors = ["r", "g", "b", "k"]
    for w in range(len(wavs_compare)):

        fig = plt.figure(1)
        ax1 = plt.subplot(111)

        for s in range(len(lums)):

            legend = ""

            for c in range(len(cavity_theta)):

                if (c == 0):
                    legend = r"Cavity-only ($\theta = 10^{\circ}$)"
                    marker = "^"
                    fc = "none"
                elif (c == 1):
                    legend = r"$\theta = 5^{\circ}$"
                    marker = ">"
                    fc = "none"
                elif (c == 2):
                    legend = r"$\theta = 10^{\circ}$"
                    marker = "<"
                    fc = "none"
#
                for i in range(len(inclins)):
#
                    if (s == 0) and (i == 3):
                        plt.scatter(lums_label[s], flam_ratio[cut][s][c][0][i][w], \
                         color=colors[i], label = legend, \
                         marker = marker, s = 20*(s+1), facecolors = fc )
                    else:
                        plt.scatter(lums_label[s], flam_ratio[cut][s][c][0][i][w], \
                         color=colors[i], \
                         marker = marker, s = 20*(s+1), facecolors = fc )

#
            for r in range(len(mult_tags)*2 + 1):

                if (r ==  0):
                    continue
                elif (r == 1):
                    legend = r"Disc x0.1"
                    marker = "1"
                    fc = "none"
                elif (r == 2):
                    legend = r"Disc x2.0"
                    marker = "2"
                    fc = "none"
                elif (r == 3):
                    legend = r"Envelope x0.1"
                    marker = "+"
                    fc = "none"
                elif (r == 4):
                    legend = r"Envelope x2.0"
                    marker = "x"
                    fc = "none"
#
                for i in range(len(inclins)):
#
                    if (r == 0):
                        continue

                    if (s == 0) and (i == 3):
                        plt.scatter(lums_label[s], flam_ratio[cut][s][0][r][i][w], \
                         color=colors[i], label = legend, \
                         marker = marker, s = 20*(s+1), facecolors = fc )
                    else:
                        plt.scatter(lums_label[s], flam_ratio[cut][s][0][r][i][w], \
                         color=colors[i], \
                         marker = marker, s = 20*(s+1), facecolors = fc )

#

        plt.xlabel(r"L$_{*, o}$ / L$_{*, q}$", fontsize = cs.fontsize, \
         labelpad=0.5)
        plt.ylabel(r"F$_{\lambda, o}$ / F$_{\lambda, q}$", \
         fontsize = cs.fontsize, labelpad=0.5)
#    ax1.set_ylim(0.8*flam_ratio[0][0][0][3][w], 1.2*flam_ratio[3][0][3][0][w])
        plt.text(0.8, 0.1,r'$\lambda$ = '+str(wavs_compare[w])+r'$\mu$m', \
         ha='center', va='center', transform=ax1.transAxes, fontsize = cs.fontsize)
        if (wavs_compare[w] == 70) or (wavs_compare[w] == 850):
            plt.legend(loc = "upper left", fontsize=cs.leg_fontsize)

        plt.tight_layout()
        plt.savefig(plt_dir+"/"+r_cores[cut]+"_"+ \
         str(wavs_compare[w])+"micron_param_ratios."+plt_form, format = plt_form)
        plt.clf()
#
#
    colors = ["r", "g", "b", "k"]
    fig = plt.figure(1)
    ax1 = plt.subplot(111)

    for s in range(len(lums) + 1):

        legend = ""

        for c in range(len(cavity_theta)):

            if (c == 0):
                legend = r"Cavity-only ($\theta = 10^{\circ}$)"
                marker = "^"
                fc = "none"
            elif (c == 1):
                legend = r"$\theta = 5^{\circ}$"
                marker = ">"
                fc = "none"
            elif (c == 2):
                legend = r"$\theta = 10^{\circ}$"
                marker = "<"
                fc = "none"
#
            for i in range(len(inclins)):
#
                if (s == 0) and (i == 3):
                    plt.scatter(T_bol[cut][s][c][0][i][0], L_bol[cut][s][c][0][i][0], \
                     color=colors[i], label = legend, \
                     marker = marker, s = 20*(s+1), facecolors = fc )
                else:
                    plt.scatter(T_bol[cut][s][c][0][i][0], L_bol[cut][s][c][0][i][0], \
                     color=colors[i], \
                     marker = marker, s = 20*(s+1), facecolors = fc )

#
            for r in range(len(mult_tags)*2 + 1):

                if (r ==  0):
                    continue
                elif (r == 1):
                    legend = r"Disc x0.1"
                    marker = "1"
                    fc = "none"
                elif (r == 2):
                    legend = r"Disc x2.0"
                    marker = "2"
                    fc = "none"
                elif (r == 3):
                    legend = r"Envelope x0.1"
                    marker = "+"
                    fc = "none"
                elif (r == 4):
                    legend = r"Envelope x2.0"
                    marker = "x"
                    fc = "none"
#
                for i in range(len(inclins)):
#
                    if (r == 0):
                        continue

                    if (s == 0) and (c == 0) and (i == 3):
                        plt.scatter(T_bol[cut][s][c][r][i][0], L_bol[cut][s][c][r][i][0], \
                         color=colors[i], label = legend, \
                         marker = marker, s = 20*(s+1), facecolors = fc )
                    else:
                        plt.scatter(T_bol[cut][s][c][r][i][0], L_bol[cut][s][c][r][i][0], \
                         color=colors[i], \
                         marker = marker, s = 20*(s+1), facecolors = fc )
#
    for e in range(len(sim_tags)):

        if (e == 0):
            marker = "d"
            legend = r"E1"
        elif (e == 1):
            marker = "D"
            legend = r"E2"

        for t in range(len(sim_tags[0])):

            for i in range(len(inclins)):

                if (t == 0) and (i == 3):
                    plt.scatter(sim_T_bol[e][t][i][0], sim_L_bol[e][t][i][0], \
                      color = colors[i], label = legend, marker = marker, \
                      s = 30)
                else:
                    plt.scatter(sim_T_bol[e][t][i][0], sim_L_bol[e][t][i][0], \
                      color = colors[i], marker = marker, \
                      s = 30)

#
    plt.ylabel(r"L$_{BOL}$ (L$_{\odot}$)", fontsize = cs.fontsize, labelpad=0.5)
    ax1.set_xscale("log")
    ax1.set_xlim(1000,10)
    plt.xlabel(r"T$_{BOL}$ (K)", fontsize = cs.fontsize, \
     labelpad=0.5)
    ax1.set_yscale("log")
    ax1.set_ylim(0.5,2500)
    plt.axvline(x=70.00, linestyle="dashed", color="k")
    legend1= plt.legend(loc = "upper right", fontsize=cs.leg_fontsize-2)
    plt.gca().add_artist(legend1)

    i1, = plt.plot([1e-99,1e-98], [1e-99,1e-98], linestyle = "-", color = "r")
    i2, = plt.plot([1e-99,1e-98], [1e-99,1e-98], linestyle = "-", color = "g")
    i3, = plt.plot([1e-99,1e-98], [1e-99,1e-98], linestyle = "-", color = "b")
    i4, = plt.plot([1e-99,1e-98], [1e-99,1e-98], linestyle = "-", color = "k")
    leg_add = [] ; leg_add.append([i1,i2,i3,i4])
    legend2 = plt.legend( leg_add[0], \
     ["i = 0 deg.", "i = 30 deg.", "i = 60 deg.", "i = 90 deg."], \
     fontsize=cs.leg_fontsize-2, loc="upper left")
    plt.gca().add_artist(legend2)

    plt.tight_layout()
    plt.savefig(plt_dir+"/Lbol_Tbol_"+r_cores[cut]+"."+plt_form, \
     format = plt_form)
    plt.clf()

    # Compute, and plot difference of bolometric properties for all models,
    # between 10 000 AU and 1000 AU core radius for SED integration

T_bol = np.array(T_bol) ; L_bol = np.array(L_bol)
Tbol_diff = T_bol[0] - T_bol[1] ; Lbol_diff = L_bol[0] - L_bol[1]

colors = ["r", "g", "b", "k"]
fig = plt.figure(1)
ax1 = plt.subplot(111)

for s in range(len(lums) + 1):

    legend = ""

    for c in range(len(cavity_theta)):

        if (c == 0):
            legend = r"Cavity-only ($\theta = 10^{\circ}$)"
            marker = "^"
            fc = "none"
        elif (c == 1):
            legend = r"$\theta = 5^{\circ}$"
            marker = ">"
            fc = "none"
        elif (c == 2):
            legend = r"$\theta = 10^{\circ}$"
            marker = "<"
            fc = "none"
#
        for i in range(len(inclins)):
#
            if (s == 0) and (i == 3):
                plt.scatter(Tbol_diff[s][c][0][i][0], Lbol_diff[s][c][0][i][0], \
                 color=colors[i], label = legend, \
                 marker = marker, s = 20*(s+1), facecolors = fc )
            else:
                plt.scatter(Tbol_diff[s][c][0][i][0], Lbol_diff[s][c][0][i][0], \
                 color=colors[i], \
                 marker = marker, s = 20*(s+1), facecolors = fc )
#
        for r in range(len(mult_tags)*2 + 1):

            if (r ==  0):
                continue
            elif (r == 1):
                legend = r"Disc x0.1"
                marker = "1"
                fc = "none"
            elif (r == 2):
                legend = r"Disc x2.0"
                marker = "2"
                fc = "none"
            elif (r == 3):
                legend = r"Envelope x0.1"
                marker = "+"
                fc = "none"
            elif (r == 4):
                legend = r"Envelope x2.0"
                marker = "x"
                fc = "none"
#
            for i in range(len(inclins)):
#
                if (r == 0):
                    continue

                if (s == 0) and (c == 0) and (i == 3):
                    plt.scatter(Tbol_diff[s][c][r][i][0], Lbol_diff[s][c][r][i][0], \
                     color=colors[i], label = legend, \
                     marker = marker, s = 20*(s+1), facecolors = fc )
                else:
                    plt.scatter(Tbol_diff[s][c][r][i][0], Lbol_diff[s][c][r][i][0], \
                     color=colors[i], \
                     marker = marker, s = 20*(s+1), facecolors = fc )
#
plt.ylabel(r"$\Delta$ L$_{BOL}$ (L$_{\odot}$)", fontsize = cs.fontsize, labelpad=0.5)
ax1.set_xlim(150,-20)
plt.xlabel(r"$\Delta$ T$_{BOL}$ (K)", fontsize = cs.fontsize, \
 labelpad=0.5)
ax1.set_ylim(-40,5)
legend1= plt.legend(loc = "lower left", fontsize=cs.leg_fontsize-2)
plt.gca().add_artist(legend1)

i1, = plt.plot([1e-99,1e-98], [1e-99,1e-98], linestyle = "-", color = "r")
i2, = plt.plot([1e-99,1e-98], [1e-99,1e-98], linestyle = "-", color = "g")
i3, = plt.plot([1e-99,1e-98], [1e-99,1e-98], linestyle = "-", color = "b")
i4, = plt.plot([1e-99,1e-98], [1e-99,1e-98], linestyle = "-", color = "k")
leg_add = [] ; leg_add.append([i1,i2,i3,i4])
legend2 = plt.legend( leg_add[0], \
 ["i = 0 deg.", "i = 30 deg.", "i = 60 deg.", "i = 90 deg."], \
 fontsize=cs.leg_fontsize-2, loc="center left")
plt.gca().add_artist(legend2)

plt.tight_layout()
plt.savefig(plt_dir+"/Lbol_Tbol_diff."+plt_form, format = plt_form)
plt.clf()

    # Same process, but for flux ratios.

flam_ratio = np.array(flam_ratio)
flam_ratio_diff = flam_ratio[0] - flam_ratio[1]

colors = ["r", "g", "b", "k"]
for w in range(len(wavs_compare)):

    fig = plt.figure(1)
    ax1 = plt.subplot(111)

    for s in range(len(lums)):

        legend = ""

        for c in range(len(cavity_theta)):

            if (c == 0):
                legend = r"Cavity-only ($\theta = 10^{\circ}$)"
                marker = "^"
                fc = "none"
            elif (c == 1):
                legend = r"$\theta = 5^{\circ}$"
                marker = ">"
                fc = "none"
            elif (c == 2):
                legend = r"$\theta = 10^{\circ}$"
                marker = "<"
                fc = "none"
#
            for i in range(len(inclins)):
#
                if (s == 0) and (i == 3):
                    plt.scatter(lums_label[s], flam_ratio_diff[s][c][0][i][w], \
                     color=colors[i], label = legend, \
                     marker = marker, s = 20*(s+1), facecolors = fc )
                else:
                    plt.scatter(lums_label[s], flam_ratio_diff[s][c][0][i][w], \
                     color=colors[i], \
                     marker = marker, s = 20*(s+1), facecolors = fc )

#
        for r in range(len(mult_tags)*2 + 1):

            if (r ==  0):
                continue
            elif (r == 1):
                legend = r"Disc x0.1"
                marker = "1"
                fc = "none"
            elif (r == 2):
                legend = r"Disc x2.0"
                marker = "2"
                fc = "none"
            elif (r == 3):
                legend = r"Envelope x0.1"
                marker = "+"
                fc = "none"
            elif (r == 4):
                legend = r"Envelope x2.0"
                marker = "x"
                fc = "none"
#
            for i in range(len(inclins)):
#
                if (r == 0):
                    continue

                if (s == 0) and (i == 3):
                    plt.scatter(lums_label[s], flam_ratio_diff[s][0][r][i][w], \
                     color=colors[i], label = legend, \
                     marker = marker, s = 20*(s+1), facecolors = fc )
                else:
                    plt.scatter(lums_label[s], flam_ratio_diff[s][0][r][i][w], \
                     color=colors[i], \
                     marker = marker, s = 20*(s+1), facecolors = fc )

#

    plt.xlabel(r"L$_{*, o}$ / L$_{*, q}$", fontsize = cs.fontsize, \
     labelpad=0.5)
    plt.ylabel(r"$\Delta$ (F$_{\lambda, o}$ / F$_{\lambda, q}$)", \
     fontsize = cs.fontsize, labelpad=0.5)

    if (wavs_compare[w] == 70):
        plt.text(0.15, 0.4,r'$\lambda$ = '+str(wavs_compare[w])+r'$\mu$m', \
         ha='center', va='center', transform=ax1.transAxes, \
         fontsize = cs.fontsize)
    elif (wavs_compare[w] == 250):
        plt.text(0.15, 0.1,r'$\lambda$ = '+str(wavs_compare[w])+r'$\mu$m', \
         ha='center', va='center', transform=ax1.transAxes, \
         fontsize = cs.fontsize)
    elif (wavs_compare[w] == 850):
        plt.text(0.15, 0.6,r'$\lambda$ = '+str(wavs_compare[w])+r'$\mu$m', \
         ha='center', va='center', transform=ax1.transAxes, \
         fontsize = cs.fontsize)
    elif (wavs_compare[w] == 1300):
        plt.text(0.15, 0.9,r'$\lambda$ = '+str(wavs_compare[w])+r'$\mu$m', \
         ha='center', va='center', transform=ax1.transAxes, \
         fontsize = cs.fontsize)

    if (wavs_compare[w] == 70):
        plt.legend(loc = "lower left", fontsize=cs.leg_fontsize)
    if (wavs_compare[w] == 850):
        plt.legend(loc = "upper left", fontsize=cs.leg_fontsize)

    plt.tight_layout()
    plt.savefig(plt_dir+"/"+str(wavs_compare[w])+ \
     "micron_param_ratios_diff."+plt_form, format = plt_form)
    plt.clf()
