'''

outburst.py

Program to read in SED data for OUTBURST models, and tabulate [350,450,850,1300]
micron flux ratios for snapshots between all models

and plot SEDs to display 100 - 10000 micron flux for


Author: Benjamin MacFarlane
Date: 02/07/2018
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
# User modules
import constants as cs
from functions import adjustFigAspect
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - VARIABLE DEFINITIONS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
snaps = ["0901","0951","1611","1686"]   # Str. Arr.: DE05 tags
#   Corresponding times: ["87.00 kyr","87.50 kyr","94.10 kyr","94.85 kyr"]
#
isrf = ["N","Y"]         # Str. Arr.: ISRF usage (yes (Y) and/or no (N))
cut = ["NOT","POST"]     # Str. Arr.: Restriction usage (NOT and/or POST)
inclins = [0,60,90]      # Int. Arr.: Inclinations used for snapshots
#
sed_norm = False         # Bool: Normalise to stellar flux, or to distance?
dist = 140.0            # Float: Distance to source (pc)
#
wavs_compare = [70, 250, 350, 450, 850, 1300]   # Int. Arr.: Wavelength values at which
#                                       flux ratios are compared
#
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
arch_dir = cwd + "/../.."
dat_dir = arch_dir + "/runs/OUTBURST"
plt_dir = dat_dir + "/plots_analysis"

if not os.path.isdir(plt_dir):
    os.system("mkdir {0}/".format(plt_dir) )
#
t_star = [] ; r_star = [] ; l_star = []
#
# flam[SNAPS][CUT][ISRF][INCLINS], to be appended
# to flam[SNAPS][CUT][ISRF][INCLINS][llambda] (equivalent for wavs)
#
wavs = [ [ [ [ [] for i in range(len(inclins)) ]
  for r in range(len(isrf)) ]
  for c in range(len(cut)) ]
  for s in range(len(snaps)) ]
flam = [ [ [ [ [] for i in range(len(inclins)) ]
  for r in range(len(isrf)) ]
  for c in range(len(cut)) ]
  for s in range(len(snaps)) ]
#
# flam_vals[SNAP][CUT][ISRF][INCLINS][WAVS_COMPARE]
#
flam_vals = [ [ [ [ [ [0.0] for l in range(len(wavs_compare)) ]
  for i in range(len(inclins)) ]
  for r in range(len(isrf)) ]
  for c in range(len(cut)) ]
  for s in range(len(snaps)) ]
#
wavs_s = [ [] for s in range(len(snaps))]
flam_s = [ [] for s in range(len(snaps))]
#
# rloc[SNAP][CUT], to be appended to
# rloc[SNAP][CUT][rloc] (rloc_xy with equivalent structure)
#
rloc = [ [ [] for c in range(len(cut)) ] for s in range(len(snaps)) ]
rloc_xy = [ [ [] for c in range(len(cut)) ] for s in range(len(snaps)) ]
#
# temp[SNAP][CUT][ISRF], to be appended to
# temp[SNAP][CUT][ISRF][rloc]
#
temp = [ [ [ [] for r in range(len(isrf)) ]
  for c in range(len(cut)) ]
  for s in range(len(snaps)) ]
#

#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
f_flam = os.getcwd()+"/../flam_values.dat"
f_flam = open(f_flam, "w")
f_flam.write("#snapshot\tcut\tisrf\tinclin\twav\tflam\n")
print("\nFlux values for all selected wavelengths are "+ \
  "written in /analysis/flam_values.dat\n")
#
for s in range(len(snaps)):
#
    # Extract protostellar properties for snapshot
#
    f_star = dat_dir+"/"+snaps[s]+"/"+cut[0]+"/dat/stars.inp"
    f_star = open(f_star, "r")
    trash = f_star.readline() ; trash = f_star.readline()
    r_star.append(float(f_star.readline().split(" ")[0]) * cs.rsol_per_cm)
    t_star.append(-float(f_star.readlines()[-1]))
#
    l_star.append(round((r_star[s])**2.0 * (t_star[s] / 5780.)**4.0, 3) )
    round_ind = 2 - int(math.floor(math.log10(l_star[s])))
    l_star[s] = round(l_star[s], round_ind)
#
    f_star.close()
#
    # Append wavs and flam for protostellar data, normalising as necessary
#
    f_psed = arch_dir+"/protostar/spectrum"+str(int(t_star[s]))+".out"
    if not os.path.isfile(f_psed):
        print("Protostellar SED template not found, exiting...")
        exit()
    f_psed = open(f_psed,"r")
    for l in range(3):
        trash = f_psed.readline()
    for lines in f_psed:
        lines = lines.strip() ; columns = lines.split()
        wavs_s[s].append(float(columns[0]))
        flam_s[s].append( float(columns[1]) * \
          ( cs.c_cgs / (float(columns[0]) * cs.cm_per_micron) ) )
    f_psed.close()
#
    if sed_norm:
        fnorm = ( (r_star[s] * cs.cm_per_rsol)**2. * cs.sigma_sb_cgs * \
           (t_star[s])**4. ) / (cs.cm_per_pc)**2.
        for l in range(len(flam_s[s]) ):
            flam_s[s][l] /= fnorm
    elif not sed_norm:
        fnorm = dist**2.0
        for l in range(len(flam_s[s])):
            flam_s[s][l] /= fnorm
#
    # Now loop over cut, isrf and inclination regimes, reading in and appending
    # SED data to wavs and flam arrays, and normalising as necessary
#
    for c in range(len(cut)):
#
        f_loc = dat_dir+"/"+snaps[s]+"/"+cut[c]+"/dat/loc_grid.dat"
        f_loc = open(f_loc,"r")
        for lines in f_loc:
            lines = lines.strip() ; columns = lines.split()
            rloc[s][c].append(float(columns[3]) * cs.au_per_cm)
            rloc_xy[s][c].append( np.sqrt(float(columns[0])**2.0 + \
              float(columns[1])**2.0) )
        f_loc.close()
#
        for r in range(len(isrf)):
#
            if (isrf[r] == "N"):
                isrf_ext = ""
            elif (isrf[r] == "Y"):
                isrf_ext = "_ISRF"
#
            f_temp = dat_dir+"/"+snaps[s]+"/"+cut[c]+ \
              "/dat/dust_temperature"+str(int(t_star[s]))+isrf_ext+".dat"
            f_temp = open(f_temp,"r")
            for l in range(3):
                trash = f_temp.readline()
            for lines in f_temp:
                lines = lines.strip() ; columns = lines.split()
                temp[s][c][r].append( float( columns[0] ) )
            f_temp.close()
#
            for i in range(len(inclins)):
#
    # Test whether rays interesecting bin boundaries has impact of inclinations
    # depedent SED
#
#                if (inclins[i] == 0) and (cut[c] == "POST") and \
#                  ((snaps[s] == "1611") or (snaps[s] == "1686")):
#                    f_sed = dat_dir+"/"+snaps[s]+"/"+cut[c]+ \
#                      "/dat/spectrum"+str(int(t_star[s]))+isrf_ext+ \
#                      "_1i.out"
#                else:
#                    f_sed = dat_dir+"/"+snaps[s]+"/"+cut[c]+ \
#                      "/dat/spectrum"+str(int(t_star[s]))+isrf_ext+ \
#                      "_"+str(int(inclins[i]))+"i.out"
#
                f_sed = dat_dir+"/"+snaps[s]+"/"+cut[c]+ \
                  "/dat/spectrum"+str(int(t_star[s]))+isrf_ext+ \
                  "_"+str(int(inclins[i]))+"i.out"
                f_sed = open(f_sed, "r")
                for l in range(3):
                    trash = f_sed.readline()
                for lines in f_sed:
                    lines = lines.strip() ; columns = lines.split()
                    wavs[s][c][r][i].append(float(columns[0]))
                    flam[s][c][r][i].append( float(columns[1]) * \
                      ( cs.c_cgs / (float(columns[0])*cs.cm_per_micron) ) )
                f_sed.close()
#
                if sed_norm:
                    fnorm = ( (r_star[s] * cs.cm_per_rsol)**2. \
                      * cs.sigma_sb_cgs * (t_star[s])**4. ) / \
                      (cs.cm_per_pc)**2.
                    for l in range(len(flam[s][c][r][i])):
                        flam[s][c][r][i][l] /= fnorm
                elif not sed_norm:
                    fnorm = dist**2.0
                    for l in range(len(flam[s][c][r][i])):
                        flam[s][c][r][i][l] /= fnorm
#
    # Interpolate flam values, then find flam for specified wavelengths and
    # write to file (and terminal)
#
                flam_int = \
                  interp1d(wavs[s][c][r][i],flam[s][c][r][i], kind="linear")
#
                for w in range(len(wavs_compare)):
                    flam_vals[s][c][r][i][w] = \
                      flam_int(wavs_compare[w])
#
                    verbose = False
                    if verbose:
                        print('''
                        Flux of {0} snapshot,
                        {1} cut,
                        with {2} ISRF and
                        inclination of {3}
                        at {4} microns
                        is {5} Jy'''.format(snaps[s], cut[c], isrf[r], \
                          inclins[i], wavs_compare[w], \
                          flam_vals[s][c][r][i][w]) )
#
                    f_flam.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format( \
                      snaps[s], cut[c], isrf[r], inclins[i], \
                      wavs_compare[w], \
                      flam_vals[s][c][r][i][w]) )
f_flam.close()
#
#
    # Read in sink data and plot the temporal evolution of accretion rate
    # and protostellar luminosity
#
sim_sfile = arch_dir+"/ics/OUTBURST/DE05.sink1"
#
fsink = open(sim_sfile,"r")
s_time = [] ; s_mass = [] ; s_lum = [] ; s_temp = [] ; s_macc = []

for line in fsink:
    line = line.strip()
    column = line.split()
    s_time.append(float(column[1]) * 1000.)
    s_mass.append(float(column[22]))
    s_lum.append(float(column[20]))
    s_temp.append(float(column[21]))
    s_macc.append(float(column[24]))
fsink.close()
#
fig = plt.figure(1)
ax1 = plt.subplot(311)
plt.plot(s_time, s_macc, color = "k")
plt.xlim(85.00, 96.85)
ax1.xaxis.set_major_formatter(plt.NullFormatter())
plt.ylim(7.e-7, ax1.get_ylim()[1])
ax1.set_yscale("log")
plt.axvline(x=87.00, linestyle="dashed", color="r")
plt.axvline(x=87.50, linestyle="dashed", color ="b")
plt.axvline(x=94.10, linestyle="dashed", color="r")
plt.axvline(x=94.85, linestyle="dashed", color ="b")
plt.ylabel(r"$\dot{M_{*}}$ (M$_{\odot}$ yr$^{-1}$)", fontsize = cs.fontsize, \
  labelpad=0.5)

fig = plt.figure(1)
ax1 = plt.subplot(312)
plt.plot(s_time, s_mass, color = "k")
plt.xlim(85.00, 96.85)
ax1.set_yticks([0.1, 0.2, 0.3, 0.4, 0.5])
ax1.xaxis.set_major_formatter(plt.NullFormatter())
plt.ylim(0.1, 0.5)
plt.axvline(x=87.00, linestyle="dashed", color="r")
plt.axvline(x=87.50, linestyle="dashed", color ="b")
plt.axvline(x=94.10, linestyle="dashed", color="r")
plt.axvline(x=94.85, linestyle="dashed", color ="b")
plt.ylabel(r"$M_{*}}$ (M$_{\odot}$)", fontsize = cs.fontsize, labelpad=0.5)

ax1 = plt.subplot(313)
plt.plot(s_time, s_lum, color = "k")
plt.xlim(85.00, 96.85)
plt.xlabel("Time (kyr)", fontsize = cs.fontsize, labelpad=0.5)
plt.ylim(5.e-1, ax1.get_ylim()[1])
ax1.set_yscale("log")
plt.axvline(x=87.00, linestyle="dashed", color="r")
plt.axvline(x=87.50, linestyle="dashed", color ="b")
plt.axvline(x=94.10, linestyle="dashed", color="r")
plt.axvline(x=94.85, linestyle="dashed", color ="b")
plt.ylabel(r"$L_{*}}$ (L$_{\odot}$)", fontsize = cs.fontsize, labelpad=0.5)

plt.tight_layout()

plt.savefig(plt_dir + "/pstar_properties."+plt_form, format = plt_form)
plt.clf()

#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
colors_i = ['r','g','b']            # Plot options for inclinations - colors
#
colors_t = ["r","b","gold","cyan"]  # Plot options for temperature - colors
#
leg_ext = [""," - w/ ISRF"]     # Plot options for ISRF - legend
isrf_ext = ["","_ISRF"]         # Plot options for ISRF, plot extension
linestyle = ["-","--",]         # Plot options for ISRF - data plot
#
linewidth = [0.5, 1.5]            # Plot options for accretion state - data plot
name_ea = ["Quiescent","Outbursting"] # Plot option for accretion state - legend
#
name_event = ["ACC1","ACC2"]    # Plot options for accretion tag - filename
#
    # Temperature plots for single snapshot, single cut, multiple ISRF selection
#
print("Plotting single event, single cut, multiple ISRF temperature plots\n")
for c in range(len(cut)):
    for s in range(len(snaps)/2):
        s1 = s*2 ; s2 = (s*2)+1
#
        rloc_plt = np.array(rloc[s2][c]) ; temp_plt = np.array(temp[s2][c][0])
        restr = np.where(temp_plt != 0.)
        rloc_plt = rloc_plt[restr] ; temp_plt = temp_plt[restr]
        rloc_min = min(rloc_plt) ; rloc_max = max(rloc_plt)
        temp_min = 3.0 ; temp_max = max(temp_plt)
#
        fig = plt.figure(1)
        ax1 = plt.subplot(111)
        for r in range(len(isrf)):
            plt.scatter(rloc[s1][c], temp[s1][c][r], color = colors_t[r*2], \
              label = name_ea[0]+": "+str(l_star[s1])+r" $L_{\odot}$"+ \
              leg_ext[r], s = 0.1 )
            plt.scatter(rloc[s2][c], temp[s2][c][r], color = colors_t[(r*2)+1], \
              label = name_ea[1]+": "+str(l_star[s2])+ \
              r" $L_{\odot}$"+leg_ext[r], s = 0.1 )
        plt.legend(loc = "upper right", fontsize=cs.leg_fontsize, \
          scatterpoints = 20)
        plt.xlabel("Spherical Radius (AU)", fontsize = cs.fontsize, \
          labelpad=0.5)
        ax1.set_xscale("log")
        ax1.set_xlim(rloc_min, rloc_max)
        ax1.axvline(10000.00, color = "k", linestyle = "dashed")
        plt.ylabel("Temperature (K)", fontsize = cs.fontsize, labelpad=0.5)
        ax1.set_yscale("log") ; ax1.set_ylim(temp_min+0.1, temp_max*2.0)
        plt.tight_layout()
        plt.savefig(plt_dir+"/T_r_sph_"+name_event[s]+"_"+str(cut[c])+ \
          "."+plt_form, format = plt_form)
        plt.clf()
#
    # Single event, single cut, single ISRF selection, multiple inclinations
#
print("Plotting single event, single cut, single ISRF selection, "+ \
  "with multiple inclinations\n")
for r in range(len(isrf)):
    for c in range(len(cut)):
        for s in range(len(snaps)/2):
            s1 = s*2 ; s2 = (s*2)+1
#
            fig = plt.figure(1)
            ax1 = plt.subplot(111)
#
            for i in range(len(inclins)):
                plt.plot( wavs[s1][c][r][i], flam[s1][c][r][i], \
                  label = name_ea[0]+": "+str(inclins[i])+" deg.", \
                  color=colors_i[i], linestyle = linestyle[1] )
                plt.plot( wavs[s2][c][0][i], flam[s2][c][0][i], \
                  label = name_ea[1]+": "+str(inclins[i])+" deg.", \
                  color=colors_i[i], linestyle = linestyle[0] )
                if sed_norm:
                    ymax = 50.0 ; ymin = 1.e-12
                    plt.ylabel(r"($\lambda$ F$_{\lambda}$) / F$_{*}$", \
                      fontsize = cs.fontsize, labelpad=0.5)
                else:
                    ymax = max(flam[-1][0][0][0])*10. ; ymin = ymax * 1.e-9
                    plt.ylabel(r"$\lambda$ F$_{\lambda}$ (erg cm$^{-2}$ s$^{-1}$)", \
                      fontsize = cs.fontsize, labelpad=0.5)
            plt.xlabel("Wavelength ("+(r"$\mu$m")+")", \
              fontsize = cs.fontsize, labelpad=0.5)
            ax1.set_xscale("log")
            ax1.set_xlim(2.0,10000.)
            ax1.set_yscale("log") ; ax1.set_ylim( ymin, ymax )
            plt.legend(loc = "upper right", fontsize=cs.leg_fontsize, \
              scatterpoints = 20)
            plt.tight_layout()
            plt.savefig(plt_dir + "/SED_"+str(name_event[s])+"_"+str(cut[c])+ \
              isrf_ext[r]+"."+plt_form, format = plt_form)
            plt.clf()
#
    # Single event, single cut, single inclination, multiple ISRF selections
#
print("Plotting single event, single cut, single inclination, "+ \
  "with multiple ISRF selections\n")
for i in range(len(inclins)):
    for c in range(len(cut)):
        for s in range(len(snaps)/2):
            s1 = s*2 ; s2 = (s*2)+1
#
            fig = plt.figure(1)
            ax1 = plt.subplot(111)
#
            for r in range(len(isrf)):
                plt.plot( wavs[s1][c][r][i], flam[s1][c][r][i], \
                  label = name_ea[0]+": "+str(l_star[s1])+ \
                  r" $L_{\odot}$"+leg_ext[r], color=colors_t[0], \
                  linestyle = linestyle[r], linewidth = linewidth[0])
                plt.plot( wavs[s2][c][r][i], flam[s2][c][r][i], \
                  label = name_ea[1]+": "+str(l_star[s2])+ \
                  r" $L_{\odot}$"+leg_ext[r], color=colors_t[1], \
                  linestyle = linestyle[r], linewidth = linewidth[1])
                if sed_norm:
                    ymax = 50.0 ; ymin = 1.e-12
                    plt.ylabel(r"($\lambda$ F$_{\lambda}$) / F$_{*}$", \
                      fontsize = cs.fontsize, labelpad=0.5)
                else:
                    ymax = max(flam[-1][0][0][0])*10. ; ymin = ymax * 1.e-9
                    plt.ylabel(r"$\lambda$ F$_{\lambda}$ (erg cm$^{-2}$ s$^{-1}$)", \
                      fontsize = cs.fontsize, labelpad=0.5)
            plt.xlabel("Wavelength ("+(r"$\mu$m")+")", \
              fontsize = cs.fontsize, labelpad=0.5)
            ax1.set_xscale("log")
#            plt.axvline(x=1300.0, linestyle="dashed", color="k")
#            plt.axvline(x=850.0, linestyle="dashed", color ="k")
#            plt.axvline(x=450.0, linestyle="dashed", color="k")
#            plt.axvline(x=350.0, linestyle="dashed", color ="k")
            ax1.set_xlim(2.0,10000.)
            ax1.set_yscale("log") ; ax1.set_ylim( ymin, ymax )
#            plt.legend(loc = "lower center", fontsize=cs.leg_fontsize, \
#             scatterpoints = 20)
            plt.tight_layout()
            plt.savefig(plt_dir + "/SED_"+str(name_event[s])+"_"+str(cut[c])+ \
              "_"+str(int(inclins[i]))+"i."+plt_form, format = plt_form)
            plt.clf()
#
            fig = plt.figure(1)
            ax1 = plt.subplot(111)
#
            for r in range(len(isrf)):
                plt.plot( wavs[s1][c][r][i], flam[s1][c][r][i], \
                  label = name_ea[0]+": "+str(l_star[s1])+ \
                  r" $L_{\odot}$"+leg_ext[r], color=colors_t[0], \
                  linestyle = linestyle[r], linewidth = linewidth[0])
                plt.plot( wavs[s2][c][r][i], flam[s2][c][r][i], \
                  label = name_ea[1]+": "+str(l_star[s2])+ \
                  r" $L_{\odot}$"+leg_ext[r], color=colors_t[1], \
                  linestyle = linestyle[r], linewidth = linewidth[1])
                if sed_norm:
                    ymax = 50.0 ; ymin = 1.e-12
                    plt.ylabel(r"($\lambda$ F$_{\lambda}$) / F$_{*}$", \
                      fontsize = cs.fontsize, labelpad=0.5)
                else:
                    ymax = 1.e-7 ; ymin = 1.e-12
                    plt.ylabel(r"$\lambda$ F$_{\lambda}$ (erg cm$^{-2}$ s$^{-1}$)", \
                      fontsize = cs.fontsize, labelpad=0.5)
            plt.xlabel("Wavelength ("+(r"$\mu$m")+")", \
              fontsize = cs.fontsize, labelpad=0.5)
            ax1.set_xscale("log")
#            plt.axvline(x=1300.0, linestyle="dashed", color="k")
#            plt.axvline(x=850.0, linestyle="dashed", color ="k")
#            plt.axvline(x=450.0, linestyle="dashed", color="k")
#            plt.axvline(x=350.0, linestyle="dashed", color ="k")
#            plt.axvline(x=250.0, linestyle="dashed", color ="k")
#            plt.axvline(x=70.0, linestyle="dashed", color ="k")
            ax1.set_xlim(250.0,1500.)
            ax1.set_yscale("log") ; ax1.set_ylim( ymin, ymax )
            plt.legend(loc = "upper right", fontsize=cs.leg_fontsize, \
              scatterpoints = 20)
            plt.tight_layout()
            plt.savefig(plt_dir + "/SED_zoom_"+str(name_event[s])+"_"+ \
              str(cut[c])+"_"+str(int(inclins[i]))+"i."+plt_form, \
              format = plt_form)
            plt.clf()
#
    # Print flux ratios at defined wavelengths for single events, for models
    # with ISRF, at different inclinations
#
print("Flux ratios at specified wavelengths for each event, for models with" + \
  " ISRF, for singular cut regime, and singular inclinations are written" + \
  " in /analysis/flam_ratios_ISRF.dat\n")
f_ratios = os.getcwd()+"/../flam_ratios_ISRF.dat"
f_ratios = open(f_ratios, "w")
f_ratios.write("#name_event\tcut\tinclin\twav\tflam_ratio\n")
#
for s in range(len(snaps)/2):
    s1 = s*2 ; s2 = (s*2)+1
#
    for c in range(len(cut)):
#
        for i in range(len(inclins)):
#
            for w in range(len(wavs_compare)):
                flam_1 = float(flam_vals[s1][c][1][i][w])
                flam_2 = float(flam_vals[s2][c][1][i][w])
                flam_ratio = flam_2 / flam_1
#
                print(''' For {0} ({1} cut),
                    ratio of {2} micron flux
                    between i={3} degree
                    "quiescent" and "accreting"
                    snapshots for models
                    with ISRF is: {4}\n'''.format( \
                    name_event[s], cut[c], wavs_compare[w], inclins[i], \
                    flam_ratio ) )
                f_ratios.write('''{0}\t{1}\t{2}\t{3}\t{4}\n'''.format( \
                  name_event[s], cut[c], inclins[i], wavs_compare[w], \
                  flam_ratio ) )
f_ratios.close()
#
print("Flux ratios at specified wavelengths for each event, for models without" + \
  " ISRF, for singular cut regime, and singular inclinations are written" + \
  " in /analysis/flam_ratios.dat\n")
f_ratios = os.getcwd()+"/../flam_ratios.dat"
f_ratios = open(f_ratios, "w")
f_ratios.write("#name_event\tcut\tinclin\twav\tflam_ratio\n")
#
for s in range(len(snaps)/2):
    s1 = s*2 ; s2 = (s*2)+1
#
    for c in range(len(cut)):
#
        for i in range(len(inclins)):
#
            for w in range(len(wavs_compare)):
                flam_1 = float(flam_vals[s1][c][0][i][w])
                flam_2 = float(flam_vals[s2][c][0][i][w])
                flam_ratio = flam_2 / flam_1
#
                print(''' For {0} ({1} cut),
                    ratio of {2} micron flux
                    between i={3} degree
                    "quiescent" and "accreting"
                    snapshots for models
                    with ISRF is: {4}\n'''.format( \
                    name_event[s], cut[c], wavs_compare[w], inclins[i], \
                    flam_ratio ) )
                f_ratios.write('''{0}\t{1}\t{2}\t{3}\t{4}\n'''.format( \
                  name_event[s], cut[c], inclins[i], wavs_compare[w], \
                  flam_ratio ) )
f_ratios.close()
#
    # Compare the beam-dependent flux for E1 snapshots with ISRF heating
#
r_beam = [[] for i in range(2)]
lam_flam = [ [[] for i in range(len(wavs_compare))] for i in range(2) ]

for s in range(len(snaps) / 2):
    f_beam = dat_dir+"/"+snaps[s]+"/MULTI/lam_flam_multi_ISRF.dat"
    f_beam = open(f_beam,"r")
    for lines in f_beam:
        lines = lines.strip() ; columns = lines.split()
        r_beam[s].append(float(columns[0]))
        for w in range(len(wavs_compare)):
            lam_flam[s][w].append(float(columns[1+w]))
    f_beam.close()

cols = ["g","r","b","k"]
lines = ["-","--"]
phase = ["Quiescent","Outbursting"]

fig = plt.figure(1)
ax1 = plt.subplot(111)
for s in range(len(snaps) / 2):
    for w in range(len(wavs_compare)):
        if wavs_compare[w] == 350 or wavs_compare[w] == 1300:
            plt.plot(r_beam[s], lam_flam[s][w], \
              color=cols[w], linestyle = lines[s], \
              label = str(phase[s])+": "+str(int(wavs_compare[w]))+r" $\mu$m" )
plt.xlabel("Radius (AU)", fontsize = cs.fontsize, labelpad=0.5)
ax1.set_xlim(0, max(r_beam[0]))
plt.ylabel(r"$\lambda$ F$_{\lambda}$ (erg cm$^{-2}$ s$^{-1}$)", \
  fontsize = cs.fontsize, labelpad=0.5)
ax1.set_yscale("log")
plt.legend(loc = "center right", fontsize=cs.leg_fontsize)
plt.tight_layout()
plt.savefig(plt_dir+"/lamflam_beam"+"."+plt_form, format = plt_form)
plt.clf()



exit()
