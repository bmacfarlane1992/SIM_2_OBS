'''

envelope.py

Routine to check both SED and temperature profiles of RT model computed for:

(a) Envelope with central protostar RT
(b) Envelope with ISRF RT
(c) Envelope with both central protostar and ISRF RT

Protostar is 1 R_sol, 6000K, and envelope is 2 M_sol, extending
between 100-10000 AU. RT adopts ISRF of Andre et al. (2003), and opacity
tables based on OH5 (Ossenkopf & Henning, 1994) column of thin ice mantles
at gas density of 10^6 g/cm^3.


Author: Benjamin MacFarlane
Date: 29/01/2018
Contact: bmacfarlane@uclan.ac.uk

'''
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - MODULE IMPORTS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Standard Python modules
import os
import matplotlib.pyplot as plt
#
from constants import *
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - VARIABLE DEFINITIONS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
names = ["PROTOSTAR","ISRF","BOTH"]
color = ["b","r","c"]
#
plt_form = "png"    # Str.: Output plot fomat ["png","eps"]
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
		# # # - - - MAIN PROGRAM - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
#
cwd = os.getcwd()

arch_dir = cwd + "/../.."	# Project architecture directory
run_dir = arch_dir + "/runs/ENVELOPE"
plt_dir = run_dir+"/plots_analysis"
#
    # Initialise arrays for reading in data
#
wav = [[] for i in range(len(names))]
flam = [[] for i in range(len(names))]
rloc = [[] for i in range(len(names))]
temp = [[] for i in range(len(names))]
#
for i in range(len(names)):
    path = run_dir + "/ENVELOPE_"+names[i]+"/dat"
#
    # Check run directories actually exist. If not, exit.
#
    if not (os.path.isdir(path)):
        print("Run _{0} does not exist. Code exiting...".format(names[i]) )
        exit()
#
    fsed = path + "/spectrum.out"
    floc = path + "/loc_grid.dat"
    ftemp = path + "/dust_temperature.dat"
#
    # Begin reading in data (f1 = SED ; f2 = radius ; f3 = temperature)
#
    f1 = open(fsed,"r") ; f2 = open(floc,"r") ; f3 = open(ftemp,"r")
#
    for j in range(0, 3):
        header = f1.readline()
    for lines in f1:
        lines = lines.strip() ; columns = lines.split()
        wav[i].append(float(columns[0]))
        flam[i].append( float(columns[1]) * \
          ( c_cgs / (float(columns[0])*cm_per_micron) ) )
    f1.close()
#
    for lines in f2:
        lines = lines.strip() ; columns = lines.split()
        rloc[i].append(float(columns[3]) / cm_per_au)
    f2.close()
#
    for j in range(3):
        trash = f3.readline()
    for lines in f3:
        temp[i].append(float(lines))
    f3.close()
#
    # Now plot SED then T(r), looping over $names$ to plot in $color$,
    # consistent with DS tests
#
    # Before plotting SED, read in protostellar spectrum
#
f_star = arch_dir+"/protostar/spectrum6000.out"
if not(os.path.isfile(f_star)):
    print("Protostellar data not present, so can't plot :(")
else:
    wav_star = [] ; flam_star = [] ; color_star = "g"
    f1 = open(f_star, "r")
    for j in range(0, 3):
        header = f1.readline()
    for lines in f1:
        lines = lines.strip() ; columns = lines.split()
        wav_star.append(float(columns[0]))
        flam_star.append( float(columns[1]) * \
          ( c_cgs / (float(columns[0])*cm_per_micron) ) )
    f1.close()
#
    # Now plot SED, with or without protostellar data
#
fig = plt.figure(1)
ax1 = plt.subplot(111)
for i in range(len(names)):
    plt.plot( wav[i], flam[i], label = names[i], color=color[i] )
ymax = max(flam[2])*10. ; ymin = ymax / 5.e6
if os.path.isfile(f_star):
    plt.plot( wav_star, flam_star, label = "Protostar", color = color_star )
plt.xlabel("Wavelength ("+(r"$\mu$m")+")", fontsize = 18, labelpad=0.5)
plt.xticks(fontsize = 15) ;   ax1.set_xscale("log") ; ax1.set_xlim(1.e-2,2500.)
plt.ylabel(r"$\lambda$ F$_{\lambda}$", fontsize = 18, labelpad=0.5)
ax1.set_yscale("log") ; ax1.set_ylim( ymin, ymax )
plt.legend(loc = "upper left", fontsize=10)
plt.tight_layout()
plt.savefig(plt_dir+"/SED."+plt_form, format = plt_form)
plt.clf()
#
    # Now plot radial temperature profiles
#
fig = plt.figure(1)
ax1 = plt.subplot(111)
for i in range(len(names)):
    plt.scatter(rloc[i], temp[i], color=color[i])
plt.xlabel("Radius (AU)", fontsize = 18, labelpad=0.5)
plt.xticks(fontsize = 15) ; ax1.set_xlim(0, 1e4)
plt.ylabel("Temperature (K)", fontsize = 18, labelpad=0.5)
ax1.set_ylim(0., 40.0)
plt.tight_layout()
plt.savefig(plt_dir+"/T_r."+plt_form, format = plt_form)
plt.clf()
#
#
exit()
