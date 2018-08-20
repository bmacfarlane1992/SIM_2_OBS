'''

evolution.py

Program to read in SED data for EVOLUTION models, and compare between
    - Single snapshots for multiple inclinations
    - Multiple snapshots, for inclination limits

Author: Benjamin MacFarlane
Date: 19/01/2018
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
# User modules
from constants import *
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
dat_dir = os.getcwd()+"/../../../runs/EVOLUTION/"
plt_dir = os.getcwd()+"/../plots/"
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
    fname = dat_dir+snaps[i]+"/stars.inp"
    f = open(fname, 'r')
    trash = f.readline() ; trash = f.readline()
    r_star.append(float(f.readline().split(' ')[0]) * rsol_per_cm)
    t_star.append(-float(f.readlines()[-1]))
    f.close()
#
    fname = dat_dir+snaps[i]+"/spectrum"+str(int(t_star[i]))+".out"
    f = open(fname, 'r')
    for l in range(3):
        trash = f.readline()
    for lines in f:
        lines = lines.strip() ; columns = lines.split()
        wavs_s[i].append(float(columns[0]))
        flam_s[i].append( float(columns[1]) * \
           ( c_cgs / (float(columns[0])*cm_per_micron) ) )
    f.close()
    if sed_norm:
        fnorm = ( (r_star[i] * cm_per_rsol)**2. * sigma_sb_cgs * \
           (t_star[i])**4. ) / (cm_per_pc)**2.
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
        fname = dat_dir+snaps[i]+"/spectrum_"+str(inclins[j])+"i.out"
        f = open(fname, 'r')
        for l in range(3):
            trash = f.readline()
        for lines in f:
            lines = lines.strip() ; columns = lines.split()
            if (j == 0):
                wavs[i].append(float(columns[0]))
            flam[i][j].append( float(columns[1]) * \
               ( c_cgs / (float(columns[0])*cm_per_micron) ) )
        f.close()
#
        if sed_norm:
            fnorm = ( (r_star[i] * cm_per_rsol)**2. * sigma_sb_cgs * \
               (t_star[i])**4. ) / (cm_per_pc)**2.
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
    plt.savefig(plt_dir + "SED_"+snaps[i]+".png") ; plt.clf()
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
plt.savefig(plt_dir + "SED_compare.png") ; plt.clf()
