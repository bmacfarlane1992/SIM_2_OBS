'''

opac_gen.py

Module to generate opacity tables, in a variety of ways. Uses wavelength range
and binning as per inp_gen.py module. For generation of dustkappa_XXXX.inp and
dustopac.inp files, first option: use a single opacity table for the entire
model. Second option: Split into components of protostellar system as"

   - Disc defined as H < H_crit = h(r) / r, with currently only one
      absorbing/scattering dust population
   - Envelope, (H > H_crit = h(r) / r) with:
     --- Inner sphere, in which no ice mantles reside
     --- Intermediate sphere, with ice mantles
     --- Outer sphere, with thick ice mantles
Disc opacities read from /opac directory from file named dustkappa_XXXX,
where as envelope opacities used from Ossenkopf & Henning (1994) Jena databases.
Currently species are spatially distinct, i.e. there are no thin/thick ice
mantle regions.

Last Modified: 21/08/2018

'''

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - MODULE IMPORTS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

import numpy as np
import math
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

import constants as cs
import params_run as ps

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - MAIN PROGRAM - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

def gen(globals, rt_params):

    str_id = ps.opacity_inp.split("_")[0]
    if (str_id == "grain") or (str_id == "ice") or (str_id == "thickice"):
        dir_prefix = "/opac/oh94"
    else:
        dir_prefix = "/opac"

    wavs_opac = [] ; kappa_abs_opac = [] ; kappa_scat_opac = []
    kappa_tot = []

    fr = open(globals.arch_dir+dir_prefix+ \
     "/dustkappa_"+ps.opacity_inp+".inp","r")

    iformat = int(fr.readline()) ; n_wavs_opac = int(fr.readline())
    for i in range(n_wavs_opac):
        line = fr.readline() ; column = line.split()

# Read opacity data, converting to dust opacities for OH data

        if (ps.opacity_inp != "OH5") and (dir_prefix == "/opac/oh94"):
            wavs_opac.append(float(column[0]))
            kappa_abs_opac.append(float(column[1]) / cs.dust_to_gas)
            kappa_tot.append(kappa_abs_opac[i])
        elif (ps.opacity_inp != "OH5") and (dir_prefix == "/opac"):
            wavs_opac.append(float(column[0]))
            kappa_abs_opac.append(float(column[1]))
            kappa_scat_opac.append(float(column[2]))
            kappa_tot.append(kappa_abs_opac[i] + kappa_scat_opac[i])
        elif (ps.opacity_inp == "OH5"):
            wavs_opac.append(float(column[0]))
            kappa_abs_opac.append(float(column[1]) / cs.dust_to_gas)
            kappa_scat_opac.append(float(column[2]) / cs.dust_to_gas)
            kappa_tot.append(kappa_abs_opac[i] + kappa_scat_opac[i])

# Interpolate, and write abs (and scat) opacities to dustkappa_ file, then
# general information of opacities in dustopac.inp

    kappa_abs_int = interp1d(wavs_opac, kappa_abs_opac,kind="linear")
    if ((ps.opacity_inp != "OH5") and (dir_prefix == "/opac")) or \
      (ps.opacity_inp == "OH5") :
        kappa_scat_int = interp1d(wavs_opac, kappa_scat_opac,kind="linear")

    fw = open(globals.dat_dir+"/dustkappa_"+ps.opacity_inp+".inp","w")
    fw.write("{0}\n".format(iformat))
    fw.write("{0}\n".format(rt_params.n_wavs))
    if (dir_prefix == "/opac/oh94"):
        for i in range(rt_params.n_wavs):
            fw.write("{0} {1}\n".format(rt_params.wavs[i], \
              kappa_abs_int(rt_params.wavs[i]) ) )
    else:
        for i in range(rt_params.n_wavs):
            fw.write("{0} {1} {2}\n".format(rt_params.wavs[i], \
             kappa_abs_int(rt_params.wavs[i]), \
             kappa_scat_int(rt_params.wavs[i]) ) )
    fr.close() ; fw.close()

    f1 = open(globals.dat_dir+"/dustopac.inp","w")
    f1.write("2               Format number of this file\n")
    f1.write("1               Nr of dust species\n")
    f1.write("==========================================================\n")
    f1.write("1               Way in which this dust species is read\n")
    f1.write("0               0=Thermal grain, 1=Quantum heated\n")
    f1.write('''{0}
        Extension of name of dustkappa_XXXX file\n'''.format(ps.opacity_inp) )
    f1.write("----------------------------------------------------------\n")
    f1.close()

    if ps.opac_plot:

        fig = plt.figure(1)
        ax1 = plt.subplot(111)
        plt.plot(wavs_opac, kappa_abs_opac, color = "k", \
          linestyle = "-", label = r"Absorption, $\kappa_\lambda$")
        if (dir_prefix != "/opac/oh94"):
            plt.plot(wavs_opac, kappa_scat_opac, color = "k", \
              linestyle = "--", label = r"Scattering, $k_\lambda$")
        plt.xlabel("Wavelength ("+(r"$\mu$m")+")", fontsize = cs.fontsize, \
           labelpad=0.5)
        ax1.set_xscale("log"); ax1.set_yscale("log")
        plt.ylabel(r"Dust Opacity (cm$^{2}$ g$^{-1}$)", fontsize = cs.fontsize)
        ax1.set_xlim(7.e-2,1.e4)
        plt.xticks(fontsize = cs.fontsize)
        plt.yticks(fontsize = cs.fontsize)
        plt.legend(loc = "upper right")
        plt.tight_layout()
        plt.savefig(globals.plt_dir+"/Opacity."+ps.plt_form, \
         format=str(ps.plt_form))
        plt.clf()
