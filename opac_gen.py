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

Last Modified: 18/07/2018

'''

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - MODULE IMPORTS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

import numpy as np
import math
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

import constants as cs
from params_run import *

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - MAIN PROGRAM - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

def gen(arch_dir, dat_dir, plt_dir, llambda, wav_bins):

    # For single opacity selection:

    if not multicomp:

        str_id = opacity_inp.split("_")[0]
        if (str_id == "grain") or (str_id == "ice") or (str_id == "thickice"):
            dir_prefix = "/opac/oh94"
        else:
            dir_prefix = "/opac"

        llambda_opac = [] ; kappa_abs_opac = [] ; kappa_scat_opac = []
        kappa_tot = []

        fr = open(arch_dir+dir_prefix+"/dustkappa_"+opacity_inp+".inp","r")

        iformat = int(fr.readline()) ; wav_bins_opac = int(fr.readline())
        for i in range(wav_bins_opac):
            line = fr.readline() ; column = line.split()

    # Read opacity data, converting to dust opacities for OH data

            if (opacity_inp != "OH5") and (dir_prefix == "/opac/oh94"):
                llambda_opac.append(float(column[0]))
                kappa_abs_opac.append(float(column[1]) / cs.dust_to_gas)
                kappa_tot.append(kappa_abs_opac[i])
            elif (opacity_inp != "OH5") and (dir_prefix == "/opac"):
                llambda_opac.append(float(column[0]))
                kappa_abs_opac.append(float(column[1]))
                kappa_scat_opac.append(float(column[2]))
                kappa_tot.append(kappa_abs_opac[i] + kappa_scat_opac[i])
            elif (opacity_inp == "OH5"):
                llambda_opac.append(float(column[0]))
                kappa_abs_opac.append(float(column[1]) / cs.dust_to_gas)
                kappa_scat_opac.append(float(column[2]) / cs.dust_to_gas)
                kappa_tot.append(kappa_abs_opac[i] + kappa_scat_opac[i])

    # Interpolate, and write abs (and scat) opacities to dustkappa_ file, then
    # general information of opacities in dustopac.inp

        kappa_abs_int = interp1d(llambda_opac,kappa_abs_opac,kind="linear")
        if ((opacity_inp != "OH5") and (dir_prefix == "/opac")) or \
          (opacity_inp == "OH5") :
            kappa_scat_int = \
              interp1d(llambda_opac,kappa_scat_opac,kind="linear")

        fw = open(dat_dir+"/dustkappa_"+opacity_inp+".inp","w")
        fw.write("{0}\n".format(iformat))
        fw.write("{0}\n".format(wav_bins))
        if (dir_prefix == "/opac/oh94"):
            for i in range(wav_bins):
                fw.write("{0} {1}\n".format(llambda[i], \
                  kappa_abs_int(llambda[i]) ) )
        else:
            for i in range(wav_bins):
                fw.write("{0} {1} {2}\n".format(llambda[i], \
                  kappa_abs_int(llambda[i]), kappa_scat_int(llambda[i]) ) )
        fr.close() ; fw.close()

        f1 = open(dat_dir+"/dustopac.inp","w")
        f1.write("2               Format number of this file\n")
        f1.write("1               Nr of dust species\n")
        f1.write("==========================================================\n")
        f1.write("1               Way in which this dust species is read\n")
        f1.write("0               0=Thermal grain, 1=Quantum heated\n")
        f1.write('''{0}
            Extension of name of dustkappa_XXXX file\n'''.format(opacity_inp) )
        f1.write("----------------------------------------------------------\n")
        f1.close()

        if opac_plot:
            fig = plt.figure(1)
            ax1 = plt.subplot(111)
            plt.plot(llambda_opac, kappa_abs_opac, color = "k", \
              linestyle = "-", label = r"Absorption, $\kappa_\lambda$")
            if (dir_prefix != "/opac/oh94"):
                plt.plot(llambda_opac, kappa_scat_opac, color = "k", \
                  linestyle = "--", label = r"Scattering, $k_\lambda$")
            plt.xlabel("Wavelength ("+(r"$\mu$m")+")", fontsize = cs.fontsize, \
               labelpad=0.5)
            ax1.set_xscale("log"); ax1.set_yscale("log")
            plt.ylabel(r"Dust Opacity (cm$^{2}$ g$^{-1}$)", \
              fontsize = cs.fontsize)
            ax1.set_xlim(7.e-2,1.e4)
            plt.xticks(fontsize = cs.fontsize)
            plt.yticks(fontsize = cs.fontsize)
            plt.legend(loc = "upper right")
            plt.tight_layout()
            plt.savefig(plt_dir+"/Opacity."+plt_form, format=str(plt_form))
            plt.clf()


        # For opacities selected in specific regions (see module notes)

    elif multicomp:

        for i in range(4):

    # Disc dust population

            if (i == 0):
                llambda_opac = [] ; kappa_abs_opac = [] ; kappa_scat_opac = []
                fr = open(arch_dir+"/opac/dustkappa_"+opacity_multi[0]+".inp", \
                  "r")
                fw = open(dat_dir+"/dustkappa_"+opacity_multi[0]+".inp","w")

                iformat = int(fr.readline())
                wav_bins_opac = int(fr.readline())
                for j in range(wav_bins_opac):
                    line = fr.readline() ; column = line.split()
                    llambda_opac.append(float(column[0]))
                    kappa_abs_opac.append(float(column[1]))
                    kappa_scat_opac.append(float(column[2]))

                kappa_abs_int = \
                   interp1d(llambda_opac,kappa_abs_opac,kind="linear")
                kappa_scat_int = \
                   interp1d(llambda_opac,kappa_scat_opac,kind="linear")

                fw.write(str(iformat)+"\n")
                fw.write(str(wav_bins)+"\n")
                for i in range(wav_bins):
                    fw.write("{0} {1} {2}\n".format(llambda[i], \
                      kappa_abs_int(llambda[i]), \
                      kappa_scat_int(llambda[i]) ) )

    # Envelope dust populations

            else:
                llambda_opac = [] ; kappa_abs_opac = []
                fr = open(arch_dir+"/opac/oh94/"+"dustkappa_oh94_"+ \
                   opacity_multi[i]+".inp","r")
                fw = open(dat_dir+"/dustkappa_"+opacity_multi[i]+".inp","w")

                iformat = int(fr.readline())
                wav_bins_opac = int(fr.readline())
                for j in range(wav_bins_opac):
                    line = fr.readline() ; column = line.split()
                    llambda_opac.append(float(column[0]) / cs.dust_to_gas)
                    kappa_abs_opac.append(float(column[1]) / cs.dust_to_gas)

                kappa_abs_int = \
                   interp1d(llambda_opac,kappa_abs_opac,kind="linear")

                fw.write("{0}\n".format(iformat) )
                fw.write("{0}\n".format(wav_bins) )
                for i in range(wav_bins):
                    fw.write("{0} {1}\n".format(llambda[i], \
                      kappa_abs_int(llambda[i]) ) )

    # Read in location and density data from globally gridded data, converting
    # to AU for envelope boundary and computing H value

        xloc, yloc, zloc, rloc = np.loadtxt(dat_dir+"/loc_grid.dat",unpack=True)
        xloc = xloc * cs.au_per_cm ; yloc = yloc * cs.au_per_cm
        zloc = zloc * cs.au_per_cm ; rloc = rloc * cs.au_per_cm
        r_xyloc = [] ; H = []
        for i in range(len(rloc)):
            r_xyloc.append( math.hypot(xloc[i],yloc[i]) )
            H.append(math.fabs(zloc[i] / r_xyloc[i]))

        rhobin = np.loadtxt(dat_dir+"/dust_density.inp", skiprows=3, \
          unpack=True)

    # Split dust distribution into components, in order of
    # 1: Disc, 2: Inner envelope, 3: Intermediate envelope, 4: Outer envelope
    # then re-write dust_density.inp and write dustopac.inp

        rho_disc = [0.0 for i in range(len(rhobin))]
        rho_env1 = [0.0 for i in range(len(rhobin))]
        rho_env2 = [0.0 for i in range(len(rhobin))]
        rho_env3 = [0.0 for i in range(len(rhobin))]

        for i in range(len(rhobin)):

            if (r_xyloc[i] < disc_lim) and (H[i] < H_crit):
                rho_disc[i] = rhobin[i]
            elif ((rloc[i] < env_bound[0]) and (H[i] > H_crit)) \
               or ((rloc[i] < env_bound[0]) and (H[i] < H_crit) \
               and (r_xyloc[i] > disc_lim)):
                rho_env1[i] = rhobin[i]
            elif ((rloc[i] > env_bound[0]) and (rloc[i] < env_bound[1]) \
               and (H[i] > H_crit)) or ((rloc[i] > env_bound[0]) \
               and (rloc[i] < env_bound[1]) and (H[i] < H_crit) \
               and (r_xyloc[i] > disc_lim)):
                rho_env2[i] = rhobin[i]
            elif ((rloc[i] > env_bound[1]) and (H[i] > H_crit)) \
               or ((rloc[i] > env_bound[1]) and (H[i] < H_crit) \
               and (r_xyloc[i] > disc_lim)):
                rho_env3[i] = rhobin[i]

        f = open(dat_dir+"/dust_density.inp","w")
        f.write("1\n")
        f.write("{0}\n".format(len(rhobin) ) )
        f.write("4\n")
        for i in range(len(rhobin)):
            f.write("{0}\n".format(rho_disc[i]) )
        for i in range(len(rhobin)):
            f.write("{0}\n".format(rho_env1[i]) )
        for i in range(len(rhobin)):
            f.write("{0}\n".format(rho_env2[i]) )
        for i in range(len(rhobin)):
            f.write("{0}\n".format(rho_env3[i]) )
        f.close()

        f1 = open(dat_dir+"/dustopac.inp","w")
        f1.write("2\n")
        f1.write("4\n")
        f1.write("==========================================================\n")
        for i in range(4):
            f1.write("1\n")
            f1.write("0\n")
            f1.write("{0}\n".format(opacity_multi[i]) )
            f1.write("------------------------------------------------------\n")
        f1.close()

    return
