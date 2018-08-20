'''

inp_gen.py

Module to generate most .inp files as required for RADMC-3D, including:
       - wavelength_micron.inp
       - external_source.inp
       - stars.inp (with/without seeded planet)
       - radmc3d.inp

Last Modified: 29/06/2018

'''

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - MODULE IMPORTS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

import os
import shutil
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

import constants as cs
from params_run import *

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - MAIN PROGRAM - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

def gen(arch_dir, dat_dir, plt_dir):

    global t_star, r_star, l_star, t_planet, r_planet, isrf_mult

    wav_bins = 200

    # First, write wavelength_micron.inp file, computing wavelengths
    # (in microns) dependent on user defined range

    llambda = []
    f1 = open("{0}/wavelength_micron.inp".format(dat_dir), "w")
    f1.write("{0}\n".format(wav_bins))
    d0 = math.log(wav_lims[0],10.0) ; d1 = math.log(wav_lims[1],10.0)
    dinc = (d1 - d0) / float(wav_bins)

    for i in range(wav_bins):
        llambda.append( 10.0**( d0 + dinc*float(i) ) )
        f1.write("{0}\n".format(llambda[i]) )
    f1.close()

### ------------------------------------------------------------------------ ###

    # If called for, generate the ISRF file external_source.inp, interpolating
    # Mathis+ (1983) data with required unit conversion user defined wavelength
    # grid as required. If not used, remove from directory

    if not incl_isrf and os.path.isfile(dat_dir+"/external_source.inp"):
        os.remove(dat_dir+"/external_source.inp")

    elif incl_isrf:

        f1 = open(dat_dir+"/external_source.inp","w")
        wav, fnu = np.loadtxt(arch_dir+"/isrf/andre03.txt", unpack=True)

        for i in range(len(fnu)):
            if ((wav[i] * 1000.0) >= 500):
                continue
            else:
                fnu[i] *= isrf_mult

    # Interpolate onto specified wavelength grid, write to file and plot if
    # called. Flux assumed to be zero if ISRF wavelength out of wav_lims range

        fnu_int = interp1d(wav,fnu,kind="linear")
        f1.write("2\n")
        f1.write("{0}\n".format(wav_bins) )
        for i in range(wav_bins):
            f1.write("{0}\n".format(llambda[i]) )
        for i in range(wav_bins):
            if (llambda[i] < min(wav)) or (llambda[i] > max(wav)):
                f1.write("0\n")
            else:
                f1.write("{0}\n".format(fnu_int(llambda[i])) )
        f1.close()

        lam_flam_isrf = []
        for w in range(len(wav)):
            lam_flam_isrf.append( fnu[w] * \
              ( cs.c_cgs / (wav[w] * cs.cm_per_micron)  ) )
#            lam_flam_isrf.append(fnu[w] / (cs.Jy_cgs * 1.0e6))
        if plot_isrf:
            fig = plt.figure(1)
            ax1 = plt.subplot(111)
            plt.plot(wav, lam_flam_isrf, color = "k")
            plt.xlabel("Wavelength ("+(r"$\mu$m")+")", fontsize = cs.fontsize, \
               labelpad=0.5)
            plt.xticks(fontsize = 15) ;   ax1.set_xscale("log")
            ax1.set_xlim(7.0e-2,max(wav))
            ax1.set_yscale("log")
            ax1.set_ylim( 3.e-5, max(lam_flam_isrf)*1.2)
#            ax1.set_ylim( 1.0e-2, 1.0e3)
            plt.ylabel(r"$\lambda I_\lambda$ (erg cm$^{-2}$ s$^{-1}$ sr$^{-1}$)", \
              fontsize = cs.fontsize)
#            plt.ylabel(r"$I_\nu$ (MJy sr$^{-1}$)", fontsize = cs.fontsize)
            plt.xticks(fontsize = cs.fontsize)
            plt.yticks(fontsize = cs.fontsize)
            plt.tight_layout()
            plt.savefig(plt_dir+"/ISRF."+plt_form, format=str(plt_form))
            plt.clf()

### ------------------------------------------------------------------------ ###

    # If protostellar properties are file read, identify snapshot time with
    # rdisc.1, read in sink properties vs. time, and identify protostellar
    # radius and temperature for snapshot time

    if f_star and sim_inp:

        frdisc = open(sim_r1file,"r")
        for line in frdisc:
            line = line.strip()
            column = line.split()
            if ( str(column[13]).split(".")[2] == \
              sim_cfile.split("/")[-1].split(".")[2]):
                snaptime = float(column[0]) / 1000.0
                break
            else:
                continue
        frdisc.close()

        fsink = open(sim_sfile,"r")
        s_time = [] ; s_radius = [] ; s_lum = [] ; s_temp = [] ; s_macc = []

        for line in fsink:
            line = line.strip()
            column = line.split()
            s_time.append(float(column[1]))
            s_radius.append(float(column[19]))
            s_lum.append(float(column[20]))
            s_temp.append(float(column[21]))
            s_macc.append(float(column[25]))
            if ( round(float(column[1]) * 1000.0, 3) == round(snaptime, 3) ):
                r_star = float(column[19]) * cs.rsol_per_pc
                t_star = round(float(column[21]))
                l_star = float(column[20])
            else:
                continue
        fsink.close()

    f = open(dat_dir+"/stars.inp","w")
    f.write("2 \n")
    if (run_tag == "EC53"):
        print("Luminosity of central protostar is: {0} L_sol".format(l_star))
    else:
        print("Luminosity of central protostar is: "+ \
         str(round((r_star)**2.0 * (t_star / 5780.)**4.0, 3) )+" L_sol\n")

    if incl_planet:

        if not incl_star:
            print("Stars must be included if planets are, exiting...\n")
            exit()

        print("Luminosity of seeded planet is: "+ \
        str(round((r_planet)**2.0 * (t_planet / 5780.)**4.0, 3) )+" L_sol\n")

        f.write("2 {0}\n".format(wav_bins) )
        f.write("{0} 1.998e+33 0. 0. 0.\n".format(r_star * cs.cm_per_rsol) )
        f.write("{0} 1.998+30 {1} {2} {3}\n".format(r_planet * cs.cm_per_rsol, \
          loc_planet[0] * cs.cm_per_au, loc_planet[1] * cs.cm_per_au, \
          loc_planet[2] * cs.cm_per_au,))
        for i in range(wav_bins):
            f.write("{0}\n".format(llambda[i]))
        f.write("-{0}\n".format(t_star))
        f.write("-{0}\n".format(t_planet))
        f.close()
    else:
        f.write("1 {0}\n".format(wav_bins) )
        f.write("{0} 1.998e+33 0. 0. 0.\n".format(r_star * cs.cm_per_rsol) )
        for i in range(wav_bins):
            f.write("{0}\n".format(llambda[i]))
        f.write("-{0}\n".format(int(t_star)))
        f.close()

    # From protostellar temperatures (and isrf multipliers where relevant),
    # make data pointers

    if incl_isrf and (isrf_mult <= 0.0):
        print "\nVariable isrf_mult must be a positive float. Setting to 1.0\n"

    if (isrf_mult < 1.0):
        isrf_ext = str(isrf_mult).split(".")[0]+str(isrf_mult).split(".")[1]
        isrf_ext += "ISRF"
    elif (isrf_mult == 1.0):
        isrf_ext = "ISRF"
    elif (isrf_mult > 1.0):
        isrf_ext = str(int(isrf_mult))+"ISRF"

    if incl_star and incl_isrf:

        if (run_tag == "EC53"):
            plt_ext = str(int(l_star))+"L_"+isrf_ext
            tdat_ext = plt_ext+".dat"
            pdat_ext = plt_ext+".out"
        else:
            plt_ext = str(int(t_star))+"_"+isrf_ext
            tdat_ext = plt_ext+".dat"
            pdat_ext = plt_ext+".out"

    elif incl_star and not incl_isrf:

        if (run_tag == "EC53"):
            plt_ext = str(int(l_star))+"L"
            tdat_ext = plt_ext+".dat"
            pdat_ext = plt_ext+".out"

        else:
            plt_ext = str(int(t_star))
            tdat_ext = plt_ext+".dat"
            pdat_ext = plt_ext+".out"

    elif not incl_star and incl_isrf:
        plt_ext = "ISRF"
        tdat_ext = plt_ext+".dat"
        pdat_ext = plt_ext+".out"

    if not incl_star:
        if not os.path.isdir(dat_dir+"/bin"):
            os.makedirs(dat_dir+"/bin")
        shutil.move(dat_dir+"/stars.inp", dat_dir+"/bin/stars.inp")

### ------------------------------------------------------------------------ ###

    f = open(dat_dir+"/radmc3d.inp","w")
    f.write("nphot = {0}\n".format(nphot) )
    f.write("lines_mode = 1\n")
    if mod_rw:
        f.write("modified_random_walk = 1 \n")
    f.write("tgas_eq_tdust = 1 \n")
    f.write("debug_write_stats = 1\n")
    if (threads > 0):
        f.write("setthreads = {0}\n".format(threads))
    if incl_star:
        f.write("scattering_mode_max = 1\n")
        f.write("camera_incl_stars = 1\n")
    elif not incl_star:
        f.write("scattering_mode_max = 0\n")
        f.write("camera_incl_stars = 0\n")
    f.close()

### ------------------------------------------------------------------------ ###

    # Check if protostellar template exists in /protostar/ directory. If not,
    # generate SED with MCRT and RRT over dummy input files.

    f_pstar = arch_dir+"/protostar/spectrum"+str(int(t_star))+".out"
    path_pstar = dat_dir+"/protostar"
    if not (os.path.isfile(f_pstar)):
        import grid_gen
        os.makedirs(path_pstar)
        os.chdir(path_pstar+"/")
        shutil.copy2(dat_dir+"/radmc3d.inp", path_pstar+"/")
        shutil.copy2(dat_dir+"/wavelength_micron.inp", path_pstar+"/")

        grid_gen.protostar(path_pstar, wav_bins, llambda, t_star, r_star)

        os.system("{0}radmc3d sed incl 0".format(exec_loc) )
        shutil.copy2(path_pstar+"/spectrum.out", \
          arch_dir+"/protostar/spectrum"+str(int(t_star))+".out")
        os.chdir(arch_dir+"/src/")
        os.system("rm -rf {0}/".format(path_pstar))

    # Now do the same, but for any planets.

    if incl_planet:

        f_pstar = arch_dir+"/protostar/spectrum"+str(int(t_planet))+".out"
        if not (os.path.isfile(f_pstar)):
            import grid_gen
            os.makedirs(path_pstar)
            os.chdir(path_pstar+"/")
            shutil.copy2(dat_dir+"/radmc3d.inp", path_pstar+"/")
            shutil.copy2(dat_dir+"/wavelength_micron.inp", path_pstar+"/")

            grid_gen.protostar(path_pstar, wav_bins, llambda, \
              t_planet, r_planet)

            os.system("{0}radmc3d sed incl 0".format(exec_loc) )
            shutil.copy2(path_pstar+"/spectrum.out", \
              arch_dir+"/protostar/spectrum"+int(t_planet)+".out")
            os.chdir(arch_dir+"/src/")
            os.system("rm -rf {0}/".format(path_pstar))

    else:
        t_planet = 0 ; r_planet = 0

    return plt_ext, tdat_ext, pdat_ext, t_star, r_star, t_planet, \
      r_planet, llambda, wav_bins
