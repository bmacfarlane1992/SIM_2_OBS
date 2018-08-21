'''

inp_gen.py

Module to generate most .inp files as required for RADMC-3D, including:
       - wavelength_micron.inp
       - external_source.inp
       - stars.inp (with/without seeded planet)
       - radmc3d.inp

Last Modified: 21/08/2018

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
import params_run as ps

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - MAIN PROGRAM - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

def gen(globals, rt_params):

    rt_params.n_wavs = 200

    # First, write wavelength_micron.inp file, computing wavelengths
    # (in microns) dependent on user defined range

    f1 = open("{0}/wavelength_micron.inp".format(globals.dat_dir), "w")
    f1.write("{0}\n".format(rt_params.n_wavs))
    d0 = math.log(ps.wav_lims[0],10.0)
    d1 = math.log(ps.wav_lims[1],10.0)
    dinc = (d1 - d0) / float(rt_params.n_wavs)

    for i in range(rt_params.n_wavs):
        rt_params.wavs.append( 10.0**( d0 + dinc*float(i) ) )
        f1.write("{0}\n".format( rt_params.wavs[i]) )
    f1.close()

### ------------------------------------------------------------------------ ###

    # If called for, generate the ISRF file external_source.inp, interpolating
    # Mathis+ (1983) data with required unit conversion user defined wavelength
    # grid as required. If not used, remove from directory

    if not ps.incl_isrf and \
     os.path.isfile(globals.dat_dir+"/external_source.inp"):

        os.remove(globals.dat_dir+"/external_source.inp")

    elif ps.incl_isrf:

        f1 = open(globals.dat_dir+"/external_source.inp","w")
        isrf_wavs, isrf_fnu = \
         np.loadtxt(globals.arch_dir+"/isrf/andre03.txt", unpack=True)

        for i in range(len(isrf_fnu)):

            if ((isrf_wavs[i] * 1000.0) >= 500):
                continue
            else:
                isrf_fnu[i] *= ps.isrf_mult

    # Interpolate onto specified wavelength grid, write to file and plot if
    # called. Flux assumed to be zero if ISRF wavelength out of wav_lims range

        isrf_fnu_int = interp1d(isrf_wavs, isrf_fnu, kind="linear")
        f1.write("2\n")
        f1.write("{0}\n".format(rt_params.n_wavs) )
        for i in range(rt_params.n_wavs):
            f1.write("{0}\n".format(rt_params.wavs[i]) )
        for i in range(rt_params.n_wavs):
            if (rt_params.wavs[i] < min(isrf_wavs)) or \
             (rt_params.wavs[i] > max(isrf_wavs)):
                f1.write("0\n")
            else:
                f1.write("{0}\n".format(isrf_fnu_int(rt_params.wavs[i])) )
        f1.close()

        if ps.plot_isrf:

            Jy = False

            isrf_lam_flam = []
            for w in range(len(isrf_wavs)):

                if Jy:
                    lam_flam_isrf.append(isrf_fnu[w] / (cs.Jy_cgs * 1.0e6))
                else:
                    lam_flam_isrf.append( isrf_fnu[w] * \
                     ( cs.c_cgs / (isrf_wavs[w] * cs.cm_per_micron)  ) )

            fig = plt.figure(1)
            ax1 = plt.subplot(111)
            plt.plot(isrf_wavs, isrf_lam_flam, color = "k")
            plt.xlabel("Wavelength ("+(r"$\mu$m")+")", fontsize = cs.fontsize, \
               labelpad=0.5)
            plt.xticks(fontsize = 15) ;   ax1.set_xscale("log")
            ax1.set_xlim(7.0e-2,max(isrf_wavs))
            ax1.set_yscale("log")

            if Jy:
                ax1.set_ylim( 1.0e-2, 1.0e3)
                plt.ylabel(r"$I_\nu$ (MJy sr$^{-1}$)", fontsize = cs.fontsize)
            else:
                ax1.set_ylim( 3.e-5, max(isrf_lam_flam)*1.2)
                plt.ylabel(r"$\lambda I_\lambda$ (erg cm$^{-2}$ s$^{-1}$ sr$^{-1}$)", \
                  fontsize = cs.fontsize)

            plt.xticks(fontsize = cs.fontsize)
            plt.yticks(fontsize = cs.fontsize)
            plt.tight_layout()
            plt.savefig(globals.plt_dir+"/ISRF."+ps.plt_form, \
             format=str(ps.plt_form))
            plt.clf()

### ------------------------------------------------------------------------ ###

    # If protostellar properties are file read, identify snapshot time with
    # rdisc.1, read in sink properties vs. time, and identify protostellar
    # radius and temperature for snapshot time

    if ps.f_star and ps.sim_inp:

        frdisc = open(ps.sim_r1file,"r")
        for line in frdisc:
            line = line.strip()
            column = line.split()
            if ( str(column[13]).split(".")[2] == \
              ps.sim_cfile.split("/")[-1].split(".")[2]):
                snaptime = float(column[0]) / 1000.0
                break
            else:
                continue
        frdisc.close()

        fsink = open(ps.sim_sfile,"r")
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
                ps.r_star = float(column[19]) * cs.rsol_per_pc
                ps.t_star = round(float(column[21]))
                ps.l_star = float(column[20])
            else:
                continue
        fsink.close()

    # Set rt_params instances to (potentially) modified luminosity source variables

    rt_params.t_star = ps.t_star ; rt_params.r_star = ps.r_star
    rt_params.l_star = ps.l_star
    rt_params.t_planet = ps.t_planet ; rt_params.r_planet = ps.r_planet

    # Write to stars.inp file, with relevant parameters for included source

    f = open(globals.dat_dir+"/stars.inp","w")
    f.write("2 \n")
    if (ps.run_tag == "EC53"):
        print("Luminosity of central protostar is: {0} L_sol".format(rt_params.l_star))
    else:
        print("Luminosity of central protostar is: "+ \
         str(round((rt_params.r_star)**2.0 * \
          (rt_params.t_star / 5780.)**4.0, 3) )+" L_sol\n")

    if ps.incl_planet:

        if not ps.incl_star:
            print("Stars must be included if planets are, exiting...\n")
            exit()

        print("Luminosity of seeded planet is: "+ \
        str(round((rt_params.r_planet)**2.0 * \
         (rt_params.t_planet / 5780.)**4.0, 3) )+" L_sol\n")

        f.write("2 {0}\n".format(rt_params.n_wavs) )
        f.write("{0} 1.998e+33 0. 0. 0.\n".format(rt_params.r_star * cs.cm_per_rsol) )
        f.write("{0} 1.998+30 {1} {2} {3}\n".format(rt_params.r_planet * cs.cm_per_rsol, \
          ps.loc_planet[0] * cs.cm_per_au, ps.loc_planet[1] * cs.cm_per_au, \
          ps.loc_planet[2] * cs.cm_per_au,))
        for i in range(rt_params.n_wavs):
            f.write("{0}\n".format(rt_params.wavs[i]))
        f.write("-{0}\n".format(rt_params.t_star))
        f.write("-{0}\n".format(rt_params.t_planet))
        f.close()
    else:
        f.write("1 {0}\n".format(rt_params.n_wavs) )
        f.write("{0} 1.998e+33 0. 0. 0.\n".format(rt_params.r_star * cs.cm_per_rsol) )
        for i in range(rt_params.n_wavs):
            f.write("{0}\n".format(rt_params.wavs[i]))
        f.write("-{0}\n".format(int(rt_params.t_star)))
        f.close()

    # From protostellar temperatures (and isrf multipliers where relevant),
    # make data pointers

    if ps.incl_isrf and (ps.isrf_mult <= 0.0):
        print "\nVariable isrf_mult must be a positive float. Setting to 1.0\n"

    if (ps.isrf_mult < 1.0):
        isrf_ext = str(ps.isrf_mult).split(".")[0]+str(ps.isrf_mult).split(".")[1]
        isrf_ext += "ISRF"
    elif (ps.isrf_mult == 1.0):
        isrf_ext = "ISRF"
    elif (ps.isrf_mult > 1.0):
        isrf_ext = str(int(ps.isrf_mult))+"ISRF"

    if ps.incl_star and ps.incl_isrf:

        if (ps.run_tag == "EC53"):
            globals.plt_ext = str(int(ps.l_star))+"L_"+isrf_ext
            globals.tdat_ext = globals.plt_ext+".dat"
            globals.pdat_ext = globals.plt_ext+".out"
        else:
            globals.plt_ext = str(int(ps.t_star))+"_"+isrf_ext
            globals.tdat_ext = globals.plt_ext+".dat"
            globals.pdat_ext = globals.plt_ext+".out"

    elif ps.incl_star and not ps.incl_isrf:

        if (ps.run_tag == "EC53"):
            globals.plt_ext = str(int(ps.l_star))+"L"
            globals.tdat_ext = globals.plt_ext+".dat"
            globals.pdat_ext = globals.plt_ext+".out"

        else:
            globals.plt_ext = str(int(ps.t_star))
            globals.tdat_ext = globals.plt_ext+".dat"
            globals.pdat_ext = globals.plt_ext+".out"

    elif not ps.incl_star and ps.incl_isrf:
        globals.plt_ext = "ISRF"
        globals.tdat_ext = globals.plt_ext+".dat"
        globals.pdat_ext = globals.plt_ext+".out"

    if not ps.incl_star:
        if not os.path.isdir(globals.dat_dir+"/bin"):
            os.makedirs(globals.dat_dir+"/bin")
        shutil.move(globals.dat_dir+"/stars.inp", globals.dat_dir+"/bin/stars.inp")

### ------------------------------------------------------------------------ ###

    f = open(globals.dat_dir+"/radmc3d.inp","w")
    f.write("nphot = {0}\n".format(ps.nphot) )
    f.write("lines_mode = 1\n")
    if ps.mod_rw:
        f.write("modified_random_walk = 1 \n")
    f.write("tgas_eq_tdust = 1 \n")
    f.write("debug_write_stats = 1\n")
    if (ps.threads > 0):
        f.write("setthreads = {0}\n".format(ps.threads))
    if ps.incl_star:
        f.write("scattering_mode_max = 1\n")
        f.write("camera_incl_stars = 1\n")
    elif not ps.incl_star:
        f.write("scattering_mode_max = 0\n")
        f.write("camera_incl_stars = 0\n")
    f.close()

### ------------------------------------------------------------------------ ###

    # Check if protostellar template exists in /protostar/ directory. If not,
    # generate SED with MCRT and RRT over dummy input files.

    f_pstar = globals.arch_dir+"/protostar/spectrum"+str(int(ps.t_star))+".out"
    path_pstar = globals.dat_dir+"/protostar"
    if not (os.path.isfile(f_pstar)):
        import grid_gen
        os.makedirs(path_pstar)
        os.chdir(path_pstar+"/")
        shutil.copy2(globals.dat_dir+"/radmc3d.inp", path_pstar+"/")
        shutil.copy2(globals.dat_dir+"/wavelength_micron.inp", path_pstar+"/")

        grid_gen.protostar(path_pstar, rt_params)

        os.system("{0}radmc3d sed incl 0".format(ps.exec_loc) )
        shutil.copy2(path_pstar+"/spectrum.out", \
          globals.arch_dir+"/protostar/spectrum"+str(int(rt_params.t_star))+".out")
        os.chdir(globals.arch_dir+"/src/")
        os.system("rm -rf {0}/".format(path_pstar))

    # Now do the same, but for any planets.

    if ps.incl_planet:

        f_pstar = globals.arch_dir+"/protostar/spectrum"+str(int(rt_params.t_planet))+".out"
        if not (os.path.isfile(f_pstar)):
            import grid_gen
            os.makedirs(path_pstar)
            os.chdir(path_pstar+"/")
            shutil.copy2(globals.dat_dir+"/radmc3d.inp", path_pstar+"/")
            shutil.copy2(globals.dat_dir+"/wavelength_micron.inp", path_pstar+"/")

            grid_gen.protostar(path_pstar, rt_params.n_wavs, rt_params.wavs, \
              rt_params.t_planet, rt_params.r_planet)

            os.system("{0}radmc3d sed incl 0".format(ps.exec_loc) )
            shutil.copy2(path_pstar+"/spectrum.out", \
              globals.arch_dir+"/protostar/spectrum"+int(rt_params.t_planet)+".out")
            os.chdir(globals.arch_dir+"/src/")
            os.system("rm -rf {0}/".format(path_pstar))

    return globals, rt_params
