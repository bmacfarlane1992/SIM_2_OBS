'''

multi_beam.py

Module to run multiple RRT analyses, for decreasing coverage over entire RT
domain. Used to evaluate flux levels attained whilst mimicking smaller beams.

Last Modified: 29/06/2018

'''

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - MODULE IMPORTS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

import numpy as np
import math
import shutil
import os
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

import constants as cs
from params_run import *

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - MAIN PROGRAM - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

def run(dat_dir, tdat_ext, ext):

    # Firstly, create MULTI directory above run_dir, and copy all NOT data here

    f_opac = []

    if run_tag == "EC53":
        not_dir = dat_dir
        multi_dir = dat_dir+"/MULTI"
        f_opac.append(not_dir+"/dustkappa_kmhnew_extrap_dust.inp")
        f_opac.append(not_dir+"/dustkappa_ww04_dust.inp")
        f_opac.append(not_dir+"/dustkappa_www003_extrap_dust.inp")
        f_opac.append(not_dir+"/dustkappa_oh5_dust.inp")

    else:
        multi_dir = dat_dir+"/../../MULTI"
        not_dir = dat_dir+"/../../NOT/dat"
        f_opac.append(not_dir+"/dustkappa_"+opacity_inp+".inp")

    if not os.path.isdir(multi_dir):
        os.makedirs(multi_dir)

    for f in range(len(f_opac)):
        shutil.copy2(f_opac, multi_dir)

    shutil.copy2(not_dir+"/amr_grid.inp", multi_dir)
    shutil.copy2(not_dir+"/dust_density.inp", multi_dir)
    shutil.copy2(not_dir+"/wavelength_micron.inp", multi_dir)
    shutil.copy2(not_dir+"/dustopac.inp", multi_dir)
    shutil.copy2(not_dir+"/radmc3d.inp", multi_dir)
    shutil.copy2(not_dir+"/stars.inp", multi_dir)
    shutil.copy2(not_dir+"/loc_grid.dat", multi_dir)

    # Now read in data

    xloc, yloc, zloc, rloc, xwid = \
       np.loadtxt(multi_dir+"/loc_grid.dat",unpack=True)
    r_max = max(rloc)

    f = open(multi_dir+"/dust_density.inp","r")
    trash = f.readline()
    ngrid = int(f.readline())
    nrspecies = int(f.readline())
    rho = [[] for s in range(nrspecies)]
    for s in range(nrspecies):
        for j in range(ngrid):
            rho[s].append(float(f.readline()))
    f.close()

    for i in range(len(inclins)):

        shutil.copy2(not_dir+"/dust_temperature"+tdat_ext, \
         multi_dir+"/dust_temperature.dat")

        f = open(multi_dir+"/dust_temperature.dat","r")
        trash = f.readline() ; trash = f.readline() ; trash = f.readline()
        temp = [[] for s in range(nrspecies)]
        for s in range(nrspecies):
            for j in range(ngrid):
                temp[s].append(float(f.readline()))
        f.close()

    # Iterate over decreasing radial cuts, and re-write MULTI dataset

        r_cut = []
        lam_flam_cut = [ [] for w in range(len(multi_rrt_wav))]

        f_multi = multi_dir+"/lam_flam_multi"+ext+"_"+str(inclins[i])+"i.dat"
        f_multi = open(f_multi, "w")

        for d in range(n_multi_rrt):

            RRT_cut = r_max - ( d * ( r_max / float(n_multi_rrt) ) )

            r_cut.append(RRT_cut * cs.au_per_cm)
            f_multi.write("{0}\t".format(r_cut[d]))

            f = open(multi_dir+"/dust_temperature.dat","w")
            f.write("1\n")
            f.write("{0}\n".format(ngrid))
            f.write("{0}\n".format(nrspecies))
            for s in range(nrspecies):
                for j in range(ngrid):
                    if (RRT_cut > 0.) and (rloc[j] < RRT_cut):
                        f.write("{0}\n".format(temp[s][j]) )
                    elif (RRT_cut < 0.) and (rloc[j] > RRT_cut):
                        f.write("{0}\n".format(temp[s][j]) )
                    else:
                        f.write("0.0\n")
            f.close()

    # Now run the RRT over the refined dataset, for only the long-wavelengths

            os.chdir(multi_dir+"/")
            if incl_star:
                os.system('''{0}radmc3d sed incl {1}
                 lambdarange 300 1500 nlam 10'''.format(exec_loc, inclins[i]) )
            elif not incl_star:
                os.system('''{0}radmc3d sed incl {1} nostar
                 lambdarange 300 1500 nlam 10'''.format(exec_loc, inclins[i]))

    # With output spectra, interpolate, find required flux values,
    # then write outputs to file

            f = open(multi_dir+"/spectrum.out")
            wav = [] ; fnu = [] ; flam = [] ; lam_flam = []
            for j in range(3):
                header = f.readline()
            for lines in f:
                lines = lines.strip() ; columns = lines.split()
                wav.append(float(columns[0]))
                fnu.append(float(columns[1]))
                flam.append(float(columns[1]) * \
                 ( cs.c_cgs / (float(columns[0]) * cs.cm_per_micron)**2.0 ) )
                lam_flam.append(float(columns[1]) * \
                 ( cs.c_cgs / (float(columns[0]) * cs.cm_per_micron) ) / \
                  dist**(2.0))
            f.close()

            lam_flam_int = interp1d(wav,lam_flam,kind="linear")

            for w in range(len(multi_rrt_wav)):
                lam_flam_cut[w].append(lam_flam_int(float(multi_rrt_wav[w])) )
                f_multi.write("{0}\t".format(lam_flam_cut[w][d]))
            f_multi.write("\n")

    # Finally, plot the radius/flux values for each wavelengths

        cols = ["r","g","b","k","c"]

        fig = plt.figure(1)
        ax1 = plt.subplot(111)
        for w in range(len(multi_rrt_wav)):
            plt.plot(r_cut, lam_flam_cut[w], color=cols[w], \
             label = str(int(multi_rrt_wav[w]) )+r" $\mu$m" )
        plt.xlabel("Radius (AU)", fontsize = 18, labelpad=0.5)
        plt.xticks(fontsize = 15)
        ax1.set_xlim(0, max(r_cut))
        plt.ylabel(r"$\lambda$ F$_{\lambda}$", fontsize = 18, labelpad=0.5)
        ax1.set_yscale("log")
        plt.legend(loc = "upper right", fontsize=10)
        plt.tight_layout()
        plt.savefig(multi_dir+"/mutli_rrt_lamflam"+ext+"_"+str(inclins[i])+"i." \
         +plt_form, format = plt_form)
        plt.clf()

    return
