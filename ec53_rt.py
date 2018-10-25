''''

ec53_rt.py

Module to rescale simulation mass for EC53 study to 5.80 M_sol, then allocate
different dust species to cavity, envelope, disc and disc atmosphere. Writes
to dust_density.inp, dustkappa_XXXX.inp and dustopac.inp files.

Setup adopts approach/parameters of Baek, G. 2D modelling, with opacities
called for relevant structural components. Briefly,

CAVITY_theta = 20

R_star1 (disc scale) = 2.09 R_sol
T_star = 4000 K
L_star    = !!!FREE PARAMETER!!! [L_sol]
R_star2 (input) = R_sol * sqrt(L_star) * (T_sol / T_star)**2
CAVITY_rho = 1.00e-17
R_disk_in = 14.25 * R_star1

Module generates and transfers:
   - stars.inp
   - wavelength_micron.inp
   - external_source.inp
   - dustopac.inp
   - dustkappa_kmhnew_extrap_dust.inp
   - dustkappa_oh5_dust.inp
   - dustkappa_ww04_dust.inp
   - dustkappa_www03_extrap_dust.inp
   - dust_density.inp

i.e. module by itself serves (partial) purposes of inp_gen.py, scale_cut.py
and opac_gen.py modules.

Last Modified: 23/05/2018


'''

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - MODULE IMPORTS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

import numpy as np
import os
import shutil
import math
import matplotlib.pyplot as plt

import constants as cs
from params_run import *

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - MAIN PROGRAM - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

def gen(arch_dir, dat_dir, plt_dir):

    global r_star,  t_star, CAVITY_rho, CAVITY_theta, CAVITY_beta, \
      CAVITY_alpha, CAVITY_ralpha, FULL_mass_scale

    FULL_mass_scale = 5.80  # As per Enoch+ (2009) and Dunham+ (2015) ; (M_sol)

    hfact1 = 1.0            # Float: Scale height multiplier for disc midplane
    hfact2 = 3.0            # Float: As above, for disc atmosphere

    scale_yso = False       # Bool.: Is sim. extent rescaled? If no, is cut.
    R_env = 10000.0         # Float: Radial extent of YSO envelope (AU)

    CAVITY_beta = 1.5       # Float: Cavity profile exponent
    CAVITY_theta = 20.0     # Float: Cavity opening angle (deg.)
    CAVITY_ralpha = 10000.0 # Float: R where opening angle is achieved (AU)
    CAVITY_rho = 1.e-19     # Float: Paramaters as per G. Baek 2D modelling

    flat_cavity = False     # Bool.: Is cavity density constant? (If no, taper)
    R_taper = 100.0         # Float: Radius (AU) where taper occurs, if called


### ------------------------------------------------------------------------ ###


    # Read wavelength_micron.inp

    shutil.copy2(arch_dir+"/isrf/dust_opacity_and_ISRF/wavelength_micron.inp", \
      dat_dir)

    f = dat_dir+"/wavelength_micron.inp"
    f = open(f,"r")
    wavs = []
    n_wavs = int(f.readline())
    for i in range(n_wavs):
        wavs.append(float(f.readline()))
    f.close()


### ------------------------------------------------------------------------ ###


    # Transfer ISRF data, apply relevant flux rescaling, and plot if called for.
    # Delete from dat_dir if run doesn't include ISRF

    if incl_isrf:
        shutil.copy2(arch_dir+ \
          "/isrf/dust_opacity_and_ISRF/external_source.inp", dat_dir)

        isrf_wavs = [] ; isrf_fnu = []
        f = open(dat_dir+"/external_source.inp","r")
        trash = f.readline()
        isrf_n_wavs = int(f.readline())
        for i in range(isrf_n_wavs):
            isrf_wavs.append(float(f.readline()))
        for i in range(isrf_n_wavs):
            isrf_fnu.append(float(f.readline()))
        f.close()

        for i in range(len(isrf_fnu)):
            if ((isrf_wavs[i] * 1000.0) >= 500):
                continue
            else:
                isrf_fnu[i] *= isrf_mult

        f = open(dat_dir+"/external_source.inp","w")
        f.write("2\n")
        f.write("{0}\n".format(isrf_n_wavs))
        for i in range(isrf_n_wavs):
            f.write("{0}\n".format(isrf_wavs[i]))
        for i in range(isrf_n_wavs):
            f.write("{0}\n".format(isrf_fnu[i]))
        f.close()

        isrf_lam_flam = []
        for w in range(isrf_n_wavs):
            isrf_lam_flam.append( isrf_fnu[w] * \
              ( cs.c_cgs / (isrf_wavs[w] * cs.cm_per_micron)  ) )

        if plot_isrf:

            fig = plt.figure(1)
            ax1 = plt.subplot(111)
            plt.plot(isrf_wavs, isrf_lam_flam, color = "k")
            plt.xlabel("Wavelength ("+(r"$\mu$m")+")", fontsize = cs.fontsize, \
               labelpad=0.5)
            plt.xticks(fontsize = 15) ;   ax1.set_xscale("log")
            ax1.set_xlim(0.1,10000)
            ax1.set_yscale("log")
            ax1.set_ylim( 3.e-5, max(lam_flam_isrf)*1.2)
            plt.ylabel(r"$\lambda I_\lambda$ (erg cm$^{-2}$ s$^{-1}$ sr$^{-1}$)", \
              fontsize = cs.fontsize)
            plt.xticks(fontsize = cs.fontsize)
            plt.yticks(fontsize = cs.fontsize)
            plt.tight_layout()
            plt.savefig(plt_dir+"/ISRF."+plt_form, format=str(plt_form))
            plt.clf()

    elif not incl_isrf and os.path.isfile(dat_dir+"/external_source.inp"):
        os.remove(dat_dir+"/external_source.inp")


### ------------------------------------------------------------------------ ###


    # Transfer all opacity files, and plot if called for.

    opac_files = ["kmhnew_extrap_dust","oh5_dust","ww04_dust", \
     "www003_extrap_dust"]

    shutil.copy2(arch_dir+"/isrf/dust_opacity_and_ISRF/"+ \
      "dustopac.inp", dat_dir)
    for i in range(len(opac_files)):
        shutil.copy2(arch_dir+"/isrf/dust_opacity_and_ISRF/"+ \
         "dustkappa_"+opac_files[i]+".inp", dat_dir)

    wavs_opac = [[] for i in range(4)]
    kappa_abs_opac = [[] for i in range(4)]
    kappa_scat_opac = [[] for i in range(4)]
    kappa_tot = [[] for i in range(4)]

    for i in range(len(opac_files)):

        f = open(arch_dir+"/isrf/dust_opacity_and_ISRF/dustkappa_"+ \
         opac_files[i]+".inp", "r")
        trash = f.readline() ; trash = f.readline()
        for lines in f:
            lines = lines.strip() ; columns = lines.split()
            wavs_opac[i].append(float(columns[0]))
            kappa_abs_opac[i].append(float(columns[1]))
            kappa_scat_opac[i].append(float(columns[2]))
            kappa_tot[i].append(float(columns[1])+float(columns[2]))
        f.close()

    if opac_plot:

        print "here!!"
        colors = ["k","b","g","r"]
        leg_tag = ["Cavity","Envelope","Disc Midplane","Disc Atmosphere"]
        fig = plt.figure(1)
        ax1 = plt.subplot(111)
        for i in range(4):
            plt.plot(wavs_opac[i], kappa_tot[i], color = colors[i], \
             linestyle = "-", label = leg_tag[i])
        plt.xlabel("Wavelength ("+(r"$\mu$m")+")", fontsize = cs.fontsize, \
           labelpad=0.5)
        ax1.set_xscale("log"); ax1.set_yscale("log")
        plt.ylabel(r"Total Dust Opacity, $\kappa_\nu$ (cm$^{2}$ g$^{-1}$)", \
          fontsize = cs.fontsize)
        ax1.set_ylim(1.e-2, 6.e5)
        ax1.set_xlim(1.e-1,1.e3)
        plt.xticks(fontsize = cs.fontsize)
        plt.yticks(fontsize = cs.fontsize)
        legend1 = plt.legend(loc = "upper right", fontsize = cs.leg_fontsize)
        plt.gca().add_artist(legend1)

        plt.tight_layout()
        plt.savefig(plt_dir+"/Opacity."+plt_form, format=str(plt_form))
        plt.clf()


### ------------------------------------------------------------------------ ###


    # Generate stars.inp, computing protostar radius from luminosity.
    # From computed protostellar properties, compute dust destruction radius.

    f = open(dat_dir+"/stars.inp", "w")
    r_star = np.sqrt(l_star) * (5780.00 / 4000.0)**2.00
    f.write("2\n")
    f.write("1 {0}\n".format(n_wavs) )
    f.write("{0} 1.998e+33 0. 0. 0.\n".format(r_star * cs.cm_per_rsol) )
    for i in range(n_wavs):
        f.write("{0}\n".format(wavs[i]))
    f.write("-{0}\n".format(int(t_star)))
    f.close()


    R_DESTR = (r_star / 2.0) * (t_star / 1000.0 )**(6.0/2.0)
    R_DESTR *= (cs.cm_per_rsol / cs.cm_per_au)

    if (R_DESTR < 1.0):
        R_DESTR = 1.0


### ------------------------------------------------------------------------ ###

    # Set how YSO extent is adjusted to R_env. Either scale entire simulation
    # by ratio of simulation to parameter extent, or simply cut simulation
    # envelope extent.

    xloc, yloc, zloc, rloc, xwid = \
       np.loadtxt(dat_dir+"/loc_grid.dat",unpack=True)
    rloc = rloc * cs.au_per_cm
    ngrid = len(rloc)

    if scale_yso:

        # Rescale amr_grid.inp limits, if entire snpashot scaled

        f = open(dat_dir+"/amr_grid.inp","r")
        lines = f.readlines()
        lim = float(lines[7].split()[1])

        scale = R_env / (lim * cs.au_per_cm)

        for i in range(3):
            lines[i+7] = "-"+str(lim * scale)+" "+str(lim * scale)+"\n"
        f.close()

        f = open(dat_dir+"/amr_grid.inp","w")
        f.writelines(lines)
        f.close()

        # Apply rescale to both loc_grid.dat, and dust_density.inp

        xloc = xloc * scale
        yloc = yloc * scale
        zloc = zloc * scale
        rloc = (rloc * cs.cm_per_au) * scale
        xwid = xwid * scale
        f = open(dat_dir+"/loc_grid.dat","w")
        for i in range(ngrid):
            f.write("{0} {1} {2} {3} {4}\n".format(xloc[i], yloc[i], zloc[i], \
             rloc[i], xwid[i]))
        f.close()

        rho = np.loadtxt(dat_dir+"/dust_density.inp", skiprows=3, unpack=True)
        rho = rho * scale**(-3.0)

    elif not scale_yso:

        # Simply set density exterior to R_env to zero. Re-read cut data to

        rho = np.loadtxt(dat_dir+"/dust_density.inp", skiprows=3, unpack=True)

        f = open(dat_dir+"/dust_density.inp","w")
        f.write("1\n")
        f.write("{0}\n".format(ngrid))
        f.write("1\n")
        for i in range(ngrid):

            if (rloc[i] > R_env):
                f.write("0\n")
            else:
                f.write("{0}\n".format(rho[i]))

        f.close()


### ------------------------------------------------------------------------ ###


    # Set up new dust_density.inp file

    f = open(dat_dir+"/dust_density.inp","w")
    f.write("1\n")
    f.write("{0}\n".format(ngrid) )
    f.write("4\n")

    # Re-read loc_grid.dat, to begin setting component boundaries

    xloc, yloc, zloc, rloc, xwid = \
       np.loadtxt(dat_dir+"/loc_grid.dat",unpack=True)
    zloc *= cs.au_per_cm
    rloc *= cs.au_per_cm
    rloc_xy = np.hypot(xloc, yloc) * cs.au_per_cm

    # Compute cavity boundary, writing to dust_density.inp as required

    CAVITY_theta = CAVITY_theta * (math.pi / 180.0)
    CAVITY_alpha = CAVITY_ralpha / \
     (CAVITY_ralpha * math.tan(CAVITY_theta) )**CAVITY_beta

    for i in range(ngrid):

        if (rloc[i] > R_env):
            f.write("0\n")
            continue

        zcav = CAVITY_alpha * rloc_xy[i]**(CAVITY_beta)
        if (abs(zloc[i]) > zcav):

    # Write different cavity densities, depending on tapering method

            if flat_cavity:

                if (rho[i] < CAVITY_rho):
                    f.write("{0}\n".format(rho[i]))
                else:
                    f.write("{0}\n".format(CAVITY_rho) )

            elif not flat_cavity:

                if (rloc[i] <= R_taper):

                    if (rho_tmp > rho[i]):
                        f.write("{0}\n".format(rho[i]))
                    else:
                        f.write("{0}\n".format(CAVITY_rho))

                elif (rloc[i] > R_taper):
                    rho_tmp = CAVITY_rho * (R_taper / rloc[i])**(2.0)

                    if (rho_tmp > rho[i]):
                        f.write("{0}\n".format(rho[i]))
                    else:
                        f.write("{0}\n".format(rho_tmp))

        else:
            f.write("0\n")


    # Use rdisc files to compute radial extent, and scale height of disc

    if not (os.path.isfile(sim_pfile) or os.path.isfile(sim_r1file)):
        print("rdisc.1 or pdisc.1 files can't be found, exiting...")
        exit()

    frdisc = open(sim_r1file, "r")
    for line in frdisc:
        line = line.strip()
        column = line.split()
        if ( str(column[13]).split(".")[2] == \
          sim_cfile.split("/")[-1].split(".")[2] ):
            DISC_r_cut = float(column[5])
            print("Radial restiction of disc, set by surface density "+ \
              "restriction is: {0} AU\n".format(float(column[5]) ) )
            print("Mass of simulation disc (radial infall restriction) "+ \
              "is: {0} M_sol\n".format( float(column[6]) ) )
            break
        else:
            continue
    frdisc.close()

    scale_height = 0.0
    count = 0
    fpdisc = open(sim_pfile, 'r')
    header = fpdisc.readline()
    for line in fpdisc:
        line = line.strip()
        column = line.split()
        r = float(column[1])
        if ( (int(r) != 0) and (r < DISC_r_cut) ) :
            h = np.sqrt( (cs.k_b * float(column[3]) ) / \
              (cs.mu*cs.m_p) ) / (float(column[13]) * 1000.)
            if (h != 'nan'):
                count += 1
                scale_height += h

    DISC_h_cut = hfact1 * (scale_height / float(count) )
    print('''Scale height restriction of disc midplane, from {0} * average h(r)
     value for disc < DISC_r_cut is: {1}\n'''.format(hfact1, DISC_h_cut) )

    ATMOSPHERE_h_cut = hfact2 * (scale_height / float(count) )
    print('''Scale height restriction of disc atmosphere, from {0} * average
      h(r) value for disc < DISC_r_cut is: {1}\n'''.format(hfact2, \
      ATMOSPHERE_h_cut) )

    fpdisc.close()

    # Now set envelope densities for locations not in cavity, or in either
    # disc dust Species

    for i in range(ngrid):

        if (rloc[i] > R_env):
            f.write("0\n")
            continue

        zcav = CAVITY_alpha * rloc_xy[i]**(CAVITY_beta)
        if (abs(zloc[i]) < zcav):

            if (abs(zloc[i] / rloc_xy[i]) > ATMOSPHERE_h_cut) or \
             (rloc_xy[i] > DISC_r_cut):
                f.write("{0}\n".format(rho[i]) )
            else:
                f.write("0\n")
        else:
            f.write("0\n")

    # Now for disc midplane

    for i in range(ngrid):

        if (rloc[i] > R_env):
            f.write("0\n")
            continue

        zcav = CAVITY_alpha * rloc_xy[i]**(CAVITY_beta)
        if (abs(zloc[i]) < zcav):

            if (abs(zloc[i] / rloc_xy[i]) < DISC_h_cut):

                if (rloc_xy[i] < DISC_r_cut):
                    f.write("{0}\n".format(rho[i]))
                else:
                    f.write("0\n")
            else:
                f.write("0\n")
        else:
            f.write("0\n")

    # Now for disc atmosphere

    for i in range(ngrid):

        if (rloc[i] > R_env):
            f.write("0\n")
            continue

        zcav = CAVITY_alpha * rloc_xy[i]**(CAVITY_beta)
        if (abs(zloc[i]) < zcav):

            if (abs(zloc[i] / rloc_xy[i]) > DISC_h_cut) and \
             (abs(zloc[i] / rloc_xy[i]) < ATMOSPHERE_h_cut):

                if (rloc_xy[i] < DISC_r_cut):
                    f.write("{0}\n".format(rho[i]) )
                else:
                    f.write("0\n")
            else:
                f.write("0\n")
        else:
            f.write("0\n")
    f.close()


### ------------------------------------------------------------------------ ###


    # Now read new dust_density.inp with all four dust species

    f = open(dat_dir+"/dust_density.inp","r")
    trash = f.readline()
    ngrid = int(f.readline())
    nrspecies = int(f.readline())
    rho = [[] for i in range(nrspecies)]
    for i in range(nrspecies):
        for j in range(ngrid):
            rho[i].append(float(f.readline()))
    f.close()

    # Rescale model, setting total cavity mass constant, and changing density
    # for grid bins attributed to "other" components.
    # Either way, cut 1 AU cavity from central region to remove artificially
    # high density regions around protostar

    cavity_mass = 0.00
    other_mass = 0.00

    for i in range(nrspecies):

        for j in range(ngrid):

            mass_bin = rho[i][j] * (2.0*xwid[j])**3.0

            if (i == 0):
                cavity_mass += mass_bin
            else:
                other_mass += mass_bin

    other_scale = ( (FULL_mass_scale * cs.dust_to_gas) - \
     (cavity_mass * cs.msol_per_g) ) / (other_mass * cs.msol_per_g )

# Print a few model scaling/RT diagnostics to terminal

    print("\nAdopting dust destruction radius of {0} AU\n".format(R_DESTR))

    print("Gas mass model scaled to: {0} M_sol".format(FULL_mass_scale))
    print("Non-cavity gas mass: {0} M_sol".format(other_mass * \
     cs.msol_per_g / cs.dust_to_gas))
    print("Cavity gas mass: {0} M_sol".format(cavity_mass * \
     cs.msol_per_g / cs.dust_to_gas))
    print("Non-cavity mass rescale ratio: {0}\n".format(other_scale))

    f = open(dat_dir+"/dust_density.inp","w")
    f.write("1\n")
    f.write("{0}\n".format(ngrid))
    f.write("{0}\n".format(nrspecies))
    for i in range(nrspecies):
        for j in range(ngrid):

            if (i == 0):
                f.write("{0}\n".format(rho[i][j]))
            else:
                f.write("{0}\n".format(rho[i][j] * other_scale) )
    f.close()

    return
