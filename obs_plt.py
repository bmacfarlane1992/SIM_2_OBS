'''

obs_plt.py

Module to plot synthetic observations for either:

    - Continuum dust images, with radmc3dPy modules
    - Spectral energy distributions

Last Modified: 27/02/2018

'''

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - MODULE IMPORTS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

import imp
import os
import shutil
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


try:
    imp.find_module("radmc3dPy")
    import radmc3dPy
except:
    print("radmc3dPy modules not found for image plots...\n")
    plot_image = False

import constants as cs
from params_run import *

from functions import Trapezoidal

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - MAIN PROGRAM - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

def image(plt_dir, inclin):

    if not plot_image:
        return

    for w in range(len(im_wavs)):

    # First, make image through MCRT run, and read into RADMC-3D Python package

        os.system('''{0}radmc3d image lambda {1} incl {2} posang {3} zoomau \
          -{4} {4} -{4} {4} npix 1000'''.format(exec_loc, im_wavs[w], inclin, \
          pos_ang, im_size))
        imag = radmc3dPy.image.readImage()

    # Output to .fits files, which can be processed in CASA, or plot in X window

        if fits_out:
            fits_name = plt_dir+"/IM_"+str(im_wavs[w])+"micron_"+ \
              str(inclin)+"i.fits"
            os.remove(fits_name)
            imag.writeFits(fits_name, dpc=dist, coord="03h10m05s -10d05m30s")

        elif not fits_out:
            fits_name = ""

        # If called, convolve image with Gaussian with FWHM and position angle
        # defined in variables

            if im_conv is True:
                global im_conv_FWHM, im_conv_posang
                imag = imag.imConv(fwhm=[im_conv_FWHM[0], im_conv_FWHM[1]], \
                   pa=im_conv_posang, dpc=dist)

    # Finally, plot in either angular or physical scale, in S_nu

            if im_ang:
                radmc3dPy.image.plotImage(imag, arcsec=True, \
                   dpc=dist, bunit="snu", log=True, maxlog=5, \
                   cmap = plt.cm.gist_heat)
            elif not im_ang:
                radmc3dPy.image.plotImage(imag, au=True, \
                   dpc=dist, bunit="snu", log=True, maxlog=5, \
                   cmap = plt.cm.gist_heat)

    return fits_name

### ------------------------------------------------------------------------ ###

def sed(arch_dir, dat_dir, sdat_ext, plt_dir, plt_ext, inclin, \
  t_star, r_star, t_planet, r_planet, beam_r):

    wav = [] ; nu = []
    fnu = [] ; flam = [] ; lam_flam = []
    flux_bol = 0.0 ; flux_submm = 0.0

    f = open(dat_dir+"/spectrum"+sdat_ext,"r")

    # Read in data, converting F_nu to F_lambda and lambda F_lambda.
    # lambda F_lambda normalised to observer distance, F_lambda at 1 pc

    for j in range(3):
        header = f.readline()
    for lines in f:
        lines = lines.strip() ; columns = lines.split()
        wav.append(float(columns[0]))
        fnu.append(float(columns[1]))
        nu.append(cs.c_cgs / (float(columns[0]) * cs.cm_per_micron) )
        flam.append(float(columns[1]) * \
           ( cs.c_cgs / (float(columns[0]) * cs.cm_per_micron)**2.0 ) )
        lam_flam.append(float(columns[1]) * \
           ( cs.c_cgs / (float(columns[0]) * cs.cm_per_micron) ) / dist**(2.0))
    f.close()
    nu = nu[::-1]

    if incl_star and (mod_inp != "protostar"):
        shutil.copy2(arch_dir+"/protostar/spectrum"+str(int(t_star))+".out",
          dat_dir+"/")
        wav_star = [] ; nu_star = []
        fnu_star = [] ; flam_star = [] ; lam_flam_star = []
        f = open(dat_dir+"/spectrum"+str(int(t_star))+".out","r")
        for j in range(3):
            header = f.readline()
        for lines in f:
            lines = lines.strip() ; columns = lines.split()
            wav_star.append(float(columns[0]))
            nu_star.append(cs.c_cgs / (float(columns[0]) * cs.cm_per_micron) )
            fnu_star.append(float(columns[1]))
            flam_star.append(float(columns[1]) * \
              ( cs.c_cgs / (float(columns[0]) * cs.cm_per_micron)**2.0 ) )
            lam_flam_star.append( float(columns[1]) * \
              ( cs.c_cgs / (float(columns[0]) * cs.cm_per_micron) ) \
              / dist**(2.0))
        f.close()
        nu_star = nu_star[::-1]

        if incl_planet:
            shutil.copy2(arch_dir+"/protostar/spectrum"+str(int(t_planet))+ \
              ".out", dat_dir+"/")
            wav_planet = [] ; nu_planet = []
            fnu_planet = [] ; flam_planet = [] ; lam_flam_planet = []
            f = open(dat_dir+"/spectrum"+str(int(t_planet))+".out","r")
            for j in range(3):
                header = f.readline()
            for lines in f:
                lines = lines.strip() ; columns = lines.split()
                wav_planet.append(float(columns[0]))
                nu_planet.append(cs.c_cgs / (float(columns[0]) * cs.cm_per_micron) )
                fnu.append(float(columns[1]))
                flam_planet.append(float(columns[1]) * \
                  ( cs.c_cgs / (float(columns[0]) * cs.cm_per_micron)**2.0 ) )
                lam_flam_planet.append( float(columns[1]) * \
                  (cs.c_cgs / (float(columns[0]) * cs.cm_per_micron)) / \
                  dist**(2.0))
            f.close()
            nu_planet = nu_planet[::-1]

    # Dependent on sed_norm, either set limits, then compute total flux
    # from central source to normalise SED, or set flam limits based on
    # absolute values

    if sed_norm:
        norm = ( (r_star * cs.cm_per_rsol)**2. * cs.sigma_sb_cgs * \
          (t_star)**4. ) / (cs.cm_per_pc)**2.
        for j in range(len(lam_flam)):
            lam_flam[j] = lam_flam[j] / norm
        ymax = max(lam_flam) * 10.0

        if incl_star and (mod_inp != "protostar"):
            for j in range(len(lam_flam_star)):
                lam_flam_star[j] = lam_flam_star[j] / norm

    elif not sed_norm:
        ymax = max(lam_flam) * 10.0

    ymin = ymax / 1.0e9

    fig = plt.figure(1)
    ax1 = plt.subplot(111)
    plt.plot(wav, lam_flam, color = "k", linewidth = 1, label = "Full Model")

    if incl_star and (mod_inp != "protostar"):
        plt.plot(wav_star, lam_flam_star, color = "g", linewidth = 1, \
          label="Protostar - "+str(int(t_star))+"K")
        if incl_planet:
            plt.plot(wav_planet, lam_flam_planet, color = "r", linewidth = 1, \
               label = "Planet - "+str(int(t_planet))+"K")
        plt.legend(loc = "upper left", fontsize=10)

    plt.xlabel("Wavelength ("+(r"$\mu$m")+")", fontsize = 18, labelpad=0.5)
    plt.xticks(fontsize = 15) ;   ax1.set_xscale("log")
    ax1.set_xlim(5.e-3,10000.)
    if sed_norm:
        plt.ylabel(r"$\lambda$ F$_{\lambda}$ / F", fontsize = 18, labelpad=0.5)
    elif not sed_norm:
        plt.ylabel(r"$\lambda$ F$_{\lambda}$", fontsize = 18, labelpad=0.5)
    ax1.set_yscale("log") ; ax1.set_ylim( ymin , ymax)
    plt.tight_layout()
    plt.savefig(plt_dir+"/SED_"+plt_ext+"_"+str(inclin)+"i."+plt_form, \
      format = plt_form)
    plt.clf()

    # Compute flux values for [350, 450, 850, 1300] microns

    flux_values(wav, fnu, beam_r)

    # Compute diagnostics for YSO classification from SED data

    yso_class(nu, fnu, nu_star, fnu_star, t_star, r_star)

    return

### ------------------------------------------------------------------------ ###

def flux_values(wav, fnu, beam_r):

    # First, interpolate fnu values

    fnu_int = interp1d(wav, fnu, kind = "cubic")

#    wavs = [70, 350, 450, 850, 1300]
    wavs = [450, 850]

    for i in range(len(wavs)):

        fnu_cgs = fnu_int(wavs[i])


        if (wavs[i] == 450):
            beam_r = dist * 7.9 * cs.cm_per_au
        elif (wavs[i] == 850):
            beam_r = dist * 13.0 * cs.cm_per_au

        delt_omega = (math.pi * beam_r**2.0) / (dist*cs.cm_per_pc)**2.0
        fnu_Jy = fnu_cgs * delt_omega / cs.Jy_cgs

#        delt_omega = (math.pi * beam_r**2.0) / (dist*cs.cm_per_pc)**2.0
#        fnu_Jy = fnu_cgs * delt_omega / cs.Jy_cgs

        print('''Flux density at {0} microns is {1} Jy/beam
        '''.format(wavs[i], fnu_Jy) )

    return

### ------------------------------------------------------------------------ ###

def yso_class(nu, fnu, nu_star, fnu_star, t_star, r_star):

    # Set trapezium rule for numerical integration

    int_sample = int(1e5)

    # Protostellar bolometric luminosity

    fnu_star = fnu_star[::-1]
    flux_int = interp1d(nu_star, fnu_star, kind="cubic")
    a = min(nu_star) ; b = max(nu_star)
    def g(lam):
        return flux_int(lam)
    result = Trapezoidal(g, a, b, int_sample)
    L_star = (4.0 * math.pi * (cs.cm_per_pc**2.0) * result) / cs.Lsol_cgs

    # Model bolometric luminosity

    fnu = fnu[::-1]
    flux_int = interp1d(nu, fnu, kind="cubic")
    a = min(nu) ; b = max(nu)
    def g(lam):
        return flux_int(lam)
    result = Trapezoidal(g, a, b, int_sample)
    L_bol = (4.0 * math.pi * (cs.cm_per_pc**2.0) * result) / cs.Lsol_cgs

    # Model sub-mm luminosity

    a = min(nu) ; b = cs.c / 350.e-6
    def g(lam):
        return flux_int(lam)
    result = Trapezoidal(g, a, b, int_sample)
    L_submm = (4.0 * math.pi * (cs.cm_per_pc**2.0) * result) / cs.Lsol_cgs

    # Print values to terminal

    ratio_submm_bol = (L_submm / L_bol) * 100.0
    print("\nProtostellar luminosity (SED) is: {0} L_sol\n".format(L_star))
    print("Bolometric luminosity is: {0} L_sol\n".format(L_bol))
    print("Sub-mm luminosity is: {0} L_sol\n".format(L_submm))
    print("Ratio of Bolometric and sub-mm "+ \
      "luminosity is {0} %\n".format(ratio_submm_bol) )

    # Compute bolometric temperature for model SED

    flux_int = interp1d(nu, fnu, kind="cubic")
    a = min(nu) ; b = max(nu)
    def g(lam):
        return flux_int(lam)
    result1 = Trapezoidal(g, a, b, int_sample)

    nu_fnu = []
    for w in range(len(nu)):
        nu_fnu.append(fnu[w] * nu[w])

    energy_int = interp1d(nu, nu_fnu, kind = "cubic")
    def g(lam):
        return energy_int(lam)
    result2 = Trapezoidal(g, a, b, int_sample)

    T_bol = 1.25e-11 * (result2 / result1)

    print("Bolometric temperature is: {0} K\n".format(T_bol))

    return
