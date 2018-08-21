'''

obs_plt.py

Module to plot synthetic observations for either:

    - Continuum dust images, with radmc3dPy modules
    - Spectral energy distributions

Last Modified: 20/08/2018

'''

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - MODULE IMPORTS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

    # Import system modules

import imp
import os
import shutil
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

    # (Try to) Import local modules

try:
    imp.find_module("radmc3dPy")
    import radmc3dPy
except:
    print("radmc3dPy modules not found for image plots...\n")
    plot_image = False

import constants as cs
from params_run import *
import functions as fs

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

    # Read in SED, as per sdat_ext, usinf SED_read() function

    file_read = dat_dir + "/spectrum" + sdat_ext
    wav, flam, lam_flam, nu, fnu, nu_fnu = fs.SED_read(file_read)

    # Normalise plotted lam_flam list to user defined distance to source

    lam_flam /= dist**(2.0)

    # Do the same, for protostellar SED, as comparison to RT output

    if incl_star and (mod_inp != "protostar"):

        shutil.copy2(arch_dir+"/protostar/spectrum"+str(int(t_star))+".out",
         dat_dir+"/")

        file_read = dat_dir + "/spectrum"+str(int(t_star)) + ".out"
        wav_star, flam_star, lam_flam_star, \
         nu_star, fnu_star, nu_fnu_star = fs.SED_read(file_read)

        lam_flam_star /= dist**(2.0)

    # If a planet is present in the SED calculation, apply the same read
    # and normalisation to the relevant data

        if incl_planet:

            shutil.copy2(arch_dir+"/protostar/spectrum"+str(int(t_planet))+ \
             ".out", dat_dir+"/")

            file_read = dat_dir + "/spectrum"+str(int(t_planet)) + ".out"
            wav_planet, flam_planet, lam_flam_planet, \
             nu_planet, fnu_planet, nu_fnu_planet = fs.SED_read(file_read)

            lam_flam_planet /= dist**(2.0)

    # Dependent on sed_norm, either set limits, then compute total flux
    # from central source to normalise SED, or set flam limits based on
    # absolute values

    if sed_norm:

        norm = ( (r_star * cs.cm_per_rsol)**2. * cs.sigma_sb_cgs * \
          (t_star)**4. ) / (cs.cm_per_pc)**2.

        lam_flam /= norm

        if incl_star and (mod_inp != "protostar"):
            lam_flam_star /= norm

        ymax = max(lam_flam) * 10.0

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

    yso_class(nu_star, fnu_star, nu, fnu, nu_fnu, t_star, r_star)

    return

### ------------------------------------------------------------------------ ###

def flux_values(wav, fnu, beam_r):

    # First, interpolate fnu values. fnu must be flipped, to account for
    # prior flip to account for integration

    fnu = fnu[::-1]
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

def yso_class(nu_star, fnu_star, nu, fnu, nu_fnu, t_star, r_star):

    # Protostellar bolometric luminosity

    L_star = fs.Lbol_calc(nu_star, fnu_star)

    # Model bolometric luminosity

    L_bol = fs.Lbol_calc(nu, fnu)

    # Ratio of model sub-mm to bolometric luminosity

    L_ratio = fs.L_ratio_calc(nu, fnu)

    # Compute bolometric temperature for model SED

    T_bol = fs.Tbol_calc(nu, fnu, nu_fnu)

    # Print values to terminal

    print("\nProtostellar luminosity (SED) is: {0} L_sol\n".format(L_star))
    print("Bolometric luminosity is: {0} L_sol\n".format(L_bol))
    print("Ratio of Bolometric and sub-mm "+ \
      "luminosity is {0} %\n".format(L_ratio) )
    print("Bolometric temperature is: {0} K\n".format(T_bol))

    return
