'''

diagnose.py

Module to evaluate model characteristics. At present, includes functions to
evaluate:

  - Mass within a core centred on [0,0,0]

Author: Benjamin MacFarlane
Date: 02/07/2018
Contact: bmacfarlane@uclan.ac.uk

'''

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - MODULE IMPORTS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

import numpy as np
import shutil

import constants as cs
from params_run import *

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - MAIN PROGRAM - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

def mass(dat_dir):

    global mass_rcore

    # Read in dust_density.inp

    f = open(dat_dir+"/dust_density.inp","r")
    trash = f.readline()
    ngrid = int(f.readline())
    nrspecies = int(f.readline())
    rho = [[] for s in range(nrspecies)]
    for s in range(nrspecies):
        for j in range(ngrid):
            rho[s].append(float(f.readline()))
    f.close()

    # Read in loc_grid.dat

    xloc, yloc, zloc, rloc, xwid = \
      np.loadtxt(dat_dir+"/loc_grid.dat", unpack = True)

    # If no relevant area set, simply use full model volume

    if not mass_rcore:
        mass_rcore = max(rloc) * cs.au_per_cm

    # Now loop over dust species in rho array, summating mass in relevant volume

    mass_sum = 0.0
    for s in range(nrspecies):
        for j in range(ngrid):
            if (rloc[j] * cs.au_per_cm <= mass_rcore):
                mass_sum += ( rho[s][j] * (2.0 * xwid[j])**3.0)

    mass_sum *= ( cs.msol_per_g * cs.gas_to_dust )

    print("\nGas mass in {0} AU core is {1} M_sol".format(mass_rcore, mass_sum))

### ------------------------------------------------------------------------ ###

def lum(dat_dir, tdat_ext, plt_ext):

    global lum_rcore, wav_probe

    max_r = lum_rcore

    # Estimate luminosity from emissivity of tau = 1 -> domain edge

    # First, read in temeperature data

    f = open(dat_dir+"/dust_temperature"+tdat_ext,"r")
    trash = f.readline()
    ngrid = int(f.readline())
    nrspecies = int(f.readline())
    temp = [[] for s in range(nrspecies)]
    for s in range(nrspecies):
        for j in range(ngrid):
            temp[s].append(float(f.readline()))
    f.close()
    temp = np.array(temp)

    # Now, read loc_grid.dat, finding maximum radial limits

    xloc, yloc, zloc, rloc, xwid = \
     np.loadtxt(dat_dir+"/loc_grid.dat", unpack = True)
    rloc *= cs.au_per_cm

    # Copy data file under inspection to radmc3d input, ready for tau test

    shutil.copy2(dat_dir+"/dust_temperature"+tdat_ext, \
     dat_dir+"/dust_temperature.dat")

    # Loop over inclinations

    L = [0.0 for i in range(len(inclins))]

    for i in range(len(inclins)):

        sdat_ext = plt_ext+"_"+str(inclins[i])+"i.out"

    # If no probing wavelenght set, read and identifying SED peak wavelength

        if not wav_probe:

            wav = [] ; lam_flam = []
            f = open(dat_dir+"/spectrum"+sdat_ext)
            for j in range(3):
                header = f.readline()
            for lines in f:
                lines = lines.strip() ; columns = lines.split()
                wav.append(float(columns[0]))
                lam_flam.append(  float(columns[1]) * \
                 ( cs.c_cgs / (float(columns[0])*cs.cm_per_micron) ) )
            f.close()

            wav = np.array(wav) ; lam_flam = np.array(lam_flam)

            wav_probe = np.average(wav, weights = lam_flam)

    # Now run tau test

        npix = 1000
        os.system('''{0}radmc3d tausurf 1 lambda {1} incl {2} \\
            zoomau {3} {4} {3} {4} npix {5} secondorder'''.format(exec_loc, \
            wav_probe, inclins[i], - max_r, max_r, npix) )

    # Find average radial location from resultant tausurface_3d

        r_im  = []
        f = open(dat_dir+"/image.out","r")
        for l in range(6):
            trash = f.readline()
        for y in range(npix):
            yloc = -max_r + (max_r / npix) + y*(2*max_r / npix)

            for x in range(npix):

                xloc = -max_r + (max_r / npix) + x*(2*max_r / npix)

                z_tau = float(f.readline())

                if (z_tau == -1.e91):
                    continue
                else:
                    r_im.append( np.sqrt(xloc**2.0 + yloc**2.0 + \
                     (z_tau * cs.au_per_cm)**2.0) )
        f.close()
        r_tau = np.mean(r_im)

    # Now compute average temperature of regions from
    # tau = 1 surface to domain limits

        count = 0.0 ; temp_ave = 0.0
        for s in range(nrspecies):
            for j in range(ngrid):
                if (rloc[j] > r_tau) and (rloc[j] < max_r):
                    count += 1.0
                    temp_ave += temp[s][j]
        if (count > 0):
            temp_ave /= count
        else:
            print "No temperatures in shell found.."

    # With radial limits, and average temperature, estimate
    # luminosity of emissive envelope

        L[i] = ( (max_r - r_tau) * cs.m_per_au )**2.0 * (temp_ave**4.0)
        L[i] /= ( (cs.m_per_rsol)**2.0 * (5780.0)**4.0 )


    for i in range(len(inclins)):
        print('''\n Luminosity of i={0} deg. model \
is {1} L_sol\n'''.format(inclins[i], L[i]))

    return
