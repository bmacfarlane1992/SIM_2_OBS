'''

inp_plt.py

Module to plot AMR gridded/RADMC-3D input data. Functionality includes

- Temperature/Density profiles, for:
  -- Spherical radial distribtion
  -- Azimuthally averaged distribution
- Slice surface plot

Last Modified: 17/07/2018

'''

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - MODULE IMPORTS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

import os
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib import cm
import scipy.interpolate as scipy_interp
import scipy.ndimage as ndimage

import constants as cs
from params_run import *

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - MAIN PROGRAM - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

def read(dat_dir, plt_ext, tdat_ext, pdat_ext):

    print("Reading data, and generating diagnostic plots for {0} run".format( \
      run_tag) )

    if not os.path.isfile(dat_dir+"/loc_grid.dat"):
        print("loc_grid.dat not found in run_dir, exiting...")
        exit()

    # Read in data from any found dust_density.inp, dust_temperature.dat,
    # photon_statistics.out and loc_grid.dat files, saving to arrays

    xloc, yloc, zloc, rloc, xwid = \
       np.loadtxt(dat_dir+"/loc_grid.dat",unpack=True)

    xloc = xloc * cs.au_per_cm ; yloc = yloc * cs.au_per_cm
    zloc = zloc * cs.au_per_cm ; rloc = rloc * cs.au_per_cm
    xwid = xwid * cs.au_per_cm ; rloc_xy = np.hypot(xloc,yloc)
    h = abs(zloc/rloc_xy)

    # If simulation restricted by disc, loop over grid and include/restrict
    # as necessary. Save indices to then point to array location used in plots

    ngrid = len(rloc)

    xloc_plt = [] ; yloc_plt = [] ; zloc_plt = [] ;
    rloc_plt = [] ; xwid_plt = [] ; rloc_xy_plt = []
    restr = []

    for i in range(ngrid):
        if (disc_h > 0) and \
           ( ((h[i] > disc_h) and (rloc_xy[i] < disc_r)) or \
           ((h[i] < disc_h) and (rloc_xy[i] > disc_r)) or \
           (rloc[i] > disc_r) ):
            restr.append(i)
            xloc_plt.append(xloc[i]) ; yloc_plt.append(yloc[i])
            zloc_plt.append(zloc[i]) ; rloc_plt.append(rloc[i])
            xwid_plt.append(xwid[i]) ; rloc_xy_plt.append(rloc_xy[i])
        elif (disc_h < 0) and ((h[i] < abs(disc_h)) and (rloc_xy[i] < disc_r)):
            restr.append(i)
            xloc_plt.append(xloc[i]) ; yloc_plt.append(yloc[i])
            zloc_plt.append(zloc[i]) ; rloc_plt.append(rloc[i])
            xwid_plt.append(xwid[i]) ; rloc_xy_plt.append(rloc_xy[i])
        elif (disc_h == 0):
            restr.append(i)
            xloc_plt.append(xloc[i]) ; yloc_plt.append(yloc[i])
            zloc_plt.append(zloc[i]) ; rloc_plt.append(rloc[i])
            xwid_plt.append(xwid[i]) ; rloc_xy_plt.append(rloc_xy[i])

    if os.path.isfile(dat_dir+"/dust_density.inp"):
        print("Plotting diagnostics for density\n")
        diag_rho = True

        if os.path.isfile(dat_dir+"/CAVITY_UNCUT_dust_density.inp"):
            f = open(dat_dir+"/CAVITY_UNCUT_dust_density.inp","r")
        else:
            f = open(dat_dir+"/dust_density.inp","r")
        trash = f.readline()
        trash = int(f.readline())   # Is actually ngrid, but just binned here..
        nrspecies = int(f.readline())
        rho = [[] for i in range(nrspecies)]
        rho_plt = [[] for i in range(nrspecies)]
        for i in range(nrspecies):
            for j in range(ngrid):
                rho[i].append(float(f.readline()))
            for j in range(len(restr)):
                rho_plt[i].append(rho[i][restr[j]])
        f.close()
    else:
        rho_plt = [] ; nrspecies = 0 ; diag_rho = False

    if os.path.isfile(dat_dir+"/dust_temperature"+tdat_ext):
        diag_temp = True
        print("Plotting diagnostics for temperature\n")
        f = open(dat_dir+"/dust_temperature"+tdat_ext,"r")
        trash = f.readline() ; trash = f.readline() ; trash = f.readline()
        temp = [[] for i in range(nrspecies)]
        temp_plt = [[] for i in range(nrspecies)]
        for i in range(nrspecies):
            for j in range(ngrid):
                temp[i].append(float(f.readline()))
            for j in range(len(restr)):
                temp_plt[i].append(temp[i][restr[j]])
        f.close()
    else:
        temp_plt = [] ; diag_temp = False

    if os.path.isfile(dat_dir+"/photon_statistics"+pdat_ext):
        diag_nphot = True
        print("Plotting diagnostics for photon statistics\n")
        f = open(dat_dir+"/photon_statistics"+pdat_ext,"r")
        trash = f.readline()
        nphot = [] ; nphot_plt = []
        for i in range(ngrid):
            nphot.append(float(f.readline()))
        for i in range(len(restr)):
            nphot_plt.append(nphot[restr[i]])
        f.close()
    else:
        nphot_plt = [] ; diag_nphot = False

    for i in range(len(inclins)):

        sdat_ext = plt_ext+"_"+str(inclins[i])+"i.out"

        if os.path.isfile(dat_dir+"/spectrum"+sdat_ext):
            diag_sed = True
        else:
            print("{0} SED data not found...".format(sdat_ext))
            diag_sed = False

    rloc_tmp = np.array(rloc_plt) ; rho_tmp = np.array(rho_plt[0])
    restr = np.where(rho_tmp != 0.)
    rloc_tmp = rloc_tmp[restr]
    beam_r = max(rloc_tmp) * cs.cm_per_au

    return xloc_plt, yloc_plt, zloc_plt, rloc_plt, xwid_plt, rloc_xy_plt, \
       nrspecies, rho_plt, diag_rho, temp_plt, diag_temp, \
       nphot_plt, diag_nphot, diag_sed, beam_r

### ------------------------------------------------------------------------ ###

def radial(dat_dir, plt_dir, plt_ext, rloc, rloc_xy, nrspecies, rho, diag_rho, \
  R_DESTR, temp, diag_temp):

    rloc_min = 1. ; rloc_max = 50000.
    rho_min = 1.e-23 ; rho_max = 1.e-12
    temp_min = 1.0 ; temp_max = 5000.

    print("Analysing radial profiles\n")

    color = ["k","b","g","r"]

    # Run data manipulation to generate azimuthally averaged profiles. Radial
    # slices split into 1000 logarithmically spaced bins

    xlabel = "Azimuthal Radius (AU)" ; ext = "azi"

    rmin = min(rloc) ; rmax = max(rloc) ; rinc = (rmax - rmin) / 1000.0
    rloc_azi = [0.0] * 1000

    if diag_rho:
        rho_azi = [ [0.0 for j in range(1000)] for i in range(nrspecies)]
    if diag_temp:
        temp_azi = [ [0.0 for j in range(1000)] for i in range(nrspecies)]

    for i in range(1000):
        r0 = rmin + ( rinc*float(i) ) ; r1 = rmin + (rinc*float(i+1))
        rloc_azi[i] = (r0 + r1) / 2.0
        count = 0
        for j in range(len(rloc)):
            if ( (rloc_xy[j] > r0) and (rloc_xy[j] < r1) ):
                count += 1
                for s in range(nrspecies):

                    if diag_rho:
                        rho_azi[s][i] = rho_azi[s][i] + rho[s][j]
                    if diag_temp:
                        temp_azi[s][i] = temp_azi[s][i] + temp[s][j]

        if (count != 0):
            for s in range(nrspecies):
                if diag_rho:
                    rho_azi[s][i] = rho_azi[s][i] / count
                if diag_temp:
                    temp_azi[s][i] = temp_azi[s][i] / count
        else:
            continue

    if diag_temp:

        rloc_plt = np.array(rloc_azi) ; temp_plt = np.array(temp_azi)

        fig = plt.figure(1)
        ax1 = plt.subplot(111)
        for s in range(nrspecies):
            plt.plot(rloc_plt, temp_plt[s], color=color[s], linewidth = 1)
        plt.xlabel(xlabel, fontsize = 18, labelpad=0.5)
        plt.xticks(fontsize = 15) ;   ax1.set_xscale("log")
        ax1.set_xlim(rloc_min, rloc_max)
        plt.ylabel("Temperature (K)", fontsize = 18, labelpad=0.5)
        ax1.set_yscale("log") ; ax1.set_ylim(temp_min+0.1, temp_max*2.0)
        plt.tight_layout()
        plt.savefig(plt_dir+"/T_"+ext+"_"+plt_ext+"."+plt_form, \
          format = plt_form)
        plt.clf()

    if diag_rho:

        rloc_plt = np.array(rloc_azi) ; rho_plt = np.array(rho_azi)

        fig = plt.figure(1)
        ax1 = plt.subplot(111)
        for s in range(nrspecies):
            plt.plot(rloc_plt, rho_plt[s], color=color[s], linewidth = 1)
        plt.xlabel(xlabel, fontsize = 18, labelpad=0.5)
        plt.xticks(fontsize = 15) ; ax1.set_xscale("log")
        ax1.set_xlim(rloc_min, rloc_max)
        plt.ylabel("Volume Density (g/cm"+r"$^{3}$)", fontsize = 18, \
          labelpad=0.5)
        ax1.set_yscale("log")
        ax1.set_ylim(rho_min/2.0, rho_max*2.0)
        plt.tight_layout()
        plt.savefig(plt_dir+"/rho_"+ext+"."+plt_form, format = plt_form)
        plt.clf()

    # Spherical radius...

    xlabel = "Spherical Radius (AU)" ; ext = "spher"

    if diag_temp:

        rloc_plt = np.array(rloc) ; temp_plt = np.array(temp)

        fig = plt.figure(1)
        ax1 = plt.subplot(111)
        for s in range(nrspecies):
            plt.scatter(rloc_plt, temp_plt[s], color=color[s], s=0.5)
        plt.xlabel(xlabel, fontsize = 18, labelpad=0.5)
        plt.xticks(fontsize = 15) ;   ax1.set_xscale("log")
        ax1.set_xlim(rloc_min, rloc_max)
        plt.ylabel("Temperature (K)", fontsize = 18, labelpad=0.5)
        ax1.set_yscale("log") ; ax1.set_ylim(temp_min, temp_max*2.0)
        plt.tight_layout()
        plt.savefig(plt_dir+"/T_"+ext+"_"+plt_ext+"."+plt_form, \
          format = plt_form)
        plt.clf()

    if diag_rho:

        rloc_plt = np.array(rloc) ; rho_plt = np.array(rho)

        fig = plt.figure(1)
        ax1 = plt.subplot(111)

        leg_tag = ["Cavity","Envelope","Disc Midplane","Disc Atmosphere"]
        for s in range(nrspecies):
            if (run_tag == "EC53"):
                plt.scatter(rloc_plt, rho_plt[s], color=color[s], s=0.5, \
                 label = leg_tag[s])
            else:
                plt.scatter(rloc_plt, rho_plt[s], color=color[s], s=0.5)

        plt.xlabel(xlabel, fontsize = 18, labelpad=0.5)
        plt.xticks(fontsize = 15) ; ax1.set_xscale("log")

        if os.path.isfile(dat_dir+"/CAVITY_UNCUT_dust_density.inp"):
            plt.axvline(x = R_DESTR, linestyle=":", color = "k")

        ax1.set_xlim(rloc_min, rloc_max)
        plt.ylabel("Volume Density (g/cm"+r"$^{3}$)", fontsize = 18, \
          labelpad=0.5)
        ax1.set_yscale("log")
        ax1.set_ylim(rho_min/2.0, rho_max*2.0)

        if (run_tag == "EC53"):
            plt.legend(loc = "upper right", scatterpoints = 20, \
             fontsize = cs.leg_fontsize)

        plt.tight_layout()
        plt.savefig(plt_dir+"/rho_"+ext+"."+plt_form, format = plt_form)
        plt.clf()

    return

### ------------------------------------------------------------------------ ###

def run_slice(plt_dir, plt_ext, xloc, yloc, zloc, xwid, nrspecies, \
  rho, diag_rho, temp, diag_temp, nphot, diag_nphot):

    print("Now generating slice plots across x-, y- and z- axis\n")

    # Loop over three dimensions for XY // XZ // YZ slice, setting relevant
    # plot/axis variables as required

    for d in range(3):

        if (d == 0):
            dimpoint0 = 2 ; dimpoint1 = 0 ; dimpoint2 = 1
            pos0 = zloc ; pos1 = xloc ; pos2 = yloc ; xwid = xwid
            plt_x = "X (AU)" ; plt_y = "Y (AU)"
            temp_plt_name = "temp_xy"+plt_ext+"."+plt_form
            rho_plt_name = "rho_xy."+plt_form
            nphot_plt_name = "nphot_xy"+plt_ext+"."+plt_form
            print("Analysing x-y slice")
        elif (d == 1):
            dimpoint0 = 1 ; dimpoint1 = 0 ; dimpoint2 = 2
            pos0 = yloc ; pos1 = xloc ; pos2 = zloc ; xwid = xwid
            plt_x = "X (AU)" ; plt_y = "Z (AU)"
            temp_plt_name = "temp_xz"+plt_ext+"."+plt_form
            rho_plt_name = "rho_xz."+plt_form
            nphot_plt_name = "nphot_xz"+plt_ext+"."+plt_form
            print("Analysing x-z slice")
        elif (d == 2):
            dimpoint0 = 0 ; dimpoint1 = 1 ; dimpoint2 = 2
            pos0 = xloc ; pos1 = yloc ; pos2 = zloc ; xwid = xwid
            plt_x = "Y (AU)" ; plt_y = "Z (AU)"
            temp_plt_name = "temp_yz"+plt_ext+"."+plt_form
            rho_plt_name = "rho_yz."+plt_form
            nphot_plt_name = "nphot_yz"+plt_ext+"."+plt_form
            print("Analysing y-z slice\n")

    # Now constrain grid bins that fall within slice domain, and
    # compute base grid from minimum width of that subset

        sub_pos1 = [] ; sub_pos2  = [] ; sub_xwid = []
        sub_temp = [] ; sub_rho = [] ; sub_nphot = []
        for i in range(len(pos0)):
            if (pos1[i] >= -slice_hw) and (pos1[i] <= slice_hw) and \
               (pos2[i] >= -slice_hw) and (pos2[i] <= slice_hw) and \
               ( (abs(pos0[i]) - xwid[i]) <= 0.0) and \
               (xwid[i] != 0.0) :

                sub_pos1.append(pos1[i])
                sub_pos2.append(pos2[i])
                sub_xwid.append(xwid[i])

                if diag_temp:
                    temp_pix = 0.0
                    for s in range(nrspecies):
                        temp_pix += temp[s][i]
                    sub_temp.append(temp_pix)
                if diag_rho:
                    rho_pix = 0.0
                    for s in range(nrspecies):
                        rho_pix += rho[s][i]
                    sub_rho.append(rho_pix)
                if diag_nphot:
                    sub_nphot.append(nphot[i])
            else:
                continue
        min_xwid = min(sub_xwid)

    # Now define the number of pixels in each axis, and then initialise
    # arrays to be filled with pixel positions, and temperature/density
    #
    # If image has more than 1024 pixels in x- axis, artificially coarsen

        npix1 = int( ( 2.0 * slice_hw) / min_xwid ) ; npix2 = npix1

        while (npix1 > 2048):
            npix1 /= 1.01
            npix2 /= 1.01
            min_xwid *= 1.01
        npix1 = int(npix1) ; npix2 = int(npix2)

        pix_loc = [[ [0.0 for i in range(npix2)] for i in range(npix1)] \
           for i in range(2)]
        count = [ [0 for i in range(npix2)] for i in range(npix1)]

        if diag_temp:
            temp_slice = [ [0.0 for i in range(npix2)] for i in range(npix1)]
            temp_slice = np.array(temp_slice)
        if diag_rho:
            rho_slice = [ [0.0 for i in range(npix2)] for i in range(npix1)]
            rho_slice = np.array(rho_slice)
        if diag_nphot:
            nphot_slice = [ [0.0 for i in range(npix2)] for i in range(npix1)]
            nphot_slice = np.array(nphot_slice)

        pix_loc = np.array(pix_loc) ; count = np.array(count)

    # Now loop over npix for each axis, and set location (bin centre)

        for ii in range(npix1):
            for jj in range(npix2):
                pix_loc[0][ii][jj] = -slice_hw + (0.5 * min_xwid) + \
                  ( min_xwid * ii )
                pix_loc[1][ii][jj] = -slice_hw + (0.5 * min_xwid) + \
                  ( min_xwid * jj )

    # Now loop over subset of grid points, computing the pix_loc indices.
    # Using sub_xwid, add density/temperature to running total for both
    # grid centre, and respective bin volume

        for i in range(len(sub_pos1)):

            x_ind = int((sub_pos1[i] + slice_hw ) / min_xwid)
            y_ind = int((sub_pos2[i] + slice_hw ) / min_xwid)
            bin_inf = int(math.ceil(sub_xwid[i] / min_xwid))

    # While adding to running totals, if grid volume exceeds image domain,
    # simply discard.

            for ii in range(-bin_inf,bin_inf+1):

                if (x_ind + ii < 0) or (x_ind + ii > npix1 - 1):
                    continue

                else:
                    for jj in range(-bin_inf,bin_inf+1):

                        if (y_ind + jj < 0) or (y_ind + jj > npix2 - 1):
                            continue

                        else:

                            if diag_temp:
                                temp_slice[x_ind + ii][y_ind + jj] += \
                                  sub_temp[i]
                            if diag_rho:
                                rho_slice[x_ind + ii][y_ind + jj] += sub_rho[i]
                            if diag_nphot:
                                nphot_slice[x_ind + ii][y_ind + jj] += \
                                  sub_nphot[i]

                            count[x_ind + ii][y_ind + jj] += 1

    # Now generate plots, normalising to number of grid bins in each pixel.
    # Also define image interpolation, and rotate to preserve correct axes

        im_interp = "gaussian"
        slice_convolve = int(float(npix1) / 75.0)

    # Temperature surface plots for slice taken

        if diag_temp:

            temp_scale = "linear"

            temp_vmin = np.where(temp_slice != 0)
            if (temp_scale == "linear"):
                temp_vmin = temp_slice[temp_vmin].min()
            elif (temp_scale == "log"):
                temp_vmin = np.log10(temp_slice[temp_vmin].min())

            for ii in range(npix1):
                for jj in range(npix2):
                    if (count[ii][jj] == 0) or (temp_slice[ii][jj] == 0):
                        temp_slice[ii][jj] = temp_vmin
                    else:
                        if (temp_scale == "linear"):
                            temp_slice[ii][jj] /= count[ii][jj]
                        elif (temp_scale == "log"):
                            temp_slice[ii][jj] = np.log10(temp_slice[ii][jj])

            plt.figure(1)
            ax1 = plt.subplot(111)

            temp_slice = np.rot90(temp_slice,1)
            temp_slice = ndimage.gaussian_filter(temp_slice, \
              sigma=(slice_convolve, slice_convolve), order=0)
            temp_vmax = temp_slice.max() ; temp_vmin = temp_slice.min()+1

            plt.imshow(temp_slice, interpolation = im_interp , \
              vmin=temp_vmin, vmax=temp_vmax, origin="lower",
              extent=[pix_loc[0].min(), pix_loc[0].max(), \
              pix_loc[1].min(), pix_loc[1].max()], \
              aspect = "auto", cmap=cm.jet)
            cbar = plt.colorbar()
            cbar.ax.set_ylabel("Temperature (K)", labelpad = 2.5, \
              fontsize = cs.fontsize)
            plt.xlabel(plt_x, labelpad=0.5, fontsize = cs.fontsize)
            plt.ylabel(plt_y, labelpad=0.5, fontsize = cs.fontsize)
            plt.xlim(pix_loc[0].min(), pix_loc[0].max())
            plt.ylim(pix_loc[1].min(), pix_loc[1].max())
            plt.tight_layout()
            plt.savefig(plt_dir+"/"+temp_plt_name, format = plt_form)
            plt.clf()

    # Density surface plot for slice taken

        if diag_rho:

            rho_vmin = np.where(rho_slice != 0.)
            rho_vmin = np.log10(rho_slice[rho_vmin].min())

            for ii in range(npix1):
                for jj in range(npix2):
                    if (count[ii][jj] == 0):
                        rho_slice[ii][jj] = rho_vmin
                    else:
                        rho_slice[ii][jj] /= count[ii][jj]
                        rho_slice[ii][jj] = np.log10(rho_slice[ii][jj])

            plt.figure(1)
            ax1 = plt.subplot(111)

            rho_slice = np.rot90(rho_slice,1)
            rho_slice = ndimage.gaussian_filter(rho_slice, \
              sigma=(slice_convolve, slice_convolve), order=0)
            rho_vmax = rho_slice.max() ; rho_vmin = rho_slice.min()

            plt.imshow(rho_slice, interpolation = im_interp, \
              vmin=rho_vmin, vmax=rho_vmax, origin="lower",
              extent=[pix_loc[0].min(), pix_loc[0].max(), \
              pix_loc[1].min(), pix_loc[1].max()], \
              aspect = "auto", cmap=cm.jet)
            cbar = plt.colorbar()
            cbar.ax.set_ylabel("Density "+r"(g cm$^{-3}$)")
            plt.xlabel(plt_x, labelpad=0.5) ; plt.ylabel(plt_y, labelpad=0.5)
            plt.xlim(pix_loc[0].min(), pix_loc[0].max())
            plt.ylim(pix_loc[1].min(), pix_loc[1].max())
            plt.tight_layout()
            plt.savefig(plt_dir+"/"+rho_plt_name, format = plt_form)
            plt.clf()

    # Photon count plot for slice taken

        if diag_nphot:

            nphot_vmin = np.where(nphot_slice != 0)
            nphot_vmin = nphot_slice[nphot_vmin].min()

            for ii in range(npix1):
                for jj in range(npix2):
                    if (count[ii][jj] == 0):
                        nphot_slice[ii][jj] = nphot_vmin
                    else:
                        nphot_slice[ii][jj] /= count[ii][jj]

            plt.figure(1)
            ax1 = plt.subplot(111)

            nphot_slice = np.rot90(nphot_slice,1)
            nphot_slice = ndimage.gaussian_filter(nphot_slice, \
              sigma=(slice_convolve, slice_convolve), order=0)
            nphot_vmax = nphot_slice.max() ; nphot_vmin = nphot_slice.min()

            plt.imshow(nphot_slice, interpolation = im_interp, \
              vmin=nphot_vmin, vmax=nphot_vmax, origin="lower",
              extent=[pix_loc[0].min(), pix_loc[0].max(), \
              pix_loc[1].min(), pix_loc[1].max()], \
              aspect = "auto", cmap=cm.jet)
            cbar = plt.colorbar()
            cbar.ax.set_ylabel("Number of photons")
            plt.xlabel(plt_x, labelpad=0.5) ; plt.ylabel(plt_y, labelpad=0.5)
            plt.xlim(pix_loc[0].min(), pix_loc[0].max())
            plt.ylim(pix_loc[1].min(), pix_loc[1].max())
            plt.tight_layout()
            plt.savefig(plt_dir+"/"+nphot_plt_name, format = plt_form)
            plt.clf()

    return
