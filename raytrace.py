'''

raytrace.py

Module to compute column density and optical depth of rays from protostar
to observer location as a function of position angle and inclination

Last Modified: 21/03/2018

'''

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - MODULE IMPORTS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

import numpy as np
import math
import os
import shutil
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from matplotlib import cm

import constants as cs
from params_run import *
from functions import Trapezoidal
from functions import Planck
from functions import Column_Density

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - MAIN PROGRAM - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

def column(dat_dir, plt_dir, plt_ext, plt_form):

    renew = True    # Is data regenerated (True), or just plotted (False)?

    n_phi = 180      # Number of samplings in azimuthal angle
    n_inclin = 45      # Number of samplings in inclination

    phi_inc = (2.0 * math.pi) / float(n_phi + 1)
    inclin_inc = math.pi / float(n_inclin)

    # If plot data isn't to be regenerated, read, plot and return to make.py

    fmap = dat_dir+"/raymap.dat"
    if os.path.isfile(fmap) and not renew:

        fmap = open(fmap, "r")
        n_phi = int(fmap.readline()) ; n_inclin = int(fmap.readline())

        phi_plt = [ [ 0.0 for i in range(n_inclin) ] \
          for p in range(n_phi) ]
        inclin_plt = [ [ 0.0 for i in range(n_inclin) ] \
          for p in range(n_phi) ]
        rho_column_plt = [ [ 0.0 for i in range(n_inclin) ] \
          for p in range(n_phi) ]

        for p in range(n_phi):
            for i in range(n_inclin):
                line = fmap.readline()
                phi_plt[p][i] = float(line.split(" ")[0])
                inclin_plt[p][i] = float(line.split(" ")[1])
                rho_column_plt[p][i] = float(line.split(" ")[2])

        column_analyse(plt_dir, plt_ext, plt_form, phi_plt, inclin_plt, \
          rho_column_plt)

        return

    # Otherwise, compute column density as function of observer location

    elif renew:

        fmap = open(fmap, "w")
        fmap.write("{0}\n".format(n_phi) )
        fmap.write("{0}\n".format(n_inclin) )

    # Read in location/density data, computing phi/i value ranges for each bin.
    # Also find total ray path length

        xloc, yloc, zloc, rloc, xwid = \
          np.loadtxt(dat_dir+"/loc_grid.dat",unpack=True)
        rloc_xy = np.hypot(xloc,yloc)

        xloc = np.array(xloc) ; yloc = np.array(yloc) ; zloc = np.array(zloc)
        rloc = np.array(rloc) ; rloc_xy = np.array(rloc_xy)
        xwid = np.array(xwid)

        ngrid = len(xloc)

#        phi = [ [] for i in range(ngrid) ]
#        inclin = [ [] for i in range(ngrid) ]

#        for i in range(ngrid):

#            d_phi = 2.0 * math.atan2((np.sqrt(8.0) * xwid[i]) , rloc_xy[i])

#            d_inclin = 2.0 * math.atan2((np.sqrt(8.0) * xwid[i]) , rloc[i])

#            if (yloc[i] >= 0.0):
#                phi[i].append( math.atan2(yloc[i] , xloc[i]) )
#                phi[i].append( d_phi )
#            elif (yloc[i] < 0.0):
#                phi[i].append( (2. * math.pi) - \
#                  math.atan2(abs(yloc[i]) , xloc[i]) )
#                phi[i].append( d_phi )
#            elif ((yloc[i] == 0.0) and (xloc[i] == 0.0)):
#                phi[i].append(0)
#                phi[i].append(0)

#            if (zloc[i] != 0.):
#                inclin[i].append( math.atan2(rloc_xy[i] , zloc[i]) )
#                inclin[i].append( d_inclin )
#            elif (zloc[i] == 0.0):
#                inclin[i].append( 0 )
#                inclin[i].append( 0 )

        f = open(dat_dir+"/dust_density.inp","r")
        trash = f.readline() ; trash = int(f.readline())
        nrspecies = int(f.readline())
        rho = [ [] for i in range(nrspecies)]
        for i in range(nrspecies):
            for j in range(ngrid):
                rho[i].append(float(f.readline()))
        f.close()
        rho = np.array(rho[0])

#        ray_lim = max(rloc) ; not_zero = np.where(xwid != 0)
#        ray_n = int(max(rloc) / min(xwid[not_zero]))

    # Loop over phi/i values, setting to non-zero to prevent missing bins
    # due to ray passing through two bin edges

        phi_plt = [ [ 0.0 for i in range(n_inclin) ] \
          for p in range(n_phi) ]
        inclin_plt = [ [ 0.0 for i in range(n_inclin) ] \
          for p in range(n_phi) ]
        rho_column_plt = [ [ 0.0 for i in range(n_inclin) ] \
          for p in range(n_phi) ]

        phi_ray = 0.0
        for p in range(n_phi):

            inclin_ray = 0.0
            for i in range(n_inclin):

                phi_plt[p][i] = phi_ray
                inclin_plt[p][i] = inclin_ray

    # Refine density/location bin data based on phi/i value being used

#                bin_dens = []
#                bin_xloc = [] ; bin_yloc = [] ; bin_zloc = [] ; bin_xwid = []

#                for b in range(ngrid):

#                    phi_diff = abs(phi_ray - phi[b][0] )
#                    inclin_diff = abs(inclin_ray - inclin[b][0])

#                    if (phi_diff <= phi[b][1] ) and \
#                      (inclin_diff <= inclin[b][1]):

#                        bin_dens.append(rho[b])
#                        bin_xloc.append(xloc[b])
#                        bin_yloc.append(yloc[b])
#                        bin_zloc.append(zloc[b])
#                        bin_xwid.append(xwid[b])

#                    else:
#                        continue

    # Loop over ray path length, adding location density/radius to array
    # (N.B., ray path from centre of envelope to envelope edge)

#                r_ray = 0.0 ; r_path = [] ; rho_path = []
#                for l in range(ray_n):

#                    r_ray += (ray_lim / ray_n)

#                    z_ray = r_ray * math.cos(inclin_ray)
#                    rxy_ray = abs(r_ray) * math.sin(inclin_ray)

#                    x_ray = abs(rxy_ray) * math.cos(phi_ray)
#                    y_ray = abs(rxy_ray) * math.sin(phi_ray)

#                    dens_ray = 0.0
#                    for b in range(len(bin_dens)):

#                        if (x_ray >= (bin_xloc[b] - bin_xwid[b]) ) and \
#                          (x_ray <= (bin_xloc[b] + bin_xwid[b]) ) and \
#                          (y_ray >= (bin_yloc[b] - bin_xwid[b]) ) and \
#                          (y_ray <= (bin_yloc[b] + bin_xwid[b]) ) and \
#                          (z_ray >= (bin_zloc[b] - bin_xwid[b]) ) and \
#                          (z_ray <= (bin_zloc[b] + bin_xwid[b]) ):

#                            dens_ray += bin_dens[b]

#                            break

#                        else:
#                            continue

#                    r_path.append(r_ray) ; rho_path.append(dens_ray)

    # Sanity check to ensure density along first ray is sufficiently resolved

#            if p == 0 and i == 0:
#                plt.scatter(r_path, np.log10(rho_path))
#                plt.show()

    # Integrate densiy profile for ray to compute column density through ray

#                rho_int = interp1d(r_path, rho_path)
#                def density_r(r):
#                    return rho_int(r)

#                rho_column_plt[p][i] = Trapezoidal(density_r, \
#                  min(r_path), max(r_path), ray_n)

                rho_column_plt[p][i] = Column_Density(inclin_ray, phi_ray, \
                  xloc, yloc, zloc, xwid, rho)
#
                fmap.write("{0} {1} {2}\n".format(phi_ray, inclin_ray, \
                  rho_column_plt[p][i] ) )

                print phi_ray, inclin_ray, rho_column_plt[p][i]

    # Increment in observer location (phi//inclin)

                inclin_ray += inclin_inc
            phi_ray += phi_inc
        fmap.close()

    # Now plot 1d and 2d column density distributions

        column_analyse(plt_dir, plt_ext, plt_form, phi_plt, inclin_plt, \
          rho_column_plt)

        return

### ------------------------------------------------------------------------ ###

def column_analyse(plt_dir, plt_ext, plt_form, phi, inclin, rho_column):

    # Translate to RADMC-3D angular system, and convert to degrees

    for p in range(len(phi)):
        for i in range(len(phi[0])):

            phi[p][i] = (2.0 * math.pi) - phi[p][i]
            if (phi[p][i] <= (math.pi / 2.0) ):
                phi[p][i] += (3.0/2.0) * math.pi
            else:
                phi[p][i] -= (math.pi / 2.0)

            phi[p][i] = math.degrees(phi[p][i])
            inclin[p][i] = math.degrees(inclin[p][i])

    # Now rotate all arrays to plot phi (x-axis) vs. inclination (y-axis)

    rho_column = np.log10(np.rot90(rho_column,-1))
    inclin = np.rot90(inclin,-1)
    phi = np.rot90(phi,-1)

    # Surface plot for phi/i depedent column density

    plt.figure(1) ; ax1 = plt.subplot(111)

    plt.imshow(rho_column, interpolation = "none", \
      vmin=rho_column.min(), vmax=rho_column.max(), origin="lower",
      extent=[np.amin(phi), np.amax(phi), np.amin(inclin), np.amax(inclin)], \
      aspect = "auto", cmap=cm.jet)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel("Column Density "+r"(g cm$^{-2}$)")
    plt.xlabel("Azimuthal angle (deg.)", labelpad=0.5)
    plt.ylabel("Inclination angle (deg.)", labelpad=0.5)
    plt.xlim(np.amin(phi), np.amax(phi))
    plt.ylim(np.amin(inclin), np.amax(inclin))
    plt.tight_layout()
    plt.savefig(plt_dir+"/rho_map."+plt_form, format = plt_form) ; plt.clf()

    # Azimuthally averaged column density profile

    mean =  [] ; std_dev = [] ; inclin_plt = []
    for i in range(len(phi)):
        rho_ref = []
        for p in range(len(phi[0])):
            rho_ref.append(rho_column[i][p])

        rho_ref = np.array(rho_ref)
        inclin_plt.append(inclin[i][0])
        mean.append(np.log10(np.mean(rho_ref)))
        std_dev.append(0.434 * np.std(rho_ref) / np.mean(rho_ref))

    plt.figure(1) ; ax1 = plt.subplot(111)

    ax1.errorbar(inclin_plt, mean, yerr=std_dev)
    plt.ylabel(r"log Mean Column Density (g cm$^{-2}$)", labelpad=0.5)
    plt.xlabel("Inclination angle (deg.)", labelpad=0.5)
    plt.tight_layout()
    plt.savefig(plt_dir+"/meanrho."+plt_form, format = plt_form) ; plt.clf()

    # phi = 0 column density profile

    rho_plt = [] ; inclin_plt = []
    for i in range(len(phi)):
        inclin_plt.append(inclin[i][0])
        rho_plt.append(np.log10(rho_column[i][0]))

    plt.figure(1) ; ax1 = plt.subplot(111)

    plt.plot(inclin_plt, rho_plt)

    plt.ylabel("Column Density "+r"(g cm$^{-2}$)", labelpad=0.5)
    plt.xlabel("Inclination angle (deg.)", labelpad=0.5)
    ax1.set_xlim(0, 180) ; ax1.set_ylim(min(rho_plt), 1.1*max(rho_plt))
    plt.tight_layout()
    plt.savefig(plt_dir+"/inclin_rho_0phi."+plt_form, format = plt_form)
    plt.clf()

### ------------------------------------------------------------------------ ###

def tau1d(dat_dir, plt_dir, plt_form, tdat_ext):

    global wavelength_val

    renew = False    # Is data regenerated (True), or just plotted (False)?

    phi_ray = 1.e-10
    inclin_ray = [1.e-10, (math.pi / 3.0) + 1.e-10, (math.pi / 2.0) + 1.e-10]

    ray_n_coarse = 2.5  # Resolution factor decrease to sample rho along full ray
    cols_i = ["r","g","b"]   # Plot colors for inclinations

    f_tauI = dat_dir+"/tau_I_"+str(int(wavelength_val))+"micron.dat"

    # If data isn't renewed and file is found, read data then go to plotting

    if not renew and os.path.isfile(f_tauI):

        f_tauI = open(f_tauI, "r")
        n_inclin = int(f_tauI.readline())
        n_path = []
        for i in range(n_inclin):
            n_path.append( int( f_tauI.readline() ) - 1)
        r_path = [ [] for i in range(n_inclin) ]
        tau_path = [ [] for i in range(n_inclin) ]
        intensity_path = [ [] for i in range(n_inclin) ]
        sum_intensity_path = [ [] for i in range(n_inclin) ]
        for i in range(n_inclin):
            for p in range(n_path[i]):
                line = f_tauI.readline() ; columns = line.split()
                r_path[i].append(float(columns[0]))
                tau_path[i].append(float(columns[1]))
                intensity_path[i].append(float(columns[2]))
                sum_intensity_path[i].append(float(columns[3]))
        f_tauI.close()

    # Otherwise, generate data

    elif renew or not os.path.isfile(f_tauI):

    # Read in AMR grid, and allocate inclination/azimuthal angle values

        xloc, yloc, zloc, rloc, xwid = \
          np.loadtxt(dat_dir+"/loc_grid.dat",unpack=True)
        rloc_xy = np.hypot(xloc,yloc)

        ngrid = len(xloc)

        phi = [ [] for i in range(ngrid) ]
        inclin = [ [] for i in range(ngrid) ]
        for i in range(ngrid):

            d_phi = 2.0 * math.atan2((np.sqrt(8.0) * xwid[i]) , rloc_xy[i])

            d_inclin = 2.0 * math.atan2((np.sqrt(8.0) * xwid[i]) , rloc[i])

            if (yloc[i] >= 0.0):
                phi[i].append( math.atan2(yloc[i] , xloc[i]) )
                phi[i].append( d_phi )
            elif (yloc[i] < 0.0):
                phi[i].append( (2. * math.pi) - \
                  math.atan2(abs(yloc[i]) , xloc[i]) )
                phi[i].append( d_phi )
            elif ((yloc[i] == 0.0) and (xloc[i] == 0.0)):
                phi[i].append(0)
                phi[i].append(0)

            if (zloc[i] != 0.):
                inclin[i].append( math.atan2(rloc_xy[i] , zloc[i]) )
                inclin[i].append( d_inclin )
            elif (zloc[i] == 0.0):
                inclin[i].append( 0 )
                inclin[i].append( 0 )

    # Now read in both density and temperature data to compute tau and intensity

        f_rho = open(dat_dir+"/dust_density.inp","r")
        trash = f_rho.readline() ; trash = int(f_rho.readline())
        nrspecies = int(f_rho.readline())
        rho = [ [] for i in range(nrspecies)]
        for i in range(nrspecies):
            for j in range(ngrid):
                rho[i].append(float(f_rho.readline()))
        f_rho.close()

        f_temp = open(dat_dir+"/dust_temperature"+tdat_ext,"r")
        trash = f_temp.readline() ; trash = int(f_temp.readline())
        trash = int(f_temp.readline())
        temp = [ [] for i in range(nrspecies)]
        for i in range(nrspecies):
            for j in range(ngrid):
                temp[i].append(float(f_temp.readline()))
        f_temp.close()

        rho = np.array(rho[0]) ; temp = np.array(temp[0])

    # Read, and interpolate for opacity for selected wavelength

        llambda_opac = [] ; opac = []

        f_opac = open(dat_dir+"/dustkappa_"+opacity_inp+".inp","r")
        iformat = f_opac.readline()
        wav_bins_opac = int(f_opac.readline())
        for i in range(wav_bins_opac):
            line = f_opac.readline() ; column = line.split()
            llambda_opac.append(float(column[0]))
            opac.append( float(column[1]) ) # Use only abs. opacity
        f_opac.close()

        opac_int = interp1d(llambda_opac,opac,kind="linear")

        opacity_val = opac_int(wavelength_val)

    # Compute ray-specific variables

        ray_lim = max(rloc) ; not_zero = np.where(xwid != 0)
        ray_n = int(max(rloc) / ( ray_n_coarse * min(xwid[not_zero]) ) )

    # Arrays for refined bin data based on inclination value being used

        bin_dens = [ [] for i in range(len(inclin_ray))]
        bin_temp = [ [] for i in range(len(inclin_ray))]
        bin_xloc = [ [] for i in range(len(inclin_ray))]
        bin_yloc = [ [] for i in range(len(inclin_ray))]
        bin_zloc = [ [] for i in range(len(inclin_ray))]
        bin_xwid = [ [] for i in range(len(inclin_ray))]

    # Arrays for location/density/optical depth/intensity along ray

        r_path = [ [] for i in range(len(inclin_ray)) ]
        rho_path = [ [] for i in range(len(inclin_ray)) ]
        tau_path = [ [] for i in range(len(inclin_ray)) ]
        intensity_path = [ [] for i in range(len(inclin_ray)) ]
        sum_intensity_path = [ [] for i in range(len(inclin_ray)) ]

    # Loop over inclinations, loop for ray path length and evaluate properties
    # (N.B., ray path from envelope edge to centre)

        for i in range(len(inclin_ray)):

            for b in range(ngrid):

                phi_diff = abs(phi_ray - phi[b][0] )
                inclin_diff = abs(inclin_ray[i] - inclin[b][0])

                if (phi_diff <= phi[b][1] ) and (inclin_diff <= inclin[b][1]):

                    bin_dens[i].append(rho[b])
                    bin_temp[i].append(temp[b])
                    bin_xloc[i].append(xloc[b])
                    bin_yloc[i].append(yloc[b])
                    bin_zloc[i].append(zloc[b])
                    bin_xwid[i].append(xwid[b])

                else:
                    continue

            r_ray = ray_lim ; l_ray = 0.0 ; path_count = 0
            tau_sum = 0.0
            print('''\nDistance along ray (AU) ||| Optical depth |||\
            Sum (Intesity)''')

            min_rinc = 0.01 * min(xwid[not_zero])
            max_rinc = 1000 * min_rinc
            ray_n = int( (2.0 * ray_lim) / (min_rinc + max_rinc) )

            for l in range(ray_n):

                if (tau_sum > 2.5):
                    break

                rinc = max_rinc - ( (float(l) / float(ray_n)) * \
                  (max_rinc - min_rinc) )
                r_ray -= rinc
                l_ray += rinc
#                r_ray -= (ray_lim / ray_n)  # Radial loc. of ray w.r.t. [0,0,0]
#                l_ray += (ray_lim / ray_n)  # Radial distance along ray

                z_ray = r_ray * math.cos(inclin_ray[i])
                rxy_ray = abs(r_ray) * math.sin(inclin_ray[i])
                x_ray = abs(rxy_ray) * math.cos(phi_ray)
                y_ray = abs(rxy_ray) * math.sin(phi_ray)

    # Evaluate optical depth, and emissivity along ray path

                for b in range(len(bin_dens[i])):

                    if (x_ray >= (bin_xloc[i][b] - bin_xwid[i][b]) ) and \
                      (x_ray <= (bin_xloc[i][b] + bin_xwid[i][b]) ) and \
                      (y_ray >= (bin_yloc[i][b] - bin_xwid[i][b]) ) and \
                      (y_ray <= (bin_yloc[i][b] + bin_xwid[i][b]) ) and \
                      (z_ray >= (bin_zloc[i][b] - bin_xwid[i][b]) ) and \
                      (z_ray <= (bin_zloc[i][b] + bin_xwid[i][b]) ):

                        r_path[i].append( l_ray )

                        rho_path[i].append( bin_dens[i][b] )

                        tau_incr = rho_path[i][path_count] * \
                          (ray_lim / ray_n) * opacity_val

                        path_count += 1
                        if (path_count < 2):
                            tau_path[i].append(0)
                            intensity_path[i].append(0)
                            sum_intensity_path[i].append(0)
                            break

                        rho_int = interp1d(r_path[i], rho_path[i])
                        def density_r(r):
                            return rho_int(r)

                        int_n = path_count
                        rho_column_sum = Trapezoidal(density_r, \
                          min(r_path[i]), max(r_path[i]), int_n)
                        tau_sum = rho_column_sum * opacity_val

                        tau_path[i].append( tau_sum )

                        wavlength_B = wavelength_val * 1.00e-6
                        B_nu = Planck(wavelength_val, bin_temp[i][b])
                        intensity = B_nu * math.exp( - tau_sum ) * tau_incr

                        intensity_path[i].append( intensity )

                        sum = 0.0
                        for s in range(len(intensity_path[i])):
                            sum += intensity_path[i][s]

                        if (path_count % 100 == 0):
                            print l_ray*cs.au_per_cm, "\t", tau_sum, "\t", sum

                        sum_intensity_path[i].append(sum)

                        break

                    else:
                        continue

    # If data has been renewed, write to file

        f_tauI = open(f_tauI, "w")
        f_tauI.write("{0}\n".format(len(inclin_ray)) )

        for i in range(len(inclin_ray)):
            f_tauI.write("{0}\n".format( len(r_path[i]) ) )

        for i in range(len(inclin_ray)):
            for p in range(len(r_path[i])):
                if p == 0:
                    continue
                f_tauI.write('''{0} {1} {2} {3}\n'''.format(r_path[i][p], \
                  tau_path[i][p], intensity_path[i][p], \
                  sum_intensity_path[i][p]) )
        f_tauI.close()

    r_path = np.array(r_path)
    tau_path = np.array(tau_path)
    intensity_path = np.array(intensity_path)
    sum_intensity_path = np.array(sum_intensity_path)

    max_r = max(r_path[0])
    for i in range(len(inclin_ray)):
        for p in range(len(r_path[i])):
#            r_path[i][p] = max_r - r_path[i][p]
            r_path[i][p] *= cs.au_per_cm
            r_path[i][p] = 10000 - r_path[i][p]

    # Plot intensity as function of optical depth and inclination

    plt.figure(1) ; ax1 = plt.subplot(111)

    for i in range(len(inclin_ray)):
        plt.plot(tau_path[i], intensity_path[i], color = cols_i[i], \
          label = "i = "+str(int(math.degrees(inclin_ray[i]) ) )+" deg." )
    plt.legend(loc = "upper right", fontsize=10)

    plt.xlabel("Optical Depth at "+str(int(wavelength_val))+" microns", \
      labelpad=0.5)
    ax1.set_xlim(0, 3)
    ax1.axvline(1, ymin = min(intensity_path[0]), \
      ymax = max(intensity_path[2])*100, linestyle = ":" )
    ax1.set_ylim(min(intensity_path[0]), max(intensity_path[2])*100 )
    plt.tight_layout()
    plt.savefig(plt_dir+"/tau1d_I_"+str(int(wavelength_val))+"micron."+ \
      plt_form, format = plt_form)
    plt.clf()

    # Plot radius as function of optical depth and inclination

    plt.figure(1) ; ax1 = plt.subplot(111)

    for i in range(len(inclin_ray)):
        plt.plot(r_path[i][::20], tau_path[i][::20], color = cols_i[i], \
          label = "i = "+str(int(math.degrees(inclin_ray[i]) ) )+" deg." )
    plt.legend(loc = "upper right", fontsize=cs.leg_fontsize)

    plt.ylabel(r"$\tau_\lambda$ ($\lambda$ = "+str(int(wavelength_val))+r" $\mu$m)", \
      fontsize = cs.fontsize, labelpad=0.5)
    ax1.set_ylim(0., 1.5)
#    ax1.axhline(1.0, xmin = min(r_path[2]), \
#      xmax = max(r_path[2]), linestyle = ":", color = "k" )
#    ax1.axhline(1.0, xmin = 0, \
#      xmax = 100, linestyle = ":", color = "k" )
    plt.xlabel(r"R (AU)", labelpad=0.5, fontsize = cs.fontsize)
#    ax1.set_xlim(min(r_path[0]), max(r_path[2]))
    ax1.set_xlim(0, 100 )
    plt.tight_layout()
    plt.savefig(plt_dir+"/tau1d_R_"+str(int(wavelength_val))+"micron."+ \
      plt_form, format = plt_form)
    plt.clf()

    # Plot summation of intensity as function of optical depth and inclination

    plt.figure(1) ; ax1 = plt.subplot(111)

    for i in range(len(inclin_ray)):
        plt.plot(tau_path[i], sum_intensity_path[i], color = cols_i[i], \
          label = "i = "+str(int(math.degrees(inclin_ray[i]) ) )+" deg." )

    plt.legend(loc = "lower right", fontsize=10)
    plt.xlabel("Optical Depth at "+str(wavelength_val)+" microns", labelpad=0.5)
    ax1.set_xlim(0.0, 3)
    ax1.axvline(1, ymin = min(sum_intensity_path[0]), \
      ymax = max(sum_intensity_path[2])*100, linestyle = ":" )
    plt.ylabel(r"$\Sigma I_{\nu}$", labelpad=0.5) ; ax1.set_yscale("log")
    ax1.set_ylim(min(sum_intensity_path[0]), max(sum_intensity_path[2])*100 )
    plt.tight_layout()
    plt.savefig(plt_dir+"/tau1d_sumI_"+str(int(wavelength_val))+"micron."+ \
      plt_form, format = plt_form)
    plt.clf()

    return

### ------------------------------------------------------------------------ ###

def tau3d(dat_dir, plt_dir, plt_form, tdat_ext):

    npix = 1000             # Number of pixels in each image axis for RRT
    inclins = [0, 60, 90]   # Inclinations (deg.) to be analysed

    xloc, yloc, zloc, rloc, xwid = \
       np.loadtxt(dat_dir+"/loc_grid.dat",unpack=True)
    rloc_xy = np.hypot(xloc,yloc)

    max_r = max(rloc) * cs.au_per_cm

    z_ave = []

    shutil.copy2(dat_dir+"/dust_temperature"+tdat_ext, \
      dat_dir+"/dust_temperature.dat" )
    for i in range(len(inclins)):

        os.system('''{0}radmc3d tausurf {1} lambda {2} incl {3} \\
        zoomau {4} {5} {4} {5} npix {6} secondorder'''.format(exec_loc, \
        tausurf_val, wavelength_val, inclins[i], - max_r, max_r, npix) )

        x_im = [] ; y_im = [] ; z_im = []
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
                    x_im.append( xloc )
                    y_im.append( yloc )
                    z_im.append( z_tau * cs.au_per_cm )

                z_ave.append(np.mean(z_im))
        f.close()

    # Scatter for phi/i, with color scaled by radius

        plt.figure(1) ; ax1 = plt.subplot(111)
        vmax = max(z_im) ; vmin = min(z_im)
        plt.scatter(x_im, y_im, c = z_im, edgecolors = "none", \
          vmin=vmin, vmax=vmax, cmap=cm.jet)

        cbar = plt.colorbar()
        cbar.ax.set_ylabel(r"$Z_{proj}$ ($\tau$ = "+str(int(tausurf_val))+ \
          ") at "+ str(wavelength_val)+" micron (AU)")

        plt.xlabel(r"$X_{proj}$ (AU)", labelpad=0.5)
        plt.xlim(min(x_im), max(x_im) )
        plt.ylabel(r"$Y_{proj}$ (AU)", labelpad=0.5)
        plt.ylim(min(y_im), max(y_im) )
        plt.tight_layout()
        plt.savefig(plt_dir+"/r_tau"+str(int(tausurf_val))+"_"+ \
          str(int(wavelength_val))+"micron_"+ \
          str(int(inclins[i]))+"i."+plt_form, format = plt_form)
        plt.clf()

    print("Mean Z_proj values are:")
    for i in range(len(inclins)):
        print("i = {0} \t mean(Z_proj) = {1} AU".format(inclins[i], \
          z_ave[i]))


    return
