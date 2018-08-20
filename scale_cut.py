''''

scale_cut.py

Module to scale mass, and apply radial restrictions to models prior to MCRT and/
or prior to RRT in RADMC-3D.

Module also modifies the dust_density.inp file for artificially including a
cavity, dependent on location and profile as set by CAVITY_ parameter variables

Prior to any rescaling/cuts, a copy of the uncut/unscaled inputs are copied.

Last Modified: 18/07/2018

'''

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - MODULE IMPORTS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

import numpy as np
import os
import shutil
import math

import constants as cs
from params_run import *

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - MAIN PROGRAM - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

def pre_full_scale(dat_dir):

    # Read in density/bin data, and add to running mass total with bin mass.
    # Then compute mass scaling factor and re-write dust_density.inp

    xloc, yloc, zloc, rloc, xwid = \
       np.loadtxt(dat_dir+"/loc_grid.dat",unpack=True)
    rloc = rloc * cs.au_per_cm

    if not os.path.isfile(dat_dir+"/MCRT_UNSCALED_dust_density.inp"):
        shutil.copy2(dat_dir+"/dust_density.inp", \
          dat_dir+"/MCRT_UNSCALED_dust_density.inp")

    f = open(dat_dir+"/dust_density.inp","r")
    trash = f.readline()
    ngrid = int(f.readline())
    nrspecies = int(f.readline())

    rho = [[] for i in range(nrspecies)]
    amr_mass = 0.00
    for i in range(nrspecies):
        for j in range(ngrid):
            rho[i].append(float(f.readline()))
            amr_mass += rho[i][j] * (2.0*xwid[j])**3.0
    f.close()

    f_scale = (FULL_mass_scale * cs.dust_to_gas) / \
      (amr_mass * cs.msol_per_g)

    f = open(dat_dir+"/dust_density.inp","w")
    f.write("1\n")
    f.write("{0}\n".format(ngrid))
    f.write("{0}\n".format(nrspecies))
    for i in range(nrspecies):
        for j in range(ngrid):
            f.write("{0}\n".format(rho[i][j]*f_scale) )
    f.close()

    return

### ------------------------------------------------------------------------ ###

def pre_comp_scale(dat_dir):

    # With simulation/user defined component domains, run as per
    # pre_full_scale(), isolating disc and envelope components

    global DISC_h_cut, DISC_r_cut

    if f_DISC_params:

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
                  "restriction is: {0} AU\n".format(float(column[11]) ) )
                print("Mass of simulation disc (radial infall restriction) "+ \
                  "is: {0} M_sol\n".format( float(column[6]) ) )
                break
            else:
                continue
        frdisc.close()

        DISC_h_cut = 0.0
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
                    DISC_h_cut += h
        hfact = 3.0
        DISC_h_cut = hfact * (DISC_h_cut / float(count) )
        print('''Scale height restriction of disc, set by {0} * average h(r)
         value for disc < DISC_r_cut is: {1}\n'''.format(hfact, DISC_h_cut) )
        fpdisc.close()

    xloc, yloc, zloc, rloc, xwid = \
       np.loadtxt(dat_dir+"/loc_grid.dat",unpack=True)
    rloc *= cs.au_per_cm
    rloc_xy = np.sqrt( xloc**2.0 + yloc**2.0 ) * cs.au_per_cm
    hloc = ( zloc * cs.au_per_cm ) / rloc_xy

    if not os.path.isfile(dat_dir+"/MCRT_UNSCALED_dust_density.inp"):
        shutil.copy2(dat_dir+"/dust_density.inp", \
          dat_dir+"/MCRT_UNSCALED_dust_density.inp")

    f = open(dat_dir+"/dust_density.inp","r")
    trash = f.readline()
    ngrid = int(f.readline())
    nrspecies = int(f.readline())

    rho = [[] for i in range(nrspecies)]
    disc_mass = 0.00 ; envelope_mass = 0.00

    for i in range(nrspecies):
        for j in range(ngrid):
            rho[i].append(float(f.readline()))
            if (abs(hloc[j]) < DISC_h_cut) and (rloc_xy[j] < DISC_r_cut):
                disc_mass += rho[i][j] * (2.0*xwid[j])**3.0
            else:
                envelope_mass += rho[i][j] * (2.0*xwid[j])**3.0
    f.close()

    f_scale_disc = (DISC_mass_scale * cs.dust_to_gas) / \
      (disc_mass * cs.msol_per_g)
    f_scale_env = (ENV_mass_scale * cs.dust_to_gas) / \
      (envelope_mass * cs.msol_per_g)

    f = open(dat_dir+"/dust_density.inp","w")
    f.write("1\n")
    f.write("{0}\n".format(ngrid))
    f.write("{0}\n".format(nrspecies))
    for i in range(nrspecies):
        for j in range(ngrid):

            if (abs(hloc[j]) < DISC_h_cut) and (rloc_xy[j] < DISC_r_cut):
                f.write("{0}\n".format(rho[i][j] * f_scale_disc) )
            else:
                f.write("{0}\n".format(rho[i][j] * f_scale_env) )
    f.close()

    return

### ------------------------------------------------------------------------ ###

def pre_comp_mult(dat_dir):

    # With simulation/user defined component domains, run as per
    # pre_full_scale(), isolating disc and envelope components

    global DISC_h_cut, DISC_r_cut

    if f_DISC_params:

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
                  "restriction is: {0} AU\n".format(float(column[11]) ) )
                print("Mass of simulation disc (radial infall restriction) "+ \
                  "is: {0} M_sol\n".format( float(column[6]) ) )
                break
            else:
                continue
        frdisc.close()

        DISC_h_cut = 0.0
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
                    DISC_h_cut += h
        hfact = 3.0
        DISC_h_cut = hfact * (DISC_h_cut / float(count) )
        print('''Scale height restriction of disc, set by {0} * average h(r)
         value for disc < DISC_r_cut is: {1}\n'''.format(hfact, DISC_h_cut) )
        fpdisc.close()

    xloc, yloc, zloc, rloc, xwid = \
       np.loadtxt(dat_dir+"/loc_grid.dat",unpack=True)
    rloc *= cs.au_per_cm
    rloc_xy = np.sqrt( xloc**2.0 + yloc**2.0 ) * cs.au_per_cm
    hloc = ( zloc * cs.au_per_cm ) / rloc_xy

    if not os.path.isfile(dat_dir+"/MCRT_UNSCALED_dust_density.inp"):
        shutil.copy2(dat_dir+"/dust_density.inp", \
          dat_dir+"/MCRT_UNSCALED_dust_density.inp")

    f = open(dat_dir+"/dust_density.inp","r")
    trash = f.readline()
    ngrid = int(f.readline())
    nrspecies = int(f.readline())

    rho = [[] for i in range(nrspecies)]

    for i in range(nrspecies):
        for j in range(ngrid):
            rho[i].append(float(f.readline()))
    f.close()

    f = open(dat_dir+"/dust_density.inp","w")
    f.write("1\n")
    f.write("{0}\n".format(ngrid))
    f.write("{0}\n".format(nrspecies))
    for i in range(nrspecies):
        for j in range(ngrid):

            if (abs(hloc[j]) < DISC_h_cut) and (rloc_xy[j] < DISC_r_cut):
                f.write("{0}\n".format(rho[i][j] * DISC_mass_mult) )
            else:
                f.write("{0}\n".format(rho[i][j] * ENV_mass_mult) )
    f.close()

    return

### ------------------------------------------------------------------------ ###

def pre_cut(dat_dir):

    # Loop over grid, and set density/location value to zero if grid bin
    # is not within spatial restriction

    if not os.path.isfile(dat_dir+"/MCRT_UNCUT_dust_density.inp"):
        shutil.copy2(dat_dir+"/dust_density.inp", \
          dat_dir+"/MCRT_UNCUT_dust_density.inp")
    if not os.path.isfile(dat_dir+"/MCRT_UNCUT_loc_grid.dat"):
        shutil.copy2(dat_dir+"/loc_grid.dat", \
          dat_dir+"/MCRT_UNCUT_loc_grid.dat")

    xloc, yloc, zloc, rloc, xwid = \
       np.loadtxt(dat_dir+"/loc_grid.dat",unpack=True)
    rloc = rloc * cs.au_per_cm

    f = open(dat_dir+"/dust_density.inp","r")
    trash = f.readline()
    ngrid = int(f.readline())
    nrspecies = int(f.readline())
    rho = [[] for i in range(nrspecies)]
    for i in range(nrspecies):
        for j in range(ngrid):
            rho[i].append(float(f.readline()))
    f.close()

    f = open(dat_dir+"/dust_density.inp","w")
    f.write("1\n")
    f.write("{0}\n".format(ngrid) )
    f.write("{0}\n".format(nrspecies) )
    for i in range(nrspecies):
        for j in range(ngrid):

            if (MCRT_cut > 0.0) and (rloc[j] < MCRT_cut):
                f.write("{0}\n".format(rho[i][j]) )
            elif (MCRT_cut < 0.0) and (rloc[j] > abs(MCRT_cut)):
                f.write("{0}\n".format(rho[i][j]) )
            else:
                f.write("0.0\n")
    f.close()

    f = open(dat_dir+"/loc_grid.dat","w")
    for i in range(ngrid):

        if (MCRT_cut > 0.0) and (rloc[i] < MCRT_cut):
            f.write("{0} {1} {2} {3} {4}\n".format(xloc[i], yloc[i], zloc[i], \
              rloc[i] / cs.au_per_cm, xwid[i]) )
        elif (MCRT_cut < 0.0) and (rloc[i] > abs(MCRT_cut)):
            f.write("{0} {1} {2} {3} {4}\n".format(xloc[i], yloc[i], zloc[i], \
              rloc[i] / cs.au_per_cm, xwid[i]) )
        else:
            f.write("0.0 0.0 0.0 0.0 0.0\n")
    f.close()

    return

### ------------------------------------------------------------------------ ###

def inner_cavity(dat_dir, t_star, r_star):

    global R_DESTR

    # Loop over grid, and set density/location value to zero if grid bin
    # is within DUST_DESTR

    if not os.path.isfile(dat_dir+"/CAVITY_UNCUT_dust_density.inp"):
        shutil.copy2(dat_dir+"/dust_density.inp", \
          dat_dir+"/CAVITY_UNCUT_dust_density.inp")

    xloc, yloc, zloc, rloc, xwid = \
       np.loadtxt(dat_dir+"/loc_grid.dat",unpack=True)
    rloc *= cs.au_per_cm
    xwid *= cs.au_per_cm

    f = open(dat_dir+"/dust_density.inp","r")
    trash = f.readline()
    ngrid = int(f.readline())
    nrspecies = int(f.readline())
    rho = [[] for i in range(nrspecies)]
    for i in range(nrspecies):
        for j in range(ngrid):
            rho[i].append(float(f.readline()))
    f.close()

    # Either adopt user variable, or compute dust destruction radius. If smaller
    # than 1 AU from protostellar properties, hard set to 1 AU.

    if isinstance(R_DESTR, bool):
        R_DESTR = (r_star / 2.0) * (t_star / 1000.0 )**(6.0/2.0)
        R_DESTR *= (cs.cm_per_rsol / cs.cm_per_au)

        if (R_DESTR < 1.0):
            R_DESTR = 1.0

        print("Dust destruction radius is {0} AU".format(R_DESTR))
    elif not isinstance(R_DESTR, float):
        print("\nData type for variable R_DESTR is incorect. Either set to\n")
        print("float is hard set, or to bool if dustr destruction radius\n")
        print("is computed by merit of protostellar temperature/radius\n")

    f = open(dat_dir+"/dust_density.inp","w")
    f.write("1\n")
    f.write("{0}\n".format(ngrid) )
    f.write("{0}\n".format(nrspecies) )
    for i in range(nrspecies):
        for j in range(ngrid):

            if ( (rloc[j] - xwid[j]) < R_DESTR):
                f.write("0.0\n")
            else:
                f.write("{0}\n".format(rho[i][j]) )
    f.close()

    return R_DESTR

### ------------------------------------------------------------------------ ###

def pre_cut_scale(dat_dir):

    # As per pre_full_scale(), but for spatially restricted region

    xloc, yloc, zloc, rloc, xwid = \
       np.loadtxt(dat_dir+"/loc_grid.dat",unpack=True)
    rloc = rloc * cs.au_per_cm

    if not os.path.isfile(dat_dir+"/MCRT_CUT_UNSCALED_dust_density.inp"):
        shutil.copy2(dat_dir+"/dust_density.inp", \
          dat_dir+"/MCRT_CUT_UNSCALED_dust_density.inp")

    f = open(dat_dir+"/dust_density.inp","r")
    trash = f.readline()
    ngrid = int(f.readline())
    nrspecies = int(f.readline())

    rho = [[] for i in range(nrspecies)]
    amr_mass = 0.00
    for i in range(nrspecies):
        for j in range(ngrid):
            rho[i].append(float(f.readline()))
            amr_mass += rho[i][j] * (2.0*xwid[j])**3.0
    f.close()

    f_scale = (MCRT_cut_scale * cs.dust_to_gas) / (amr_mass * cs.msol_per_g)

    f = open(dat_dir+"/dust_density.inp","w")
    f.write("1\n")
    f.write("{0}\n".format(ngrid))
    f.write("{0}\n".format(nrspecies))
    for i in range(nrspecies):
        for j in range(ngrid):
            f.write("{0}\n".format(rho[i][j]*f_scale) )
    f.close()

    return

### ------------------------------------------------------------------------ ###

def envelope_gen(dat_dir):

    # Set an envelope on constant density based on component domain limits
    # as per simulation or user defined values

    global DISC_h_cut, DISC_r_cut, ENV_rho

    ENV_rho = ENV_rho * cs.g_per_H2 * cs.dust_to_gas

    if f_DISC_params:

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
                  "restriction is: {0} AU\n".format(float(column[11]) ) )
                print("Mass of simulation disc (radial infall restriction) "+ \
                  "is: {0} M_sol\n".format( float(column[6]) ) )
                break
            else:
                continue
        frdisc.close()

        DISC_h_cut = 0.0
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
                    DISC_h_cut += h
        hfact = 3.0
        DISC_h_cut = hfact * (DISC_h_cut / float(count) )
        print('''Scale height restriction of disc, set by {0} * average h(r)
         value for disc < DISC_r_cut is: {1}\n'''.format(hfact, DISC_h_cut) )
        fpdisc.close()

    xloc, yloc, zloc, rloc, xwid = \
       np.loadtxt(dat_dir+"/loc_grid.dat",unpack=True)
    rloc *= cs.au_per_cm
    rloc_xy = np.sqrt( xloc**2.0 + yloc**2.0 ) * cs.au_per_cm
    hloc = ( zloc * cs.au_per_cm ) / rloc_xy

    if not os.path.isfile(dat_dir+"/MCRT_UNSCALED_dust_density.inp"):
        shutil.copy2(dat_dir+"/dust_density.inp", \
          dat_dir+"/MCRT_UNSCALED_dust_density.inp")

    f = open(dat_dir+"/dust_density.inp","r")
    trash = f.readline()
    ngrid = int(f.readline())
    nrspecies = int(f.readline())

    rho = [[] for i in range(nrspecies)]
    for i in range(nrspecies):
        for j in range(ngrid):
            rho[i].append(float(f.readline()))
    f.close()

    f = open(dat_dir+"/dust_density.inp","w")
    f.write("1\n")
    f.write("{0}\n".format(ngrid))
    f.write("{0}\n".format(nrspecies))
    for i in range(nrspecies):
        for j in range(ngrid):

            if (abs(hloc[j]) < DISC_h_cut) and (rloc_xy[j] < DISC_r_cut):
                f.write("{0}\n".format(rho[i][j]) )
            else:
                f.write("{0}\n".format(ENV_rho) )
    f.close()

    return

### ------------------------------------------------------------------------ ###

def cavity_gen(dat_dir):

    # As per envelope_gen(), for bipolar cavity using domain set by user

    global CAVITY_rho, CAVITY_theta, CAVITY_beta, CAVITY_ralpha

    CAVITY_rho = CAVITY_rho * cs.g_per_H2 * cs.dust_to_gas
    CAVITY_theta = CAVITY_theta * (math.pi / 180.0)
    CAVITY_alpha = CAVITY_ralpha /  \
      (CAVITY_ralpha * math.tan(CAVITY_theta) )**CAVITY_beta

    xloc, yloc, zloc, rloc, xwid = \
       np.loadtxt(dat_dir+"/loc_grid.dat",unpack=True)
    rloc = rloc * cs.au_per_cm ; zloc = zloc * cs.au_per_cm
    xwid = xwid * cs.au_per_cm
    rloc_xy = np.sqrt(xloc**2.0 + yloc**2.0) * cs.au_per_cm

    if not os.path.isfile(dat_dir+"/MCRT_UNSCALED_dust_density.inp"):
        shutil.copy2(dat_dir+"/dust_density.inp", \
          dat_dir+"/MCRT_UNSCALED_dust_density.inp")

    f = open(dat_dir+"/dust_density.inp","r")
    iformat = int(f.readline())
    ngrid = int(f.readline())
    nrspecies = int(f.readline())
    rho = [[] for i in range(nrspecies)]
    for i in range(nrspecies):
        for j in range(ngrid):
            rho[i].append(float(f.readline()))
    f.close()

    f = open(dat_dir+"/dust_density.inp",'w')
    f.write("{0}\n".format(iformat) )
    f.write("{0}\n".format(ngrid) )
    f.write("{0}\n".format(nrspecies) )
    for i in range(nrspecies):
        for j in range(ngrid):
            zcav = CAVITY_alpha * rloc_xy[j]**(CAVITY_beta)
            if (abs(zloc[j])-xwid[j] > zcav):
                if (rho[i][j] < CAVITY_rho):
                    f.write("{0}\n".format(rho[i][j]))
                else:
                    f.write("{0}\n".format(CAVITY_rho) )
            else:
                f.write("{0}\n".format(rho[i][j]) )
    f.close()

    return

### ------------------------------------------------------------------------ ###

def post_cut(dat_dir, tdat_ext):

    # As per pre_cut(), for post-MCRT // pre-RRT (also set temperature to 0.0)

    if not os.path.isfile(dat_dir+"/RRT_UNCUT_dust_temperature"+tdat_ext):
        shutil.copy2(dat_dir+"/dust_temperature"+tdat_ext, \
          dat_dir+"/RRT_UNCUT_dust_temperature"+tdat_ext)

    xloc, yloc, zloc, rloc, xwid = \
       np.loadtxt(dat_dir+"/loc_grid.dat",unpack=True)
    rloc = rloc * cs.au_per_cm

    f = open(dat_dir+"/dust_density.inp","r")
    trash = f.readline()
    ngrid = int(f.readline())
    nrspecies = int(f.readline())
    rho = [[] for i in range(nrspecies)]
    for i in range(nrspecies):
        for j in range(ngrid):
            rho[i].append(float(f.readline()))
    f.close()

    f = open(dat_dir+"/dust_temperature"+tdat_ext,"r")
    trash = f.readline() ; trash = f.readline() ; trash = f.readline()
    temp = [[] for i in range(nrspecies)]
    for i in range(nrspecies):
        for j in range(ngrid):
            temp[i].append(float(f.readline()))
    f.close()

    f = open(dat_dir+"/dust_temperature"+tdat_ext,"w")
    f.write("1\n")
    f.write("{0}\n".format(ngrid))
    f.write("{0}\n".format(nrspecies))
    for i in range(nrspecies):
        for j in range(ngrid):
            if (RRT_cut > 0.) and (rloc[j] < RRT_cut):
                f.write("{0}\n".format(temp[i][j]) )
            elif (RRT_cut < 0.) and (rloc[j] > RRT_cut):
                f.write("{0}\n".format(temp[i][j]) )
            else:
                f.write("0.0\n")
    f.close()

    return
