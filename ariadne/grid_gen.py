''''

grid_gen.py

Module to translate density distribution from SPH particle -> AMR grid. This
uses HYPERION modules, refining AMR grid based on either the maximum refinement
level, or maximum number of SPH particles allowable in each grid point volume.

Module generates dust_density.inp, amr_grid.inp and loc_grid.dat to be used in
MCRT/RRT or diagnostic analyses.

Last Modified: 05/02/2018

'''
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - MODULE IMPORTS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
import numpy as np
import os
import shutil
import imp
import math
#
import constants as cs
from params_run import *
try:
    imp.find_module("sphtool")
    sphtool_find = True
    import sphtool
except ImportError:
    sphtool_find = False
try:
    imp.find_module("hyperion")
    from hyperion.importers.sph import construct_octree
    hyperion_find = True
except ImportError:
    hyperion_find = False
# To open the hyperion from the source, open a shell then inputs the following:
# >>> import os
# >>> import hyperion.importers.sph
# >>> os.system("atom "+hyperion.importers.sph.__file__)
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - MAIN PROGRAM - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
def translate_hyp(dat_dir, pos, rho, h, m_part, r, extras):
#
    if not hyperion_find:
        print("HYPERION modules not found, exiting...\n")
        exit()
#
    plims = [-np.amax(pos), np.amax(pos)]
    WORLD_SIZE = (plims[1]-plims[0])
    dx = 0.5*WORLD_SIZE ; dy = dx ; dz = dx
    n_octbase = 2 ; octbase = dx ; base0 = 0. - dx
    m_part = np.array( [m_part]*len(pos[0]) )
#
    # Run Hyperion octree code, using outputs to find nleafsmax and nbranchmax
#
    def stop(x, y, z, dx, dy, dz, px, py, pz, sigma):
        return len(px) <= partmax
#
    o = construct_octree(0., 0., 0., dx, dy, dz, \
       np.array(pos[0]), np.array(pos[1]), np.array(pos[2]), \
       np.array(h), m_part, n_levels = levelmax, stopping_criterion=stop)
#
    xcen = np.array(o["xcen"][0].array) ; ycen = np.array(o["ycen"][0].array)
    zcen = np.array(o["zcen"][0].array) ; xwid = np.array(o["xwid"][0].array)
    density = np.array(o["density"][0].array)
#
    n = len(o.refined)
    nleafsmax = 0 ; nbranchmax = 0
    for i in range(0, n):
        if o.refined[i] == True:
            nbranchmax += 1
        elif o.refined[i] == False:
            nleafsmax += 1
    ncells = nleafsmax
#
    # Use results to ouput grid data into amr_grid.inp for RADMC-3D
#
    print("Writing to amr_grid.inp\n")
    f = open(dat_dir+"/amr_grid.inp","w")
    f.write("1\n")
    f.write("1\n")
    f.write("0\n")
    f.write("0\n")
    f.write("1 1 1\n")
    f.write("1 1 1\n ")
    f.write("{0} {1} {1}\n".format(levelmax, n) )
    f.write("{0} {1}\n".format(base0, -base0) )
    f.write("{0} {1}\n".format(base0, -base0) )
    f.write("{0} {1}\n".format(base0, -base0) )
    for i in range(0, n):
        if o.refined[i] == True:
            f.write("1\n")
        else:
            f.write("0\n")
    f.close()
#
    # Use results to ouput grid to indicate radial locations of leaf nodes, then
    # write to file for any desired I/O analyses.
#
    print("Writing to loc_grid.inp\n")
#
    rgrid = [0. * n]
    if not sim_inp and (mod_inp == "ivezic_env"):
        r_dust = extras[1]
        for i in range(n):
            if not o.refined[i]:
                xcen[i] = xcen[i] / r_dust
                ycen[i] = ycen[i] / r_dust
                zcen[i] = zcen[i] / r_dust
                rgrid = np.sqrt( xcen**2. + ycen**2. + zcen**2. )
                xwid[i] = xwid[i] / r_dust
    else:
        rgrid = np.sqrt( xcen**2. + ycen**2. + zcen**2. )
#
    f = open(dat_dir+"/loc_grid.dat","w")
    for i in range(0, n):
        if not o.refined[i]:
            f.write("{0} {1} {2} {3} {4}\n".format(xcen[i], ycen[i], zcen[i], \
              rgrid[i], 0.5 * xwid[i]) )
        else:
            continue
    f.close()
#
    # With output grid structure from HYPERION modules, write AMR grid density
    # structure to file as required for RADMC-3D
#
    print("Writing to dust_density.inp\n")
    ngrid = ncells
    f = open(dat_dir+"/dust_density.inp","w")
    f.write("1\n")
    f.write("{0}\n".format(ngrid) )
    f.write("1\n")
    for i in range(n):
        if not o.refined[i]:
            f.write("{0}\n".format(density[i]) )
        else:
            continue
    f.close()
#
    return
#
### ------------------------------------------------------------------------ ###
#
def translate_spht(dat_dir, pos, rho, h, m_part, r, extras):
#
    if not sphtool_find:
        print("SPHTOOL modules not found, exiting...\n")
        exit()
#
    m_part = np.array([m_part for i in range(len(pos[0]))] )
    pos = np.array(pos) ; rho = np.array(rho)
#
    tool = sphtool.SPHTool(nThreads=threads)
    tool.setSPHData(x=pos[0], y=pos[1], z=pos[2], h=h, rho=rho, pmass=m_part)
    tool.init()
    scalarFloor = np.zeros(1, dtype=float)+1e-30
#
    if rhograd:
        tool.regridToAMRGrid(refineDensityGradient=rhograd, \
           nDensTrialPoints=n_rhograd, densVarMax = rhograd_var, \
           maxAMRParticlesPerCell=partmax, maxAMRTreeDepth=levelmax, \
           scalarFloor=scalarFloor)
    else:
        tool.regridToAMRGrid(maxAMRParticlesPerCell=partmax, \
           maxAMRTreeDepth=levelmax, scalarFloor=scalarFloor)
#
    tool.writeGrid()
    tool.writeScalar(fname="dust_density.inp", ivar=0, nscal=1, binary=False)
#
    shutil.move(os.getcwd()+"/amr_grid.inp", dat_dir+"/amr_grid.inp")
    shutil.move(os.getcwd()+"/dust_density.inp", dat_dir+"/dust_density.inp")
#
    return
#
### ------------------------------------------------------------------------ ###
#
def protostar(dat_dir):
#
    # Simply set the amr_grid.inp based on a single 1AU**3 box, with zero
    # density set for grid point in dust_density.
#
    f = open(dat_dir+"/amr_grid.inp","w")
    f.write("1\n")
    f.write("0\n")
    f.write("0\n")
    f.write("0\n")
    f.write("1 1 1\n")
    f.write("1 1 1\n")
    f.write("-1.496e13 1.496e13\n")
    f.write("-1.496e13 1.496e13\n")
    f.write("-1.496e13 1.496e13\n")
    f.close()
#
    f = open(dat_dir+"/dust_density.inp","w")
    f.write("1\n")
    f.write("1\n")
    f.write("1\n")
    f.write("0\n")
    f.close()
#
    f = open(dat_dir+"/dust_temperature.dat","w")
    f.write("1\n")
    f.write("1\n")
    f.write("1\n")
    f.write("0\n")
    f.close()
#
    return
#
### ------------------------------------------------------------------------ ###
#
def yso(dat_dir):
#
    # Firstly, set up base grid for core, and bipolar outflow
#
    rho_base = [ [[0.0] for i in range(ngrid)] for i in range(3)]
    loc_base = [ [[0.0] for i in range(ngrid+1)] for i in range(3)]
#
    for d in range(3):
        for i in range(ngrid+1):
            loc_base[d][i] = -r_out4 + ( i * ( (2.0 * r_out4) / ngrid) )
            loc_base[d][i] *= cs.cm_per_au
#
    print("Not yet functioning, exiting...\n")
    exit()
#
    return
#
### ------------------------------------------------------------------------ ###
#
def pre_cut(dat_dir):
#
    # Duplicate both density and temperature files, producing uncut files
#
    if not os.path.isfile(dat_dir+"/UNCUT_dust_density.inp"):
        shutil.copy2(dat_dir+"/dust_density.inp", \
          dat_dir+"/UNCUT_dust_density.inp")
    if not os.path.isfile(dat_dir+"/UNCUT_loc_grid.dat"):
        shutil.copy2(dat_dir+"/loc_grid.dat", dat_dir+"/UNCUT_loc_grid.dat")
#
    # Now read in all location/density/temperature data
#
    xloc, yloc, zloc, rloc, xwid = \
       np.loadtxt(dat_dir+"/loc_grid.dat",unpack=True)
    rloc = rloc * cs.au_per_cm
#
    f = open(dat_dir+"/dust_density.inp","r")
    trash = f.readline()
    ngrid = int(f.readline())
    nrspecies = int(f.readline())
    rho = [[] for i in range(nrspecies)]
    for i in range(nrspecies):
        for j in range(ngrid):
            rho[i].append(float(f.readline()))
    f.close()
#
    # Now reloop over entire grid, re-writing density/temperature/loc_grid files
    # based on whether bins fall in RRT_cut
#
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
#
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
#
    return
#
### ------------------------------------------------------------------------ ###
#
def pre_full_scale(arch_dir, dat_dir):
#
    # Read in AMR grid data
#
    xloc, yloc, zloc, rloc, xwid = \
       np.loadtxt(dat_dir+"/loc_grid.dat",unpack=True)
    rloc = rloc * cs.au_per_cm
#
    # Copy non-scaled dataset here, for "clean" density distributions
#
    if not os.path.isfile(dat_dir+"/UNSCALED_dust_density.inp"):
        shutil.copy2(dat_dir+"/dust_density.inp", \
          dat_dir+"/UNSCALED_dust_density.inp")
    if not os.path.isfile(dat_dir+"/UNSCALED_loc_grid.dat"):
        shutil.copy2(dat_dir+"/loc_grid.dat", dat_dir+"/UNSCALED_loc_grid.dat")
#
    f = open(dat_dir+"/dust_density.inp","r")
    trash = f.readline()
    ngrid = int(f.readline())
    nrspecies = int(f.readline())
    rho = [[] for i in range(nrspecies)]
    for i in range(nrspecies):
        for j in range(ngrid):
            rho[i].append(float(f.readline()))
    f.close()
#
    # For SPH particle data, calculate total mass of matter in $MCRT_cut$
#
    os.system("{0} to ascii {1}".format(xsplash, sim_cfile) )
    inp_file = sim_cfile+".ascii"
#
    sim_mass = 0.00
    f = open(inp_file, "r")
    for i in range(12):
        trash = f.readline()
    for lines in f:
        line = lines.strip() ; columns = line.split()
        if (float(columns[10]) == 1):
            x = float(columns[0])
            y = float(columns[1])
            z = float(columns[2])
            mpart = float(columns[9])
            r = ( np.sqrt(x**2.0 + y**2.0 + z**2.0) * cs.au_per_cm )
#
    # Add mass for particles in radial (un)cut region, to rescale
    #
    # Note: MCRT_cut_scale (Bool.) is used to ensure that MCRT_cut region
    # density is consistent with uncut simulation mass scaled to FULL_mass_scale
#
            if not MCRT_cut:
                sim_mass += mpart
            elif MCRT_cut and not MCRT_cut_scale:
                sim_mass += mpart
            elif MCRT_cut and MCRT_cut_scale \
              and (MCRT_cut > 0.0) and (r < MCRT_cut):
                sim_mass += mpart
            elif MCRT_cut and MCRT_cut_scale \
              and (MCRT_cut < 0.0) and (r < abs(MCRT_cut)):
                sim_mass += mpart
#
            else:
                continue
        else:
            continue
    f.close()
    os.remove(inp_file)
#
    # Set scaling factor based on above, and rewrite dust_density.inp
#
    f_scale = (FULL_mass_scale / sim_mass)
#
    f = open(dat_dir+"/dust_density.inp","w")
    f.write("1\n")
    f.write("{0}\n".format(ngrid))
    f.write("{0}\n".format(nrspecies))
    for i in range(nrspecies):
        for j in range(ngrid):
            f.write("{0}\n".format(rho[i][j]*f_scale) )
#
    f.close()
#
    return
#
#
### ------------------------------------------------------------------------ ###
#
def pre_comp_scale(arch_dir, dat_dir):
#
    global DISC_h_cut, DISC_r_cut
#
    if f_DISC_params:
#
    # Define/read files to be read for sim specific disc parameters, then
    # replace DISC_h_cut and DISC_r_cut parameters
#
        if not (os.path.isfile(sim_pfile) or os.path.isfile(sim_r1file)):
            print("rdisc.1 or pdisc.1 files can't be found, exiting...")
            exit()
#
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
                print("Mass of simulation envelope (10 000 AU) "+ \
                  "is: {0} M_sol\n".format( float(column[19]) ) )
                break
            else:
                continue
        frdisc.close()
#
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
        hfact = 1.0
        DISC_h_cut = hfact * (DISC_h_cut / float(count) )
        print("Scale height restriction of disc, set by {0} * average h(r) "+ \
          "value for disc < DISC_r_cut is: {1}\n".format(hfact, DISC_h_cut) )
        fpdisc.close()
#
    # Read in AMR grid data
#
    xloc, yloc, zloc, rloc, xwid = \
       np.loadtxt(dat_dir+"/loc_grid.dat",unpack=True)
    rloc *= cs.au_per_cm
    rloc_xy = np.sqrt( xloc**2.0 + yloc**2.0 ) * cs.au_per_cm
    hloc = ( zloc * cs.au_per_cm ) / rloc_xy
#
    # Copy non-scaled dataset here, for "clean" density distributions
#
    if not os.path.isfile(dat_dir+"/UNSCALED_dust_density.inp"):
        shutil.copy2(dat_dir+"/dust_density.inp", \
          dat_dir+"/UNSCALED_dust_density.inp")
    if not os.path.isfile(dat_dir+"/UNSCALED_loc_grid.dat"):
        shutil.copy2(dat_dir+"/loc_grid.dat", dat_dir+"/UNSCALED_loc_grid.dat")
#
    f = open(dat_dir+"/dust_density.inp","r")
    trash = f.readline()
    ngrid = int(f.readline())
    nrspecies = int(f.readline())
    rho = [[] for i in range(nrspecies)]
    mass = 0
    for i in range(nrspecies):
        for j in range(ngrid):
            rho[i].append(float(f.readline()))
            mass += rho[i][j] * (2.0*xwid[j])**3.0
    f.close()
    mass /= cs.g_per_msol * cs.dust_to_gas ; print mass
#
    os.system("{0} to ascii {1}".format(xsplash, sim_cfile) )
    inp_file = sim_cfile+".ascii"
#
    # For SPH particle data, calculate total mass of matter in $MCRT_cut$
#
    envelope_mass = 0.00 ; disc_mass = 0.0
    f = open(inp_file, "r")
    for i in range(12):
        trash = f.readline()
    for lines in f:
        line = lines.strip() ; columns = line.split()
        if (float(columns[10]) == 1):
            x = float(columns[0])
            y = float(columns[1])
            z = float(columns[2])
            mpart = float(columns[9])
            r = np.sqrt(x**2.0 + y**2.0 + z**2.0) * cs.au_per_cm
            r_xy = np.sqrt(x**2.0 + y**2.0) * cs.au_per_cm
            h = (z * cs.au_per_cm) / r_xy
#
            if MCRT_cut and (MCRT_cut > 0.0) and (r < MCRT_cut):
                if (abs(h) < DISC_h_cut) and (r_xy < DISC_r_cut):
                    disc_mass += mpart
                else:
                    envelope_mass += mpart
#
            elif MCRT_cut and (MCRT_cut < 0.0) and (r > MCRT_cut):
                if (abs(h) < DISC_h_cut) and (r_xy < DISC_r_cut):
                    disc_mass += mpart
                else:
                    envelope_mass += mpart
#
            elif not MCRT_cut:
                if (abs(h) < DISC_h_cut) and (r_xy < DISC_r_cut):
                    disc_mass += mpart
                else:
                    envelope_mass += mpart
#
            else:
                continue
        else:
            continue
    f.close()
    os.remove(inp_file)
#
    # Set scaling factor based on above, and rewrite dust_density.inp
#
    f_scale_disc = (DISC_mass_scale / disc_mass)
    f_scale_env = (ENV_mass_scale / envelope_mass)
    print disc_mass, envelope_mass
    print f_scale_disc, f_scale_env
#
    f = open(dat_dir+"/dust_density.inp","w")
    f.write("1\n")
    f.write("{0}\n".format(ngrid))
    f.write("{0}\n".format(nrspecies))
    for i in range(nrspecies):
        for j in range(ngrid):
#
            if MCRT_cut and (MCRT_cut > 0.0) and (rloc[j] < MCRT_cut):
                if (abs(hloc[j]) < DISC_h_cut) and (rloc_xy[j] < DISC_r_cut):
                    f.write("{0}\n".format(rho[i][j] * f_scale_disc) )
                else:
                    f.write("{0}\n".format(rho[i][j] * f_scale_env) )
#
            elif MCRT_cut and (MCRT_cut < 0.0) and (rloc[j] > MCRT_cut):
                if (abs(hloc[j]) < DISC_h_cut) and (rloc_xy[j] < DISC_r_cut):
                    f.write("{0}\n".format(rho[i][j] * f_scale_disc) )
                else:
                    f.write("{0}\n".format(rho[i][j] * f_scale_env) )
#
            elif not MCRT_cut:
                if (abs(hloc[j]) < DISC_h_cut) and (rloc_xy[j] < DISC_r_cut):
                    f.write("{0}\n".format(rho[i][j] * f_scale_disc) )
                else:
                    f.write("{0}\n".format(rho[i][j] * f_scale_env) )
            else:
                f.write("0.0\n")
#
    f.close()
#
    return
#
### ------------------------------------------------------------------------ ###
#
def cavity_gen(dat_dir):
#
    global CAVITY_rho, CAVITY_theta, CAVITY_beta, CAVITY_alpha, CAVITY_ralpha
#
    CAVITY_rho = CAVITY_rho * cs.g_per_H2 * cs.dust_to_gas
    CAVITY_theta = CAVITY_theta * (math.pi / 180.0)
    CAVITY_alpha = (CAVITY_ralpha**(1.0 - CAVITY_beta)) / math.tan(CAVITY_theta)
#
    xloc, yloc, zloc, rloc, xwid = \
       np.loadtxt(dat_dir+"/loc_grid.dat",unpack=True)
#
    rloc = rloc * cs.au_per_cm ; zloc = zloc * cs.au_per_cm
    xwid = xwid * cs.au_per_cm
    rloc_xy = np.sqrt(xloc**2.0 + yloc**2.0) * cs.au_per_cm
#
    # Copy non-scaled dataset here, for "clean" density distributions
#
    if not os.path.isfile(dat_dir+"/UNSCALED_dust_density.inp"):
        shutil.copy2(dat_dir+"/dust_density.inp", \
          dat_dir+"/UNSCALED_dust_density.inp")
    if not os.path.isfile(dat_dir+"/UNSCALED_loc_grid.dat"):
        shutil.copy2(dat_dir+"/loc_grid.dat", dat_dir+"/UNSCALED_loc_grid.dat")
#
    f = open(dat_dir+"/dust_density.inp","r")
    iformat = int(f.readline())
    ngrid = int(f.readline())
    nrspecies = int(f.readline())
    rho = [[] for i in range(nrspecies)]
    for i in range(nrspecies):
        for j in range(ngrid):
            rho[i].append(float(f.readline()))
    f.close()
#
    f = open(dat_dir+"/dust_density.inp",'w')
    f.write("{0}\n".format(iformat) )
    f.write("{0}\n".format(ngrid) )
    f.write("{0}\n".format(nrspecies) )
    for i in range(nrspecies):
        for j in range(ngrid):
            zcav = CAVITY_alpha * rloc_xy[j]**(CAVITY_beta)
            if (abs(zloc[j])-xwid[j] < zcav):
                f.write("{0}\n".format(rho[i][j]) )
            elif (abs(zloc[j])-xwid[j] > zcav):
                f.write("{0}\n".format(CAVITY_rho) )
#
    f.close()
#
    return
#
### ------------------------------------------------------------------------ ###
#
def post_cut(dat_dir, tdat_ext):
#
    # Duplicate both density and temperature files, producing uncut files
#
    if not os.path.isfile(dat_dir+"/UNCUT_dust_density.inp"):
        shutil.copy2(dat_dir+"/dust_density.inp", \
          dat_dir+"/UNCUT_dust_density.inp")
    if not os.path.isfile(dat_dir+"/UNCUT_dust_temperature"+tdat_ext):
        shutil.copy2(dat_dir+"/dust_temperature"+tdat_ext, \
          dat_dir+"/UNCUT_dust_temperature"+tdat_ext)
    if not os.path.isfile(dat_dir+"/UNCUT_loc_grid.dat"):
        shutil.copy2(dat_dir+"/loc_grid.dat", dat_dir+"/UNCUT_loc_grid.dat")
#
    # Now read in all location/density/temperature data
#
    xloc, yloc, zloc, rloc, xwid = \
       np.loadtxt(dat_dir+"/loc_grid.dat",unpack=True)
    rloc = rloc * cs.au_per_cm
#
    f = open(dat_dir+"/dust_density.inp","r")
    trash = f.readline()
    ngrid = int(f.readline())
    nrspecies = int(f.readline())
    rho = [[] for i in range(nrspecies)]
    for i in range(nrspecies):
        for j in range(ngrid):
            rho[i].append(float(f.readline()))
    f.close()
#
    f = open(dat_dir+"/dust_temperature"+tdat_ext,"r")
    trash = f.readline() ; trash = f.readline() ; trash = f.readline()
    temp = [[] for i in range(nrspecies)]
    for i in range(nrspecies):
        for j in range(ngrid):
            temp[i].append(float(f.readline()))
    f.close()
#
    # Now reloop over entire grid, re-writing density/temperature/loc_grid files
    # based on whether bins fall in RRT_cut.
    # If called for, rescale post-cut mass in this loop too.
#
    f = open(dat_dir+"/dust_density.inp","w")
    f.write("1\n")
    f.write("{0}\n".format(ngrid))
    f.write("{0}\n".format(nrspecies))
    for i in range(nrspecies):
        for j in range(ngrid):
            if (rloc[j] < RRT_cut) and RRT_cut_scale:
                f.write("{0}\n".format(rho[i][j] * RRT_cut_scale) )
            elif (rloc[j] < RRT_cut) and not RRT_cut_scale:
                f.write("{0}\n".format(rho[i][j]) )
            else:
                f.write("0.0\n")
    f.close()
#
    f = open(dat_dir+"/dust_temperature"+tdat_ext,"w")
    f.write("1\n")
    f.write("{0}\n".format(ngrid))
    f.write("{0}\n".format(nrspecies))
    for i in range(nrspecies):
        for j in range(ngrid):
            if (rloc[j] < RRT_cut):
                f.write("{0}\n".format(temp[i][j]) )
            else:
                f.write("0.0\n")
    f.close()
#
    f = open(dat_dir+"/loc_grid.dat","w")
    for i in range(ngrid):
        if (rloc[i] < RRT_cut):
            f.write("{0} {1} {2} {3} {4}\n".format(xloc[i], yloc[i], zloc[i], \
              rloc[i] / cs.au_per_cm, xwid[i]) )
        else:
            f.write("0.0 0.0 0.0 0.0 0.0\n")
    f.close()
#
    return
