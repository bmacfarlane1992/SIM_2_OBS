''''

grid_gen.py

Module to translate density distribution from SPH particle -> AMR grid. This
uses HYPERION modules, refining AMR grid based on either the maximum refinement
level, or maximum number of SPH particles allowable in each grid point volume.

Module generates dust_density.inp, amr_grid.inp and loc_grid.dat to be used in
MCRT/RRT or diagnostic analyses.

Module can also generate required files for MCRT and RRT in protostar and
YSO modelling, where SPH treatment is not needed/overkill.

Last Modified: 21/08/2018

'''

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - MODULE IMPORTS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

import numpy as np
import os
import sys
import shutil
import imp
import math
import random
import matplotlib.pyplot as plt

import constants as cs
import params_run as ps
sys.path.append(os.getcwd()+"/oct/")
from gen_oct import construct_octree, compute_octree_geometry

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - MAIN PROGRAM - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

def translate(globals, inp_model):

    # Translate SPH -> AMR density, writing outputs to RADMC-3D input files

    dx = np.amax(inp_model.pos) ; dy = dx ; dz = dx
    m_part = np.array( [inp_model.m_part]*len(inp_model.pos[0]) )

    def stop(x, y, z, dx, dy, dz, px, py, pz, sigma):
        return len(px) <= partmax

    o = construct_octree(0.0, 0.0, 0.0, dx, dy, dz, \
     np.array(inp_model.pos[0]), np.array(inp_model.pos[1]), \
     np.array(inp_model.pos[2]), np.array(inp_model.h), np.array(m_part), \
     n_levels = levelmax, stopping_criterion=stop)

    n_amr = len(o.refined) ; n_leaf = len(np.where(o.refined == False)[0])

    xcen = np.array(o["xcen"][0].array) ; ycen = np.array(o["ycen"][0].array)
    zcen = np.array(o["zcen"][0].array) ; xwid = np.array(o["xwid"][0].array)
    density = np.array(o["density"][0].array)

    print("Writing to amr_grid.inp\n")
    f = open(globals.dat_dir+"/amr_grid.inp","w")
    f.write("1\n")
    f.write("1\n")
    f.write("0\n")
    f.write("0\n")
    f.write("1 1 1\n")
    f.write("1 1 1\n ")
    f.write("{0} {1} {1}\n".format(levelmax, n_amr) )
    f.write("{0} {1}\n".format(-dx, dx) )
    f.write("{0} {1}\n".format(-dx, dx) )
    f.write("{0} {1}\n".format(-dx, dx) )
    for i in range(n_amr):
        if o.refined[i]:
            f.write("1\n")
        else:
            f.write("0\n")
    f.close()

    print("Writing to loc_grid.inp\n")
    rgrid = np.sqrt( xcen**2. + ycen**2. + zcen**2. )
    f = open(globals.dat_dir+"/loc_grid.dat","w")
    for i in range(n_amr):
        if o.refined[i]:
            continue
        else:
            f.write("{0} {1} {2} {3} {4}\n".format(xcen[i], ycen[i], zcen[i], \
              rgrid[i], 0.5 * xwid[i]) )
    f.close()

    print("Writing to dust_density.inp\n")
    f = open(globals.dat_dir+"/dust_density.inp","w")
    f.write("1\n")
    f.write("{0}\n".format(n_leaf) )
    f.write("1\n")
    for i in range(n_amr):
        if o.refined[i]:
            continue
        else:
            f.write("{0}\n".format(density[i]) )
    f.close()

    return

### ------------------------------------------------------------------------ ###

def protostar(path, rt_params):

    # Generate dummy RADMC-3D input files

    f = open(path+"/amr_grid.inp","w")
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

    f = open(path+"/dust_density.inp","w")
    f.write("1\n")
    f.write("1\n")
    f.write("1\n")
    f.write("0\n")
    f.close()

    f = open(path+"/dust_temperature.dat","w")
    f.write("1\n")
    f.write("1\n")
    f.write("1\n")
    f.write("0\n")
    f.close()

    f = open(path+"/dustopac.inp","w")
    f.write("2 \n")
    f.write("1 \n")
    f.write("==========================================================\n")
    f.write("1 \n")
    f.write("0 \n")
    f.write("dum\n")
    f.write("------------------------------------------------------\n")
    f.close()

    f = open(path+"/dustkappa_dum.inp","w")
    f.write("2\n")
    f.write("{0}\n".format(rt_params.n_wavs) )
    for i in range(rt_params.n_wavs):
        f.write("{0} 1 1\n".format(rt_params.wavs[i]) )
    f.close()

    f = open(path+"/stars.inp","w")
    f.write("2 \n")
    f.write("1 {0}\n".format(rt_params.n_wavs) )
    f.write("{0} 1.998e+33 0. 0. 0.\n".format(rt_params.r_star*cs.cm_per_rsol) )
    for i in range(rt_params.n_wavs):
        f.write("{0}\n".format(rt_params.wavs[i]) )
    f.write("-{0}\n".format(rt_params.t_star) )
    f.close()

    return

### ------------------------------------------------------------------------ ###

def yso(globals, rt_model):

    print("Creating YSO disc/envelope\n")

    pos = [ [ [0.0] for i in range(2*ps.nsph)] for i in range(3) ]

    # Determine inner/outer limits for disc/envelope

    r_0_d = ( (rt_model.r_star / 2.0) * (rt_model.t_star / ps.t_dust)**3.00 ) \
     * cs.cm_per_rsol
    r_in_d = r_0_d
    r_out_d = ps.r_out_d * cs.cm_per_au
    w_in_d = (r_in_d**2.0) / (r_out_d**2.0)
    w_out_d = (r_out_d**2.0) / (r_0_d**2.0)

    r_in_e = ps.r_in_e * cs.cm_per_au
    r_out_e = ps.r_out_e * cs.cm_per_au
    exp = 1.0 - (float(ps.p_e) / 3.0)
    w_e = r_out_e**(3.0) / r_in_e**(3.0)

    for i in range(ps.nsph):

    # Generate particle locations for disc

        Rs1 = random.random() ; Rs2 = random.random() ; Rs3 = random.random()
        w = ( (1.0 + w_in_d)**(1.0 - (ps.p_d / 2.0) ) + \
         Rs1 * ( (1.0 + w_out_d)**(1.0 - (ps.p_d / 2.0) ) - \
         (1.0 + w_in_d)**(1.0 - (ps.p_d / 2.0) ) ))**( 2.0 / (2.0 - ps.p_d) )-1.0

        r = r_0_d * w**(0.5)
        phi = 2.0 * math.pi * Rs2

        pos[0][i] = r * math.cos( phi )
        pos[1][i] = r * math.sin( phi )
        z_0 = 3.0 * ps.alpha_d * r
        pos[2][i] = z_0 * (2.0 / math.pi) * math.asin( 2.0 * Rs3 - 1.0)

    # Generate particle locations for envelope (with new random seeds)

        Rs1 = random.random() ; Rs2 = random.random() ; Rs3 = random.random()
        w = ( ( ( w_e**exp - 1.0) * Rs1 ) + 1.0 )**(1.0/exp)
        r = r_in_e * ( w**(1.0/3.0) )
        theta = math.acos(1.0 - (2.0 * Rs2) )
        phi = 2.0 * math.pi * Rs3
        pos[0][i+nsph] = ( r * (math.sin(theta) * math.cos( phi ) ) )
        pos[1][i+nsph] = ( r * ( math.sin(theta) * math.sin( phi ) ) )
        pos[2][i+nsph] = ( r * math.cos(theta))

    # Generate octree geometry grid particles, with dummy mass//sigma arrays

    dx = np.amax(pos) ; dy = dx ; dz = dx
    mass = [ [0.0] for i in range(len(pos[0]) ) ] ; sigma = mass

    def stop(x, y, z, dx, dy, dz, px, py, pz, sigma):
        return len(px) <= partmax

    o = compute_octree_geometry(0.0, 0.0, 0.0, dx, dy, dz, \
      np.array(pos[0]), np.array(pos[1]), np.array(pos[2]), \
      np.array(sigma), np.array(mass), \
      n_levels = levelmax, stopping_criterion=stop)

    n_amr = len(o.refined) ; n_leaf = len(np.where(o.refined == False)[0])

    xcen = np.array(o["xcen"][0].array) ; ycen = np.array(o["ycen"][0].array)
    zcen = np.array(o["zcen"][0].array) ; xwid = np.array(o["xwid"][0].array)

    # Write to loc_grid.dat and amr_grid.inp with grid octree information

    print("Writing to amr_grid.inp\n")
    f = open(globals.dat_dir+"/amr_grid.inp","w")
    f.write("1\n")
    f.write("1\n")
    f.write("0\n")
    f.write("0\n")
    f.write("1 1 1\n")
    f.write("1 1 1\n ")
    f.write("{0} {1} {1}\n".format(levelmax, n_amr ) )
    f.write("{0} {1}\n".format(-dx, dx) )
    f.write("{0} {1}\n".format(-dx, dx) )
    f.write("{0} {1}\n".format(-dx, dx) )
    for i in range(n_amr):
        if o.refined[i]:
            f.write("1\n")
        else:
            f.write("0\n")
    f.close()

    print("Writing to loc_grid.inp\n")
    rgrid = np.sqrt( xcen**2. + ycen**2. + zcen**2. )
    f = open(globals.dat_dir+"/loc_grid.dat","w")
    for i in range(n_amr):
        if o.refined[i]:
            continue
        else:
            f.write("{0} {1} {2} {3} {4}\n".format(xcen[i], ycen[i], zcen[i], \
              rgrid[i], 0.5 * xwid[i]) )
    f.close()

    # Allocate density contributions at each location, writing to file

    print("Writing to dust_density.inp\n")
    f = open(globals.dat_dir+"/dust_density.inp","w")
    f.write("1\n")
    f.write("{0}\n".format( n_leaf ) )
    f.write("1\n")

    # Disc scaling properties

    m_d = ps.m_d * (cs.g_per_msol * cs.dust_to_gas)

    SD_0_d = ( (m_d * (2.0-ps.p_d)) / (2.0 * math.pi*r_0_d**(2.0) ) ) * \
     ( ((r_0_d**2.0 + r_out_d**2.0)/ r_0_d**2.0)**(ps.p_d/2.0) - \
     ((r_0_d**2.0 + r_in_d**2.0)/ r_0_d**2.0)**(ps.p_d/2.0) )**(-1.0)

    # Envelope scaling properties

    m_e =  ps.m_e * (cs.g_per_msol * cs.dust_to_gas)
    rho_0_e = ( (m_e  * (3.0 * exp) ) / (4.0 * math.pi * r_in_e**(3.0)) ) * \
      (w_e**(exp) - 1.0)**(-1.0)
#    rho_0_e = ( m_e / (4.0 * math.pi * r_0_d**2.0) ) * \
#      ( ((r_0_d**2.0 + r_out_d**2.0)/ r_0_d**2.0)**(0.5) - 1.0 )**(-1.0)

    # Cavity scaling properties

    CAVITY_rho = ps.CAVITY_rho * cs.g_per_H2 * cs.dust_to_gas
    CAVITY_theta = ps.CAVITY_theta * (math.pi / 180.0)
    CAVITY_alpha = ((ps.CAVITY_ralpha*cs.cm_per_au)**(1.0 - ps.CAVITY_beta)) / \
      math.tan(CAVITY_theta)

    for i in range(len(xcen)):
        if not o.refined[i]:
            rho = 0.0
            rgrid_xy = np.sqrt(xcen[i]**2.0 + ycen[i])

    # Envelope contributions

            rho += rho_0_e * (r_in_e / rgrid[i])**ps.p_e

    # Disc contribution

            if ( rgrid_xy < r_out_d ) and \
              ( abs(zcen[i]) / rgrid_xy  <  3.0 * ps.alpha_d ):

                rho_0_d = (math.pi * SD_0_d) / \
                  (4.0 * (3.0 * ps.alpha_d * rgrid_xy) ) * \
                  (r_0_d / (rgrid[i] + r_0_d) )**(ps.p_d / 2.0)

                rho += rho_0_d * (rgrid_xy / r_0_d)**(-2.0) * \
                  math.exp((-1.0/2.0) * (zcen[i]/(ps.alpha_d * rgrid_xy) )**2.0 )

    # Cavity (if CAVITY_rho != 0.0)

            if CAVITY_rho:
                zcav = CAVITY_alpha * \
                  np.sqrt(xcen[i]**2.0 + ycen[i]**2.0)**(CAVITY_beta)
                if (abs(zcen[i]) - xwid[i] < zcav):
                    rho = rho
                elif (abs(zcen[i]) - xwid[i] > zcav):
                    rho = CAVITY_rho

            f.write("{0}\n".format(rho) )

    f.close()

    return
