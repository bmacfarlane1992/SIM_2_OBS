'''

dens_gen.py

Module to generate/read the density structure of an SPH model/simulation. For
simulations, use SPLASH to convert cdisc file to ascii format and read. For
models, generate SPH distribution. With respective SPH distributions, write
into arrays for SPH particle -> AMR grid translation.

Density distribution is for dust mass as per requirements of RADMC-3D.

Last Modified: 21/08/2018

'''

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - MODULE IMPORTS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

import os
import numpy as np
import math
import random

import constants as cs
import params_run as ps

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - MAIN PROGRAM - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

# Function that reads in simulation, and outputs SPH gas data to arrays for
# further analyses in RT modelling.

def sim():

    os.system("{0} to ascii {1}".format(ps.xsplash, ps.sim_cfile) )
    inp_file = ps.sim_cfile+".ascii"

    f = open(inp_file, "r")
    for i in range(12):
        trash = f.readline()
    for lines in f:
        line = lines.strip() ; columns = line.split()

        if (float(columns[10]) == 1):
            inp_model.pos[0].append(float(columns[0]) * cs.cm_per_pc)
            inp_model.pos[1].append(float(columns[1]) * cs.cm_per_pc)
            inp_model.pos[2].append(float(columns[2]) * cs.cm_per_pc)
            inp_model.m_part = float(columns[9]) * cs.g_per_msol * cs.dust_to_gas
            inp_model.h.append(float(columns[7]) * cs.cm_per_pc)
            inp_model.rho.append(float(columns[8]))
            inp_model.r.append( np.sqrt(float(columns[0])**2.0 + \
               float(columns[1])**2.0 + float(columns[2])**2.0  ) )
        else:
            continue
    f.close()
    os.remove(inp_file)

    return inp_model

### ------------------------------------------------------------------------ ###

    # Main modelling function, which calls mod_inp dependent functions.

def model():

    N_neigh = 50.0         # Float: No. of neighbours density is smoothed over

    if (mod_inp == "disc"):
        inp_model = disc(N_neigh)
    elif (mod_inp == "env"):
        inp_model = envelope(N_neigh)
    elif (mod_inp == "plum"):
        inp_model = plum_sphere(N_neigh)
    else:
        print("Incorrect mod_inp selection made, exiting program\n")
        exit()

    return inp_model

### ------------------------------------------------------------------------ ###

def disc(N_neigh):

    print("Creating disc\n")

    pos = [ [[0.0] for i in range(ps.nsph)] for i in range(3)]
    r = [0.0] * ps.nsph ; rho = [0.0] * ps.nsph ; h = [0.0] * ps.nsph

    # Compute particle *gas* mass, then convert mass/distance values to c.g.s

    m_d = ps.m_d * cs.g_per_msol
    m_part = m_d / float(ps.nsph)
    m_star = ps.m_star * cs.g_per_msol
    r_0 = ps.r_0_d * cs.cm_per_au
    r_in = ps.r_in_d * cs.cm_per_au
    r_out = ps.r_out_d * cs.cm_per_au

    # Define w values from which radial profiles can be attained

    w_in=(r_in**2.0) / (r_out**2.0)
    w_out=(r_out**2.0) / (r_0**2.0)

    for i in range(ps.nsph):

        Rs1 = random.random() ; Rs2 = random.random() ; Rs3 = random.random()

        w_r = ( (1.0 + w_in)**(1.0 - (ps.p_d / 2.0) ) + \
         Rs1 * ( (1.0 + w_out)**(1.0 - (ps.p_d / 2.0) ) - \
         (1.0 + w_in)**(1.0 - (ps.p_d / 2.0) ) ))**( 2.0 / (2.0 - ps.p_d) ) -1.0

        r[i] = r_0 * w_r**(1.0/2.0)
        phi = 2.0 * math.pi * Rs2
        pos[0][i] = r[i] * math.cos( phi ) ; pos[1][i] = r[i] * math.sin( phi )

    # Compute temperature structure

        t = (ps.t_0**2.0 * ( (r[i]**2.0 + r_in**2.0) * \
          cs.au_per_cm**2.0 )**(-ps.q_d) + ps.t_inf**2.0 ) ** (0.5)
        u = ( 1.5 * cs.k_b * t ) / (m_part * cs.kg_per_g)
        c_s2 = (cs.k_b * t) / (cs.mu * cs.m_p)

    # Compute disc thickness from surface density/temperature to find z location

        SD_0=( (ps.m_d*(2.0-ps.p_d)) / (2.0*math.pi*r_0**2.0) ) * \
           ( ((r_0**2.0 + r_out**2.0)/ r_0**2.0)**(ps.p_d/2.0) - \
           ((r_0**2.0 + r_in**2.0)/ r_0**2.0)**(ps.p_d/2.0) )**(-1.0)
        SD_R = SD_0 * ( r_0**2.0 / (r_0**2.0 + r[i]**2.0) )**(ps.p_d/2.0)

        z_0 = -((math.pi * SD_R * (r[i]**3.0)) / (2.0 * m_star)) + \
          ( ( (math.pi * SD_R * (r[i]**3.0) ) / (2.0 * m_star) )**2.0 + \
          ( (r[i]**3.0) / (cs.G_grav * (cs.cm_per_m**3.0) \
          * (cs.kg_per_g) * m_star) ) * c_s2 )**(0.5)
        pos[2][i] = z_0 * (2.0/math.pi) * math.asin( 2.0 * Rs3 - 1.0)

        rho[i] = ((math.pi * SD_0) / (4.0 * z_0 )) * \
           ( r_in**2.0 / (r_in**2.0 + r[i]**2.0 ))**(ps.p_d/2.0) * \
           math.cos((math.pi/2.0) * (pos[2][i] / z_0 ))
        h[i] = ( (3.0 * N_neigh * m_part) / (32.0 * math.pi * rho[i]))**(0.333)

    m_part = m_part * cs.dust_to_gas

    # Allocate lists to inp_model object

    inp_model.pos = pos
    inp_model.rho = rho
    inp_model.h = h
    inp_model.m_part = m_part
    inp_model.r = r

    return inp_model

### ------------------------------------------------------------------------ ###

def envelope(N_neigh):

    print("Creating envelope\n")

    pos = [ [[0.0] for i in range(ps.nsph)] for i in range(3)]
    r = [0.0] * ps.nsph ; rho = [0.0] * ps.nsph ; h = [0.0] * ps.nsph

    r_in = ps.r_in_e * cs.cm_per_au
    r_out = ps.r_out_e * cs.cm_per_au
    m_e = ps.m_e * cs.g_per_msol * cs.dust_to_gas
    m_part = m_e / float(ps.nsph)

    exp = 1.0 - (float(ps.p_e) / 3.0)
    w_env = r_out**(3.0) / r_in**(3.0)

    rho_0 = ( (m_e  * (3.0 * exp) ) / (4.0 * math.pi * r_in**(3.0)) ) * \
       (w_env**(exp) - 1.0)**(-1.0)

    for i in range(ps.nsph):

        Rs1 = random.random() ; Rs2 = random.random() ; Rs3 = random.random()

        w_r = ( ( ( w_env**exp - 1.0) * Rs1 ) + 1.0 )**(1.0/exp)

        r[i]= r_in * ( (w_r)**(1.0/3.0) )
        theta = math.acos(1.0 - (2.0 * Rs2) )
        phi = 2.0 * math.pi * Rs3

        pos[0][i] = ( r[i] * (math.sin(theta) * math.cos( phi ) ) )
        pos[1][i] = ( r[i] * ( math.sin(theta) * math.sin( phi ) ) )
        pos[2][i] = ( r[i] * math.cos(theta))

        rho[i] = rho_0 * (r_in / r[i])**(float(ps.p_e))
        h[i] = ( (3.0 * N_neigh * m_part) / (32.0 * math.pi * rho[i]))**(0.333)


    # Allocate lists to inp_model object

    inp_model.pos = pos
    inp_model.rho = rho
    inp_model.h = h
    inp_model.m_part = m_part
    inp_model.r = r
    inp_model.extras.append(rho_0)

    return inp_model

### ------------------------------------------------------------------------ ###

def plum_sphere(N_neigh):

    print("Creating Plummer sphere\n")

    pos = [ [[0.0] for i in range(ps.nsph)] for i in range(3)]
    r = [0.0] * ps.nsph ; rho = [0.0] * ps.nsph ; h = [0.0] * ps.nsph

    m_e = ps.m_e * cs.g_per_msol * cs.dust_to_gas
    m_part = m_e / ps.nsph
    a_plum = ps.a_plum * cs.cm_per_au
    r_env = ps.a_plum * ps.n_aplum

    for i in range(ps.nsph):

        Rs1 = random.random() ; Rs2 = random.random() ; Rs3 = random.random()

        alpha = ( r_env**(3.0) / ( r_env**(2.0) + a_plum**(2.0) )**(1.5) ) \
            * Rs1

        r[i]= ( a_plum * alpha**(0.333) ) / ( 1.0 - alpha**(0.666) )**(0.5)
        theta = math.acos(1.0 - (2.0 * Rs2) )
        phi = 2.0 * math.pi * Rs3

        pos[0][i] = ( r[i] * (math.sin(theta) * math.cos( phi ) ) )
        pos[1][i] = ( r[i] * ( math.sin(theta) * math.sin( phi ) ) )
        pos[2][i] = ( r[i] * math.cos(theta))

        rho[i] = ( (3.0 * m_e) / (4.0 * math.pi * a_plum**(3.0) ) ) / \
           ( 1.0 + (r[i] / a_plum)**(2.0) )**(2.5)
        h[i] = ( (3.0 * N_neigh * m_part) / (32.0 * math.pi * rho[i]))**(0.333)


    # Allocate lists to inp_model object

    inp_model.pos = pos
    inp_model.rho = rho
    inp_model.h = h
    inp_model.m_part = m_part
    inp_model.r = r

    return inp_model
