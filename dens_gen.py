'''

dens_gen.py

Module to generate/read the density structure of an SPH model/simulation. For
simulations, use SPLASH to convert cdisc file to ascii format and read. For
models, generate SPH distribution. With respective SPH distributions, write
into arrays for SPH particle -> AMR grid translation.

Density distribution is for dust mass as per requirements of RADMC-3D.

Last Modified: 31/01/2018

'''

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - MODULE IMPORTS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

import os
import numpy as np
import math
import random

import constants as cs
from params_run import *

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - MAIN PROGRAM - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

# Function that reads in simulation, and outputs SPH gas data to arrays for
# further analyses in RT modelling.

def sim(arch_dir):

    os.system("{0} to ascii {1}".format(xsplash, sim_cfile) )
    inp_file = sim_cfile+".ascii"

    pos = [[] for i in range(3)] ; h = [] ; rho = [] ; r = []
    f = open(inp_file, "r")
    for i in range(12):
        trash = f.readline()
    for lines in f:
        line = lines.strip() ; columns = line.split()

        if (float(columns[10]) == 1):
            pos[0].append(float(columns[0]) * cs.cm_per_pc)
            pos[1].append(float(columns[1]) * cs.cm_per_pc)
            pos[2].append(float(columns[2]) * cs.cm_per_pc)
            m_part = float(columns[9]) * cs.g_per_msol * cs.dust_to_gas
            h.append(float(columns[7]) * cs.cm_per_pc)
            rho.append(float(columns[8]))
            r.append( np.sqrt(float(columns[0])**2.0 + \
               float(columns[1])**2.0 + float(columns[2])**2.0  ) )
        else:
            continue
    f.close()
    os.remove(inp_file)

    extras = []
    return pos, rho, h, m_part, r, extras

### ------------------------------------------------------------------------ ###

    # Main modelling function, which calls mod_inp dependent functions.

def model():

    N_neigh = 50.0         # Float: No. of neighbours density is smoothed over

    if (mod_inp == "disc"):
        pos, rho, h, m_part, r, extras = disc(N_neigh)
    elif (mod_inp == "env"):
        pos, rho, h, m_part, r, extras = envelope(N_neigh)
    elif (mod_inp == "plum"):
        pos, rho, h, m_part, r, extras = plum_sphere(N_neigh)
    else:
        print("Incorrect mod_inp selection made, exiting program\n")
        exit()

    return pos, rho, h, m_part, r, extras

### ------------------------------------------------------------------------ ###

def disc(N_neigh):

    global m_star, m_d, r_0_d, r_in_d, r_out_d

    print("Creating disc\n")

    pos = [ [[0.0] for i in range(nsph)] for i in range(3)]
    r = [0.0] * nsph ; rho = [0.0] * nsph ; h = [0.0] * nsph

    # Compute particle *gas* mass, then convert mass/distance values to c.g.s

    m_part = (m_d * cs.g_per_msol) / float(nsph)
    m_d = m_d * cs.g_per_msol ; m_star = m_star * cs.g_per_msol
    r_0 = r_0_d * cs.cm_per_au
    r_in = r_in_d * cs.cm_per_au ; r_out = r_out_d * cs.cm_per_au

    # Define w values from which radial profiles can be attained

    w_in=(r_in**2.0) / (r_out**2.0) ; w_out=(r_out**2.0) / (r_0**2.0)

    for i in range(nsph):

        Rs1 = random.random() ; Rs2 = random.random() ; Rs3 = random.random()

        w_r = ( (1.0 + w_in)**(1.0 - (p_d / 2.0) ) + \
           Rs1 * ( (1.0 + w_out)**(1.0 - (p_d / 2.0) ) - \
           (1.0 + w_in)**(1.0 - (p_d / 2.0) ) )) ** ( 2.0 / (2.0 - p_d) ) - 1.0

        r[i] = r_0 * w_r**(1.0/2.0)
        phi= 2.0 * math.pi * Rs2
        pos[0][i] = r[i] * math.cos( phi ) ; pos[1][i] = r[i] * math.sin( phi )

    # Compute temperature structure

        t = (t_0**2.0 * ( (r[i]**2.0 + r_in**2.0) * \
          cs.au_per_cm**2.0 )**(-q_d) + t_inf**2.0 ) ** (0.5)
        u = ( 1.5 * cs.k_b * t ) / (m_part * cs.kg_per_g)
        c_s2 = (cs.k_b * t) / (cs.mu * cs.m_p)

    # Compute disc thickness from surface density/temperature to find z location

        SD_0=( (m_d*(2.0-p_d)) / (2.0*math.pi*r_0**2.0) ) * \
           ( ((r_0**2.0 + r_out**2.0)/ r_0**2.0)**(p_d/2.0) - \
           ((r_0**2.0 + r_in**2.0)/ r_0**2.0)**(p_d/2.0) )**(-1.0)
        SD_R = SD_0 * ( r_0**2.0 / (r_0**2.0 + r[i]**2.0) )**(p_d/2.0)

        z_0 = -((math.pi * SD_R * (r[i]**3.0)) / (2.0 * m_star)) + \
          ( ( (math.pi * SD_R * (r[i]**3.0) ) / (2.0 * m_star) )**2.0 + \
          ( (r[i]**3.0) / (cs.G_grav * (cs.cm_per_m**3.0) \
          * (cs.kg_per_g) * m_star) ) * c_s2 )**(0.5)
        pos[2][i] = z_0 * (2.0/math.pi) * math.asin( 2.0 * Rs3 - 1.0)

        rho[i] = ((math.pi * SD_0) / (4.0 * z_0 )) * \
           ( r_in**2.0 / (r_in**2.0 + r[i]**2.0 ))**(p_d/2.0) * \
           math.cos((math.pi/2.0) * (pos[2][i] / z_0 ))
        h[i] = ( (3.0 * N_neigh * m_part) / (32.0 * math.pi * rho[i]))**(0.333)

    m_part = m_part * cs.dust_to_gas

    extras = []
    return pos, rho, h, m_part, r, extras

### ------------------------------------------------------------------------ ###

def envelope(N_neigh):

    global r_in_e, r_out_e, m_e

    print("Creating envelope\n")

    pos = [ [[0.0] for i in range(nsph)] for i in range(3)]
    r = [0.0] * nsph ; rho = [0.0] * nsph ; h = [0.0] * nsph

    r_in = r_in_e * cs.cm_per_au ; r_out = r_out_e * cs.cm_per_au
    m_e = m_e * cs.g_per_msol * cs.dust_to_gas
    m_part = m_e / float(nsph)

    exp = 1.0 - (float(p_e) / 3.0)
    w_env = r_out**(3.0) / r_in**(3.0)

    rho_0 = ( (m_e  * (3.0 * exp) ) / (4.0 * math.pi * r_in**(3.0)) ) * \
       (w_env**(exp) - 1.0)**(-1.0)

    for i in range(nsph):

        Rs1 = random.random() ; Rs2 = random.random() ; Rs3 = random.random()

        w_r = ( ( ( w_env**exp - 1.0) * Rs1 ) + 1.0 )**(1.0/exp)

        r[i]= r_in * ( (w_r)**(1.0/3.0) )
        theta = math.acos(1.0 - (2.0 * Rs2) )
        phi = 2.0 * math.pi * Rs3

        pos[0][i] = ( r[i] * (math.sin(theta) * math.cos( phi ) ) )
        pos[1][i] = ( r[i] * ( math.sin(theta) * math.sin( phi ) ) )
        pos[2][i] = ( r[i] * math.cos(theta))

        rho[i] = rho_0 * (r_in / r[i])**(float(p_e))
        h[i] = ( (3.0 * N_neigh * m_part) / (32.0 * math.pi * rho[i]))**(0.333)

    extras = []
    extras.append(rho_0)
    extras = np.array(extras)

    return pos, rho, h, m_part, r, extras

### ------------------------------------------------------------------------ ###

def plum_sphere(N_neigh):

    global m_e, a_plum, n_aplum

    print("Creating Plummer sphere\n")

    pos = [ [[0.0] for i in range(nsph)] for i in range(3)]
    r = [0.0] * nsph ; rho = [0.0] * nsph ; h = [0.0] * nsph

    m_e = m_e * cs.g_per_msol * cs.dust_to_gas ; m_part = m_e / nsph
    a_plum = a_plum * cs.cm_per_au ; r_env = a_plum * n_aplum

    for i in range(nsph):

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

    extras = []
    return pos, rho, h, m_part, r, extras
