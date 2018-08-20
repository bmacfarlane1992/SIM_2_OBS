'''

functions.py

Common functions to be called by any Python procedures.

Current list of functions and brief outline:

    - Trapezoidal() - Integrate function numerically, using trapezium rule

    - Quad_int() - Intgrate function numerically, using Fortra lib QUADPACK

    - SED_read() - Read SED data, calculating fluxes in terms of lambda, and nu

    - Lbol_calc() - Use integrated SED data, to compute bolometric luminosity

    - L_ratio_calc() - Compute sub-mm to bolometric luminosity ratio

    - Tbol_calc() - Use integrated SED data, to compute bolometric temperature

    - Column_Density() - Estimate column density of line-of-sight to protostar

Author: Benjamin MacFarlane
Date: 20/03/2018
Contact: bmacfarlane@uclan.ac.uk

'''

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - MODULE IMPORTS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

import numpy as np
import math
from scipy.interpolate import interp1d
import scipy.integrate as integrate

import constants as cs
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - VARIABLE DEFINITIONS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
int_scheme = "quad"     # Integration scheme: ["trapezium","quad"]
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - MAIN PROGRAM - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

### ------------------------------------------------------------------------ ###
    # Function that numerically integrates a function between limits start #
    # and end, by splitting into n_increments and using trapezium          #
    # run_slice to calculate the area under the functions                  #
### ------------------------------------------------------------------------ ###


def Trapezoidal(x, y, start, end):

    interp = interp1d(x, y, kind="cubic")
    function = lambda x: interp(x)

    n_increments = int(1e5)

    width = (end - start) / float(n_increments)
    height = 0.5 * (function(start) + function(end))
    for i in range(1, n_increments, 1):
        height = height + function(start + i * width)

    return width * height



### ------------------------------------------------------------------------ ###
    # Function that numerically integrates a function between limits start #
    # and end, by adopting Fortran QUADPACK integration methods            #
### ------------------------------------------------------------------------ ###


def Quad_int(x, y, start, end):

    interp = interp1d(x, y, kind="cubic")
    function = lambda x: interp(x)

    soln = integrate.quad(function, start, end, limit=200)

    return soln[0]


### ------------------------------------------------------------------------ ###
    # Function that read an SED file, and return all data as required for  #
    # SED plots/analysis.                                                  #
    # Note: All fluxes are normalised to 1 pc as per RADMC-3D              #
### ------------------------------------------------------------------------ ###


def SED_read(file):

    # Initialise lists to be appended to

    wav = [] ; flam = [] ; lam_flam = [] ; nu = [] ; fnu = [] ; nu_fnu = []

    # Read file, and append to relevant lists - deriving values are needed

    f = open(file, "r")
    for j in range(3):
        header = f.readline()
    for lines in f:
        lines = lines.strip() ; columns = lines.split()
        wav.append(float(columns[0]))
        flam.append(float(columns[1]) * \
           ( cs.c_cgs / (float(columns[0]) * cs.cm_per_micron)**2.0 ) )
        lam_flam.append(float(columns[1]) * \
           ( cs.c_cgs / (float(columns[0]) * cs.cm_per_micron) ) )
        nu.append(cs.c_cgs / (float(columns[0]) * cs.cm_per_micron) )
        fnu.append(float(columns[1]))
        nu_fnu.append( (cs.c_cgs / (float(columns[0]) * cs.cm_per_micron) ) * \
         float(columns[1]) )
    f.close()

    # Reverse frequency based lists, to ensure correct SED integration

    nu = nu[::-1] ; fnu = fnu[::-1] ; nu_fnu = nu_fnu[::-1]

    wav = np.array(wav) ; flam = np.array(flam) ; lam_flam = np.array(lam_flam)
    nu = np.array(nu) ; fnu = np.array(fnu) ; nu_fnu = np.array(nu_fnu)

    return wav, flam, lam_flam, nu, fnu, nu_fnu

### ------------------------------------------------------------------------ ###
    # Function that uses SED data, to integrate and compute L_bol          #
### ------------------------------------------------------------------------ ###


def Lbol_calc(nu, fnu):

    min_nu = min(nu) ; max_nu = max(nu)

    if (int_scheme == "trapezium"):
        result = Trapezoidal(nu, fnu, min_nu, max_nu)
    elif (int_scheme == "quad"):
        result = Quad_int(nu, fnu, min_nu, max_nu)
    else:
        print("Incorrect integration scheme selected in functions.py")
        exit()

    L_bol = (4.0 * math.pi * (cs.cm_per_pc**2.0) * result) / cs.Lsol_cgs

    return L_bol


### ------------------------------------------------------------------------ ###
    # Function that computes ratio of sub-mm to bolometric luminosity      #
### ------------------------------------------------------------------------ ###


def L_ratio_calc(nu, fnu):

    L_bol = Lbol_calc(nu, fnu)

    min_nu = min(nu) ; max_nu = cs.c / 350.e-6

    if (int_scheme == "trapezium"):
        result = Trapezoidal(nu, fnu, min_nu, max_nu)
    elif (int_scheme == "quad"):
        result = Quad_int(nu, fnu, min_nu, max_nu)
    else:
        print("Incorrect integration scheme selected in functions.py")
        exit()

    L_submm = (4.0 * math.pi * (cs.cm_per_pc**2.0) * result) / cs.Lsol_cgs

    L_ratio = (L_submm / L_bol) * 100.0

    return L_ratio

### ------------------------------------------------------------------------ ###
    # Function that uses SED data, to integrate and compute T_bol          #
### ------------------------------------------------------------------------ ###


def Tbol_calc(nu, fnu, nu_fnu):

    min_nu= min(nu) ; max_nu = max(nu)

    if (int_scheme == "trapezium"):
        result1 = Trapezoidal(nu, fnu, min_nu, max_nu)
        result2 = Trapezoidal(nu, nu_fnu, min_nu, max_nu)
    elif (int_scheme == "quad"):
        result1 = Quad_int(nu, fnu, min_nu, max_nu)
        result2 = Quad_int(nu, nu_fnu, min_nu, max_nu)
    else:
        print("Incorrect integration scheme selected in functions.py")

    T_bol = 1.25e-11 * (result2 / result1)

    return T_bol

### ------------------------------------------------------------------------ ###
    # Calculate the Planck function variable B_nu, using input of          #
    # wavelength (m) and temperature (K)                                   #
### ------------------------------------------------------------------------ ###

def Planck(wavelength, temperature):

    nu = cs.c / wavelength

    a = ( 2.0 * cs.h * nu**(3.0) ) / ( cs.c**(2.0) )

    val = ( cs.h * nu ) / ( cs.k_b * temperature )
    try:
        exp = math.exp( val )
        b = 1.0 / (exp - 1.0)
    except OverflowError:
        b = 0.0

    B_nu = a * b

    return B_nu

### ------------------------------------------------------------------------ ###
    # For an AMR grid centred on [0,0,0], compute the column density       #
    # for a ray from [0,0,0] to the observer at infinity, parameterised    #
    # by the inclination and azimuthal angle                               #
### ------------------------------------------------------------------------ ###

def Column_Density(inclin, phi, x, y, z, xw, rho):

    eps = 1.e-1 # math.pi / 90

    # First, identify quadrant that is valid for inclin/phi pair

    quad = 4 * math.floor(inclin / (math.pi / 2.0) )
    quad += math.floor(phi / (math.pi / 2.0) )
    quad = int(quad)

    # Refine AMR grid arrays to only this quadrant, then redefine arrays/angles
    # to positive x-y-z quadrant for quadrant-independent analysis

    if (quad == 0):
        refined = np.where( (z > 0.0) & (x > 0.0) & (y > 0.0) )
        phi_offset = 0
    elif (quad == 1):
        refined = np.where( (z > 0.0) & (x < 0.0) & (y > 0.0) )
        phi_offset = 1
    elif (quad == 2):
        refined = np.where( (z > 0.0) & (x < 0.0) & (y < 0.0) )
        phi_offset = 2
    elif (quad == 3):
        refined = np.where( (z > 0.0) & (x > 0.0) & (y < 0.0) )
        phi_offset = 3
    elif (quad == 4):
        refined = np.where( (z < 0.0) & (x > 0.0) & (y > 0.0) )
        phi_offset = 0
    elif (quad == 5):
        refined = np.where( (z < 0.0) & (x < 0.0) & (y > 0.0) )
        phi_offset = 1
    elif (quad == 6):
        refined = np.where( (z < 0.0) & (x < 0.0) & (y < 0.0) )
        phi_offset = 2
    elif (quad == 7):
        refined = np.where( (z < 0.0) & (x > 0.0) & (y < 0.0) )
        phi_offset = 3

    x = x[refined] ; y = y[refined] ; z = z[refined]
    xw = xw[refined] ; rho = rho[refined]
    r_xy = np.sqrt(x**2.0 + y**2.0)
    r = np.sqrt(x**2.0 + y**2.0 + z**2.0)

    # Further refine by bins that only the rays pass through for inclin and phi

    phi_bin = [] ; inclin_bin = []
    d_phi1 = [] ; d_phi2 = [] ; d_inclin1 = [] ; d_inclin2 = []

    for b in range(len(x)):

        # For azimuthal angle

        if (y[b] >= 0.0):
            phi_bin.append( math.atan2(y[b], x[b]) )
        elif (y[b] < 0.0):
            phi_bin.append( (2.0 * math.pi) - math.atan2(abs(y[b]), x[b]) )

        r_e1 = np.sqrt( (abs(x[b]) + xw[b])**2.0 + (abs(y[b]) - xw[b])**2.0 )
        frac = ( r_e1**2.0 - 2.0 * xw[b]**2.0 + r_xy[b]**2.0 ) / \
          (2.0 * r_e1 * r_xy[b] )
        d_phi1.append( math.acos(frac))

        r_e2 = np.sqrt( (abs(x[b]) - xw[b])**2.0 + (abs(y[b]) + xw[b])**2.0 )
        frac = ( r_e2**2.0 - 2.0 * xw[b]**2.0 + r_xy[b]**2.0 ) / \
          (2.0 * r_e2 * r_xy[b] )
        d_phi2.append( math.acos(frac) )

        # Now inclination

        inclin_bin.append( math.atan2(r_xy[b], z[b]) )

        r_e1 = np.sqrt( (r_xy[b] - xw[b])**2.0 + (abs(z[b]) + xw[b])**2.0 )
        frac = ( r_e1**2.0 - 2.0*xw[b]**2.0 + r[b]**2.0 ) / \
          ( 2.0 * r_e1 * r[b] )
        d_inclin1.append( math.acos(frac))

        r_e2 = np.sqrt( (r_xy[b] + xw[b])**2.0 + (abs(z[b]) - xw[b])**2.0 )
        frac = ( r_e2**2.0 - 2.0*xw[b]**2.0 + r[b]**2.0 ) / \
          (2.0 *  r_e2 * r[b] )
        d_inclin2.append( math.acos(frac) )

    phi_lim1 = np.array(phi_bin) - np.array(d_phi1)
    phi_lim2 = np.array(phi_bin) + np.array(d_phi2)

    inclin_lim1 = np.array(inclin_bin) - np.array(d_inclin1)
    inclin_lim2 = np.array(inclin_bin) + np.array(d_inclin2)

    refined = np.where( \
      ( phi_lim1 < (phi + eps) ) & ( phi_lim2 > (phi + eps) ) & \
      ( inclin_lim1 < (inclin + eps) ) & ( inclin_lim2 > (inclin + eps) ) )

    x = abs(x[refined]) ; y = abs(y[refined]) ; z = abs(z[refined])
    xw = xw[refined] ; r_xy = r_xy[refined] ; rho = rho[refined]

    if (quad > 3):
        inclin -= math.pi / 2.0
    phi -= ( (math.pi / 2.0) * phi_offset )

    # Now for refined array, compute where ray is at limiting vertices of bin

    x0 = x - xw ; y0 = y - xw ; r_xy0 = r_xy - xw ; z0 = z - xw
    x1 = x + xw ; y1 = y + xw ; r_xy1 = r_xy + xw ; z1 = z + xw

    xs = ( (math.tan(phi+eps))**(-1.0) * y0)
    xe = ((math.tan(phi+eps))**(-1.0) * y1)

    ys = ( math.tan(phi+eps) * x0)
    ye = ( math.tan(phi+eps) * x1)

    zs = ( (math.tan(inclin+eps))**(-1.0) * r_xy0)
    ze = ( (math.tan(inclin+eps))**(-1.0) * r_xy1)

    # Now compute the length over which the ray intersects the bin

    l = np.sqrt( (xe - xs)**2.0 + (ye - ys)**2.0 + (ze - zs)**2.0 )

    # Finally, compute the column Column_Density

    rho_column = np.sum(l * rho)

    return rho_column
