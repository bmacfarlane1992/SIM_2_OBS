'''

functions.py

Common functions to be called by any Python procedures.

Author: Benjamin MacFarlane
Date: 20/03/2018
Contact: bmacfarlane@uclan.ac.uk

'''

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - MODULE IMPORTS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

import numpy as np
import math

import constants as cs

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - MAIN PROGRAM - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

### ------------------------------------------------------------------------ ###
    # Function that numerically integrates a function between limits start #
    # and end, by splitting into n_increments and using trapezium          #
    # run_slice to calculate the area under the functions                  #
### ------------------------------------------------------------------------ ###


def Trapezoidal(function, start, end, n_increments):

    width = (end - start) / float(n_increments)
    height = 0.5 * (function(start) + function(end))
    for i in range(1, n_increments, 1):
        height = height + function(start + i * width)

    return width * height


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

### ------------------------------------------------------------------------ ###
    #                               #
### ------------------------------------------------------------------------ ###


def adjustFigAspect(fig,aspect=1):
    '''
    Adjust the subplot parameters so that the figure has the correct
    aspect ratio.
    '''
    xsize,ysize = fig.get_size_inches()
    minsize = min(xsize,ysize)
    xlim = .4*minsize/xsize
    ylim = .4*minsize/ysize
    if aspect < 1:
        xlim *= aspect
    else:
        ylim /= aspect
    fig.subplots_adjust(left=.5-xlim,
                        right=.5+xlim,
                        bottom=.5-ylim,
                        top=.5+ylim)
