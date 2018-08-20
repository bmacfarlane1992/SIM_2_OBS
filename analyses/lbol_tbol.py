'''

lbol_tbol.py

Diagnostic procedure to check, for a single SED file, the bolometric
luminosity

Author: Benjamin bmacfarlane
Date: 20/08/2018
Contact: bmacfarlane@uclan.ac.uk

'''

    # Import system modules

import math
import numpy as np
import os
from scipy.interpolate import interp1d

    # Append location of constants and functions modules into pythonpath

import sys
sys.path.insert(0,"./../")

    # Import local modules

import constants as cs
import functions as fs


# ------------------------------------------------------------------------------


arch_dir = os.getcwd()+"/../../runs"

run_dir = arch_dir+"/OUTBURST/1686/POST"

file_read = run_dir+"/dat/spectrum19248_ISRF_0i.out"

# ------------------------------------------------------------------------------

    # Read SED, using relevant function

wav, flam, lam_flam, nu, fnu, nu_fnu = fs.SED_read(file_read)

    # From SED data, call relevant function to derive bolometric properties

L_bol = fs.Lbol_calc(nu, fnu)

print("\nBolometric Luminosity is {0} L_sol\n".format(L_bol))

T_bol = fs.Tbol_calc(nu, fnu, nu_fnu)

print("\nBolometric Temperature is {0} K\n".format(T_bol))

os.remove("./../constants.pyc")
os.remove("./../functions.pyc")
