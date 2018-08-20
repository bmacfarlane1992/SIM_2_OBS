'''

lbol_tbol.py

Diagnostic procedure to check, for a single SED file, the bolometric
luminosity

Author: Benjamin bmacfarlane
Date: 11/07/2018
Contact: bmacfarlane@uclan.ac.uk

'''

import math
import numpy as np
import os
from scipy.interpolate import interp1d

import constants as cs


# ------------------------------------------------------------------------------


arch_dir = os.getcwd()+"/../runs"

run_dir = arch_dir+"/ENVELOPE/ENVELOPE_ISRF"

file_read = run_dir+"/dat/spectrumISRF_0i.out"


# ------------------------------------------------------------------------------


n = int(1e5)
def Trapezoidal(g, a, b, n):
    h = (b - a) / float(n)
    s = 0.5 * (g(a) + g(b))
    for i in range(1,n,1):
        s = s + g(a + i*h)
    return h*s

wavs = [] ; flam = []
f = open(file_read, "r")
for j in range(3):
    header = f.readline()
for lines in f:
    lines = lines.strip() ; columns = lines.split()
    wavs.append(float(columns[0]))
    flam.append(float(columns[1]) * \
      ( cs.c_cgs / (float(columns[0]) * cs.cm_per_micron)**2.0 ) )
f.close()

for i in range(len(wavs)):
    wavs[i] *= cs.cm_per_micron

flux_int = interp1d(wavs, flam, kind="cubic")
a = min(wavs) ; b = max(wavs)
def g(lam):
    return flux_int(lam)
result = Trapezoidal(g, a, b, n)
L_bol = (4.0 * math.pi * (cs.cm_per_pc**2.0) * result) / cs.Lsol_cgs

print("\nBolometric Luminosity is {0} L_sol\n".format(L_bol))
