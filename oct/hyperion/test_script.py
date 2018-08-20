'''

test_script.py

Test script to import all functions for SPH - AMR translation and construct
octree over a small SPH test array

Author: Benjamin MacFarlane
Date: 09/02/2018

'''
# - - -
#
import os
import sys
import numpy as np
import random
#
from sph import construct_octree
#
# - - -
#
    # Set up octree world values, and SPH particle arrays
#
x = 0.0 ; y = 0.0 ; z = 0.0
dx = 20.0 ; dy = 20.0 ; dz = 20.0
#
px = [] ; py = [] ; pz = []
for i in range(100):
    px.append(20.0 * random.random() )
    py.append(20.0 * random.random() )
    pz.append(20.0 * random.random() )
px = np.array(px) ; py = np.array(py) ; pz = np.array(pz)
#
mpart = np.array([ 1.0 for i in range(100)])
sigma = np.array([ 2.0 for i in range(100)])
#
    # Set stopping criterion of octree based on max. number of particles in bin
#
def stop(x, y, z, dx, dy, dz, px, py, pz, sigma):
    return len(px) <= 2
#
    # Now run construct_octree function
#
o = construct_octree(x, y, z, dx, dy, dz, px, py, pz, sigma, mpart, \
  n_levels = 20, stopping_criterion=stop)
#
dir = os.getcwd()
files = os.listdir(dir)
for file in files:
    if file.endswith(".pyc"):
        os.remove(os.path.join(dir,file))
exit()
