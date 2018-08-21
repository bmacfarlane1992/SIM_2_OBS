'''

make.py

Python program to either generate models/read DRAGON or SEREN data, then
convert SPH particle -> AMR grid density distribution.
Code then generates RADMC-3D ICs, and runs MCRT and RRT to compute/generate

    - Temperature distribution (MCRT)
	- Spectral Energy Distribtions (RRT)
    - Diagnostics (Radial/Slice density/temperature plots, and Continuum
       images, using RRT and radmc3dPy modules)

Author: Benjamin MacFarlane
Date: 20/08/2018
Contact: bmacfarlane@uclan.ac.uk

'''

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
		# # # - - - MODULE IMPORTS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

# System module imports

import os
import imp
import shutil
import sys

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
		# # # - - - MAIN PROGRAM - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

    # Intitialise objects to be appended to within SIM_2_OBS

class globals:

    arch_dir = os.getcwd()+"/.."
    run_dir = ""
    dat_dir = ""
    plt_dir = ""

    plt_ext = ""
    tdat_ext = ""
    pdat_ext = ""
    sdat_ext = ""

class rt_params:

    t_star = 0.0
    r_star = 0.0
    l_star = 0.0

    t_planet = 0.0
    r_planet = 0.0

    wavs = []
    n_wavs = 0

class inp_model:

    pos = [[] for i in range(3)]
    rho = []
    h = []
    m_part = []
    r = []
    extras = []

class rt_model:

    xloc = []
    yloc = []
    zloc = []
    rloc = []
    rloc_xy = []
    xwid = []
    rho = []

### ------------------------------------------------------------------------ ###

    # Find parameters files for batch SIM_2_OBS analyses in /params/pending, or
    # simply adopt parameters.py in /params. If multiple params used, ensure
    # completed runs are transported to params/done directory

param_dir = globals.arch_dir + "/params/pending"
params_done_dir = globals.arch_dir + "/params/done"
f_params = [] ; f_ext = []
if os.path.isdir(param_dir):
    if not (os.path.isdir(params_done_dir)):
        os.makedirs(params_done_dir)
    for root, dirs, files in os.walk(param_dir):

    # Append to f_params list

        for f in files:
            if (f.split(".")[1] == "py"):
                f_params.append(param_dir+"/"+f)
                f_ext.append(f)
            else:
                continue

if (len(f_params) == 0):
    f_params.append(globals.arch_dir+"/params/parameters.py")

for files in range(len(f_params)):
    shutil.copy2( f_params[files], os.getcwd()+"/params_run.py" )
    print("Parameters being used from {0}\n".format(f_params[files]) )

    # Now /src/ modules, to ensure that "from params_run import *" can be
    # used at modular level

    import params_run as ps
    import inp_gen
    import dens_gen
    import grid_gen
    import opac_gen
    import ec53_rt
    import scale_cut
    import multi_beam
    import inp_plt
    import obs_plt
    import raytrace
    import diagnose

### ------------------------------------------------------------------------ ###

    # Define, then generate run/dat/plt directoried, where data/plots are stored

    globals.run_dir = globals.arch_dir + "/runs/" + ps.run_tag
    globals.dat_dir = globals.run_dir+"/dat"
    globals.plt_dir = globals.run_dir+"/plots"

    if not os.path.isdir(globals.run_dir):
        os.makedirs(globals.run_dir)
    if not os.path.isdir(globals.dat_dir):
        os.makedirs(globals.dat_dir)
    if not os.path.isdir(globals.plt_dir):
        os.makedirs(globals.plt_dir)

### ------------------------------------------------------------------------ ###

	# Make most .inp files. If source SEDs not present in /protostar/, generate

    globals, rt_params = inp_gen.gen(globals, rt_params)

### ------------------------------------------------------------------------ ###

    # If called for, make AMR grid from constructed SPH particle distribution
    # (simulation or model) with HYPERION, then generate opacity tables

    if ps.RT_RUN:

        if ps.GRID and (ps.mod_inp != "yso"):

            if ps.sim_inp:
                inp_model = dens_gen.sim()

            elif not ps.sim_inp:
                inp_model = dens_gen.model()

            grid_gen.translate(globals, inp_model)

        elif ps.GRID and (ps.mod_inp == "yso"):
            grid_gen.yso(globals)

        if (ps.run_tag != "EC53"):
            opac_gen.gen(globals, rt_params)


### ------------------------------------------------------------------------ ###

    # Apply cuts and/or mass rescaling for full/components of simulation

    if ps.sim_inp:

        if (ps.run_tag == "EC53") and ps.RT_RUN and ps.MCRT:
            ec53_rt.gen(globals)

        if ps.FULL_mass_scale and not (ps.DISC_mass_scale or ps.ENV_mass_scale):
            scale_cut.pre_full_scale(globals)

        if (ps.DISC_mass_scale and ps.ENV_mass_scale):
            scale_cut.pre_comp_scale(globals)

        if (ps.DISC_mass_mult and ps.ENV_mass_mult):
            scale_cut.pre_comp_mult(globals)

        if ps.MCRT_cut:
            scale_cut.pre_cut(globals)

        if ps.MCRT_cut_scale:
            scale_cut.pre_cut_scale(globals)

        if ps.CAVITY_rho:
            scale_cut.cavity_gen(globals)

        if ps.ENV_rho:
            scale_cut.envelope_gen(globals)

        if ps.R_DESTR:
            ps.R_DESTR = scale_cut.inner_cavity(globals, rt_params)

### ------------------------------------------------------------------------ ###

    # Run MCRT to generate temperature distribution

    if ps.RT_RUN and ps.MCRT:

        os.chdir(globals.dat_dir+"/")
        os.system("{0}radmc3d mctherm".format(ps.exec_loc) )

        shutil.copy2( globals.dat_dir+"/dust_temperature.dat", \
          globals.dat_dir+"/dust_temperature"+globals.tdat_ext )
        shutil.copy2( globals.dat_dir+"/photon_statistics.out", \
          globals.dat_dir+"/photon_statistics"+globals.pdat_ext )

### ------------------------------------------------------------------------ ###

    # If model is spatially restricted or mass scaled before RRT, do here

    if ps.RRT_cut:
        scale_cut.post_cut(globals)

### ------------------------------------------------------------------------ ###

    # Run RRT to generate SED. Also move .dat/.inp files dependent on which
    # radiation sources are included in raytracing

    if ps.RT_RUN and ps.RRT:

        os.chdir(globals.dat_dir+"/")
        print os.getcwd(), globals.dat_dir, globals.tdat_ext
        shutil.copy2(globals.dat_dir+"/dust_temperature"+globals.tdat_ext, \
         globals.dat_dir+"/dust_temperature.dat" )

        if not os.path.isdir(globals.dat_dir+"/isrf/"):
            os.makedirs(globals.dat_dir+"/isrf")

        if ps.incl_isrf and not ps.incl_star:
            shutil.move(globals.dat_dir+"/external_source.inp", \
             globals.dat_dir+"/isrf/external_source.inp")
            shutil.move(globals.dat_dir+"/bin/stars.inp", \
             globals.dat_dir+"/stars.inp")
        elif ps.incl_isrf and ps.incl_star:
            shutil.move(globals.dat_dir+"/external_source.inp", \
             globals.dat_dir+"/isrf/external_source.inp")

        for i in range(len(ps.inclins)):

            if ps.incl_star:
                os.system("{0}radmc3d sed incl {1}".format(ps.exec_loc, \
                 ps.inclins[i]) )
            elif not ps.incl_star:
                os.system("{0}radmc3d sed incl {1} nostar".format(ps.exec_loc, \
                 ps.inclins[i]) )

            globals.sdat_ext = globals.plt_ext+"_"+str(ps.inclins[i])+"i.out"

            shutil.copy2(globals.dat_dir+"/spectrum.out", \
             globals.dat_dir+"/spectrum"+globals.sdat_ext)

            if (ps.mod_inp == "protostar"):
                shutil.copy2(globals.dat_dir+"/spectrum"+int(rt_params.t_star)+".out", \
                  globals.arch_dir+"/protostar/spectrum"+int(rt_params.t_star)+".out")

### ------------------------------------------------------------------------ ###

    # Run RRT to generate SED, for multipled RRT_cut values. Used to quickly
    # evaluate impact of changing effective beam size on flux

    if ps.MULTI_RRT:

        multi_beam.run(globals)

### ------------------------------------------------------------------------ ###

    # Generate diagnostic plots, and SED. Run RRT for images if called for. Also
    # generate location dependent column density/optical depth properties

    if not ps.RT_RUN:

        os.chdir(globals.dat_dir)

        if ( ps.plot_radial or ps.plot_slice or ps.plot_sed):
            xloc, yloc, zloc, rloc, xwid, rloc_xy, nrspecies, rho, diag_rho, \
              temp, diag_temp, nphot, diag_nphot, diag_sed, beam_r \
              = inp_plt.read(globals)

        if ps.plot_radial and ( diag_rho or diag_temp ):
            inp_plt.radial(globals, rloc, rloc_xy, \
              nrspecies, rho, diag_rho, R_DESTR, temp, diag_temp)

        if ps.plot_slice and (diag_rho or diag_temp or diag_nphot):
            inp_plt.run_slice(globals, xloc, yloc, zloc, xwid, \
              nrspecies, rho, diag_rho, temp, diag_temp, nphot, diag_nphot)

        if ps.plot_sed and diag_sed:

            for i in range(len(inclins)):

                globals.sdat_ext = globals.plt_ext+"_"+str(ps.inclins[i])+"i.out"

                obs_plt.sed(globals, ps.inclins[i], rt_model, beam_r)

        if ps.plot_image:
            for i in range(len(ps.inclins)):
                obs_plt.image(globals, ps.inclins[i])

        if ps.plot_rays:

            # First, generate column density map

            raytrace.column(globals)

            # Then generate locations at which optical depth reaches prescribed
            # value, to diagnose as a function of observer location

            raytrace.tau3d(globals)

            # Now compute, for specified observer what the cloud emission
            # is as a function of optical depth

            raytrace.tau1d(globals)

        if mass_eval:

            diagnose.mass(globals)

        if lum_eval:

            diagnose.lum(globals)


### ------------------------------------------------------------------------ ###

    # Remove modules from system list, so new params_run file can be loaded

    # Also transfer params file to params_done_dir, to track which runs complete

    os.chdir("{0}/src/".format(globals.arch_dir))
    del sys.modules["params_run"]
    os.remove(os.getcwd()+"/params_run.py")
    del sys.modules["inp_gen"]
    del sys.modules["dens_gen"]
    del sys.modules["grid_gen"]
    del sys.modules["opac_gen"]
    del sys.modules["ec53_rt"]
    del sys.modules["scale_cut"]
    del sys.modules["raytrace"]
    del sys.modules["inp_plt"]
    del sys.modules["obs_plt"]
    del sys.modules["diagnose"]
    print("\n\n\n")

    if (len(f_params) > 1):
        shutil.move( f_params[files], \
          params_done_dir+"/"+f_ext[files])

### ------------------------------------------------------------------------ ###

    # Clean directories of garbage, incuding .pyc, ~, and SPLASH .ascii files

if os.path.isfile(globals.arch_dir+"/src/radmc3d.out"):
    os.remove(globals.arch_dir+"/src/radmc3d.out")

dirs = [globals.arch_dir+"/ics/", globals.arch_dir+"/src/", globals.arch_dir+"/"]
exts = [".pyc", ".ascii"]
for d in range(len(dirs)):
    dir = dirs[d]
    files = os.listdir(dir)
    for file in files:
        for e in range(len(exts)):
            if file.endswith(exts[e]):
                os.remove(os.path.join(dir,file))

exit()
