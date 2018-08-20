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

arch_dir = os.getcwd()+"/.."

    # Find parameters files for batch SIM_2_OBS analyses in /params/pending, or
    # simply adopt parameters.py in /params. If multiple params used, ensure
    # completed runs are transported to params/done directory

param_dir = arch_dir + "/params/pending"
params_done_dir = arch_dir + "/params/done"
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
    f_params.append(arch_dir+"/params/parameters.py")

for files in range(len(f_params)):
    shutil.copy2( f_params[files], os.getcwd()+"/params_run.py" )
    print("Parameters being used from {0}\n".format(f_params[files]) )

    # Now /src/ modules, to ensure that "from params_run import *" can be
    # used at modular level

    from params_run import *
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

    # Generate run/dat/plt directory, in which data/plots are stored

    run_dir = arch_dir+"/runs/"+run_tag
    dat_dir = run_dir+"/dat" ; plt_dir = run_dir+"/plots"

    if not (os.path.isdir(run_dir)):
        os.makedirs(run_dir)
    if not (os.path.isdir(dat_dir)):
        os.makedirs(dat_dir)
    if not (os.path.isdir(plt_dir)):
        os.makedirs(plt_dir)

### ------------------------------------------------------------------------ ###

	# Make most .inp files. If source SEDs not present in /protostar/, generate

    plt_ext, tdat_ext, pdat_ext, t_star, r_star, t_planet, r_planet, \
      llambda, wav_bins = inp_gen.gen(arch_dir, dat_dir, plt_dir)

### ------------------------------------------------------------------------ ###

    # If called for, make AMR grid from constructed SPH particle distribution
    # (simulation or model) with HYPERION, then generate opacity tables

    if RT_RUN:

        if GRID and (mod_inp != "yso"):
            if sim_inp:
                pos, rho, h, m_part, r, extras = dens_gen.sim(arch_dir)
            elif not sim_inp:
                pos, rho, h, m_part, r, extras = dens_gen.model()

            grid_gen.translate(dat_dir, pos, rho, h, m_part, r, extras)

        elif GRID and (mod_inp == "yso"):
            grid_gen.yso(dat_dir)

        if (run_tag != "EC53"):
            opac_gen.gen(arch_dir, dat_dir, plt_dir, llambda, wav_bins)

### ------------------------------------------------------------------------ ###

    # Apply cuts and/or mass rescaling for full/components of simulation

    if sim_inp:

        if (run_tag == "EC53") and RT_RUN and MCRT:
            ec53_rt.gen(arch_dir, dat_dir, plt_dir)

        if FULL_mass_scale and not (DISC_mass_scale or ENV_mass_scale):
            scale_cut.pre_full_scale(dat_dir)

        if (DISC_mass_scale and ENV_mass_scale):
            scale_cut.pre_comp_scale(dat_dir)

        if (DISC_mass_mult and ENV_mass_mult):
            scale_cut.pre_comp_mult(dat_dir)

        if MCRT_cut:
            scale_cut.pre_cut(dat_dir)

        if MCRT_cut_scale:
            scale_cut.pre_cut_scale(dat_dir)

        if CAVITY_rho:
            scale_cut.cavity_gen(dat_dir)

        if ENV_rho:
            scale_cut.envelope_gen(dat_dir)

        if R_DESTR:
            R_DESTR = scale_cut.inner_cavity(dat_dir, t_star, r_star)

### ------------------------------------------------------------------------ ###

    # Run MCRT to generate temperature distribution

    if RT_RUN and MCRT:

        os.chdir(dat_dir+"/")
        os.system("{0}radmc3d mctherm".format(exec_loc) )

        shutil.copy2( dat_dir+"/dust_temperature.dat", \
          dat_dir+"/dust_temperature"+tdat_ext )
        shutil.copy2( dat_dir+"/photon_statistics.out", \
          dat_dir+"/photon_statistics"+pdat_ext )

### ------------------------------------------------------------------------ ###

    # If model is spatially restricted or mass scaled before RRT, do here

    if RRT_cut:
        scale_cut.post_cut(dat_dir, tdat_ext)

### ------------------------------------------------------------------------ ###

    # Run RRT to generate SED. Also move .dat/.inp files dependent on which
    # radiation sources are included in raytracing

    if RT_RUN and RRT:

        os.chdir(dat_dir+"/")
        print os.getcwd(), dat_dir, tdat_ext
        shutil.copy2(dat_dir+"/dust_temperature"+tdat_ext, \
          dat_dir+"/dust_temperature.dat" )

        if not os.path.isdir(dat_dir+"/isrf/"):
            os.makedirs(dat_dir+"/isrf")

        if incl_isrf and not incl_star:
            shutil.move(dat_dir+"/external_source.inp", \
              dat_dir+"/isrf/external_source.inp")
            shutil.move(dat_dir+"/bin/stars.inp", dat_dir+"/stars.inp")
        elif incl_isrf and incl_star:
            shutil.move(dat_dir+"/external_source.inp", \
              dat_dir+"/isrf/external_source.inp")

        for i in range(len(inclins)):

            if incl_star:
                os.system("{0}radmc3d sed incl {1}".format(exec_loc, \
                  inclins[i]) )
            elif not incl_star:
                os.system("{0}radmc3d sed incl {1} nostar".format(exec_loc, \
                  inclins[i]) )

            sdat_ext = plt_ext+"_"+str(inclins[i])+"i.out"

            shutil.copy2(dat_dir+"/spectrum.out", dat_dir+"/spectrum"+sdat_ext)

            if (mod_inp == "protostar"):
                shutil.copy2(dat_dir+"/spectrum"+int(t_star)+".out", \
                  arch_dir+"/protostar/spectrum"+int(t_star)+".out")

### ------------------------------------------------------------------------ ###

    # Run RRT to generate SED, for multipled RRT_cut values. Used to quickly
    # evaluate impact of changing effective beam size on flux

    if MULTI_RRT:

        multi_beam.run(dat_dir, tdat_ext, plt_ext)

### ------------------------------------------------------------------------ ###

    # Generate diagnostic plots, and SED. Run RRT for images if called for. Also
    # generate location dependent column density/optical depth properties

    if not RT_RUN:

        os.chdir(dat_dir)

        if ( plot_radial or plot_slice or plot_sed):
            xloc, yloc, zloc, rloc, xwid, rloc_xy, nrspecies, rho, diag_rho, \
              temp, diag_temp, nphot, diag_nphot, diag_sed, beam_r \
              = inp_plt.read(dat_dir, plt_ext, tdat_ext, pdat_ext)

        if plot_radial and ( diag_rho or diag_temp ):
            inp_plt.radial(dat_dir, plt_dir, plt_ext, rloc, rloc_xy, \
              nrspecies, rho, diag_rho, R_DESTR, temp, diag_temp)

        if plot_slice and (diag_rho or diag_temp or diag_nphot):
            inp_plt.run_slice(plt_dir, plt_ext, xloc, yloc, zloc, xwid, \
              nrspecies, rho, diag_rho, temp, diag_temp, nphot, diag_nphot)

        if plot_sed and diag_sed:

            for i in range(len(inclins)):

                sdat_ext = plt_ext+"_"+str(inclins[i])+"i.out"

                obs_plt.sed(arch_dir, dat_dir, sdat_ext, plt_dir, plt_ext, \
                  inclins[i], t_star, r_star, t_planet, r_planet, beam_r)

        if plot_image:
            for i in range(len(inclins)):
                obs_plt.image(plt_dir, inclins[i])

        if plot_rays:

            # First, generate column density map

            raytrace.column(dat_dir, plt_dir, plt_ext, plt_form)

            # Then generate locations at which optical depth reaches prescribed
            # value, to diagnose as a function of observer location

            raytrace.tau3d(dat_dir, plt_dir, plt_form, tdat_ext)

            # Now compute, for specified observer what the cloud emission
            # is as a function of optical depth

            raytrace.tau1d(dat_dir, plt_dir, plt_form, tdat_ext)

        if mass_eval:

            diagnose.mass(dat_dir)

        if lum_eval:

            diagnose.lum(dat_dir, tdat_ext, plt_ext)


### ------------------------------------------------------------------------ ###

    # Remove modules from system list, so new params_run file can be loaded

    # Also transfer params file to params_done_dir, to track which runs complete

    os.chdir("{0}/src/".format(arch_dir))
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

os.remove(arch_dir+"/src/radmc3d.out")

dirs = [arch_dir+"/ics/", arch_dir+"/src/", arch_dir+"/"]
exts = [".pyc", ".ascii"]
for d in range(len(dirs)):
    dir = dirs[d]
    files = os.listdir(dir)
    for file in files:
        for e in range(len(exts)):
            if file.endswith(exts[e]):
                os.remove(os.path.join(dir,file))

exit()
