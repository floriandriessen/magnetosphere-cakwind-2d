# Magnetosphere_CAK_2d
Make MHD model of the radiation-driven wind of a magnetic O star interacting with its magnetosphere using MPI-AMRVAC. Stellar wind is spherically symmetric and according to CAK theory including finite-disk correction. So far an isothermal MHD flow is assumed.

## Setup 

After cloning go into the local directory and fetch AMRVAC makefile (assuming you have exported the AMRVAC_DIR variable in the .bashrc (linux) or .bash_profile (macos))
```
$ setup.pl -d=2
```
and do a make (if the AMRVAC source was not compiled yet, this will be done now as well)
```
$ make -j
```

## User options

Simulations can be run using included .par files. The usual procedure to do a simulation is as follows

1. Make a relaxed CAK wind model that serves as input to the magnetosphere simulation (**usr_make_cak.par**).
2. Perform a magnetosphere simulations based on a relaxed CAK wind model.

Step 2 can be done using standard MHD techniques (**usr_magnetosphere.par**) or using Tanaka's magnetic field splitting technique (**usr_magnetosphere_tanaka.par**). Tanaka's method is particularly recommended for highly magnetic flows. 

For step 2 a restart is required from a corresponding relaxed CAK wind model. This is specified in the .par file (make sure it points to the correct file). 

It can occur that a magnetosphere simulation has to be extended. This can be done using a resume which resumes the simulation from the last generated snapshot file. The difference between a restart and a resume is thus that a restart can be done using **any** snapshot, while the resume automatically takes the last snapshot. To enable a resume the following has to be changed in the .par file

1. Comment *restart_from_file* in the *filelist*.
2. Comment *reset_time* and *time_init* in the *savelist*.
3. Ensure that *time_max* in *stoplist* is bigger than the time of the last snapshot (otherwise it cannot be extended...).

When restarting the code AMRVAC has to be told that we want to do a resume using its **-resume** flag. For example, straight from the terminal
```
mpiexec -np XXX ./amrvac -i usr_xxx.par -resume
```
but this can of course equally be done in a BATCH script on some HPC cluster.

## User parameters

The meaning of the AMRVAC parameter lists and variables in the .par file can be looked up on the AMRVAC webpage. The .par files provided are 'basic' ones that work for the problem. 

Additionally, an *star_list* is specified in the .par file containing variables relevant for our problem. The ones relevant for computations are converted in the code to unitless quantities.

+ lstar = stellar luminosity (solar units)
+ mstar = stellar mass (solar units)
+ rstar = stellar radius (solar units)
+ twind = wind temperature (Kelvin)
+ bpole = polar magnetic field strength (Gauss)
+ rhobound = boundary density (g/cm^3)
+ alpha = CAK line-force parameter (no units)
+ Qbar = Gayley's line-force parameter (no units)
+ Qmax = OCR's cut-off parameter (no units)
+ kappae = electron scattering opacity
+ beta = beta velocity law power index (no units)
+ tstat = time to start to compute statistical averages (sec)

## Notice
Tested with AMRVAC version 2.3 (Fall 2019).

## Known issues

+ using Constrained Transport for divergence B leads to strange blob developments in simulations (CT implementation excluded for now).
+ using a user specified AMR region (e.g. to resolve the magnetic equatorial plane) leads to streaks at the boundary of the refined and non-refined grid.
