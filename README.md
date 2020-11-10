
# 2D CAK Wind-Magnetosphere interaction
Make MHD model of the radiation-driven wind of a magnetic O-star interacting with its magnetosphere using MPI-AMRVAC. Stellar wind is spherically symmetric and according to CAK theory including finite-disk correction. In the current setup an isothermal MHD flow is assumed.

## Setup 

After cloning go into the local directory and fetch AMRVAC makefile (assuming you have exported the AMRVAC_DIR variable in the .bashrc (linux) or .bash_profile (macos))
```
$ setup.pl -d=2
```
and do a make (if the AMRVAC source was not compiled yet, this will be done now as well)
```
$ make -j
```

## Setup options

Simulations can be run using included .par files. The usual procedure to do a simulation is as follows

1. Make a relaxed CAK wind model that serves as input to the magnetosphere simulation (**usr_cak.par**). 
2. Perform a magnetosphere simulations based on a relaxed CAK wind model.

Step 1 creates a 2D spherically symmetric, finite-disk, CAK wind. So far only relaxed CAK models of Zeta Puppis exist, but this can be easily changed to another star. Notice that only CAK models made here can be directly used for a magnetosphere and the MHD equations are solved (the magnetic field is set to zero). This is because AMRVAC stores the physics option in the .dat file and between restarts this has to be compatible.

Step 2 can be done using standard MHD techniques (**usr_magnetosphere.par**) or using Tanaka's magnetic field splitting technique (**usr_magnetosphere_tanaka.par**). Tanaka's method is particularly recommended for highly magnetic flows (magnetic confinement > 1). Whenever magnetic confinement < 1 the Tanaka method is disabled (see Issues below). To do step 2 a restart is required from a corresponding relaxed CAK wind model. This is specified in the .par file (make sure it points to the correct file). Also ensure that the mesh settings are exactly the same as used in the CAK model.

## Extending a simulation

It can occur that a magnetosphere simulation has to be extended. This is accomplished using the resume option, which resumes the simulation from the last generated snapshot file. The difference between a restart and a resume is thus that a restart can be done using **any** snapshot, while the resume automatically takes the last snapshot. To enable a resume the following has to be commented in the .par file
```
&filelist
 !restart_from_file = string

&savelist
 !reset_time = value
 !time_init  = value
```
also ensure that *time_max* in *stoplist* is bigger than the time of the last snapshot (otherwise it cannot be extended...).

When resuming the AMRVAC call requires the additional flag **-resume**. For example,
```
mpiexec -np XXX ./amrvac -i usr_xxx.par -resume
```

## Additional user parameters

The meaning of the AMRVAC parameter lists and variables in the .par file can be looked up on the AMRVAC webpage. The .par files provided are 'basic' ones that work for the problem. 

Additionally, a *star_list* is specified in the .par file containing variables specific for our problem. The ones relevant for computations are converted in the code to unitless quantities.

+ lstar = stellar luminosity (solar units)
+ mstar = stellar mass (solar units)
+ rstar = stellar radius (solar units)
+ twind = wind temperature (Kelvin)
+ imag = switch for choosing polar magnetic field strength (imag > 0) or wind confinement (imag < 0)
+ rhobound = boundary density (g/cm^3)
+ alpha = CAK line-force parameter (no units)
+ Qbar = Gayley's line-force parameter (no units)
+ Qmax = OCR's cut-off parameter (no units)
+ kappae = electron scattering opacity
+ beta = beta velocity law power index (no units)
+ tstat = time to start to compute statistical averages (sec)
+ Wrot = ratio of equatorial rotational velocity to critical velocity

## Notice
Tested with AMRVAC version 2.3 (Fall 2019).

## Known issues

+ using Constrained Transport for divergence B leads to strange blob developments in simulations (CT implementation excluded for now).
+ using a user specified AMR region (e.g. to resolve the magnetic equatorial plane) leads to streaks at the boundary of the refined and non-refined grid.
+ For low magnetic confinement (etastar < 1) the Tanaka field splitting technique is not recommended due to the formation of 'antennae' in magnetic equatorial plane. This is strange because the original paper claims that the method can be used for any plasma beta flow (~ 1/etastar).
