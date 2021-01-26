

# 2D CAK Wind-Magnetosphere interaction

Make MHD model of the radiation-driven wind of a magnetic O-star interacting with its magnetosphere using MPI-AMRVAC. The magnetosphere can be dynamical or centrifugal. The stellar wind is spherically symmetric and according to CAK theory including finite-disk correction. In the current setup an isothermal MHD flow is assumed.

## Setup

After cloning go into the local directory and fetch AMRVAC makefile (assuming you have exported the AMRVAC_DIR variable in the .bashrc (linux) or .bash_profile (macos))
```
$ setup.pl -d=2
```
and do a make (if the AMRVAC source was not compiled yet, this will be done now as well)
```
$ make -j
```

## How-to

### Magnetosphere type

Both Dynamical (DM) and Centrifugal Magnetosphere (CM) simulations can be done. The difference is essentially set by the relative positions of Alfven and Kepler radius. A DM has R_alf < R_kep while the CM has R_kep < R_alf.

To make a DM only the polar magnetic field strength or eta parameter have to be set. In the published works on DMs the rotation is set to zero. However, as suggested by Asif adding a little rotation can be quite beneficial if studying a DM with a strong confinement. Essentially this makes the little mass outflow in the magnetic equatorial plane more stable. There is no rule of thumb for choosing the rotation rate, but several km/s is fine.

To make a CM the star should be rotating fast and have a strong magnetic confinement. This is the easiest way to enforce the CM condition. A rotation rate of half critical rotation is a good setup. Going above *Wrot = 0.75d0* is not recommended as the oblateness of the star should then be properly taken into account, which is currently not implemented.

### Performing a simulation

Simulations can be run using included .par files. Specific stellar and wind parameters are set in the `star_list` inside the .par file. 

Several options exist for doing the MHD part. So far no Constrained Transport is implemented due to numerical issues in the magnetosphere. The code runs with other divergence cleaners and **Powell** or **GLM** are recommended. **NOTE**: if you use the GLM cleaner then the additional variable 'psi' will be added in the output. In order to not output temporary arrays that are used for statistical purpose you should increase the index of `w_write` in the `filelist` with +1. These indices index the position of a variable in the AMRVAC w-array starting to count from the density index.

It might be beneficial in some user cases to turn on Tanaka's magnetic field splitting. Tanaka's method is particularly recommended for highly magnetic flows (magnetic confinement > 1) and can be switched on in `mhd_list` via
```
B0field = .true.
```

and in the code the `make_dipoleboy` subroutine is called instead to add a global, time-independent background dipole field. 

## Extending a simulation

To extend a simulation there are two options in AMRVAC: resume and restart. Both cases are mainly useful for magnetosphere simulations.

### Resuming

It can occur that a magnetosphere simulation has to be extended because, for example, it has not yet reached its intended run time. A situation like this often occurs for very long computations running on clusters where a wall-time is acting. In order to continue our simulation use the resume option, which resumes the simulation from the **last snapshot file**. 

To enable a resume the following has to be activated and commented in the magnetosphere .par file
```
&savelist
 !itsave(1,1) = 0
 !itsave(1,2) = 0

&stoplist
 reset_it = .true.
```
also ensure that *time_max* in `stoplist` is bigger than the time of the last snapshot (otherwise it cannot be extended...). Additionally commented in the `savelist` is the generation of simulation log + simulation output at startup as this is generally unwanted. Also because we want to resume we do not want to reset the time or initialise the time in `stoplist`. If the *base_filename* is pointed correctly in `filelist`, AMRVAC will just continue writing output to this directory.

When resuming the AMRVAC call requires the additional flag **-resume**. For example,
```
mpiexec -np XXX ./amrvac -i usr_xxx.par -resume
```

### Restarting

It might also be of interest to resume a simulation not from the last snapshot, but from a **different snapshot**. For this purpose the restart option has to be used. For example, this might be useful if somewhere the code crashed and you want to restart from a file before the crash with different options.

If you want to restart from some existing magnetosphere output file this is possible by
```
&filelist
 restart_from_file = string_of_path_to_datfile

&stoplist
 reset_it = .true.
```
In this case the snapshot count starts from (snapshot number +1) from where you started (if you also reset time, the snapshot count would start at 0 and start to overwrite your original data!!!, likewise we do not want to initialise the time again). Contrary to a resume, you can let AMRVAC write the data to another directory by specifying a different *base_filename* in `filelist`. This could be useful if you do not want to overwrite existing snapshots following your snapshot count.

Finally, to run the simulation you can just do the same executable call as you would normally do. However, it is also possible to specify a restart directly from the command line. In order to do so an additional flag **-if path_to_datfile** has to be included (and no `restart_from_file` has to be specified in the .par file). For example,
```
mpiexec -np XXX ./amrvac -i usr_xxx.par -if path_to_datfile_to_start
```

## Additional user parameters

The meaning of the AMRVAC parameter lists and variables in the .par file can be looked up on the AMRVAC webpage. The .par files provided are 'basic' ones that work for the problem. 

Additionally, a `star_list` is specified in the .par file containing variables specific for our problem. The ones relevant for computations are converted in the code to unitless quantities within the `make_dimless_vars` routine.

| Parameter| Explanation                                                       |
|----------|-------------------------------------------------------------------|
| lstar    | stellar luminosity (solar units)                                  |
| mstar    | stellar mass (solar units)                                        |
| rstar    | stellar radius (solar units)                                      |
| twind    | wind temperature (Kelvin)                                         |
| imag     | switch for polar magnetic field strength or wind confinement 
|          | if `imag > 0`: `bpole = imag`
|          | if `imag < 0`: `etastar = -imag`                                                       | rhobound | boundary density (g/cm^3)                                         |
| alpha    | CAK line-force parameter (no units)                               |
| Qbar     | Gayley's line-force parameter (no units)                          |
| tstat    | time to start to compute statistical averages (sec)               |
| Wrot     | ratio of equatorial rotational velocity to orbital (critical)  velocity |

## Notice

Tested with AMRVAC version 2.3 (Fall 2019).

## Known issues

+ using Constrained Transport for divergence B leads to strange blob developments in simulations (CT implementation excluded for now).
+ using a user specified AMR region (e.g. to resolve the magnetic equatorial plane) leads to streaks at the boundary of the refined and non-refined grid.
