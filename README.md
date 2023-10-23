

# 2.5D CAK Wind-Magnetosphere interaction

Make MHD model of the radiation-driven wind of a magnetic OB-star interacting with its magnetosphere using MPI-AMRVAC. The magnetosphere can be dynamical or centrifugal. The stellar wind is modelled with the CAK force parametrisation and the force can be purely radial radiation or include non-radial radiation components. The CAK wind uses the CAK radiation force module of MPI-AMRVAC. Both isothermal and non-isothermal simulations are supported.

*A simplified problem is available in the MPI-AMRVAC magnetohydrodynamic test problem folder.*

## Setup

After cloning go into the local directory and fetch AMRVAC makefile (assuming you have exported the `AMRVAC_DIR` variable in the `.bashrc` on linux or `.bash_profile` on macos)
```
$ setup.pl -d=2
```
and do a make (if the AMRVAC source was not compiled yet, this will be done now as well)
```
$ make -j
```

## How-to

### Magnetosphere type

Both Dynamical (DM) and Centrifugal Magnetosphere (CM) simulations can be done. The difference is essentially set by the relative positions of Alfven and Kepler radius. A DM has $R_{Alf} < R_{Kep}$ while the CM has $R_{Kep} < R_{Alf}$. To make a magnetosphere only the polar magnetic field strength or magnetic confinement parameter have to be set.

Formally a DM can be run in the non-rotating star case---in the published works on DMs the rotation is set to zero. However, adding a little rotation can be quite beneficial if studying a DM with a strong magnetic confinement. Essentially this makes the little mass outflow in the magnetic equatorial plane more stable. There is no rule of thumb for choosing the rotation rate, but several km/s is fine.

To make a CM the star should be rotating fast and have a strong magnetic confinement. This is the easiest way to enforce the CM condition. A rotation rate of half critical rotation is a good setup. Going above `Wrot = 0.75d0` is not recommended as the oblateness of the star should then be properly taken into account, which is currently not implemented.

### Performing a simulation

Simulations can be run using included .par files. Specific stellar parameters are set in the `star_list` inside the .par file. The CAK wind parametrisation is handled with the `cak_list`.

Several options exist for doing the MHD part. So far no Constrained Transport is implemented due to numerical issues in the magnetosphere. The code runs with other divergence cleaners and **Powell** or **GLM** are recommended. **NOTE**: if you use the GLM cleaner then the additional variable 'psi' will be added in the output.

It might be beneficial in some user cases to turn on Tanaka's magnetic field splitting (added stability + in adiabatic MHD reduced probability of getting negative pressures). Tanaka's method is particularly recommended for highly magnetic flows (magnetic confinement >> 1) and can be switched on in `mhd_list` via
```
B0field = .true.
```

and in the code the `make_dipoleboy` subroutine is called instead to add a global, time-independent background dipole field. 

## Additional physics options

### Adiabatic MHD with radiative cooling

The possibility exists to take into account the full MHD energy equation including radiative cooling in the wind-magnetosphere dynamics. A separate .par file exists for this option. Essentially in the `mhd_list` additional specifications are made
```
mhd_energy = .true.
mhd_radiative_cooling = .true.
```

It is important to stress that both switches have to be logically true. Therefore, pure adiabatic MHD is **not supported**. For the radiative cooling an additional `rc_list` is included specifying some parameters related to the cooling table. I suggest to always take the 'SPEX' cooling curve since this is presently the best on the market for our type of temperature ranges. An alternative can be the Bailey & MacDonald curve 'MB' which has been used in previous magnetic CAK wind models in the literature including energy transfer.

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
also ensure that `time_max` in `stoplist` is bigger than the time of the last snapshot (otherwise it cannot be extended...). Additionally commented in the `savelist` is the generation of simulation log + simulation output at startup as this is generally unwanted. Also because we want to resume we do not want to reset the time or initialise the time in `stoplist`. If the `base_filename` is pointed correctly in `filelist`, AMRVAC will just continue writing output to this directory.

When resuming the AMRVAC call requires the additional flag **-resume**. For example,
```
mpiexec -np <NUM> ./amrvac -i <PARFILE> -resume
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
In this case the snapshot count starts with *snapshot number +1* from where you started (if you also reset time, the snapshot count will start at 0 and start to overwrite your original data!!!, likewise we do not want to initialise the time again). Contrary to a resume, you can let AMRVAC write the data to another directory by specifying a different `base_filename` in `filelist`. This could be useful if you do not want to overwrite existing snapshots following your snapshot count.

Finally, to run the simulation you can just do the same executable call as you would normally do. However, it is also possible to specify a restart directly from the command line. In order to do so an additional flag **-if path_to_datfile** has to be included (and no `restart_from_file` has to be specified in the .par file). For example,
```
mpiexec -np <NUM> ./amrvac -i <PARFILE> -if path_to_datfile_to_start
```

## Additional user parameters

The meaning of the AMRVAC parameter lists and variables in the .par file can be looked up on the AMRVAC webpage. The .par files provided are 'basic' ones that work for the problem. 

Additionally, a `star_list` is specified in the .par file containing variables specific for our problem. The ones relevant for computations are converted in the code to unitless quantities within the `make_dimless_vars` routine.

| Parameter| Explanation                                                       |
|----------|-------------------------------------------------------------------|
| mstar    | stellar mass (solar units)                                        |
| rstar    | stellar radius (solar units)                                      |
| twind    | wind temperature (Kelvin)                                         |
| bpole, etastar | polar magnetic field strength (G) or wind magnetic confinement; one of these parameters has to be specified |
| rhobound | stellar boundary density (g/cm^3)                                 |
| tstat    | time to start to compute statistical averages (sec)               |
| Wrot     | ratio of equatorial rotational velocity to orbital (critical)  velocity |

### More info on the statistics

The `tstat` should generally be set to a time where the wind reaches a time and state where initial perturbations have left the grid. If not, the statistical averaging will contain remnants from these perturbations. 

The computation of an averaged value $\langle X \rangle$ is done in the `compute_stats` routine. Since Spring 2021 it does a simplified averaging relying only on the current hydro state and not on the previous anymore. Since averages are computed every iteration variations will smoothen out quickly, as such an additional weighting with previous iteration hydro state is obsolete.

For proper behaviour during the simulation; after each computation of $\langle X \rangle$ a normalisation occurs with the current time-weight `tnormc` (required for correct output) but in the next iteration we have to 'un-normalise' it back in time `tnormp` to compute the new $`\langle X \rangle`$ again.

The routine is called at the **end** of the iteration (i.e. after the advection has been performed) and the w-array is thus evolved. However, the *nwextra* variables (defined in the code with `var_set_extravar()`) are never evolved and so they are still at the state before advection (i.e. previous timestep). If wished for, it is straightforward to also include here new quantities of interest. **Note** that although the routine is called each iteration, the actual values are only printed every `dtsave_dat`.

## Known issues (from AMRVAC version 2.2)

+ using Constrained Transport for divergence B leads to strange blob developments in simulations (CT implementation excluded for now).
+ using a user specified AMR region (e.g. to resolve the magnetic equatorial plane) leads to streaks at the boundary of the refined and non-refined grid.
+ fast rotating stars have difficulty maintaining the CM disc, which starts to oscillate, and turns into DM behaviour.

