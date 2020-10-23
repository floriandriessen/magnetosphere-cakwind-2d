!===============================================================================
! 2D wind module for launching a CAK line-driven wind from the stellar surface
! into the magnetosphere of a magnetic massive star
!
! Setup mimicked as in ud-Doula & Owocki (2002). Testbench, similar separate
! program will be used for incorporating the LDI into the magnetosphere.
!
! Coded up by Flo for his KU Leuven PhD thesis 2018-2022
!-------------------------------------------------------------------------------
! Options to be specified in usr.par file:
!   iprob = 0  : non-magnetic, finite disk, Sobolev model (CAK)
!   iprob = 1  : magnetosphere using a relaxed non-CT/CT CAK model as input
!
! Separate .par file exist containing required setups to start
!
! If iprob=1 a restart is required (specify in .par file_list) from the
! .dat file containing the hydro variables of a relaxed 2D CAK model (iprob=0):
!     $ {executable call} -i usr_magnetosphere.par
!
! (October 2019) -- Flo
!   > setup of problem with 1D CAK wind starting from beta law in 2D geometry
!   > implementation special refinement routine near magnetic equatorial plane
!
! (November 2019) -- Flo
!   > included AMRVAC 'iprob' mesh_list parameter to switch between setups
!   > testing shows that Dedner divergence cleaners lead to negative pressures
!     near the pole, build a MPI stop to prevent using this cleaner family
!   > implemented computation and storage auxiliary vars (Alfven speed + div B)
!     using usr_aux_output
!
! (December 2019) -- Flo
!   > extension code to allow for new Constrained Transport method for B-field
!     (needs 2.5D for vector potential, but phi components of vectors are zero)
!     (can actually allow for rotation now in 2D)
!   > typedivbfix together with iprob=0 in .par file now determines whether
!     CAK simulation will be staggered (for CT) or not (for non-CT)
!     does not matter for CAK, but it is saved in AMRVAC restart configs and
!     required for restart with magnetosphere
!
! (May 2020) -- Flo
!   > renaming of unit_mass/lum/ggrav in order to compile without conflict as
!     AMRVAC v2.3 has a new particle module (with protected variable unit_mass)
!   > removed finite disk factor from output (only useful for testing)
!   > implemented routine for computing time-averaged statistical variables
!   > for post-simulation convenience at startup a log file gets created that
!     contains an overview of all relevant parameters
!
! (June 2020) -- Flo
!   > implementation to do Bfield splitting via Tanaka method by calling the
!     usr_set_B0 routine containing a stellar dipole
!     NOTICE: in AMRVAC there is build-in dipole that can be activated with
!             'Bdip' in .par file, however, this formula leads to an eta_star
!             4x bigger (AMRVAC has different way of expressing dipole)
!   > change to better velocity gradient stencil to deal with non-uniform grids
!   > removed Constrained Transport implementation, extensive testing has shown
!     that unphysical results and are unreliable (for this application here)
!
!===============================================================================

module mod_usr

  ! Include a physics module
  use mod_mhd

  implicit none

  ! The usual suspects
  real(8) :: msun=1.989d33, lsun=3.827d33, rsun=6.96d10, Ggrav=6.67d-8

  ! Stellar parameters: luminosity, mass, radius, polar magnetic field,
  !                     surface density, eff. temp., Alfven radius
  real(8) :: lstar, mstar, rstar, bpole, rhobound, twind, ralf

  ! Unit quantities that are handy: gravitational constant, luminosity, mass
  real(8) :: my_unit_ggrav, my_unit_lum, my_unit_mass

  ! log(g), eff. log(g) + scale height, year [s], mean mol. weight
  real(8) :: logg, logge, heff, year=3.15567d7, mumol

  !
  ! Wind parameters: CAK alpha, Gayley Qbar + Qmax, opacity electron scattering,
  !                  beta power velocity law, Eddington gamma, escape speed,
  !                  CAK + fd mass-loss rate, terminal wind speed, sound speed
  !                  equatorial magnetic field, wind magnetic confinement
  real(8) :: alpha, Qbar, Qmax, kappae, beta, gammae, vesc, mdot, mdotfd, &
             vinf, asound, etastar

  ! Time-step accounting radiation force, time to start statistical computation
  real(8) :: new_timestep, tstat

  ! Dimensionless variables of relevant variables above
  real(8) :: dlstar, dmstar, drstar, dbpole, drhobound, dtwind, dkappae, &
             dvesc, dvinf, dmdot, dasound, dclight, dGgrav, dgammae, detaconf, &
             dtstat

  ! Additional names for wind variables
  integer :: my_gcak

  ! Additional names for statistical wind variables
  integer :: my_rhoav, my_rho2av, my_vrav, my_vr2av, my_rhovrav, my_rho2vrav,&
             my_vpolav, my_vpol2av, my_rho2vpolav
  integer :: my_tmp1, my_tmp2, my_tmp3, my_tmp4, my_tmp5,my_tmp6, my_tmp7, &
             my_tmp8, my_tmp9

  character(len=8)  :: todayis
  character(len=99) :: inputlog

contains

  !> This routine should set user methods, and activate the physics module
  subroutine usr_init()
    use mod_global_parameters
    use mod_usr_methods

    ! Choose coordinate system: 2.5D spherical
    call set_coordinate_system("spherical_2.5D")

    ! Read in the user-defined star_list (not read in by AMRVAC)
    call usr_params_read(par_files)

    !
    ! Choose independent normalization units, only 3 have to be specified:
    !     (length,temp,ndens) or (length,vel,ndens)
    ! Numberdensity chosen such that unit density becomes boundary density
    ! IMPORTANT: AMRVAC cannot be run in CGS units <=> unit quantities = 1
    !
    unit_length        = rstar                                      ! cm
    unit_temperature   = twind                                      ! K
    unit_numberdensity = rhobound/((1.d0+4.d0*He_abundance)*mp_cgs) ! g cm^-3

    ! Invoke the physics module
    call MHD_activate()

    ! Initialize user (global) quantities after initializing AMRVAC
    usr_set_parameters => initglobaldata_usr

    ! Set up initial conditions
    usr_init_one_grid => initial_conditions

    ! Special boundary conditions
    usr_special_bc => special_bound

    ! CAK line-force computation
    usr_source => CAK_source

    ! Adjusted timestep to account for total radiation force
    usr_get_dt => special_dt

    ! Every iteration retrieve global grid info and perform operations on it
    ! to make time-averaged statistical quantities of wind
    usr_process_grid => compute_stats

    ! Compute and include Alfven speed and divergence of Bfield as output
    ! (only stored when doing magnetosphere models)
    usr_aux_output    => set_extravar_output
    usr_add_aux_names => set_extravarnames_output

    ! Background dipole field if using Tanaka field splitting
    usr_set_B0 => make_dipoleboy

    ! User source variables to store in output (only temporary in source comp.)
    my_gcak = var_set_extravar("gcak", "gcak")

    ! Statistical quantities and temporary storage variables
    my_rhoav      = var_set_extravar("rho_av", "rho_av")
    my_rho2av     = var_set_extravar("rho2_av", "rho2_av")
    my_vrav       = var_set_extravar("vrad_av", "vrad_av")
    my_vpolav     = var_set_extravar("vtheta_av", "vtheta_av")
    my_vr2av      = var_set_extravar("vrad2_av", "vrad2_av")
    my_vpol2av    = var_set_extravar("vtheta2_av", "vtheta2_av")
    my_rhovrav    = var_set_extravar("rho_vrad_av", "rho_vrad_av")
    my_rho2vrav   = var_set_extravar("rho2_vrad_av", "rho2_vrad_av")
    my_rho2vpolav = var_set_extravar("rho2_vtheta_av", "rho2_vtheta_av")

    my_tmp1 = var_set_extravar("tmp1","tmp1")
    my_tmp2 = var_set_extravar("tmp2","tmp2")
    my_tmp3 = var_set_extravar("tmp3","tmp3")
    my_tmp4 = var_set_extravar("tmp4","tmp4")
    my_tmp5 = var_set_extravar("tmp5","tmp5")
    my_tmp6 = var_set_extravar("tmp6","tmp6")
    my_tmp7 = var_set_extravar("tmp7","tmp7")
    my_tmp8 = var_set_extravar("tmp8","tmp8")
    my_tmp9 = var_set_extravar("tmp9","tmp9")

  end subroutine usr_init

!===============================================================================

  subroutine usr_params_read(files)
    !
    ! Read in the usr.par file with the problem specific list
    !
    use mod_global_parameters, only: unitpar, stagger_grid
    use mod_constants

    character(len=*), intent(in) :: files(:)
    integer :: n

    namelist /star_list/ mstar, lstar, rstar, twind, bpole, rhobound, alpha, &
                          Qbar, Qmax, kappae, beta, tstat

    do n = 1,size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, star_list, end=111)
       111 close(unitpar)
    end do

    ! Scale to cgs units
    lstar = lstar * lsun
    mstar = mstar * msun
    rstar = rstar * rsun

    if (typedivbfix == 'ct') then
      call mpistop('CT disabled. Gives strange results for this problem.')
    endif

  end subroutine usr_params_read

!===============================================================================

  subroutine initglobaldata_usr
    !
    ! Compute some quantities of interest (in CGS) before making unitless
    !
    use mod_global_parameters
    use mod_constants

    ! Stellar structure
    gammae = kappae * lstar/(4.d0*dpi * Ggrav * mstar * const_c)
    logg   = log10(Ggrav * mstar/rstar**2.0d0)
    logge  = logg * (1.0d0 - gammae)**0.5d0
    mumol  = (1.0d0 + 4.0d0*He_abundance)/(2.0d0 + 3.0d0*He_abundance)
    asound = dsqrt(twind * kB_cgs/(mumol * mp_cgs))
    heff   = asound**2.0d0 / 10.0d0**logge

    ! Wind quantities in CAK theory
    Qmax    = Qmax * Qbar
    vesc    = (2.0d0 * Ggrav * mstar*(1.0d0 - gammae)/rstar)**0.5d0
    vinf    = vesc * (alpha/(1.0d0 - alpha))**0.5d0
    mdot    = lstar/const_c**2.0d0 * alpha/(1.0d0 - alpha) &
               * (Qbar * gammae/(1.0d0 - gammae))**((1.0d0 - alpha)/alpha)
    mdotfd  = mdot/(1.0d0 + alpha)**(1.0d0/alpha)
    etastar = ((bpole/2.0d0)**2.0d0 * rstar**2.0d0)/(mdot * vinf)
    ralf    = 1.0d0 + (etastar + 0.25d0)**0.25d0 - 0.25d0**0.25d0

    ! Make all relevant variables dimensionless
    call make_dimless_vars()

    if (mype == 0) then
      ! Store overview in a log file for easy reference
      inputlog = trim(base_filename) // '_param_overview.log'
      open(unit=94,file=inputlog)

      call date_and_time(todayis)
      write(94,*) 'MPI-AMRVAC simulation ran on ', &
                    todayis(7:8), '/', todayis(5:6), '/', todayis(1:4)
      write(94,*)
      write(94,*) '======================'
      write(94,*) '   Unity quantities   '
      write(94,*) '======================'
      write(94,*) 'unit length        = ', unit_length
      write(94,*) 'unit density       = ', unit_density
      write(94,*) 'unit velocity      = ', unit_velocity
      write(94,*) 'unit numberdensity = ', unit_numberdensity
      write(94,*) 'unit pressure      = ', unit_pressure
      write(94,*) 'unit temperature   = ', unit_temperature
      if (iprob == 1) write(94,*) 'unit magneticfield = ', unit_magneticfield
      write(94,*) 'unit time          = ', unit_time
      write(94,*)
      write(94,*) '==========='
      write(94,*) '   SETUP   '
      write(94,*) '==========='
      write(94,*) 'Problem options:                                  '
      write(94,*) '  iprob = 0 : non-magnetic, finite disk, CAK model'
      write(94,*) '  iprob = 1 : CAK-magnetosphere model             '
      write(94,*) '> Chosen option = ', iprob
      if (iprob == 0) then
        write(94,*) '> Making a CAK wind to be put in a magnetosphere model'
      elseif (iprob == 1) then
        write(94,*) '> Making a magnetosphere model with cleaner:', typedivbfix
      else
        call mpistop('Choose a valid iprob {0,1}')
      endif
      write(94,*)
      write(94,*) '==============================================='
      write(94,*) '   Stellar and wind parameters in CGS units    '
      write(94,*) '==============================================='
      write(94,*) 'Lstar/Lsun             = ', lstar/lsun
      write(94,*) 'Mstar/Msun             = ', mstar/msun
      write(94,*) 'Rstar/Rsun             = ', rstar/rsun
      write(94,*) 'Twind                  = ', twind
      if (iprob == 1) write(94,*) 'Polar magnetic field   = ', bpole
      if (iprob == 1) write(94,*) 'Wind confinement eta   = ', etastar
      if (iprob == 1) write(94,*) 'Ralf/Rstar             = ', ralf
      write(94,*) 'Mean molecular weight  = ', mumol
      write(94,*) 'log(g)                 = ', logg
      write(94,*) 'eff. log(g)            = ', logge
      write(94,*) 'eff. scale height heff = ', heff
      write(94,*) 'heff/Rstar             = ', heff/rstar
      write(94,*) 'Eddington gamma        = ', gammae
      write(94,*)
      write(94,*) 'adiabatic gamma = ', mhd_gamma
      write(94,*) 'alpha           = ', alpha
      write(94,*) 'Qbar            = ', Qbar
      write(94,*) 'Qmax/Qbar       = ', Qmax/Qbar
      write(94,*) 'asound          = ', asound
      write(94,*) 'eff. vesc       = ', vesc
      write(94,*) 'vinf            = ', vinf
      write(94,*)
      write(94,*) 'surface density        = ', rhobound
      write(94,*) 'analytic Mdot CAK      = ', mdot * (year/msun)
      write(94,*) '... with FD correction = ', mdotfd * (year/msun)
      write(94,*)
      write(94,*) '========================================'
      write(94,*) '    Dimensionless AMRVAC quantities     '
      write(94,*) '========================================'
      write(94,*) 'Extra computed unit quantities:'
      write(94,*) '   unit Lum  = ', my_unit_lum
      write(94,*) '   unit Mass = ', my_unit_mass
      write(94,*) '   unit Grav = ', my_unit_ggrav
      write(94,*) 'Lstar        = ', dlstar
      write(94,*) 'Mstar        = ', dmstar
      write(94,*) 'Rstar        = ', drstar
      write(94,*) 'Twind        = ', dtwind
      if (iprob == 1) write(94,*) 'Bpole        = ', dbpole
      if (iprob == 1) write(94,*) 'Eta conf.    = ', detaconf
      write(94,*) 'Edd. gamma   = ', dgammae
      write(94,*) 'rhobound     = ', drhobound
      write(94,*) 'Mdot         = ', dmdot
      write(94,*) 'alpha        = ', alpha
      write(94,*) 'Qbar         = ', Qbar
      write(94,*) 'Qmax         = ', Qmax/Qbar
      write(94,*) 'kappae       = ', dkappae
      write(94,*) 'asound       = ', dasound
      write(94,*) 'eff. vesc    = ', dvesc
      write(94,*) 'vinf         = ', dvinf
      write(94,*) 'clight       = ', dclight
      write(94,*) 'Ggrav        = ', dGgrav
      write(94,*) 'Tstat        = ', dtstat
      write(94,*)
      close(94)
    endif

  end subroutine initglobaldata_usr

!===============================================================================

  subroutine initial_conditions(ixI^L, ixO^L, w, x)
    !
    ! Initial conditions are straight computed in dimensionless units (easy)
    !
    use mod_global_parameters

    ! Subroutine arguments
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    ! Local variables
    real(8) :: rho, vr, sfac
    integer :: i, j          ! dummy indices for radial + theta grid

    ! Make a non-magnetic CAK model to be used later in a magnetosphere
    if (iprob==0) then

      ! Small offset (asound/vinf) to avoid starting at terminal wind speed
      sfac = 1.0d0 - 1.0d-3**(1.0d0/beta)

      do j = ixImin2,ixImax2
        do i = ixImin1,ixImax1
          if (x(i,j,1) >= drstar) then

            ! Set initial radial velocity field to beta law
            vr = dvinf * (1.0d0 - sfac*(drstar/x(i,j,1)))**beta

            ! Set initial density
            rho = dmdot / (4.0d0*dpi * x(i,j,1)**2.0d0 * vr)

            w(i,j,rho_)   = rho
            w(i,j,mom(1)) = rho * vr

            ! Set no polar and azimuthal velocities
            w(i,j,mom(2)) = 0.0d0
            w(i,j,mom(3)) = 0.0d0
          else
            w(i,j,rho_)   = drhobound
            w(i,j,mom(1)) = 0.0d0
            w(i,j,mom(2)) = 0.0d0
            w(i,j,mom(3)) = 0.0d0
          endif
        enddo
      enddo

      ! Set the magnetic field components (doing HD so no magnetic field now)
      w(ixI^S,mag(1)) = 0.0d0
      w(ixI^S,mag(2)) = 0.0d0
      w(ixI^S,mag(3)) = 0.0d0
    endif

    ! Make a magnetosphere: restart is happening, only magnetic field assigned
    if (iprob==1) then
      !
      ! Setup magnetic field based on vector potential (CT), or just regular
      ! magnetic field components of dipole
      !
      if (B0field) then
        w(ixI^S,mag(:)) = 0.0d0
      else
        ! Radial magnetic field
        w(ixI^S,mag(1)) = dbpole * (drstar/x(ixI^S,1))**3.0d0 * dcos(x(ixI^S,2))

        ! Polar magnetic field
        w(ixI^S,mag(2)) = 0.5d0*dbpole * (drstar/x(ixI^S,1))**3.0d0 * dsin(x(ixI^S,2))

        ! Azimuthal magnetic field
        w(ixI^S,mag(3)) = 0.0d0
      endif

      ! Stop if using Dedner+ (2002) divergence cleaning
      if(mhd_glm) w(ixO^S,psi_)=0.d0
      !if (mhd_glm) call mpistop('Dedner divergence cleaner leads somehow to negative pressures near poles')

    endif

    ! Initialize the CAK line-force and statistical quantities
    w(ixO^S,my_gcak)       = 0.0d0
    w(ixO^S,my_rhoav)      = 0.0d0
    w(ixO^S,my_rho2av)     = 0.0d0
    w(ixO^S,my_vrav)       = 0.0d0
    w(ixO^S,my_vr2av)      = 0.0d0
    w(ixO^S,my_vpolav)     = 0.0d0
    w(ixO^S,my_vpol2av)    = 0.0d0
    w(ixO^S,my_rhovrav)    = 0.0d0
    w(ixO^S,my_rho2vrav)   = 0.0d0
    w(ixO^S,my_rho2vpolav) = 0.0d0

    w(ixO^S,my_tmp1) = 0.0d0
    w(ixO^S,my_tmp2) = 0.0d0
    w(ixO^S,my_tmp3) = 0.0d0
    w(ixO^S,my_tmp4) = 0.0d0
    w(ixO^S,my_tmp4) = 0.0d0
    w(ixO^S,my_tmp5) = 0.0d0
    w(ixO^S,my_tmp7) = 0.0d0
    w(ixO^S,my_tmp8) = 0.0d0
    w(ixO^S,my_tmp9) = 0.0d0

  end subroutine initial_conditions

!===============================================================================

  subroutine special_bound(qt,ixI^L,ixB^L,iB,w,x)
    !
    ! Modified boundary values only at lower radial boundary (star)
    !
    use mod_global_parameters

    ! Subroutine arguments
    integer, intent(in)    :: ixI^L, ixB^L, iB
    real(8), intent(in)    :: qt, x(ixI^S,1:ndim)
    real(8), intent(inout) :: w(ixI^S,1:nw)

    ! Local variables
    real(8) :: vr(ixI^S), vpol(ixI^S)
    integer :: i, j      ! dummy indices for radial + theta grid
    integer :: ixBs^L    ! Face indices cell

    vr(ixI^S) = w(ixI^S,mom(1))/w(ixI^S,rho_) ! velocity from momentum
    vpol(ixI^S) = w(ixI^S,mom(2))/w(ixI^S,rho_)

    select case (iB)

    case(1) ! Left boundary (stellar surface)

      w(ixB^S,rho_) = drhobound

      ! Radial velocity field (slope extrapolation in 1st ghost, then constant)
      do i = ixBmax1,ixBmin1,-1
        if (i==ixBmax1) then
          vr(i,:) = vr(i+1,:) - (vr(i+2,:) - vr(i+1,:)) &
                            * (x(i+1,:,1) - x(i,:,1))/(x(i+2,:,1) - x(i+1,:,1))
        else
          vr(i,:) = vr(i+1,:)
        endif

        w(i,:,mom(1)) = drhobound * vr(i,:)
      enddo

      ! Polar velocity field
      w(ixB^S,mom(2)) = 0.0d0

      ! Azimuthal velocity field
      w(ixB^S,mom(3)) = 0.0d0

      ! Non-magnetic CAK model
      if (iprob==0) then
        w(ixI^S,mag(1)) = 0.0d0
        w(ixI^S,mag(2)) = 0.0d0
        w(ixI^S,mag(3)) = 0.0d0
      endif

      ! Restart from 2D spherically symmetric non-magnetic CAK model
      if (iprob==1) then

        if (B0field) then
          w(ixB^S,mag(:)) = 0.0d0
        else
          ! Radial magnetic field component
          w(ixB^S,mag(1)) = dbpole * dcos(x(ixB^S,2))

          ! Polar magnetic field component
          if (etastar >= 1.0d0) then
            !
            ! Magnetic confinement, do linear extrapolation
            !
            do i = ixBmax1,ixBmin1,-1
              w(i,:,mag(2)) = w(i+1,:,mag(2)) &
                              - (w(i+2,:,mag(2)) - w(i+1,:,mag(2))) &
                              * (x(i+1,:,1) - x(i,:,1))/(x(i+2,:,1) - x(i+1,:,1))
            enddo
          else
            !
            ! Negligible confinement, set to zero
            !
            w(ixB^S,mag(2)) = 0.0d0
          endif

          ! Azimuthal magnetic field
          w(ixB^S,mag(3)) = 0.0d0
        endif

        if(mhd_glm) w(ixB^S,psi_)=0.d0

      endif

      !
      ! Prohibit ghosts to be supersonic, if so put on sound speed momentum
      ! Also avoid overloading too much, limit to negative sound speed momentum
      !
      do j = ixImin2,ixImax2
        do i = ixBmin1,ixBmax1
          w(i,j,mom(1)) = min(w(i,j,mom(1)), drhobound*dasound)
          w(i,j,mom(1)) = max(w(i,j,mom(1)), -drhobound*dasound)
          w(i,j,mom(2)) = min(w(i,j,mom(2)), drhobound*dasound)
          w(i,j,mom(2)) = max(w(i,j,mom(2)), -drhobound*dasound)
        enddo
      enddo

    !
    ! Right radial boundary is outflow for all quantities. In usr.par 'cont'
    !
    ! Bottom + top angular boundary are 'symm' for scalar quantities, angular
    ! vector quantities are set 'asymm'. This can be auto-achieved by setting
    ! the boundary to 'pole'
    !
    case default
      call mpistop("BC not specified")
    end select

  end subroutine special_bound

!===============================================================================

  subroutine CAK_source(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)
    !
    ! This routine does computations in CGS for ease of interpretation
    ! At start we translate from dimensionless AMRVAC to CGS and at the end the
    ! reverse is being done such that AMRVAC itself keeps working dimensionless
    !
    use mod_global_parameters

    ! Subroutine arguments
    integer, intent(in)    :: ixI^L, ixO^L, iw^LIM
    real(8), intent(in)    :: qdt, qtC, qt
    real(8), intent(in)    :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    real(8), intent(inout) :: w(ixI^S,1:nw)

    ! Local variables
    real(8) :: dvdr_up(ixI^S), dvdr_down(ixI^S), dvdr_cent(ixI^S)
    real(8) :: forw(ixI^S), backw(ixI^S), cent(ixI^S)
    real(8) :: gcak(ixI^S), geff(ixI^S)
    real(8) :: vr(ixI^S), rho(ixI^S), xc(ixI^S)
    real(8) :: beta_fd(ixI^S), fdfac(ixI^S), tdum(ixO^S)
    real(8) :: fac, fac1, fac2
    integer :: i, j, jx^L, hx^L

    !========================================================================
    ! Convert from dimensionless to CGS

    ! Define the time-centred, radial velocity from the radial momentum
    vr(ixI^S) = wCT(ixI^S,mom(1)) / wCT(ixI^S,rho_) * unit_velocity

    ! Time-centred density
    rho(ixI^S) = wCT(ixI^S,rho_) * unit_density

    ! Radial grid coordinate
    xc(ixI^S) = x(ixI^S,1) * unit_length

    !========================================================================

    !
    ! Make new indices covering whole grid by increasing +1 (j) and decreasing
    ! by -1 (h). Special, fancy syntax that AMRVAC understands
    !
    jx^L=ixO^L+kr(1,^D);
    hx^L=ixO^L-kr(1,^D);


    ! NEW TEST
    do i = ixOmin1,ixOmax1
      forw(i,:) = (xc(i,:)-xc(i-1,:))*vr(i+1,:)/((xc(i+1,:)-xc(i,:)) * (xc(i+1,:)-xc(i-1,:)))
      backw(i,:) = -(xc(i+1,:)-xc(i,:))*vr(i-1,:)/((xc(i,:)-xc(i-1,:)) * (xc(i+1,:)-xc(i-1,:)))
      cent(i,:)=(xc(i+1,:)+xc(i-1,:)-2.0d0*xc(i,:))*vr(i,:)/((xc(i,:) - xc(i-1,:))*(xc(i+1,:) - xc(i,:)))
    enddo

    ! Gradient for non-uniform grids according to Sundqvist & Veronis (1970)
    ! forward difference
    ! forw(i^%1ixO^S)  = (xc(i^%1ixO^S) - xc(hx^S)) * vr(jx^S) &
    !                     / ((xc(jx^S) - xc(i^%1ixO^S)) * (xc(jx^S) - xc(hx^S)))
    !
    ! ! backward difference
    ! backw(i^%1ixO^S) = -(xc(jx^S) - xc(i^%1ixO^S)) * vr(hx^S) &
    !                     / ((xc(i^%1ixO^S) - xc(hx^S)) * (xc(jx^S) - xc(hx^S)))
    !
    ! ! central difference
    ! cent(i^%1ixO^S)  = (xc(jx^S)+xc(hx^S)-2.0d0*xc(i^%1ixO^S)) * vr(i^%1ixO^S) &
    !                     / ((xc(i^%1ixO^S) - xc(hx^S)) * (xc(jx^S) - xc(i^%1ixO^S)))

    ! Central differenced dv/dr
    if (iprob==0) then
      ! In CAK this has to be >0, otherwise stagnant flow
      dvdr_cent(ixO^S) = abs(backw(ixO^S) + cent(ixO^S) + forw(ixO^S))
    endif

    if (iprob==1) then
      ! In magnetosphere, we actually require fallback
      dvdr_cent(ixO^S) = max(backw(ixO^S) + cent(ixO^S) + forw(ixO^S), 0.0d0)
    endif

    ! Finite disk factor parameterization (Owocki & Puls 1996)
    beta_fd(ixO^S) = (1.0d0 - vr(ixO^S)/ (xc(ixO^S) * dvdr_cent(ixO^S))) &
                      * (drstar/x(ixO^S,1))**2.0d0

    ! Check the finite disk array and determine finite disk factor
    do j = ixOmin2,ixOmax2
      do i = ixOmin1,ixOmax1

        if (beta_fd(i,j) >= 1.0d0) then
          fdfac(i,j) = 1.0d0/(1.0d0 + alpha)
        else if (beta_fd(i,j) < -1.0d10) then
          fdfac(i,j) = abs(beta_fd(i,j))**alpha / (1.0d0 + alpha)
        else if (abs(beta_fd(i,j)) > 1.0d-3) then
          fdfac(i,j) = (1.0d0 - (1.0d0 - beta_fd(i,j))**(1.0d0 + alpha)) &
                        / (beta_fd(i,j)*(1.0d0 + alpha))
        else
          fdfac(i,j) = 1.0d0 - 0.5d0*alpha*beta_fd(i,j) &
                        * (1.0d0 + 1.0d0/3.0d0 * (1.0d0 - alpha)*beta_fd(i,j))
        end if

        if (fdfac(i,j) < smalldouble) then
          fdfac(i,j) = zero
        else if (fdfac(i,j) > 5.d0) then
          fdfac(i,j) = 1.0d0
        endif

      enddo
    enddo

    ! Calculate CAK line-force
    fac1 = 1.0d0/(1.0d0 - alpha) * kappae * lstar*Qbar/(4.0d0*dpi * const_c)
    fac2 = 1.0d0/(const_c * Qbar * kappae)**alpha
    fac = fac1 * fac2

    gcak(ixO^S) = fac/xc(ixO^S)**2.0d0 * (dvdr_cent(ixO^S)/rho(ixO^S))**alpha

    ! Correct for finite extend stellar disk
    gcak(ixO^S) = gcak(ixO^S) * fdfac(ixO^S)

    ! Fill up the empty CAK array variable
    w(ixO^S,my_gcak) = gcak(ixO^S)

    !========================================================================
    ! Line-force computed and saved in CGS, now return to dimensionless form
    !
    gcak(ixO^S) = gcak(ixO^S) * unit_time**2.0d0 / unit_length

    !========================================================================

    ! Effective gravity computation
    geff(ixO^S) = - dGgrav * dmstar * (1.0d0 - dgammae)/x(ixO^S,1)**2.0d0

    ! Update conservative vars: w = w + qdt*gsource
    w(ixO^S,mom(1)) = w(ixO^S,mom(1)) &
                      + qdt * (gcak(ixO^S) + geff(ixO^S))*wCT(ixO^S,rho_)

    ! Define a new time-step corrected for continuum and line-acceleration
    tdum(ixO^S) = (x(jx^S,1) - x(ixO^S,1)) / (gcak(ixO^S) + geff(ixO^S))
    new_timestep = 0.3d0 * minval(abs(tdum(ixO^S)))**0.5d0

  end subroutine CAK_source

!===============================================================================

  subroutine special_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)
    !
    ! After first iteration assign the new time-step of computation CAK force
    !
    use mod_global_parameters

    ! Subroutine arguments
    integer, intent(in)    :: ixI^L, ixO^L
    real(8), intent(in)    :: dx^D, x(ixI^S,1:ndim)
    real(8), intent(in)    :: w(ixI^S,1:nw)
    real(8), intent(inout) :: dtnew

    if (it >= 1) then
      dtnew = new_timestep
    endif

  end subroutine special_dt

!===============================================================================

  subroutine compute_stats(igrid,level,ixI^L,ixO^L,qt,w,x)
    !
    ! Routine computes the time-averaged statistical quantity <X> via:
    !
    !  <X>_i = <X>_i-1 + 0.5*dt*( X_i + X_i-1 )
    !
    ! where <X> is the average of variable X and i','i-1' are the current and
    ! previous time step respectively
    !
    ! This routine is called after each iteration BUT does only save things
    ! every 'dtsave_dat' time specified in the .par file
    !
    ! To allow the system to relax we set a 'tstat' time in the .par file after
    ! which we start the actual computation of relevant <X> quantities, if this
    ! condition is not met we just assign 0 to the <X> quantities (initialized)
    ! The first time the condition is full-filled we save the current P state
    ! into a temporary array to be used in the next iteration for weighing
    !
    ! Due to the nature of how this routine works in AMRVAC after each
    ! computation of <X> we have to normalize with the appropriate time-weight
    ! 'tnorm' (so that it is included when saving output) but in the next
    ! iteration we have to 'un-normalize' it back in time 'tnorm - dt' to
    ! compute the current timestep <X> again
    !
    use mod_global_parameters

    ! Subroutine arguments
    integer, intent(in)             :: igrid,level,ixI^L,ixO^L
    double precision, intent(in)    :: qt,x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    ! Local variables
    logical :: first_time = .true.
    real(8) :: tnorm

    if (global_time>=dtstat) then
      !
      ! Convert conservative variables to primitives (back conversion at end)
      ! This means that mom(X) below is now the velocity in direction X
      !
      call mhd_to_primitive(ixI^L,ixO^L,w,x)

      if (first_time .and. (.not. resume_previous_run)) then
        w(ixO^S,my_tmp1) = w(ixO^S,rho_)
        w(ixO^S,my_tmp2) = w(ixO^S,rho_)**2.0d0
        w(ixO^S,my_tmp3) = w(ixO^S,mom(1))
        w(ixO^S,my_tmp4) = w(ixO^S,mom(1))**2.0d0
        w(ixO^S,my_tmp5) = w(ixO^S,rho_) * w(ixO^S,mom(1))
        w(ixO^S,my_tmp6) = w(ixO^S,rho_)**2.0d0 * w(ixO^S,mom(1))
        w(ixO^S,my_tmp7) = w(ixO^S,mom(2))
        w(ixO^S,my_tmp8) = w(ixO^S,mom(2))**2.0d0
        w(ixO^S,my_tmp9) = w(ixO^S,rho_)**2.0d0 * w(ixO^S,mom(2))
        first_time  = .false.
      else
        ! Time weight
        tnorm = global_time - dtstat

        ! Average density
        w(ixO^S,my_rhoav) = (w(ixO^S,my_rhoav)*(tnorm-dt) &
                                 + 0.5*(w(ixO^S,rho_) + w(ixO^S,my_tmp1))*dt)/tnorm
        w(ixO^S,my_tmp1)  = w(ixO^S,rho_)

        ! Average density squared
        w(ixO^S,my_rho2av) = (w(ixO^S,my_rho2av)*(tnorm-dt) &
                            + 0.5*(w(ixO^S,rho_)**2.0d0 + w(ixO^S,my_tmp2))*dt)/tnorm
        w(ixO^S,my_tmp2)   = w(ixO^S,rho_)**2.0d0

        ! Average radial velocity
        w(ixO^S,my_vrav) = (w(ixO^S,my_vrav)*(tnorm-dt) &
                           + 0.5*(w(ixO^S,mom(1)) + w(ixO^S,my_tmp3))*dt)/tnorm
        w(ixO^S,my_tmp3) = w(ixO^S,mom(1))

        ! Average radial velocity squared
        w(ixO^S,my_vr2av) = (w(ixO^S,my_vr2av)*(tnorm-dt) &
                           + 0.5*(w(ixO^S,mom(1))**2.0d0 + w(ixO^S,my_tmp4))*dt)/tnorm
        w(ixO^S,my_tmp4) = w(ixO^S,mom(1))**2.0d0

        ! Average radial momentum density (correlation density-velocity)
        w(ixO^S,my_rhovrav) = (w(ixO^S,my_rhovrav)*(tnorm-dt) &
                             + 0.5*(w(ixO^S,rho_) * w(ixO^S,mom(1)) + w(ixO^S,my_tmp5))*dt)/tnorm
        w(ixO^S,my_tmp5) = w(ixO^S,rho_) * w(ixO^S,mom(1))

        ! Average weighted clump radial velocity
        w(ixO^S,my_rho2vrav) = (w(ixO^S,my_rho2vrav)*(tnorm-dt) &
                              + 0.5*(w(ixO^S,rho_)**2.0d0 * w(ixO^S,mom(1)) + w(ixO^S,my_tmp6))*dt)/tnorm
        w(ixO^S,my_tmp6) = w(ixO^S,rho_)**2.0d0 * w(ixO^S,mom(1))

        ! Average polar velocity
        w(ixO^S,my_vpolav) = (w(ixO^S,my_vpolav)*(tnorm-dt) &
                           + 0.5*(w(ixO^S,mom(2)) + w(ixO^S,my_tmp7))*dt)/tnorm
        w(ixO^S,my_tmp7) = w(ixO^S,mom(2))

        ! Average polar velocity squared
        w(ixO^S,my_vpol2av)  = (w(ixO^S,my_vpol2av)*(tnorm-dt) &
                              + 0.5*(w(ixO^S,mom(2))**2.0d0 + w(ixO^S,my_tmp8))*dt)/tnorm
        w(ixO^S,my_tmp8) = w(ixO^S,mom(2))**2.0d0

        ! Average weighted clump polar velocity
        w(ixO^S,my_rho2vpolav) = (w(ixO^S,my_rho2vpolav)*(tnorm-dt) &
                                + 0.5*(w(ixO^S,rho_)**2.0d0 * w(ixO^S,mom(2)) + w(ixO^S,my_tmp9))*dt)/tnorm
        w(ixO^S,my_tmp9) = w(ixO^S,rho_)**2.0d0 * w(ixO^S,mom(2))

      endif

      call mhd_to_conserved(ixI^L,ixO^L,w,x)

    endif

  end subroutine compute_stats

!===============================================================================

  subroutine make_dimless_vars()
    !
    ! Normalize quantities in use to unit quantities defined and computed
    ! These quantities are actually used by AMRVAC in its computations!
    !
    use mod_global_parameters

    ! From the AMRVAC unit vars compute some extra relevant for us
    my_unit_ggrav = unit_density * unit_time**2.0d0
    my_unit_lum   = unit_density * unit_length**5.0d0 / unit_time**3.0d0
    my_unit_mass  = unit_density * unit_length**3.0d0

    drhobound = rhobound/unit_density
    dlstar    = lstar/my_unit_lum
    dmstar    = mstar/my_unit_mass
    drstar    = rstar/unit_length
    dtwind    = twind/unit_temperature
    dbpole    = bpole/unit_magneticfield
    dmdot     = mdot/(my_unit_mass/unit_time)
    dasound   = asound/unit_velocity
    dclight   = const_c/unit_velocity
    dvesc     = vesc/unit_velocity
    dvinf     = dvesc * (alpha/(1.0d0 - alpha))**0.5d0
    dkappae   = kappae * unit_density * unit_length
    dGgrav    = Ggrav * my_unit_ggrav
    dgammae   = dkappae * dlstar/(4.d0*dpi * dGgrav * dmstar * dclight)
    detaconf  = ((dbpole/2.0d0)**2.0d0 * drstar**2.0d0)/(dmdot * dvinf)
    dtstat    = tstat/unit_time

  end subroutine make_dimless_vars

!===============================================================================

  subroutine set_extravar_output(ixI^L,ixO^L,w,x,normconv)
    !
    ! Computes and stores additional variables of interest
    !
    use mod_global_parameters

    ! Subroutine arguments
    integer, intent(in) :: ixI^L,ixO^L
    real(8), intent(in) :: x(ixI^S,1:ndim)
    real(8)             :: w(ixI^S,nw+nwauxio)
    real(8)             :: normconv(0:nw+nwauxio)

    ! Local variable
    real(8) :: divbboy(ixI^S)

    ! Output the Alfven speed by summing squared Bfield in each direction
    if (B0field) then
      w(ixO^S,nw+1) = dsqrt(sum((w(ixO^S,mag(:)) + &
                            block%B0(ixO^S,:,0))**2, dim=ndim+1)/w(ixO^S,rho_))
    else
      w(ixO^S,nw+1) = dsqrt(sum(w(ixO^S,mag(:))**2, dim=ndim+1)/w(ixO^S,rho_))
    endif

    !
    ! Output divB using AMRVAC routine that computes it for us
    ! For Tanaka splitting this will be div B1 (not possible to get full)
    !
    call get_divb(w,ixI^L,ixO^L,divbboy)
    w(ixO^S,nw+2) = divbboy(ixO^S)

  end subroutine set_extravar_output

!===============================================================================

  subroutine set_extravarnames_output(varnames)
    !
    ! Newly added variables: Alfven velocity, divergence B
    !
    character(len=*) :: varnames
    varnames = 'vAlf divb'

  end subroutine set_extravarnames_output

!===============================================================================

  subroutine make_dipoleboy(ixI^L,ixO^L,x,wB0)
    !
    ! Add a steady (time-independent) potential dipole background field with
    ! polar field strength set by dimensionless 'dbpole'
    !
    use mod_global_parameters

    ! Subroutine arguments
    integer, intent(in)           :: ixI^L,ixO^L
    double precision, intent(in)  :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: wB0(ixI^S,1:ndir)

    wB0(ixI^S,1) = dbpole * (drstar/x(ixI^S,1))**3.0d0 * dcos(x(ixI^S,2))
    wB0(ixI^S,2) = 0.5d0*dbpole * (drstar/x(ixI^S,1))**3.0d0 * dsin(x(ixI^S,2))
    wB0(ixI^S,3) = 0.0d0

  end subroutine make_dipoleboy

end module mod_usr
