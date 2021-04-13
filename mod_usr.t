!===============================================================================
! 2D wind module for launching a CAK line-driven wind from the stellar surface
! into the magnetosphere of a magnetic massive star
!
! Setup mimicked as in ud-Doula & Owocki (2002). Testbench, similar separate
! program will be used for incorporating the LDI into the magnetosphere.
!
! Coded up by Flo for his KU Leuven PhD thesis 2018-2022
!-------------------------------------------------------------------------------
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
! (October 2020) -- Flo
!   > added rigid body rotation to study centrifugal magnetospheres, set in
!     .par file with dimensionless 'Wrot' parameter (ud-Doula et al. 2006)
!   > included option with parameter 'imag' to select setup in terms of field
!     strength (imag>0) or wind confinement (imag<0)
!
! (December 2020) -- Flo
!   > removed rho2vrad_av and rho2vpol_av variables from statistics
!
! (January 2021) -- Flo
!   > do not start from relaxed CAK model anymore, instead start directly from
!     beta law with Bfield (c.q. Oct 2019); adds more flexibility in multi-d
!   > removed accordingly the 'iprob' option
!
! (February 2021) -- Flo
!   > starting from beta law in high confinement leads to crash, removed mods
!     from January; instead now 1D relaxed CAK read in
!   > put effective gravity force in usr_gravity and removed from usr_source
!   > determination radiation timestep (only CAK line force now) in special_dt
!     by using gcak slot of nwextra in w-array
!   > rewriting statistics to get rid off tmp variables, AMRVAC keeps in memory
!     the w-array from previous iteration, and now called in my computation
!   > call statistics at end iteration (usr_process_adv_grid) instead of at
!     start iteration to have access to the old state 'pso'
!
! (April 2021) -- Flo
!   > no more pso-state in statistics, now only current since averaging is each
!     iteration such that weight with previous iteration unnecessary
!   > prohibit pure adiabatic cooling as it crashes, but is also unrealistic
!===============================================================================

module mod_usr

  ! Include a physics module
  use mod_mhd

  implicit none

  ! The usual suspects
  real(8), parameter :: msun=1.989d33, lsun=3.827d33, rsun=6.96d10
  real(8), parameter :: Ggrav=6.67d-8, kappae=0.34d0

  ! Unit quantities that are handy: gravitational constant, luminosity, mass
  real(8) :: my_unit_ggrav, my_unit_lum, my_unit_mass

  ! Extra input parameters:
  real(8)           :: lstar, mstar, rstar, rhobound, twind, imag, alpha, Qbar
  real(8)           :: tstat, Wrot
  character(len=99) :: cakfile

  ! Additionally useful stellar and wind parameters:
  !   Eddington gamma, escape speed, CAK + fd mass-loss rate, terminal wind
  !   speed, sound speed, log(g), eff. log(g), scale height, mean mol. weight,
  !   polar magnetic field, wind magnetic confinement parameter, Alfven +
  !   Kepler + escape radius, equatorial rotation + critical rotation speed
  real(8) :: gammae, vesc, mdot, mdotfd, vinf, asound, logg, logge, heff
  real(8) :: mumol, bpole, etastar, ralf, rkep, resc, vrot, vrotc

  ! Dimensionless variables of relevant variables above
  real(8) :: dlstar, dmstar, drstar, dbpole, drhobound, dtwind, dkappae
  real(8) :: dvesc, dvinf, dmdot, dasound, dclight, dGgrav, dgammae, detaconf
  real(8) :: dtstat, dvrot

  ! Additional names for wind and statistical variables
  integer :: my_gcak, my_rhoav, my_rho2av, my_vrav, my_vr2av, my_rhovrav
  integer :: my_vpolav, my_vpol2av, my_tav

  ! Arrays required to read and store 1D profile from file
  real(8), allocatable :: woneblock(:,:), xoneblock(:,:)

contains

  !======================================================================
  ! This routine should set user methods, and activate the physics module
  !======================================================================
  subroutine usr_init()

    call set_coordinate_system("spherical_2.5D")
    call usr_params_read(par_files)

    !
    ! Choose independent normalization units, only 3 have to be specified:
    !     (length,temp,ndens) or (length,vel,ndens)
    ! Numberdensity chosen such that unit density becomes boundary density
    !
    unit_length        = rstar                                      ! cm
    unit_temperature   = twind                                      ! K
    unit_numberdensity = rhobound/((1.d0+4.d0*He_abundance)*mp_cgs) ! g cm^-3

    call MHD_activate()

    usr_set_parameters   => initglobaldata_usr
    usr_init_one_grid    => initial_conditions
    usr_special_bc       => special_bound
    usr_gravity          => effective_gravity
    usr_source           => line_force
    usr_get_dt           => special_dt
    usr_process_adv_grid => compute_stats
    usr_aux_output       => set_extravar_output
    usr_add_aux_names    => set_extravarnames_output
    usr_set_B0           => make_dipoleboy

    my_gcak    = var_set_extravar("gcak", "gcak")
    my_rhoav   = var_set_extravar("rho_av", "rho_av")
    my_rho2av  = var_set_extravar("rho2_av", "rho2_av")
    my_vrav    = var_set_extravar("vrad_av", "vrad_av")
    my_vpolav  = var_set_extravar("vtheta_av", "vtheta_av")
    my_vr2av   = var_set_extravar("vrad2_av", "vrad2_av")
    my_vpol2av = var_set_extravar("vtheta2_av", "vtheta2_av")
    my_rhovrav = var_set_extravar("rho_vrad_av", "rho_vrad_av")

    if (mhd_energy) my_tav = var_set_extravar("twind_av", "twind_av")

  end subroutine usr_init

  !========================================================
  ! Read in the usr.par file with the problem specific list
  !========================================================
  subroutine usr_params_read(files)

    ! Subroutine argument
    character(len=*), intent(in) :: files(:)

    ! Local variable
    integer :: n

    namelist /star_list/ mstar, lstar, rstar, twind, imag, rhobound, alpha, &
                          Qbar, tstat, Wrot, cakfile

    do n = 1,size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, star_list, end=111)
       111 close(unitpar)
    end do

    ! Scale to cgs units
    lstar = lstar * lsun
    mstar = mstar * msun
    rstar = rstar * rsun

  end subroutine usr_params_read

  !====================================================================
  ! Compute some quantities of interest (in CGS) before making unitless
  !====================================================================
  subroutine initglobaldata_usr

    ! Stellar structure
    gammae = kappae * lstar/(4.d0*dpi * Ggrav * mstar * const_c)
    logg   = log10(Ggrav * mstar/rstar**2.0d0)
    logge  = logg  + log10(1.0d0 - gammae)
    mumol  = (1.0d0 + 4.0d0*He_abundance)/(2.0d0 + 3.0d0*He_abundance)
    asound = sqrt(twind * kB_cgs/(mumol * mp_cgs))
    heff   = asound**2.0d0 / 10.0d0**logge
    vrotc  = sqrt(Ggrav * mstar*(1.0d0 - gammae)/rstar)
    vrot   = vrotc * Wrot

    ! Wind quantities in CAK theory
    vesc   = sqrt(2.0d0 * Ggrav * mstar*(1.0d0 - gammae)/rstar)
    vinf   = vesc * sqrt(alpha/(1.0d0 - alpha))
    mdot   = lstar/const_c**2.0d0 * alpha/(1.0d0 - alpha) &
               * (Qbar * gammae/(1.0d0 - gammae))**((1.0d0 - alpha)/alpha)
    mdotfd = mdot/(1.0d0 + alpha)**(1.0d0/alpha)

    ! Bpole given and etastar computed or vice versa
    if (imag > 0.0d0) then
      bpole   = imag
      etastar = ((bpole/2.0d0)**2.0d0 * rstar**2.0d0)/(mdot * vinf)
    else
      etastar = -imag
      bpole   = 2.0d0 * sqrt(mdot * vinf * etastar/rstar**2.0d0)
    endif

    ! Compute Alfven, Kepler, and escape radius
    ralf = 1.0d0 + (etastar + 0.25d0)**0.25d0 - 0.25d0**0.25d0
    rkep = Wrot**(-2.0d0/3.0d0)
    resc = 2.0d0**(1.0d0/3.0d0) * rkep

    call make_dimless_and_log_vars()

    if (.not. resume_previous_run) call read_initial_oned_cak(cakfile)

    if (typedivbfix == 'ct') then
      call mpistop('CT disabled. Gives strange results for this problem.')
    endif

    if (mhd_energy .and. (.not. mhd_radiative_cooling)) then
      call mpistop('Pure adiabatic cooling not supported.')
    endif

  end subroutine initglobaldata_usr

  !==========================================================================
  ! Initial conditions start from spherically symmetric 1D relaxed CAK wind
  ! and dipole field is set here or in usr_set_B0 when doing Tanaka splitting
  !==========================================================================
  subroutine initial_conditions(ixI^L,ixO^L,w,x)

    ! Subroutine arguments
    integer, intent(in)    :: ixI^L, ixO^L
    real(8), intent(in)    :: x(ixI^S,1:ndim)
    real(8), intent(inout) :: w(ixI^S,1:nw)

    ! Local variables
    integer :: ir, irb

    ! Get density and radial velocity from broadcasted 1D CAK profile
    do ir = ixOmin1,ixOmax1
      ! Find correct block index to get correct value
      irb = minloc(abs(x(ir,nghostcells,1) - xoneblock(:,1)),1)

      w(ir,:,rho_)   = woneblock(irb,1)
      w(ir,:,mom(1)) = woneblock(irb,2)
    enddo

    ! No polar velocity, in azimuth rigid rotation for full grid (stabilizing)
    w(ixO^S,mom(2)) = 0.0d0
    w(ixO^S,mom(3)) = dvrot/drstar * x(ixO^S,1) * sin(x(ixO^S,2))

    ! Setup dipole magnetic field based on Tanaka splitting or regular
    if (B0field) then
      w(ixO^S,mag(:)) = 0.0d0
    else
      w(ixO^S,mag(1)) = dbpole * (drstar/x(ixO^S,1))**3.0d0 * cos(x(ixO^S,2))

      w(ixO^S,mag(2)) = 0.5d0 * dbpole &
                        * (drstar/x(ixO^S,1))**3.0d0 * sin(x(ixO^S,2))

      w(ixO^S,mag(3)) = 0.0d0
    endif

    if (mhd_energy) w(ixO^S,p_) = w(ixO^S,rho_)

    ! If using Dedner+(2002) divergence cleaning
    if (mhd_glm) w(ixO^S,psi_) = 0.0d0

    call mhd_to_conserved(ixI^L,ixO^L,w,x)

    ! Initialise extra vars at 0 since some compilers fill with crap
    w(ixO^S,nw-nwextra+1:nw) = 0.0d0

  end subroutine initial_conditions

  !=============================================================================
  ! Special user boundary conditions at inner (+ possibly outer) radial boundary
  !=============================================================================
  subroutine special_bound(qt,ixI^L,ixB^L,iB,w,x)

    ! Subroutine arguments
    integer, intent(in)    :: ixI^L, ixB^L, iB
    real(8), intent(in)    :: qt, x(ixI^S,1:ndim)
    real(8), intent(inout) :: w(ixI^S,1:nw)

    ! Local variable
    integer :: i

    select case (iB)
    case(1)

      call mhd_to_primitive(ixI^L,ixI^L,w,x)

      w(ixB^S,rho_) = drhobound

      ! vr (only 2nd order accurate constant slope extrapolation in 1st ghost)
      do i = ixBmax1,ixBmin1,-1
        if (i == ixBmax1) then
          w(i^%1ixB^S,mom(1)) = 1.0d0/3.0d0 * (- w(i+2^%1ixB^S,mom(1)) &
                                               + 4.0d0 * w(i+1^%1ixB^S,mom(1)))
        else
          w(i^%1ixB^S,mom(1)) = w(i+1^%1ixB^S,mom(1))
        endif
      enddo

      ! Prohibit radial velocity ghosts to be supersonic, and avoid overloading
      w(ixB^S,mom(1)) = min(w(ixB^S,mom(1)), dasound)
      w(ixB^S,mom(1)) = max(w(ixB^S,mom(1)), -dasound)

      w(ixB^S,mom(2)) = 0.0d0
      w(ixB^S,mom(3)) = dvrot * sin(x(ixB^S,2))

      !=================
      ! Tanaka splitting
      !=================
      if (B0field) then
        do i = ixBmax1,ixBmin1,-1
          ! r*r*(B0r + delta Br) = constant
          w(i^%1ixB^S,mag(1)) = (dbpole * cos(x(ixBmax1+1^%1ixB^S,2)) &
                                  + w(ixBmax1+1^%1ixB^S,mag(1))) &
                                 * (drstar / x(i^%1ixB^S,1))**2.0d0 &
                                 - block%B0(i^%1ixB^S,1,0)

          ! delta Btheta
          w(i^%1ixB^S,mag(2)) = 1.0d0/3.0d0 * (- w(i+2^%1ixB^S,mag(2)) &
                                               + 4.0d0 * w(i+1^%1ixB^S,mag(2)))

          ! delta Bphi
          w(i^%1ixB^S,mag(3)) = 1.0d0/3.0d0 * (- w(i+2^%1ixB^S,mag(3)) &
                                               + 4.0d0 * w(i+1^%1ixB^S,mag(3)))
        enddo
      !========
      ! Regular
      !========
      else
        do i = ixBmax1,ixBmin1,-1
          ! r*r*Br = constant
          w(i^%1ixB^S,mag(1)) = dbpole * cos(x(ixBmax1+1^%1ixB^S,2)) &
                                * (drstar / x(i^%1ixB^S,1))**2.0d0

          ! Btheta
          w(i^%1ixB^S,mag(2)) = 1.0d0/3.0d0 * (- w(i+2^%1ixB^S,mag(2)) &
                                               + 4.0d0 * w(i+1^%1ixB^S,mag(2)))

          ! Bphi
          w(i^%1ixB^S,mag(3)) = 1.0d0/3.0d0 * (- w(i+2^%1ixB^S,mag(3)) &
                                               + 4.0d0 * w(i+1^%1ixB^S,mag(3)))
        enddo
      endif

      ! Density is fixed, so per ideal gas law also the pressure
      if (mhd_energy) w(ixB^S,p_) = w(ixB^S,rho_)

      if (mhd_glm) w(ixB^S,psi_) = 0.0d0

      call mhd_to_conserved(ixI^L,ixI^L,w,x)

    case(2)

      do i = ixBmin1,ixBmax1
        ! r*r*rho = constant
        w(i^%1ixB^S,rho_) = w(ixBmin1-1^%1ixB^S,rho_) &
                            * (x(ixBmin1-1^%1ixB^S,1) / x(i^%1ixB^S,1))**2.0d0

        ! r*r*rho*vr = constant
        w(i^%1ixB^S,mom(1)) = w(ixBmin1-1^%1ixB^S,mom(1)) &
                              * (x(ixBmin1-1^%1ixB^S,1)/x(i^%1ixB^S,1))**2.0d0

        ! rho*vtheta = constant
        w(i^%1ixB^S,mom(2)) = w(ixBmin1-1^%1ixB^S,mom(2))

        ! r*vphi = constant
        w(i^%1ixB^S,mom(3)) = &
                (w(ixBmin1-1^%1ixB^S,mom(3)) / w(ixBmin1-1^%1ixB^S,rho_)) &
                * (x(ixBmin1-1^%1ixB^S,1) / x(i^%1ixB^S,1)) * w(i^%1ixB^S,rho_)

        ! energy = constant
        if (mhd_energy) w(i^%1ixB^S,e_) = w(ixBmin1-1^%1ixB^S,e_)

        !=================
        ! Tanaka splitting
        !=================
        if (B0field) then
          ! r*r*(B0r + delta Br) = constant
          w(i^%1ixB^S,mag(1)) = ( dbpole * cos(x(ixBmin1-1^%1ixB^S,2)) &
                                  * (drstar/x(ixBmin1-1^%1ixB^S,1))**3.0d0 &
                                  + w(ixBmin1-1^%1ixB^S,mag(1)) &
                                ) &
                                * (x(ixBmin1-1^%1ixB^S,1)/x(i^%1ixB^S,1))**2.0d0 &
                                - block%B0(i^%1ixB^S,1,0)

          ! (B0theta + delta Btheta) = constant
          w(i^%1ixB^S,mag(2)) = ( 0.5d0*dbpole * sin(x(ixBmin1-1^%1ixB^S,2)) &
                                  * (drstar/x(ixBmin1-1^%1ixB^S,1))**3.0d0 &
                                  + w(ixBmin1-1^%1ixB^S,mag(2)) &
                                ) &
                                - block%B0(i^%1ixB^S,2,0)

          ! r*(B0phi + delta Bphi) = constant
          w(i^%1ixB^S,mag(3)) = w(ixBmin1-1^%1ixB^S,mag(3)) &
                                * x(ixBmin1-1^%1ixB^S,1) / x(i^%1ixB^S,1) &
                                - block%B0(i^%1ixB^S,3,0)
        !========
        ! Regular
        !========
        else
          ! r*r*Br = constant
          w(i^%1ixB^S,mag(1)) = w(ixBmin1-1^%1ixB^S,mag(1)) &
                                * (x(ixBmin1-1^%1ixB^S,1)/x(i^%1ixB^S,1))**2.0d0

          ! Btheta = constant
          w(i^%1ixB^S,mag(2)) = w(ixBmin1-1^%1ixB^S,mag(2))

          ! r*Bphi = constant
          w(i^%1ixB^S,mag(3)) = w(ixBmin1-1^%1ixB^S,mag(3)) &
                                * x(ixBmin1-1^%1ixB^S,1) / x(i^%1ixB^S,1)
        endif
      enddo

      if (mhd_glm) w(ixB^S,psi_) = 0.0d0

    case default
      call mpistop("BC not specified")
    end select

  end subroutine special_bound

  !=======================================================================
  ! Extra source using the analytical CAK line force in Gayley's formalism
  !=======================================================================
  subroutine line_force(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)

    ! Subroutine arguments
    integer, intent(in)    :: ixI^L, ixO^L, iw^LIM
    real(8), intent(in)    :: qdt, qtC, qt
    real(8), intent(in)    :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    real(8), intent(inout) :: w(ixI^S,1:nw)

    ! Local variables
    real(8) :: vr(ixI^S), rho(ixI^S)
    real(8) :: dvdr_up(ixO^S), dvdr_down(ixO^S), dvdr_cent(ixO^S), dvdr(ixO^S)
    real(8) :: gcak(ixO^S), beta_fd(ixO^S), fdfac(ixO^S)
    real(8) :: fac, fac1, fac2
    integer :: jx^L, hx^L

    ! Define time-centred, radial velocity from the radial momentum and density
    vr(ixI^S)  = wCT(ixI^S,mom(1)) / wCT(ixI^S,rho_)
    rho(ixI^S) = wCT(ixI^S,rho_)

    ! Index +1 (j) and index -1 (h) in radial direction; kr(dir,dim)=1, dir=dim
    jx^L=ixO^L+kr(1,^D);
    hx^L=ixO^L-kr(1,^D);

    ! Get dv/dr on non-uniform grid according to Sundqvist & Veronis (1970)
    ! Forward difference
    dvdr_up(ixO^S) = (x(ixO^S,1) - x(hx^S,1)) * vr(jx^S) &
                     / ((x(jx^S,1) - x(ixO^S,1)) * (x(jx^S,1) - x(hx^S,1)))

    ! Backward difference
    dvdr_down(ixO^S) = -(x(jx^S,1) - x(ixO^S,1)) * vr(hx^S) &
                        / ((x(ixO^S,1) - x(hx^S,1)) * (x(jx^S,1) - x(hx^S,1)))

    ! Central difference
    dvdr_cent(ixO^S) = (x(jx^S,1) + x(hx^S,1) - 2.0d0*x(ixO^S,1)) * vr(ixO^S) &
                       / ((x(ixO^S,1) - x(hx^S,1)) * (x(jx^S,1) - x(ixO^S,1)))

    ! Total gradient with fallback requirement check
    dvdr(ixO^S) = dvdr_down(ixO^S) + dvdr_cent(ixO^S) + dvdr_up(ixO^S)
    dvdr(ixO^S) = max(dvdr(ixO^S), 0.0d0)

    ! Finite disk factor parameterisation (Owocki & Puls 1996)
    beta_fd(ixO^S) = ( 1.0d0 - vr(ixO^S)/(x(ixO^S,1) * dvdr(ixO^S)) ) &
                      * (drstar/x(ixO^S,1))**2.0d0

    ! Check the finite disk array and determine finite disk factor
    where (beta_fd >= 1.0d0)
      fdfac = 1.0d0/(1.0d0 + alpha)
    elsewhere (beta_fd < -1.0d10)
      fdfac = abs(beta_fd)**alpha / (1.0d0 + alpha)
    elsewhere (abs(beta_fd) > 1.0d-3)
      fdfac = (1.0d0 - (1.0d0 - beta_fd)**(1.0d0 + alpha)) &
                  / (beta_fd*(1.0d0 + alpha))
    elsewhere
      fdfac = 1.0d0 - 0.5d0*alpha*beta_fd &
                      * (1.0d0 + 1.0d0/3.0d0 * (1.0d0 - alpha)*beta_fd)
    endwhere

    where (fdfac < smalldouble)
      fdfac = 0.0d0
    elsewhere (fdfac > 5.0d0)
      fdfac = 1.0d0
    endwhere

    ! Calculate CAK line-force and correct for finite extend stellar disk
    fac1 = 1.0d0/(1.0d0 - alpha) * dkappae * dlstar*Qbar/(4.0d0*dpi * dclight)
    fac2 = 1.0d0/(dclight * Qbar * dkappae)**alpha
    fac  = fac1 * fac2

    gcak(ixO^S) = fac/x(ixO^S,1)**2.0d0 * (dvdr(ixO^S)/rho(ixO^S))**alpha
    gcak(ixO^S) = gcak(ixO^S) * fdfac(ixO^S)

    ! Fill the CAK slot variable
    w(ixO^S,my_gcak) = gcak(ixO^S)

    ! Update conservative vars: w = w + qdt*gsource
    w(ixO^S,mom(1)) = w(ixO^S,mom(1)) + qdt * gcak(ixO^S) * wCT(ixO^S,rho_)

    ! Update total energy e = e + vCT*rhoCT*qdt*gsource
    if (mhd_energy) then
      w(ixO^S,e_) = w(ixO^S,e_) + qdt * gcak(ixO^S) * wCT(ixO^S,mom(1))
    endif

  end subroutine line_force

  !========================================================================
  ! After first iteration the usr_source routine has been called, take now
  ! also timestep from CAK line force into account
  !========================================================================
  subroutine special_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)

    ! Subroutine arguments
    integer, intent(in)    :: ixI^L, ixO^L
    real(8), intent(in)    :: dx^D, x(ixI^S,1:ndim)
    real(8), intent(in)    :: w(ixI^S,1:nw)
    real(8), intent(inout) :: dtnew

    ! Local variables
    real(8) :: tdum(ixO^S), dt_cak

    ! Get dt from line force that is saved in the w-array in nwextra slot
    tdum(ixO^S) = sqrt( block%dx(ixO^S,1) / abs(w(ixO^S,my_gcak)) )
    dt_cak      = courantpar * minval(tdum(ixO^S))

    if (it >= 1) then
      dtnew = min(dtnew,dt_cak)
    endif

  end subroutine special_dt

  !===================================================================
  ! Combine stellar gravity and continuum electron scattering into an
  ! effective gravity using Eddington's gamma
  !===================================================================
  subroutine effective_gravity(ixI^L,ixO^L,wCT,x,gravity_field)

    ! Subroutine arguments
    integer, intent(in)  :: ixI^L, ixO^L
    real(8), intent(in)  :: x(ixI^S,1:ndim)
    real(8), intent(in)  :: wCT(ixI^S,1:nw)
    real(8), intent(out) :: gravity_field(ixI^S,ndim)

    gravity_field(ixO^S,:) = 0.0d0

    ! Only in radial direction
    gravity_field(ixO^S,1) = -dGgrav*dmstar * (1.0d0-dgammae)/x(ixO^S,1)**2.0d0

  end subroutine effective_gravity

  !==========================================================================
  ! Routine computes the time-averaged statistical quantity <X> via:
  !   <X>_i = <X>_i-1 + dt * X_i
  ! where <X> is the average of variable X and i','i-1' are the current and
  ! previous timestep respectively
  !==========================================================================
  subroutine compute_stats(igrid,level,ixI^L,ixO^L,qt,w,x)

    ! Subroutine arguments
    integer, intent(in)    :: igrid, level, ixI^L, ixO^L
    real(8), intent(in)    :: qt, x(ixI^S,1:ndim)
    real(8), intent(inout) :: w(ixI^S,1:nw)

    ! Local variables
    real(8) :: tnormp, tnormc

    ! Note: qt is just a placeholder for the 'global_time' variable
    if (qt > dtstat) then

      call mhd_to_primitive(ixI^L,ixO^L,w,x)

      ! Current ^(n+1) and previous ^(n) timestep normalisation weigths
      tnormc = qt + dt - dtstat
      tnormp = qt - dtstat

      ! Average density
      w(ixO^S,my_rhoav) = w(ixO^S,my_rhoav)*tnormp
      w(ixO^S,my_rhoav) = w(ixO^S,my_rhoav) + dt * w(ixO^S,rho_)
      w(ixO^S,my_rhoav) = w(ixO^S,my_rhoav)/tnormc

      ! Average density squared
      w(ixO^S,my_rho2av) = w(ixO^S,my_rho2av)*tnormp
      w(ixO^S,my_rho2av) = w(ixO^S,my_rho2av) + dt * w(ixO^S,rho_)**2.0d0
      w(ixO^S,my_rho2av) = w(ixO^S,my_rho2av)/tnormc

      ! Average radial velocity
      w(ixO^S,my_vrav) = w(ixO^S,my_vrav)*tnormp
      w(ixO^S,my_vrav) = w(ixO^S,my_vrav) + dt * w(ixO^S,mom(1))
      w(ixO^S,my_vrav) = w(ixO^S,my_vrav)/tnormc

      ! Average radial velocity squared
      w(ixO^S,my_vr2av) = w(ixO^S,my_vr2av)*tnormp
      w(ixO^S,my_vr2av) = w(ixO^S,my_vr2av) + dt * w(ixO^S,mom(1))**2.0d0
      w(ixO^S,my_vr2av) = w(ixO^S,my_vr2av)/tnormc

      ! Average radial momentum density (correlation density-velocity)
      w(ixO^S,my_rhovrav) = w(ixO^S,my_rhovrav)*tnormp
      w(ixO^S,my_rhovrav) = w(ixO^S,my_rhovrav) + dt * w(ixO^S,rho_)*w(ixO^S,mom(1))
      w(ixO^S,my_rhovrav) = w(ixO^S,my_rhovrav)/tnormc

      ! Average polar velocity
      w(ixO^S,my_vpolav) = w(ixO^S,my_vpolav)*tnormp
      w(ixO^S,my_vpolav) = w(ixO^S,my_vpolav) + dt * w(ixO^S,mom(2))
      w(ixO^S,my_vpolav) = w(ixO^S,my_vpolav)/tnormc

      ! Average polar velocity squared
      w(ixO^S,my_vpol2av) = w(ixO^S,my_vpol2av)*tnormp
      w(ixO^S,my_vpol2av) = w(ixO^S,my_vpol2av) + dt * w(ixO^S,mom(2))**2.0d0
      w(ixO^S,my_vpol2av) = w(ixO^S,my_vpol2av)/tnormc

      ! Average wind temperature
      if (mhd_energy) then
        w(ixO^S,my_tav) = w(ixO^S,my_tav)*tnormp
        w(ixO^S,my_tav) = w(ixO^S,my_tav) + dt * w(ixO^S,p_)/w(ixO^S,rho_)
        w(ixO^S,my_tav) = w(ixO^S,my_tav)/tnormc
      endif

      call mhd_to_conserved(ixI^L,ixO^L,w,x)
    endif

  end subroutine compute_stats

  !======================================================================
  ! Computes and stores additional variables of interest in convert stage
  !======================================================================
  subroutine set_extravar_output(ixI^L,ixO^L,w,x,normconv)

    ! Subroutine arguments
    integer, intent(in) :: ixI^L, ixO^L
    real(8), intent(in) :: x(ixI^S,1:ndim)
    real(8)             :: w(ixI^S,nw+nwauxio)
    real(8)             :: normconv(0:nw+nwauxio)

    ! Local variable
    real(8) :: divbboy(ixI^S)

    ! Output the Alfven speed by summing squared Bfield in each direction
    if (B0field) then
      w(ixO^S,nw+1) = sqrt(sum((w(ixO^S,mag(:)) + &
                           block%B0(ixO^S,:,0))**2, dim=ndim+1)/w(ixO^S,rho_))
    else
      w(ixO^S,nw+1) = sqrt(sum(w(ixO^S,mag(:))**2, dim=ndim+1)/w(ixO^S,rho_))
    endif

    !
    ! Output divB using AMRVAC routine that computes it for us
    ! For Tanaka splitting this will be div B1 (not possible to get full)
    !
    call get_divb(w,ixI^L,ixO^L,divbboy)
    w(ixO^S,nw+2) = divbboy(ixO^S)

  end subroutine set_extravar_output

  !==================================================================
  ! Additional auxiliary io variables: Alfven velocity, divergence B
  !==================================================================
  subroutine set_extravarnames_output(varnames)

    character(len=*) :: varnames
    varnames = 'vAlf divb'

  end subroutine set_extravarnames_output

  !=======================================================================
  ! Add a steady (time-independent) potential dipole background field with
  ! polar field strength set by dimensionless 'dbpole'
  !=======================================================================
  subroutine make_dipoleboy(ixI^L,ixO^L,x,wB0)

    ! Subroutine arguments
    integer, intent(in)    :: ixI^L, ixO^L
    real(8), intent(in)    :: x(ixI^S,1:ndim)
    real(8), intent(inout) :: wB0(ixI^S,1:ndir)

    wB0(ixI^S,1) = dbpole * (drstar/x(ixI^S,1))**3.0d0 * cos(x(ixI^S,2))
    wB0(ixI^S,2) = 0.5d0*dbpole * (drstar/x(ixI^S,1))**3.0d0 * sin(x(ixI^S,2))
    wB0(ixI^S,3) = 0.0d0

  end subroutine make_dipoleboy

  !=========================================================================
  ! Normalise relevant quantities to be used in the code + make log overview
  !=========================================================================
  subroutine make_dimless_and_log_vars()

    ! Local variables
    character(len=8)  :: todayis
    character(len=99) :: inputlog

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
    dvinf     = dvesc * sqrt(alpha/(1.0d0 - alpha))
    dkappae   = kappae * unit_density * unit_length
    dGgrav    = Ggrav * my_unit_ggrav
    dgammae   = dkappae * dlstar/(4.d0*dpi * dGgrav * dmstar * dclight)
    detaconf  = (dbpole/2.0d0)**2.0d0 * drstar**2.0d0/(dmdot * dvinf)
    dtstat    = tstat/unit_time
    dvrot     = vrot/unit_velocity

    if (mype == 0) then
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
      write(94,*) 'unit magneticfield = ', unit_magneticfield
      write(94,*) 'unit time          = ', unit_time
      write(94,*)
      write(94,*) '==============================================='
      write(94,*) '   Stellar and wind parameters in CGS units    '
      write(94,*) '==============================================='
      write(94,*) 'Lstar/Lsun             = ', lstar/lsun
      write(94,*) 'Mstar/Msun             = ', mstar/msun
      write(94,*) 'Rstar/Rsun             = ', rstar/rsun
      write(94,*) 'Twind                  = ', twind
      write(94,*) 'Polar magnetic field   = ', bpole
      write(94,*) 'Wind confinement eta   = ', etastar
      write(94,*) 'Ralf/Rstar             = ', ralf
      write(94,*) 'Rkep/Rstar             = ', rkep
      write(94,*) 'Resc/Rstar             = ', resc
      write(94,*) 'W (vrot/vrotc)         = ', Wrot
      write(94,*) 'critical vrot          = ', vrotc
      write(94,*) 'vrot                   = ', vrot
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
      write(94,*) 'asound          = ', asound
      write(94,*) 'eff. vesc       = ', vesc
      write(94,*) 'vinf            = ', vinf
      write(94,*)
      write(94,*) 'surface density        = ', rhobound
      write(94,*) 'analytic Mdot CAK      = ', mdot * (const_years/msun)
      write(94,*) '... with FD correction = ', mdotfd * (const_years/msun)
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
      write(94,*) 'Bpole        = ', dbpole
      write(94,*) 'Eta conf.    = ', detaconf
      write(94,*) 'Edd. gamma   = ', dgammae
      write(94,*) 'rhobound     = ', drhobound
      write(94,*) 'Mdot         = ', dmdot
      write(94,*) 'alpha        = ', alpha
      write(94,*) 'Qbar         = ', Qbar
      write(94,*) 'kappae       = ', dkappae
      write(94,*) 'asound       = ', dasound
      write(94,*) 'eff. vesc    = ', dvesc
      write(94,*) 'vinf         = ', dvinf
      write(94,*) 'vrot         = ', dvrot
      write(94,*) 'clight       = ', dclight
      write(94,*) 'Ggrav        = ', dGgrav
      write(94,*) 'Tstat        = ', dtstat
      write(94,*)
      close(94)
    endif

  end subroutine make_dimless_and_log_vars

  !=========================================================================
  ! Read in a relaxed 1D CAK profile stored in a .blk file that is produced
  ! with the CAKwind_1d code. Modified version of AMRVAC mod_oneblock module
  !=========================================================================
  subroutine read_initial_oned_cak(filename)

    ! Subroutine argument
    character(len=*), intent(in) :: filename

    ! Local variables
    integer :: i, ncells, unit=69
    real(8) :: tmp, tmp1
    logical :: alive

    ! Root does the reading
    if (mype == 0) then
      inquire(file=filename,exist=alive)

      if (alive) then
        open(unit,file=filename,status='unknown')
      else
        call mpistop('Input file you want to use cannot be found!')
      endif

      ! The header information:
      read(unit,*) ! skip
      read(unit,*) ncells ! no need 2nd entry, #cells ndim=1 == #cells in ndir=1
      read(unit,*) ! skip

      ! Sanity check (not so robust)
      if (ncells /= domain_nx1) then
        close(unit)
        call mpistop('Inputfile ncells /= domain_nx1. No interpolation yet.')
      endif

      ! Allocate and read the grid and variables
      allocate(xoneblock(ncells,1))
      allocate(woneblock(ncells,1:2))

      do i = 1,ncells
        read(unit,*) xoneblock(i,1),woneblock(i,:),tmp,tmp1
      enddo

      close(unit)
    endif

    call MPI_BARRIER(icomm,ierrmpi)

    ! Broadcast what mype=0 read
    if (npe > 1) then
      call MPI_BCAST(ncells,1,MPI_INTEGER,0,icomm,ierrmpi)

      if (mype /= 0) then
        allocate(xoneblock(ncells,1))
        allocate(woneblock(ncells,1:2))
      endif

      call MPI_BCAST(xoneblock,ncells*1,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
      call MPI_BCAST(woneblock,ncells*2,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
    endif

  end subroutine read_initial_oned_cak

end module mod_usr
