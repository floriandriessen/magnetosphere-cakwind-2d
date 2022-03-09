!===============================================================================
! 2D module for launching a CAK line-driven wind from the stellar surface into
! the magnetosphere (dynamical or centrifugal) of a magnetic massive star
! Setup inspired by the original work of ud-Doula & Owocki (2002)
!
! Radiation force can be pure 1-D radial (iprob=0) or full 3-D (iprob=1)
!
! Coded up by Flo for his KU Leuven PhD thesis 2018-2022
!-------------------------------------------------------------------------------
! (October 2019) -- Flo
!   > setup of problem with 1D CAK wind starting from beta law in 2D geometry
!
! (May 2020) -- Flo
!   > implemented routine for computing time-averaged statistical variables
!
! (June 2020) -- Flo
!   > implementation to do Bfield splitting via Tanaka method by calling the
!     usr_set_B0 routine containing a stellar dipole
!     NOTICE: in AMRVAC there is build-in dipole that can be activated with
!             'Bdip' in .par file, however, this formula leads to an eta_star
!             4x bigger (AMRVAC has different way of expressing dipole)
!   > change to better velocity gradient stencil to deal with non-uniform grids
!
! (October 2020) -- Flo
!   > added rigid body rotation to study centrifugal magnetospheres, set in
!     .par file with dimensionless 'Wrot' parameter (ud-Doula et al. 2006)
!   > 'imag' to select field strength (imag>0) or wind confinement (imag<0)
!
! (February 2021) -- Flo
!   > put effective gravity force in usr_gravity and removed from usr_source
!   > compute radiation timestep in special_dt by using gcak slot in w-array
!
! (April 2021) -- Flo
!   > prohibit pure adiabatic cooling as it crashes, but is also unrealistic
!   > inclusion of rotating frame with fictitious forces + update special_dt
!
! (May 2021) -- Flo
!   > modified vtheta boundary to avoid funny business at poles for etastar>100
!
! (March 2022) -- Flo
!   > inclusion of 3-D line force; now selection between forces using "iprob"
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
  real(8)           :: Qmax, tstat, Wrot
  integer           :: nthetap, nphip
  logical           :: rotframe
  character(len=99) :: cakfile

  ! Additionally useful stellar and wind parameters:
  real(8) :: gammae, vesc, mdot, mdotfd, vinf, asound, logg, logge, heff
  real(8) :: mumol, bpole, etastar, ralf, rkep, resc, vrot, vrotc

  ! Dimensionless variables of relevant variables above
  real(8) :: dlstar, dmstar, drstar, dbpole, drhobound, dtwind, dkappae
  real(8) :: dvesc, dvinf, dmdot, dasound, dclight, dGgrav, dgammae, detaconf
  real(8) :: dtstat, dvrot, dOmegarot

  ! Additional names for wind and statistical variables
  integer :: my_gcak, my_rhoav, my_rho2av, my_vrav, my_vr2av, my_rhovrav
  integer :: my_vpolav, my_vpol2av, my_tav, my_gcakr, my_gcakt, my_gcakp

  ! Arrays required to read and store 1D profile from file
  real(8), allocatable :: woneblock(:,:), xoneblock(:,:)

  ! Ray parameters
  real(8), allocatable :: ay(:), wy(:), aphi(:), wphi(:)
  real(8) :: dphi

contains

  !======================================================================
  ! This routine should set user methods, and activate the physics module
  !======================================================================
  subroutine usr_init()

    call set_coordinate_system("spherical_2.5D")
    call usr_params_read(par_files)

    ! Choose normalisation units: (length,temp,ndens) or (length,vel,ndens)
    ! numberdensity chosen such that unit density becomes boundary density
    unit_length        = rstar                                      ! cm
    unit_temperature   = twind                                      ! K
    unit_numberdensity = rhobound/((1.d0+4.d0*He_abundance)*mp_cgs) ! g cm^-3

    call MHD_activate()

    usr_set_parameters   => initglobaldata_usr
    usr_init_one_grid    => initial_conditions
    usr_special_bc       => special_bound
    usr_gravity          => effective_gravity
    usr_source           => special_source
    usr_get_dt           => special_dt
    usr_process_adv_grid => compute_stats
    usr_aux_output       => set_extravar_output
    usr_add_aux_names    => set_extravarnames_output
    usr_set_B0           => make_dipoleboy

    my_gcak    = var_set_extravar("gcak", "gcak")
    my_gcakr   = var_set_extravar("gcak_r", "gcak_r")
    my_gcakt   = var_set_extravar("gcak_theta", "gcak_theta")
    my_gcakp   = var_set_extravar("gcak_phi", "gcak_phi")
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
                          Qbar, tstat, Wrot, cakfile, Qmax, rotframe, nphip, nthetap

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
    Qmax   = Qmax * Qbar
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
    if (iprob == 1) call ray_init(nthetap,nphip)

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

    ! No polar velocity
    w(ixO^S,mom(2)) = 0.0d0

    ! In azimuth no velocity (rotating frame) or rigid rotation for full grid
    if (rotframe) then
      w(ixO^S,mom(3)) = 0.0d0
    else
      w(ixO^S,mom(3)) = dOmegarot * x(ixO^S,1) * sin(x(ixO^S,2))
    endif

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

  !=================================================================
  ! Special user boundary conditions at inner radial boundary:
  !   vr, Btheta, Bphi (extrapolated); rho, vtheta, vphi, Br (fixed)
  !=================================================================
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

      if (rotframe) then
        w(ixB^S,mom(3)) = 0.0d0
      else
        w(ixB^S,mom(3)) = dvrot * sin(x(ixB^S,2))
      endif

      ! === Tanaka splitting ===
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
      ! === Regular ===
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

      ! Enforce poloidal flow along magnetic field; outside magnetosphere from
      ! induction equation, inside magnetosphere put at zero
      if (etastar > 1.0d0) then
        if (B0field) then
          where ( abs(0.5d0*dpi - x(ixB^S,2)) >= asin(sqrt(drstar/ralf)) )
            w(ixB^S,mom(2)) = w(ixB^S,mom(1)) &
                              * (block%B0(ixB^S,2,0) + w(ixB^S,mag(2))) &
                              / (block%B0(ixB^S,1,0) + w(ixB^S,mag(1)))
          elsewhere
            w(ixB^S,mom(2)) = 0.0d0
          endwhere
        else
          where ( abs(0.5d0*dpi - x(ixB^S,2)) >= asin(sqrt(drstar/ralf)) )
            w(ixB^S,mom(2)) = w(ixB^S,mom(1)) * w(ixB^S,mag(2))/w(ixB^S,mag(1))
          elsewhere
            w(ixB^S,mom(2)) = 0.0d0
          endwhere
        endif
      else
        w(ixB^S,mom(2)) = 0.0d0
      endif

      ! Density is fixed, so per ideal gas law also the pressure
      if (mhd_energy) w(ixB^S,p_) = w(ixB^S,rho_)

      if (mhd_glm) w(ixB^S,psi_) = 0.0d0

      call mhd_to_conserved(ixI^L,ixI^L,w,x)

    case default
      call mpistop("BC not specified")
    end select

  end subroutine special_bound

  !==========================================================================
  ! Special user source term to model radiation line force and possibly added
  ! fictitious forces when rotating frame is activated
  !==========================================================================
  subroutine special_source(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)

    ! Subroutine arguments
    integer, intent(in)    :: ixI^L, ixO^L, iw^LIM
    real(8), intent(in)    :: qdt, qtC, qt
    real(8), intent(in)    :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    real(8), intent(inout) :: w(ixI^S,1:nw)

    ! Local variables
    real(8) :: gcak(ixO^S,1:3), fcent(ixO^S,1:ndir), fcor(ixO^S,1:ndir)
    real(8) :: fac, fac1, fac2

    ! Initialize forces (radial=1, polar=2, azimuthal=3)
    gcak(ixO^S,1) = 0.0d0
    gcak(ixO^S,2) = 0.0d0
    gcak(ixO^S,3) = 0.0d0

    ! Get force in each point
    select case (iprob)
      case(0)
        call line_force_radial(ixI^L,ixO^L,wCT,x,gcak)
      case(1)
        call line_force_3d(ixI^L,ixO^L,wCT,x,gcak)
      case default
        call mpistop("Select a valid iprob value!")
    end select

    ! Normalization force
    fac1 = 1.0d0/(1.0d0 - alpha) * dkappae * dlstar*Qbar/(4.0d0*dpi * dclight)
    fac2 = 1.0d0/(dclight * Qbar * dkappae)**alpha
    fac  = fac1 * fac2

    gcak(ixO^S,1:3) = gcak(ixO^S,1:3) * fac

    ! Fill the nwextra slots for output
    w(ixO^S,my_gcakr) = gcak(ixO^S,1)
    w(ixO^S,my_gcakt) = gcak(ixO^S,2)
    w(ixO^S,my_gcakp) = gcak(ixO^S,3)
    w(ixO^S,my_gcak)  = gcak(ixO^S,1) + gcak(ixO^S,2) + gcak(ixO^S,3)

    ! Update conservative vars: w = w + qdt*gsource
    w(ixO^S,mom(1)) = w(ixO^S,mom(1)) + qdt * gcak(ixO^S,1) * wCT(ixO^S,rho_)
    w(ixO^S,mom(2)) = w(ixO^S,mom(2)) + qdt * gcak(ixO^S,2) * wCT(ixO^S,rho_)
    w(ixO^S,mom(3)) = w(ixO^S,mom(3)) + qdt * gcak(ixO^S,3) * wCT(ixO^S,rho_)

    ! Update total energy e = e + vCT*rhoCT*qdt*gsource
    if (mhd_energy) then
      w(ixO^S,e_) = w(ixO^S,e_) + qdt * ( gcak(ixO^S,1) * wCT(ixO^S,mom(1)) &
                                        + gcak(ixO^S,2) * wCT(ixO^S,mom(2)) &
                                        + gcak(ixO^S,3) * wCT(ixO^S,mom(3)) )
    endif

    ! Update conservative vars if in rotating frame
    if (rotframe) then
      call get_fict_forces(ixI^L,ixO^L,wCT,x,fcent,fcor)

      w(ixO^S,mom(1)) = w(ixO^S,mom(1)) - qdt &
                         * (fcent(ixO^S,1) + fcor(ixO^S,1)) * wCT(ixO^S,rho_)

      w(ixO^S,mom(2)) = w(ixO^S,mom(2)) - qdt &
                         * (fcent(ixO^S,2) + fcor(ixO^S,2)) * wCT(ixO^S,rho_)

      w(ixO^S,mom(3)) = w(ixO^S,mom(3)) - qdt * fcor(ixO^S,3) * wCT(ixO^S,rho_)

      ! Only centrifugal force performs fictitious work
      if (mhd_energy) then
        w(ixO^S,e_) = w(ixO^S,e_) - qdt * ( fcent(ixO^S,1) * wCT(ixO^S,mom(1)) &
                                          + fcent(ixO^S,2) * wCT(ixO^S,mom(2)) )
      endif
    endif

  end subroutine special_source

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
    real(8) :: tdumr(ixO^S), tdumt(ixO^S), force(ixO^S), dt_r, dt_t
    real(8) :: fcent(ixO^S,1:ndir), fcor(ixO^S,1:ndir)

    dt_r = bigdouble
    dt_t = bigdouble

    if (rotframe) then
      call get_fict_forces(ixI^L,ixO^L,w,x,fcent,fcor)
    else
      fcent(ixO^S,:) = 0.0d0
      fcor(ixO^S,:)  = 0.0d0
    endif

    ! Get dt from force in each direction -- no phi-waves in 2.5D
    ! === Radial ===
    force(ixO^S) = w(ixO^S,my_gcakr) + (fcent(ixO^S,1) + fcor(ixO^S,1))
    tdumr(ixO^S) = sqrt( block%dx(ixO^S,1) / abs(force(ixO^S)) )
    dt_r         = courantpar * minval(tdumr(ixO^S))

    ! === Polar ===
    force(ixO^S) = w(ixO^S,my_gcakt) + (fcent(ixO^S,2) + fcor(ixO^S,2))
    tdumt(ixO^S) = sqrt( block%dx(ixO^S,1) * block%dx(ixO^S,2) / abs(force(ixO^S)) )
    dt_t         = courantpar * minval(tdumt(ixO^S))

    if (it >= 1) then
      dtnew = min(dtnew,dt_r,dt_t)
    endif

  end subroutine special_dt

  !============================================================================
  ! Combine stellar gravity and continuum electron scattering into an effective
  ! gravity using Eddington's gamma
  !============================================================================
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

  !=========================================================================
  ! Routine computes the time-averaged statistical quantity <X> via:
  !   <X>_i = <X>_i-1 + dt * X_i
  ! where <X> is the average of variable X and i','i-1' are the current and
  ! previous timestep respectively
  !=========================================================================
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

    ! Output divB, for Tanaka splitting this will be divB1
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

    ! Local variable
    character(len=8)  :: todayis

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
    dOmegarot = dvrot/drstar

    if (mype == 0) then

      call date_and_time(todayis)
      write(*,*) 'MPI-AMRVAC simulation ran on ', &
                    todayis(7:8), '/', todayis(5:6), '/', todayis(1:4)
      write(*,*)
      write(*,*) '======================'
      write(*,*) '   Unity quantities   '
      write(*,*) '======================'
      write(*,*) 'unit length        = ', unit_length
      write(*,*) 'unit density       = ', unit_density
      write(*,*) 'unit velocity      = ', unit_velocity
      write(*,*) 'unit numberdensity = ', unit_numberdensity
      write(*,*) 'unit pressure      = ', unit_pressure
      write(*,*) 'unit temperature   = ', unit_temperature
      write(*,*) 'unit magneticfield = ', unit_magneticfield
      write(*,*) 'unit time          = ', unit_time
      write(*,*)
      write(*,*) '==============================================='
      write(*,*) '   Stellar and wind parameters in CGS units    '
      write(*,*) '==============================================='
      write(*,*) 'Lstar/Lsun             = ', lstar/lsun
      write(*,*) 'Mstar/Msun             = ', mstar/msun
      write(*,*) 'Rstar/Rsun             = ', rstar/rsun
      write(*,*) 'Twind                  = ', twind
      write(*,*) 'Polar magnetic field   = ', bpole
      write(*,*) 'Wind confinement eta   = ', etastar
      write(*,*) 'Ralf/Rstar             = ', ralf
      write(*,*) 'Rkep/Rstar             = ', rkep
      write(*,*) 'Resc/Rstar             = ', resc
      write(*,*) 'Mean molecular weight  = ', mumol
      write(*,*) 'log(g)                 = ', logg
      write(*,*) 'eff. log(g)            = ', logge
      write(*,*) 'eff. scale height heff = ', heff
      write(*,*) 'heff/Rstar             = ', heff/rstar
      write(*,*) 'W (vrot/vrotc)         = ', Wrot
      write(*,*) 'critical vrot          = ', vrotc
      write(*,*) 'vrot                   = ', vrot
      write(*,*)
      write(*,*) 'adiabatic gamma = ', mhd_gamma
      write(*,*) 'Eddington gamma = ', gammae
      write(*,*) 'alpha           = ', alpha
      write(*,*) 'Qbar            = ', Qbar
      write(*,*) 'Qmax/Qbar       = ', Qmax/Qbar
      write(*,*) 'asound          = ', asound
      write(*,*) 'eff. vesc       = ', vesc
      write(*,*) 'CAK vinf        = ', vinf
      write(*,*) 'FD vinf         = ', 3.0d0 * vesc
      write(*,*)
      write(*,*) 'surface density        = ', rhobound
      write(*,*) 'analytic Mdot CAK      = ', mdot * (const_years/msun)
      write(*,*) '... with FD correction = ', mdot/(1.0d0 + alpha)**(1.0d0/alpha) * (const_years/msun)
      write(*,*)
      write(*,*) '========================================'
      write(*,*) '    Dimensionless AMRVAC quantities     '
      write(*,*) '========================================'
      write(*,*) 'Extra computed unit quantities:'
      write(*,*) '   unit Lum  = ', my_unit_lum
      write(*,*) '   unit Mass = ', my_unit_mass
      write(*,*) '   unit Grav = ', my_unit_ggrav
      write(*,*) 'Lstar        = ', dlstar
      write(*,*) 'Mstar        = ', dmstar
      write(*,*) 'Rstar        = ', drstar
      write(*,*) 'Twind        = ', dtwind
      write(*,*) 'Bpole        = ', dbpole
      write(*,*) 'Eta conf.    = ', detaconf
      write(*,*) 'rhobound     = ', drhobound
      write(*,*) 'Mdot         = ', dmdot
      write(*,*) 'kappae       = ', dkappae
      write(*,*) 'asound       = ', dasound
      write(*,*) 'eff. vesc    = ', dvesc
      write(*,*) 'vinf (CAK)   = ', dvinf
      write(*,*) 'vrot         = ', dvrot
      write(*,*) 'clight       = ', dclight
      write(*,*) 'Ggrav        = ', dGgrav
      write(*,*) 'Tstat        = ', dtstat
      write(*,*)
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

  !====================================================
  ! Analytical 3-D CAK line force in Gayley's formalism
  !====================================================
  subroutine line_force_3d(ixI^L,ixO^L,wCT,x,gcak)

    ! Subroutine arguments
    integer, intent(in)    :: ixI^L, ixO^L
    real(8), intent(in)    :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    real(8), intent(inout) :: gcak(ixO^S,1:3)

    ! Local variables
    real(8) :: vr(ixI^S), vt(ixI^S), vp(ixI^S), rho(ixI^S)
    real(8) :: vrr(ixI^S), vtr(ixI^S), vpr(ixI^S)
    real(8) :: dvrdr_up(ixO^S), dvrdr_down(ixO^S), dvrdr_cent(ixO^S), dvrdr(ixO^S)
    real(8) :: dvtdr_up(ixO^S), dvtdr_down(ixO^S), dvtdr_cent(ixO^S), dvtdr(ixO^S)
    real(8) :: dvpdr_up(ixO^S), dvpdr_down(ixO^S), dvpdr_cent(ixO^S), dvpdr(ixO^S)
    real(8) :: dvrdt_up(ixO^S), dvrdt_down(ixO^S), dvrdt_cent(ixO^S), dvrdt(ixO^S)
    real(8) :: dvtdt_up(ixO^S), dvtdt_down(ixO^S), dvtdt_cent(ixO^S), dvtdt(ixO^S)
    real(8) :: dvpdt_up(ixO^S), dvpdt_down(ixO^S), dvpdt_cent(ixO^S), dvpdt(ixO^S)
    real(8) :: dvrdp_up(ixO^S), dvrdp_down(ixO^S), dvrdp_cent(ixO^S), dvrdp(ixO^S)
    real(8) :: dvtdp_up(ixO^S), dvtdp_down(ixO^S), dvtdp_cent(ixO^S), dvtdp(ixO^S)
    real(8) :: dvpdp_up(ixO^S), dvpdp_down(ixO^S), dvpdp_cent(ixO^S), dvpdp(ixO^S)
    real(8) :: a1, a2, a3, dvdn, wyray, y, wpray, phip, wtot, mustar
    real(8) :: costp, costp2, sintp, cospp, sinpp, cott0
    integer :: ix^D, jrx^L, hrx^L, jtx^L, htx^L, itp, ipp

    ! Define time-centred velocity from the momentum and density
    vr(ixI^S)  = wCT(ixI^S,mom(1)) / wCT(ixI^S,rho_)
    vt(ixI^S)  = wCT(ixI^S,mom(2)) / wCT(ixI^S,rho_)
    vp(ixI^S)  = wCT(ixI^S,mom(3)) / wCT(ixI^S,rho_)
    rho(ixI^S) = wCT(ixI^S,rho_)
    vrr(ixI^S) = vr(ixI^S) / x(ixI^S,1)
    vtr(ixI^S) = vt(ixI^S) / x(ixI^S,1)
    vpr(ixI^S) = vp(ixI^S) / x(ixI^S,1)

    ! Index +1 (j) and index -1 (h) in radial direction; kr(dir,dim)=1, dir=dim
    jrx^L=ixO^L+kr(1,^D);
    hrx^L=ixO^L-kr(1,^D);

    ! Index +1 (j) and index -1 (h) in theta direction
    jtx^L=ixO^L+kr(2,^D);
    htx^L=ixO^L-kr(2,^D);

    ! Get dvn/dn on non-uniform grid according to Sundqvist & Veronis (1970)

    ! === RADIAL DERIVATIVES ===
    ! Forward, backward, central difference
    dvrdr_up(ixO^S) = (x(ixO^S,1) - x(hrx^S,1)) * vr(jrx^S) &
                     / ((x(jrx^S,1) - x(ixO^S,1)) * (x(jrx^S,1) - x(hrx^S,1)))

    dvrdr_down(ixO^S) = -(x(jrx^S,1) - x(ixO^S,1)) * vr(hrx^S) &
                       / ((x(ixO^S,1) - x(hrx^S,1)) * (x(jrx^S,1) - x(hrx^S,1)))

    dvrdr_cent(ixO^S) = (x(jrx^S,1) + x(hrx^S,1) - 2.0d0*x(ixO^S,1)) * vr(ixO^S) &
                      / ((x(ixO^S,1) - x(hrx^S,1)) * (x(jrx^S,1) - x(ixO^S,1)))

    dvtdr_up(ixO^S) = (x(ixO^S,1) - x(hrx^S,1)) * vt(jrx^S) &
                     / ((x(jrx^S,1) - x(ixO^S,1)) * (x(jrx^S,1) - x(hrx^S,1)))

    dvtdr_down(ixO^S) = -(x(jrx^S,1) - x(ixO^S,1)) * vt(hrx^S) &
                       / ((x(ixO^S,1) - x(hrx^S,1)) * (x(jrx^S,1) - x(hrx^S,1)))

    dvtdr_cent(ixO^S) = (x(jrx^S,1) + x(hrx^S,1) - 2.0d0*x(ixO^S,1)) * vt(ixO^S) &
                      / ((x(ixO^S,1) - x(hrx^S,1)) * (x(jrx^S,1) - x(ixO^S,1)))

    dvpdr_up(ixO^S) = (x(ixO^S,1) - x(hrx^S,1)) * vp(jrx^S) &
                    / ((x(jrx^S,1) - x(ixO^S,1)) * (x(jrx^S,1) - x(hrx^S,1)))

    dvpdr_down(ixO^S) = -(x(jrx^S,1) - x(ixO^S,1)) * vp(hrx^S) &
                        / ((x(ixO^S,1) - x(hrx^S,1)) * (x(jrx^S,1) - x(hrx^S,1)))

    dvpdr_cent(ixO^S) = (x(jrx^S,1) + x(hrx^S,1) - 2.0d0*x(ixO^S,1)) * vp(ixO^S) &
                       / ((x(ixO^S,1) - x(hrx^S,1)) * (x(jrx^S,1) - x(ixO^S,1)))

    ! === POLAR DERIVATIVES ===
    ! Forward, backward, central difference
    dvrdt_up(ixO^S) = (x(ixO^S,2) - x(htx^S,2)) * vr(jtx^S) &
                     / (x(ixO^S,1) * (x(jtx^S,2) - x(ixO^S,2)) * (x(jtx^S,2) - x(htx^S,2)))

    dvrdt_down(ixO^S) = -(x(jtx^S,2) - x(ixO^S,2)) * vr(htx^S) &
                       / ( x(ixO^S,1) * (x(ixO^S,2) - x(htx^S,2)) * (x(jtx^S,2) - x(htx^S,2)))

    dvrdt_cent(ixO^S) = (x(jtx^S,2) + x(htx^S,2) - 2.0d0*x(ixO^S,2)) * vr(ixO^S) &
                     / ( x(ixO^S,1) * (x(ixO^S,2) - x(htx^S,2)) * (x(jtx^S,2) - x(ixO^S,2)))

    dvtdt_up(ixO^S) = (x(ixO^S,2) - x(htx^S,2)) * vt(jtx^S) &
                     / (x(ixO^S,1) * (x(jtx^S,2) - x(ixO^S,2)) * (x(jtx^S,2) - x(htx^S,2)))

    dvtdt_down(ixO^S) = -(x(jtx^S,2) - x(ixO^S,2)) * vt(htx^S) &
                      / ( x(ixO^S,1) * (x(ixO^S,2) - x(htx^S,2)) * (x(jtx^S,2) - x(htx^S,2)))

    dvtdt_cent(ixO^S) = (x(jtx^S,2) + x(htx^S,2) - 2.0d0*x(ixO^S,2)) * vt(ixO^S) &
                      / ( x(ixO^S,1) * (x(ixO^S,2) - x(htx^S,2)) * (x(jtx^S,2) - x(ixO^S,2)))

    dvpdt_up(ixO^S) = (x(ixO^S,2) - x(htx^S,2)) * vp(jtx^S) &
                     / ((x(jtx^S,2) - x(ixO^S,2)) * (x(jtx^S,2) - x(htx^S,2)))

    dvpdt_down(ixO^S) = -(x(jtx^S,2) - x(ixO^S,2)) * vp(htx^S) &
                       / ((x(ixO^S,2) - x(htx^S,2)) * (x(jtx^S,2) - x(htx^S,2)))

    dvpdt_cent(ixO^S) = (x(jtx^S,2) + x(htx^S,2) - 2.0d0*x(ixO^S,2)) * vp(ixO^S) &
                     / ((x(ixO^S,2) - x(htx^S,2)) * (x(jtx^S,2) - x(ixO^S,2)))

    ! === AZIMUTHAL DERIVATIVE ===
    ! Forward, backward, central difference -- all zero in 2.5D
    dvrdp_up(ixO^S)   = 0.0d0
    dvrdp_down(ixO^S) = 0.0d0
    dvrdp_cent(ixO^S) = 0.0d0
    dvtdp_up(ixO^S)   = 0.0d0
    dvtdp_down(ixO^S) = 0.0d0
    dvtdp_cent(ixO^S) = 0.0d0
    dvpdp_up(ixO^S)   = 0.0d0
    dvpdp_down(ixO^S) = 0.0d0
    dvpdp_cent(ixO^S) = 0.0d0

    ! Total gradients in radial, polar, azimuthal direction
    dvrdr(ixO^S) = dvrdr_down(ixO^S) + dvrdr_cent(ixO^S) + dvrdr_up(ixO^S)
    dvtdr(ixO^S) = dvtdr_down(ixO^S) + dvtdr_cent(ixO^S) + dvtdr_up(ixO^S)
    dvpdr(ixO^S) = dvpdr_down(ixO^S) + dvpdr_cent(ixO^S) + dvpdr_up(ixO^S)

    dvrdt(ixO^S) = dvrdt_down(ixO^S) + dvrdt_cent(ixO^S) + dvrdt_up(ixO^S)
    dvtdt(ixO^S) = dvtdt_down(ixO^S) + dvtdt_cent(ixO^S) + dvtdt_up(ixO^S)
    dvpdt(ixO^S) = dvpdt_down(ixO^S) + dvpdt_cent(ixO^S) + dvpdt_up(ixO^S)

    dvrdp(ixO^S) = dvrdp_down(ixO^S) + dvrdp_cent(ixO^S) + dvrdp_up(ixO^S)
    dvtdp(ixO^S) = dvtdp_down(ixO^S) + dvtdp_cent(ixO^S) + dvtdp_up(ixO^S)
    dvpdp(ixO^S) = dvpdp_down(ixO^S) + dvpdp_cent(ixO^S) + dvpdp_up(ixO^S)

    ! Get total acceleration from all rays at a certain grid point
    {do ix^DB=ixOmin^DB,ixOmax^DB\}
      ! Loop over the rays; first theta then phi radiation angle
      ! Get weights from current ray and their position
      do itp = 1,nthetap
        wyray  = wy(itp)
        y      = ay(itp)

        do ipp = 1,nphip
          wpray = wphi(ipp)
          phip  = aphi(ipp)

          ! Redistribute the phi rays by a small offset
          if (mod(itp,3) == 1) then
            phip = phip + dphi/3.0d0
          elseif (mod(itp,3) == 2) then
            phip = phip - dphi/3.0d0
          endif

          ! === Geometrical factors ===
          ! Make y quadrature linear in mu, not mu**2; better for gtheta,gphi
          ! y -> mu quadrature is preserved; y=0 <=> mu=1; y=1 <=> mu=mustar
          mustar = sqrt(max(1.0d0 - (drstar/x(ix^D,1))**2.0d0, 0.0d0))
          costp  = 1.0d0 - y*(1.0d0 - mustar)
          costp2 = costp*costp
          sintp  = sqrt(max(1.0d0 - costp2, 0.0d0))
          sinpp  = sin(phip)
          cospp  = cos(phip)
          cott0  = cos(x(ix^D,2))/sin(x(ix^D,2))

          ! More weight close to star, less farther away
          wtot  = wyray * wpray * (1.0d0 - mustar)

          ! Convenients a la Cranmer & Owocki (1995)
          a1 = costp
          a2 = sintp * cospp
          a3 = sintp * sinpp

          ! Get total velocity gradient for one ray with given (theta', phi')
          dvdn = a1*a1 * dvrdr(ix^D) + a2*a2 * (dvtdt(ix^D) + vrr(ix^D))  &
                + a3*a3 * (dvpdp(ix^D) + cott0 * vtr(ix^D) + vrr(ix^D))   &
                + a1*a2 * (dvtdr(ix^D) + dvrdt(ix^D) - vtr(ix^D))         &
                + a1*a3 * (dvpdr(ix^D) + dvrdp(ix^D) - vpr(ix^D))         &
                + a2*a3 * (dvpdt(ix^D) + dvtdp(ix^D) - cott0 * vpr(ix^D))

          ! Fallback requirement check
          dvdn = max(dvdn, smalldouble)

          ! Convert gradient back from wind coordinates (r',theta',phi') to
          ! stellar coordinates (r,theta,phi)
          gcak(ix^D,1) = gcak(ix^D,1) + (dvdn/rho(ix^D))**alpha * a1 * wtot
          gcak(ix^D,2) = gcak(ix^D,2) + (dvdn/rho(ix^D))**alpha * a2 * wtot
          gcak(ix^D,3) = gcak(ix^D,3) + (dvdn/rho(ix^D))**alpha * a3 * wtot
        enddo
      enddo
    {enddo\}

    gcak(ixO^S,1:3) = gcak(ixO^S,1:3)/drstar**2.0d0

  end subroutine line_force_3d

  !==============================================================
  ! Analytical purely radial CAK line force in Gayley's formalism
  !==============================================================
  subroutine line_force_radial(ixI^L,ixO^L,wCT,x,gcak)

    ! Subroutine arguments
    integer, intent(in)    :: ixI^L, ixO^L
    real(8), intent(in)    :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    real(8), intent(inout) :: gcak(ixO^S,1:3)

    ! Local variables
    real(8) :: vr(ixI^S), rho(ixI^S)
    real(8) :: dvdr_up(ixO^S), dvdr_down(ixO^S), dvdr_cent(ixO^S), dvdr(ixO^S)
    real(8) :: beta_fd(ixO^S), fdfac(ixO^S)
    real(8) :: taum(ixO^S), taumfac(ixO^S)
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
    dvdr(ixO^S) = max(dvdr(ixO^S), smalldouble)

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

    ! Only radial force in this case
    gcak(ixO^S,1) = fdfac/x(ixO^S,1)**2.0d0 * (dvdr(ixO^S)/rho(ixO^S))**alpha

    ! Correction for opacity cutoff
    ! taum(ixO^S) = dkappae * dclight * Qmax * rho(ixO^S)/dvdr(ixO^S)
    !
    ! taumfac(ixO^S) = ((1.0d0 + taum(ixO^S))**(1.0d0 - alpha) - 1.0d0) &
    !                     / taum(ixO^S)**(1.0d0 - alpha)
    !
    ! gcak(ixO^S) = gcak(ixO^S) * taumfac(ixO^S)

  end subroutine line_force_radial

  !=============================================================================
  ! Initialize on each block (theta',phi') radiation angles coming from the star
  !=============================================================================
  subroutine ray_init(ntheta_point,nphi_point)

    ! Subroutine arguments
    integer, intent(in) :: ntheta_point, nphi_point

    ! Local variables
    real(8) :: ymin, ymax, phipmin, phipmax, dphi, adum
    integer :: ii

    ! Minimum and maximum range of theta and phi rays
    ! NOTE: theta points are cast into y-space
    ymin    = 0.0d0
    ymax    = 1.0d0
    phipmin = 0.0d0
    phipmax = 2.0d0*dpi

    ! Allocation for arrays on each block
    allocate(ay(ntheta_point))
    allocate(wy(ntheta_point))
    allocate(aphi(nphi_point))
    allocate(wphi(nphi_point))

    ! theta ray positions and weights: Gauss-Legendre
    call gauss_legendre_quadrature(ymin,ymax,ntheta_point,ay,wy)

    ! theta rays and weights: uniform
    ! dth = 1.0d0 / nthetap
    ! adum = ymin + 0.5d0*dth
    ! do ip = 1,nthetap
    !   ay(ip) = adum
    !   wy(ip) = 1.0d0/nthetap
    !   adum = adum + dth
    !   !print*,'phipoints'
    !   !print*,ip,aphi(ip),wphi(ip),dphi
    ! enddo

    ! phi ray position and weights: uniform
    dphi = (phipmax - phipmin) / nphi_point
    adum = phipmin + 0.5d0*dphi
    do ii = 1,nphi_point
      aphi(ii) = adum
      wphi(ii) = 1.0d0/nphi_point
      adum     = adum + dphi
    enddo

    if (mype == 0) then
      write(*,*) '==========================='
      write(*,*) '    Radiation ray setup    '
      write(*,*) '==========================='
      write(*,*) 'Theta ray points + weights '
      do ii = 1,ntheta_point
        write(*,*) ii,ay(ii),wy(ii)
      enddo
      write(*,*)
      write(*,*) 'Phi ray points + weights   '
      do ii = 1,nphi_point
        write(*,*) ii,aphi(ii),wphi(ii)
      enddo
      write(*,*)
    endif

  end subroutine ray_init

  !===========================================================================
  ! Given the lower and upper limits of integration xlow and xhi, and given n,
  ! this routine returns arrays x and w of length n, containing the abscissas
  ! and weights of the Gauss-Legendre N-point quadrature formula.
  ! Originated by G. Rybicki; adapted from Numerical Recipes F77, p. 145
  !===========================================================================
  subroutine gauss_legendre_quadrature(xlow,xhi,n,x,w)

    ! Subroutine arguments
    real(8), intent(in)  :: xlow, xhi
    integer, intent(in)  :: n
    real(8), intent(out) :: x(n), w(n)

    ! Local variables
    integer :: i, j, m
    real(8) :: p1, p2, p3, pp, xl, xm, z, z1
    real(8), parameter :: error=3.0d-14

    m = (n + 1)/2
    xm = 0.5d0*(xhi + xlow)
    xl = 0.5d0*(xhi - xlow)

    do i = 1,m
      z = cos( dpi * (i - 0.25d0)/(n + 0.5d0) )
      z1 = 2.0d0 * z

      do while (abs(z1 - z) > error)
        p1 = 1.0d0
        p2 = 0.0d0

        do j = 1,n
          p3 = p2
          p2 = p1
          p1 = ( (2.0d0*j - 1.0d0)*z*p2 - (j - 1.0d0)*p3 )/j
        enddo

        pp = n*(z*p1 - p2) / (z*z - 1.0d0)
        z1 = z
        z = z1 - p1/pp
      enddo

      x(i)     = xm - xl*z
      x(n+1-i) = xm + xl*z
      w(i)     = 2.0d0*xl / ((1.0d0 - z*z) * pp*pp)
      w(n+1-i) = w(i)
    enddo

  end subroutine gauss_legendre_quadrature

  !===================================================================
  ! Compute fictitious forces on spherical grid when in rotating frame
  !===================================================================
  subroutine get_fict_forces(ixI^L,ixO^L,w,x,fcent,fcor)

    ! Subroutine arguments
    integer, intent(in)  :: ixI^L, ixO^L
    real(8), intent(in)  :: x(ixI^S,1:ndim), w(ixI^S,1:nw)
    real(8), intent(out) :: fcent(ixO^S,1:ndir), fcor(ixO^S,1:ndir)

    ! Centrigufal force (radial=1, polar=2, azimuthal=3)
    fcent(ixO^S,1)   = -sin(x(ixO^S,2))**2.0d0 * x(ixO^S,1)
    fcent(ixO^S,2)   = -sin(x(ixO^S,2)) * cos(x(ixO^S,2)) * x(ixO^S,1)
    fcent(ixO^S,3)   = 0.0d0
    fcent(ixO^S,1:3) = dOmegarot**2.0d0 * fcent(ixO^S,1:3)

    ! Coriolis force
    fcor(ixO^S,1)   = -w(ixO^S,mom(3))/w(ixO^S,rho_) * sin(x(ixO^S,2))
    fcor(ixO^S,2)   = -w(ixO^S,mom(3))/w(ixO^S,rho_) * cos(x(ixO^S,2))
    fcor(ixO^S,3)   = w(ixO^S,mom(1))/w(ixO^S,rho_) * sin(x(ixO^S,2)) &
                      + w(ixO^S,mom(2))/w(ixO^S,rho_) * cos(x(ixO^S,2))
    fcor(ixO^S,1:3) = 2.0d0*dOmegarot * fcor(ixO^S,1:3)

  end subroutine get_fict_forces

end module mod_usr
