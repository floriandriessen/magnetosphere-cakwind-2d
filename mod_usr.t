!==============================================================================
! 2.5D model for simulating a CAK line-driven wind interacting with a stellar
! magnetosphere (dynamical or centrifugal) of a magnetic OB-type star.
! Setup inspired by the original work of ud-Doula & Owocki (2002), ApJ 576.
!
! This module relies on my developed mod_cak_force module in MPI-AMRVAC and a
! separate routine to read in 1-D smooth CAK wind initial conditions.
!
! Coded up by Florian Driessen: 2018-2022 (PhD), 2025
!==============================================================================
module mod_usr

  use mod_mhd
  use mod_cak_force, ONLY: set_cak_force_norm, cak_alpha, gayley_qbar
  use mod_constants, ONLY: const_c, const_G, const_LSun, const_MSun, &
       const_RSun, const_kappae, const_sigma, mp_cgs, kB_cgs, const_years

  implicit none

  ! User input parameters
  real(8) :: mstar_sol, rstar_sol, twind_cgs, rhosurf_cgs, timestat_cgs, Wrot
  real(8) :: bpole_cgs=-99.0, etastar=-99.0
  character(len=99) :: cakfile

  ! Dimensionless variables for computations
  real(8) :: mstar, rstar, bpole, rhosurf, csound, clight, vrot
  real(8) :: gmstar, timestat, ralf

  ! Additional names for extra statistical variables in output
  integer :: irhoav_, irho2av_, ivrav_, ivr2av_, irhovrav_
  integer :: ivpolav_, ivpol2av_, itav_

contains

!==============================================================================
! This routine should set user methods, and activate the physics module.
!==============================================================================
  subroutine usr_init

    call set_coordinate_system("spherical_2.5D")
    call usr_params_read(par_files)

    ! Choose normalisation units:
    unit_length      = rstar_sol * const_RSun
    unit_temperature = twind_cgs
    unit_density     = rhosurf_cgs

    usr_set_parameters   => initglobaldata_usr
    usr_init_one_grid    => initial_conditions
    usr_special_bc       => special_bound
    usr_gravity          => stellar_gravity
    usr_process_adv_grid => compute_stats
    usr_aux_output       => set_extravar_output
    usr_add_aux_names    => set_extravarnames_output
    usr_set_B0           => make_dipoleboy

    call mhd_activate()

    irhoav_   = var_set_extravar("rho_av", "rho_av")
    irho2av_  = var_set_extravar("rho2_av", "rho2_av")
    ivrav_    = var_set_extravar("vrad_av", "vrad_av")
    ivpolav_  = var_set_extravar("vtheta_av", "vtheta_av")
    ivr2av_   = var_set_extravar("vrad2_av", "vrad2_av")
    ivpol2av_ = var_set_extravar("vtheta2_av", "vtheta2_av")
    irhovrav_ = var_set_extravar("rho_vrad_av", "rho_vrad_av")

    if (mhd_energy) itav_ = var_set_extravar("twind_av", "twind_av")

  end subroutine usr_init

!==============================================================================
! Read in the usr.par file with the problem specific list.
!==============================================================================
  subroutine usr_params_read(files)

    ! Subroutine argument
    character(len=*), intent(in) :: files(:)

    ! Local variable
    integer :: n
    !--------------------------------------------------------------------------

    namelist /star_list/ mstar_sol, rstar_sol, twind_cgs, bpole_cgs, etastar, &
         rhosurf_cgs, timestat_cgs, Wrot, cakfile

    do n = 1,size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, star_list, end=111)
       111 close(unitpar)
    end do

  end subroutine usr_params_read

!==============================================================================
! Compute some quantities of interest (in CGS) before making unitless.
!==============================================================================
  subroutine initglobaldata_usr

    ! Local variables
    real(8) :: unit_ggrav, unit_lum, unit_mass
    real(8) :: lstar_cgs, mstar_cgs, rstar_cgs, vesc_cgs, mdot_cgs
    real(8) :: vinf_cgs, csound_cgs, logg_cgs, logge_cgs, heff_cgs, mumol
    real(8) :: vrot_cgs, vrotc_cgs, ralf_sol, rkep_sol, resc_sol
    real(8) :: lstar, mdot, twind, gammae, vesc, vinf, pthsurf
    !--------------------------------------------------------------------------

    mstar_cgs = mstar_sol * const_MSun
    rstar_cgs = rstar_sol * const_RSun

    ! Stellar structure
    lstar_cgs  = 4.0d0*dpi * rstar_cgs**2.0d0 * const_sigma * twind_cgs**4.0d0
    gammae     = const_kappae * lstar_cgs &
         / (4.0d0*dpi * const_G * mstar_cgs * const_c)
    logg_cgs   = log10(const_G * mstar_cgs / rstar_cgs**2.0d0)
    logge_cgs  = logg_cgs + log10(1.0d0 - gammae)
    mumol      = (1.0d0 + 4.0d0*He_abundance) / (2.0d0 + 3.0d0*He_abundance)
    csound_cgs = sqrt(twind_cgs * kB_cgs / (mumol * mp_cgs))
    heff_cgs   = csound_cgs**2.0d0 / 10.0d0**logge_cgs
    vrotc_cgs  = sqrt(const_G * mstar_cgs * (1.0d0 - gammae) / rstar_cgs)
    vrot_cgs   = vrotc_cgs * Wrot

    ! Wind quantities in CAK theory
    vesc_cgs  = sqrt( 2.0d0 * const_G * mstar_cgs * (1.0d0 - gammae) &
         / rstar_cgs )
    vinf_cgs  = vesc_cgs * sqrt(cak_alpha / (1.0d0 - cak_alpha))
    mdot_cgs  = lstar_cgs / const_c**2.0d0 * cak_alpha / (1.0d0 - cak_alpha) &
         * (gayley_qbar * gammae / (1.0d0 - gammae))**( (1.0d0 - cak_alpha) &
         / cak_alpha )

    call set_cak_force_norm(rstar_sol*const_RSun, twind_cgs)

    ! Bpole given and etastar computed or vice versa
    if (bpole_cgs > 0.0d0 .and. etastar < 0.0d0) then
      etastar = ((bpole_cgs/2.0d0)**2.0d0 * rstar_cgs**2.0d0) &
           / (mdot_cgs * vinf_cgs)
    elseif (etastar > 0.0d0 .and. bpole_cgs < 0.0d0) then
      bpole_cgs = 2.0d0 * sqrt(mdot_cgs * vinf_cgs * etastar/rstar_cgs**2.0d0)
    else
      call mpistop('Set stellar Bfield with bpole or etastar in .par file.')
    endif

    ! Compute Alfven, Kepler, and escape radius
    ralf_sol = 1.0d0 + (etastar + 0.25d0)**0.25d0 - 0.25d0**0.25d0
    rkep_sol = Wrot**(-2.0d0/3.0d0)
    resc_sol = 2.0d0**(1.0d0/3.0d0) * rkep_sol

    if (typedivbfix == 'ct') then
      call mpistop('Constrained Transport for Bfield not yet implemented.')
    endif

    if (mhd_energy .and. (.not.mhd_radiative_cooling)) then
      call mpistop('No support for adiabatic cooling. Add radiative cooling.')
    endif

    ! Code units
    unit_ggrav   = unit_density * unit_time**2.0d0
    unit_lum     = unit_density * unit_length**5.0d0 / unit_time**3.0d0
    unit_mass    = unit_density * unit_length**3.0d0

    lstar    = lstar_cgs / unit_lum
    mstar    = mstar_cgs / (unit_density * unit_length**3.0d0)
    rstar    = rstar_cgs / unit_length
    twind    = twind_cgs / unit_temperature
    bpole    = bpole_cgs / unit_magneticfield
    ralf     = ralf_sol * rstar
    rhosurf  = rhosurf_cgs / unit_density
    pthsurf  = rhosurf * twind
    mdot     = mdot_cgs / (unit_mass / unit_time)
    csound   = csound_cgs / unit_velocity
    clight   = const_c / unit_velocity
    vesc     = vesc_cgs / unit_velocity
    vinf     = vesc * sqrt(cak_alpha / (1.0d0 - cak_alpha))
    vrot     = vrot_cgs / unit_velocity
    gmstar   = const_G * unit_ggrav * mstar
    timestat = timestat_cgs / unit_time

    mhd_adiab = pthsurf / rhosurf**mhd_gamma

    if (mype == 0 .and. .not.convert) then
      print*, '======================'
      print*, '   Unity quantities   '
      print*, '======================'
      print*, 'unit length        = ', unit_length
      print*, 'unit density       = ', unit_density
      print*, 'unit velocity      = ', unit_velocity
      print*, 'unit numberdensity = ', unit_numberdensity
      print*, 'unit pressure      = ', unit_pressure
      print*, 'unit temperature   = ', unit_temperature
      print*, 'unit magneticfield = ', unit_magneticfield
      print*, 'unit time          = ', unit_time
      print*
      print*, '==============================================='
      print*, '   Stellar and wind parameters in CGS units    '
      print*, '==============================================='
      print*, 'L/Lsun                 = ', lstar_cgs / const_LSun
      print*, 'M/Msun                 = ', mstar_cgs / const_MSun
      print*, 'R/Rsun                 = ', rstar_cgs / const_RSun
      print*, 'Twind                  = ', twind_cgs
      print*, 'Polar magnetic field   = ', bpole_cgs
      print*, 'Wind confinement eta   = ', etastar
      print*, 'Ralf/Rstar             = ', ralf_sol
      print*, 'Rkep/Rstar             = ', rkep_sol
      print*, 'Resc/Rstar             = ', resc_sol
      print*, 'Mean molecular weight  = ', mumol
      print*, 'log(g)                 = ', logg_cgs
      print*, 'eff. log(g)            = ', logge_cgs
      print*, 'eff. scale height heff = ', heff_cgs
      print*, 'heff/Rstar             = ', heff_cgs / rstar_cgs
      print*, 'W (vrot/vrotc)         = ', Wrot
      print*, 'critical vrot          = ', vrotc_cgs
      print*, 'vrot                   = ', vrot_cgs
      print*, 'Eddington gamma        = ', gammae
      print*, 'adiabatic gamma        = ', mhd_gamma
      print*, 'adiabatic constant     = ', mhd_adiab
      print*, 'isothermal asound      = ', csound_cgs
      print*, 'eff. vesc              = ', vesc_cgs
      print*, 'CAK vinf               = ', vinf_cgs
      print*, 'FD vinf                = ', 3.0d0 * vesc_cgs
      print*, 'analytic Mdot CAK      = ', mdot_cgs * const_years / const_MSun
      print*, '... with FD correction = ', &
           mdot_cgs / (1.0d0 + cak_alpha)**(1.0d0/cak_alpha) &
           * const_years / const_MSun
      print*
      print*, '========================================'
      print*, '  Dimensionless computation quantities  '
      print*, '========================================'
      print*, 'Extra computed unit quantities:'
      print*, '   unit Lum     = ', unit_lum
      print*, '   unit Mass    = ', unit_mass
      print*, '   unit Grav    = ', unit_ggrav
      print*, 'Lstar        = ', lstar
      print*, 'Mstar        = ', mstar
      print*, 'Rstar        = ', rstar
      print*, 'Twind        = ', twind
      print*, 'rhosurface   = ', rhosurf
      print*, 'Mdot         = ', mdot
      print*, 'csound       = ', csound
      print*, 'eff. vesc    = ', vesc
      print*, 'vinf         = ', vinf
      print*, 'clight       = ', clight
      print*, 'timestat     = ', timestat
      print*, 'vrot         = ', vrot
      print*, 'Bpole        = ', bpole
    endif

  end subroutine initglobaldata_usr

!==============================================================================
! Initial conditions start from spherically symmetric 1-D relaxed CAK wind
! imposed onto a dipole stellar magnetic field.
!==============================================================================
  subroutine initial_conditions(ixI^L, ixO^L, w, x)

    use mod_initcak

    ! Subroutine arguments
    integer, intent(in)    :: ixI^L, ixO^L
    real(8), intent(in)    :: x(ixI^S,1:ndim)
    real(8), intent(inout) :: w(ixI^S,1:nw)

    ! Local arguments
    integer :: ir
    real(8) :: local_r, interp_rho, interp_vr
    logical :: first=.true.
    !--------------------------------------------------------------------------

    if (first) then
      first = .false.
      call read_initial_oned_cakwind(cakfile)
    endif

    do ir = ixOmin1,ixOmax1
      local_r = x(ir,ixOmin2,1)
      call interp_oned_cakwind_on_grid(local_r,interp_rho,interp_vr)

      w(ir^%1ixO^S,rho_)   = interp_rho
      w(ir^%1ixO^S,mom(1)) = interp_vr
    enddo

    ! Initial rigid body rotation
    w(ixO^S,mom(2)) = 0.0d0
    w(ixO^S,mom(3)) = vrot * sin(x(ixO^S,2)) * rstar**2.0d0 / x(ixO^S,1)

    ! Setup dipole magnetic field based on Tanaka splitting or regular
    if (B0field) then
      w(ixO^S,mag(:)) = 0.0d0
    else
      w(ixO^S,mag(1)) = bpole * cos(x(ixO^S,2)) * (rstar / x(ixO^S,1))**3.0d0
      w(ixO^S,mag(2)) = 0.5d0 * bpole * sin(x(ixO^S,2)) &
           * (rstar / x(ixO^S,1))**3.0d0
      w(ixO^S,mag(3)) = 0.0d0
    endif

    ! Initial pressure is isothermal
    if (mhd_energy) w(ixO^S,p_) = w(ixO^S,rho_)

    if (mhd_glm) w(ixO^S,psi_) = 0.0d0

    call mhd_to_conserved(ixI^L,ixO^L,w,x)

    w(ixO^S,nw-nwextra+1:nw) = 0.0d0

  end subroutine initial_conditions

!==============================================================================
! Special user boundary conditions at inner radial boundary:
!   vr, Btheta, Bphi (extrapolated); rho, vtheta, vphi, Br (fixed).
!==============================================================================
  subroutine special_bound(qt, ixI^L, ixB^L, iB, w, x)

    ! Subroutine arguments
    integer, intent(in)    :: ixI^L, ixB^L, iB
    real(8), intent(in)    :: qt, x(ixI^S,1:ndim)
    real(8), intent(inout) :: w(ixI^S,1:nw)

    ! Local variable
    integer :: ir
    !--------------------------------------------------------------------------

    select case (iB)
    case(1)

      call mhd_to_primitive(ixI^L,ixI^L,w,x)

      w(ixB^S,rho_) = rhosurf

      ! vr (2nd order accurate constant slope extrapolation in 1st ghost cell)
      do ir = ixBmax1,ixBmin1,-1
        if (ir == ixBmax1) then
          w(ir^%1ixB^S,mom(1)) = 1.0d0/3.0d0 &
               * (-w(ir+2^%1ixB^S,mom(1)) + 4.0d0*w(ir+1^%1ixB^S,mom(1)))
        else
          w(ir^%1ixB^S,mom(1)) = w(ir+1^%1ixB^S,mom(1))
        endif
      enddo

      ! Prohibit ghosts to be supersonic, also avoid overloading too much
      w(ixB^S,mom(1)) = min(w(ixB^S,mom(1)),  csound)
      w(ixB^S,mom(1)) = max(w(ixB^S,mom(1)), -csound)

      w(ixB^S,mom(3)) = vrot * sin(x(ixB^S,2))

      !******************
      ! Tanaka splitting
      !******************
      if (B0field) then
        do ir = ixBmax1,ixBmin1,-1
          ! r*r*(B0r + delta Br) = constant
          w(ir^%1ixB^S,mag(1)) = ( bpole * cos(x(ixBmax1+1^%1ixB^S,2))     &
               +                   w(ixBmax1+1^%1ixB^S,mag(1)) )           &
               * (rstar/x(ir^%1ixB^S,1))**2.0d0 - block%B0(ir^%1ixB^S,1,0)

          ! delta Btheta
          w(ir^%1ixB^S,mag(2)) = 1.0d0/3.0d0 &
               * (-w(ir+2^%1ixB^S,mag(2)) + 4.0d0*w(ir+1^%1ixB^S,mag(2)))

          ! delta Bphi
          w(ir^%1ixB^S,mag(3)) = 1.0d0/3.0d0 &
               * (-w(ir+2^%1ixB^S,mag(3)) + 4.0d0*w(ir+1^%1ixB^S,mag(3)))
        enddo
      !*********
      ! Regular
      !*********
      else
        do ir = ixBmax1,ixBmin1,-1
          ! r*r*Br = constant
          w(ir^%1ixB^S,mag(1)) = bpole * cos(x(ixBmax1+1^%1ixB^S,2)) &
               * (rstar/x(ir^%1ixB^S,1))**2.0d0

          ! Btheta
          w(ir^%1ixB^S,mag(2)) = 1.0d0/3.0d0 &
               * (-w(ir+2^%1ixB^S,mag(2)) + 4.0d0*w(ir+1^%1ixB^S,mag(2)))

          ! Bphi
          w(ir^%1ixB^S,mag(3)) = 1.0d0/3.0d0 &
               * (-w(ir+2^%1ixB^S,mag(3)) + 4.0d0*w(ir+1^%1ixB^S,mag(3)))
        enddo
      endif

      ! Enforce poloidal flow along magnetic field; outside magnetosphere from
      ! induction equation, inside magnetosphere put at zero
      ! Flo: need in strong confinement models to avoid fountain flows at pole
      if (etastar > 1.0d0) then
        if (B0field) then
          where ( abs(0.5d0*dpi - x(ixB^S,2)) >= asin(sqrt(rstar/ralf)) )
            w(ixB^S,mom(2)) = w(ixB^S,mom(1)) &
                 * (block%B0(ixB^S,2,0) + w(ixB^S,mag(2))) &
                 / (block%B0(ixB^S,1,0) + w(ixB^S,mag(1)))
          elsewhere
            w(ixB^S,mom(2)) = 0.0d0
          endwhere
        else
          where ( abs(0.5d0*dpi - x(ixB^S,2)) >= asin(sqrt(rstar/ralf)) )
            w(ixB^S,mom(2)) = w(ixB^S,mom(1)) * w(ixB^S,mag(2))/w(ixB^S,mag(1))
          elsewhere
            w(ixB^S,mom(2)) = 0.0d0
          endwhere
        endif
      else
        w(ixB^S,mom(2)) = 0.0d0
      endif

      if (mhd_energy) w(ixB^S,p_) = mhd_adiab * w(ixB^S,rho_)**mhd_gamma

      if (mhd_glm) w(ixB^S,psi_) = 0.0d0

      call mhd_to_conserved(ixI^L,ixI^L,w,x)

    case default
      call mpistop("BC not specified")
    end select

  end subroutine special_bound

!==============================================================================
! Compute stellar gravity.
!==============================================================================
  subroutine stellar_gravity(ixI^L, ixO^L, wCT, x, gravity_field)

    ! Subroutine arguments
    integer, intent(in)  :: ixI^L, ixO^L
    real(8), intent(in)  :: x(ixI^S,1:ndim), wCT(ixI^S,1:nw)
    real(8), intent(out) :: gravity_field(ixI^S,ndim)
    !--------------------------------------------------------------------------

    gravity_field(ixO^S,:) = 0.0d0

    ! Only in radial direction
    gravity_field(ixO^S,1) = -gmstar / x(ixO^S,1)**2.0d0

  end subroutine stellar_gravity

!==============================================================================
! Routine computes the time-averaged statistical quantity <X> via:
!   <X>_i = <X>_i-1 + dt * X_i
! where <X> is the average of variable X and i','i-1' are the current and
! previous timestep. NOTE: every iteration (un)normalisation has to be done.
!==============================================================================
  subroutine compute_stats(igrid, level, ixI^L, ixO^L, qt, w, x)

    ! Subroutine arguments
    integer, intent(in)    :: igrid, level, ixI^L, ixO^L
    real(8), intent(in)    :: qt, x(ixI^S,1:ndim)
    real(8), intent(inout) :: w(ixI^S,1:nw)

    ! Local variables
    real(8) :: tnormp, tnormc
    !--------------------------------------------------------------------------

    ! Note: qt is just a placeholder for the 'global_time' variable
    if (qt < timestat) RETURN

    call mhd_to_primitive(ixI^L,ixO^L,w,x)

    ! Current ^(n+1) and previous ^(n) timestep normalisation weigths
    tnormc = qt + dt - timestat
    tnormp = qt - timestat

    ! Average density
    w(ixO^S,irhoav_) = w(ixO^S,irhoav_)*tnormp
    w(ixO^S,irhoav_) = w(ixO^S,irhoav_) + dt * w(ixO^S,rho_)
    w(ixO^S,irhoav_) = w(ixO^S,irhoav_)/tnormc

    ! Average mass density squared
    w(ixO^S,irho2av_) = w(ixO^S,irho2av_)*tnormp
    w(ixO^S,irho2av_) = w(ixO^S,irho2av_) + dt * w(ixO^S,rho_)**2.0d0
    w(ixO^S,irho2av_) = w(ixO^S,irho2av_)/tnormc

    ! Average radial velocity
    w(ixO^S,ivrav_) = w(ixO^S,ivrav_)*tnormp
    w(ixO^S,ivrav_) = w(ixO^S,ivrav_) + dt * w(ixO^S,mom(1))
    w(ixO^S,ivrav_) = w(ixO^S,ivrav_)/tnormc

    ! Average radial velocity squared
    w(ixO^S,ivr2av_) = w(ixO^S,ivr2av_)*tnormp
    w(ixO^S,ivr2av_) = w(ixO^S,ivr2av_) + dt * w(ixO^S,mom(1))**2.0d0
    w(ixO^S,ivr2av_) = w(ixO^S,ivr2av_)/tnormc

    ! Average radial momentum density (correlation mass density-velocity)
    w(ixO^S,irhovrav_) = w(ixO^S,irhovrav_)*tnormp
    w(ixO^S,irhovrav_) = w(ixO^S,irhovrav_) &
         + dt * w(ixO^S,rho_)*w(ixO^S,mom(1))
    w(ixO^S,irhovrav_) = w(ixO^S,irhovrav_)/tnormc

    ! Average polar velocity
    w(ixO^S,ivpolav_) = w(ixO^S,ivpolav_)*tnormp
    w(ixO^S,ivpolav_) = w(ixO^S,ivpolav_) + dt * w(ixO^S,mom(2))
    w(ixO^S,ivpolav_) = w(ixO^S,ivpolav_)/tnormc

    ! Average polar velocity squared
    w(ixO^S,ivpol2av_) = w(ixO^S,ivpol2av_)*tnormp
    w(ixO^S,ivpol2av_) = w(ixO^S,ivpol2av_) + dt * w(ixO^S,mom(2))**2.0d0
    w(ixO^S,ivpol2av_) = w(ixO^S,ivpol2av_)/tnormc

    ! Average wind temperature
    if (mhd_energy) then
      w(ixO^S,itav_) = w(ixO^S,itav_)*tnormp
      w(ixO^S,itav_) = w(ixO^S,itav_) + dt * w(ixO^S,p_)/w(ixO^S,rho_)
      w(ixO^S,itav_) = w(ixO^S,itav_)/tnormc
    endif

    call mhd_to_conserved(ixI^L,ixO^L,w,x)

  end subroutine compute_stats

!==============================================================================
! Computes and stores additional variables of interest in convert stage.
!==============================================================================
  subroutine set_extravar_output(ixI^L, ixO^L, w, x, normconv)

    ! Subroutine arguments
    integer, intent(in) :: ixI^L, ixO^L
    real(8), intent(in) :: x(ixI^S,1:ndim)
    real(8)             :: w(ixI^S,nw+nwauxio)
    real(8)             :: normconv(0:nw+nwauxio)

    ! Local variable
    real(8) :: divbboy(ixI^S)
    !--------------------------------------------------------------------------

    ! Output the Alfven speed by summing squared Bfield in each direction
    if (B0field) then
      w(ixO^S,nw+1) = sqrt( sum((w(ixO^S,mag(:)) &
           + block%B0(ixO^S,:,0))**2, dim=ndim+1) / w(ixO^S,rho_) )
    else
      w(ixO^S,nw+1) = sqrt(sum(w(ixO^S,mag(:))**2, dim=ndim+1)/w(ixO^S,rho_))
    endif

    ! Output divB, for Tanaka splitting this will be divB1
    call get_divb(w,ixI^L,ixO^L,divbboy)
    w(ixO^S,nw+2) = divbboy(ixO^S)

  end subroutine set_extravar_output

!==============================================================================
! Additional auxiliary io variables: Alfven velocity, divergence B.
!==============================================================================
  subroutine set_extravarnames_output(varnames)

    character(len=*) :: varnames
    !--------------------------------------------------------------------------

    varnames = 'vAlf divB'

  end subroutine set_extravarnames_output

!==============================================================================
! Add a steady (time-independent) potential dipole background field with polar
! magnetic field strength set by dimensionless 'bpole' variable.
!==============================================================================
  subroutine make_dipoleboy(ixI^L, ixO^L, x, wB0)

    ! Subroutine arguments
    integer, intent(in)    :: ixI^L, ixO^L
    real(8), intent(in)    :: x(ixI^S,1:ndim)
    real(8), intent(inout) :: wB0(ixI^S,1:ndir)
    !--------------------------------------------------------------------------

    wB0(ixI^S,1) = bpole * (rstar/x(ixI^S,1))**3.0d0 * cos(x(ixI^S,2))
    wB0(ixI^S,2) = 0.5d0*bpole * (rstar/x(ixI^S,1))**3.0d0 * sin(x(ixI^S,2))
    wB0(ixI^S,3) = 0.0d0

  end subroutine make_dipoleboy

end module mod_usr
