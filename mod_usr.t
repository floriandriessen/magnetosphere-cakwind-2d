!==============================================================================
! 2.5D model for simulating a CAK line-driven wind interacting with a stellar
! magnetosphere (dynamical or centrifugal) of a magnetic OB-star.
! Setup inspired by the original work of ud-Doula & Owocki (2002), ApJ 576.
!
! This module relies on my developed mod_cak_force module in MPI-AMRVAC and a
! separate routine to read in 1-D smooth CAK wind initial conditions.
!
! Coded up by Florian Driessen as part of the KU Leuven PhD thesis 2018-2022.
!
! HISTORY
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
!
! (November 2022) -- Flo
!   > update to amrvac v.3 using build in CAK radiation force module
!==============================================================================
module mod_usr

  use mod_mhd
  use mod_cak_force, ONLY: set_cak_force_norm, cak_alpha, gayley_qbar
  use mod_constants, ONLY: const_c, const_G, const_LSun, const_MSun, &
       const_RSun, const_kappae, const_sigma, mp_cgs, kB_cgs, const_years

  implicit none

  ! User input parameters
  real(8) :: mstar, rstar, twind, rhobound, timestat, Wrot
  real(8) :: bpole=-99.0, etastar=-99.0
  character(len=99) :: cakfile

  ! Extra parameters required in computation
  real(8) :: asound, Ggrav, vrot, ralf

  ! Dimensionless variables needed throughout computations
  real(8) :: dmstar, drstar, dbpole, drhobound, dasound, dclight
  real(8) :: dGgrav, dtimestat, dvrot

  ! Additional names for extra statistical variables in output
  integer :: my_rhoav, my_rho2av, my_vrav, my_vr2av, my_rhovrav
  integer :: my_vpolav, my_vpol2av, my_tav

contains

!==============================================================================
! This routine should set user methods, and activate the physics module
!==============================================================================
  subroutine usr_init

    call set_coordinate_system("spherical_2.5D")
    call usr_params_read(par_files)

    ! Choose independent normalisation units:
    unit_length      = rstar    ! cm
    unit_temperature = twind    ! K
    unit_density     = rhobound ! g cm^-3

    usr_set_parameters   => initglobaldata_usr
    usr_init_one_grid    => initial_conditions
    usr_special_bc       => special_bound
    usr_gravity          => stellar_gravity
    usr_process_adv_grid => compute_stats
    usr_aux_output       => set_extravar_output
    usr_add_aux_names    => set_extravarnames_output
    usr_set_B0           => make_dipoleboy
    usr_init_vector_potential => initvecpot_usr

    call MHD_activate()
    call set_cak_force_norm(rstar,twind)

    my_rhoav   = var_set_extravar("rho_av", "rho_av")
    my_rho2av  = var_set_extravar("rho2_av", "rho2_av")
    my_vrav    = var_set_extravar("vrad_av", "vrad_av")
    my_vpolav  = var_set_extravar("vtheta_av", "vtheta_av")
    my_vr2av   = var_set_extravar("vrad2_av", "vrad2_av")
    my_vpol2av = var_set_extravar("vtheta2_av", "vtheta2_av")
    my_rhovrav = var_set_extravar("rho_vrad_av", "rho_vrad_av")

    if (mhd_energy) my_tav = var_set_extravar("twind_av", "twind_av")

  end subroutine usr_init

!==============================================================================
! Read in the usr.par file with the problem specific list
!==============================================================================
  subroutine usr_params_read(files)

    ! Subroutine argument
    character(len=*), intent(in) :: files(:)

    ! Local variable
    integer :: n

    namelist /star_list/ mstar, rstar, twind, bpole, etastar, rhobound, &
         timestat, Wrot, cakfile

    do n = 1,size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, star_list, end=111)
       111 close(unitpar)
    end do

    ! Scale to cgs units
    mstar = mstar * const_MSun
    rstar = rstar * const_RSun

  end subroutine usr_params_read

!==============================================================================
! Compute some quantities of interest (in CGS) before making unitless
!==============================================================================
  subroutine initglobaldata_usr

    real(8) :: lstar, mumol, vesc, gammae, logg, logge, heff, vrotc
    real(8) :: rkep, resc, vinf, mdot

    ! Stellar structure
    lstar  = 4.0d0*dpi * rstar**2.0d0 * const_sigma * twind**4.0d0
    gammae = const_kappae * lstar/(4.0d0*dpi * const_G * mstar * const_c)
    logg   = log10(const_G * mstar/rstar**2.0d0)
    logge  = logg + log10(1.0d0 - gammae)
    mumol  = (1.0d0 + 4.0d0*He_abundance)/(2.0d0 + 3.0d0*He_abundance)
    asound = sqrt(twind * kB_cgs/(mumol * mp_cgs))
    heff   = asound**2.0d0 / 10.0d0**logge
    vrotc  = sqrt(const_G * mstar*(1.0d0 - gammae)/rstar)
    vrot   = vrotc * Wrot

    ! Wind quantities in CAK theory
    vesc = sqrt(2.0d0 * const_G * mstar * (1.0d0 - gammae)/rstar)
    vinf = vesc * sqrt(cak_alpha/(1.0d0 - cak_alpha))
    mdot = lstar/const_c**2.0d0 * cak_alpha/(1.0d0 - cak_alpha) &
         * (gayley_qbar*gammae/(1.0d0 - gammae))**((1.0d0-cak_alpha)/cak_alpha)

    ! Bpole given and etastar computed or vice versa
    if (bpole > 0.0d0 .and. etastar < 0.0d0) then
      etastar = ((bpole/2.0d0)**2.0d0 * rstar**2.0d0)/(mdot * vinf)
    elseif (etastar > 0.0d0 .and. bpole < 0.0d0) then
      bpole = 2.0d0 * sqrt(mdot * vinf * etastar/rstar**2.0d0)
    else
      call mpistop('Set stellar Bfield with bpole or etastar in .par file.')
    endif

    ! Compute Alfven, Kepler, and escape radius
    ralf = 1.0d0 + (etastar + 0.25d0)**0.25d0 - 0.25d0**0.25d0
    rkep = Wrot**(-2.0d0/3.0d0)
    resc = 2.0d0**(1.0d0/3.0d0) * rkep

    if (mhd_energy .and. (.not.mhd_radiative_cooling)) then
      call mpistop('No support for adiabatic cooling. Add radiative cooling.')
    endif

    if (mype == 0 .and..not.convert) then
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
      print*, '==============================================='
      print*, '   Stellar and wind parameters in CGS units    '
      print*, '==============================================='
      print*, 'L/Lsun                 = ', lstar/const_LSun
      print*, 'M/Msun                 = ', mstar/const_MSun
      print*, 'R/Rsun                 = ', rstar/const_RSun
      print*, 'Twind                  = ', twind
      print*, 'Polar magnetic field   = ', bpole
      print*, 'Wind confinement eta   = ', etastar
      print*, 'Ralf/Rstar             = ', ralf
      print*, 'Rkep/Rstar             = ', rkep
      print*, 'Resc/Rstar             = ', resc
      print*, 'Mean molecular weight  = ', mumol
      print*, 'log(g)                 = ', logg
      print*, 'eff. log(g)            = ', logge
      print*, 'heff/Rstar             = ', heff/rstar
      print*, 'W (vrot/vrotc)         = ', Wrot
      print*, 'critical vrot          = ', vrotc
      print*, 'vrot                   = ', vrot
      print*, 'Eddington gamma        = ', gammae
      print*, 'adiabatic gamma        = ', mhd_gamma
      print*, 'isothermal asound      = ', asound
      print*, 'eff. vesc              = ', vesc
      print*, 'CAK vinf               = ', vinf
      print*, 'FD vinf                = ', 3.0d0 * vesc
      print*, 'analytic Mdot CAK      = ', mdot * const_years/const_MSun
      print*, '... with FD correction = ', &
           mdot/(1.0d0 + cak_alpha)**(1.0d0/cak_alpha) * const_years/const_MSun
    endif

    ! Make some variables needed for computations dimensionless
    dmstar    = mstar/(unit_density * unit_length**3.0d0)
    drstar    = rstar/unit_length
    dbpole    = bpole/unit_magneticfield
    dasound   = asound/unit_velocity
    drhobound = rhobound/unit_density
    dGgrav    = const_G * unit_density * unit_time**2.0d0
    dtimestat = timestat/unit_time
    dvrot     = vrot/unit_velocity

    if (mype == 0 .and..not.convert) then
      print*, '========================================'
      print*, '  Dimensionless computation quantities  '
      print*, '========================================'
      print*, 'Mstar    = ', dmstar
      print*, 'Rstar    = ', drstar
      print*, 'rhobound = ', drhobound
      print*, 'asound   = ', dasound
      print*, 'Ggrav    = ', dGgrav
      print*, 'timestat = ', dtimestat
      print*, 'vrot     = ', dvrot
      print*, 'Bpole    = ', dbpole
    endif

  end subroutine initglobaldata_usr

!==============================================================================
! Initial conditions start from spherically symmetric 1-D relaxed CAK wind
! imposed onto a dipole stellar magnetic field
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

    w(ixO^S,mom(2)) = 0.0d0
    w(ixO^S,mom(3)) = dvrot * sin(x(ixO^S,2)) * drstar**2.0d0/x(ixO^S,1)

    ! Setup dipole field using Tanaka splitting, vector potential, or regular
    if (B0field) then
      w(ixO^S,mag(:)) = 0.0d0
    elseif (stagger_grid) then
      call b_from_vector_potential(ixGs^LL,ixI^L,ixO^L,block%ws,x)

      ! Update cell-center B from cell-face B
      call mhd_face_to_center(ixO^L,block)

      w(ixO^S,mag(3)) = 0.0d0
    else
      w(ixO^S,mag(1)) = dbpole * cos(x(ixO^S,2)) * (drstar/x(ixO^S,1))**3.0d0
      w(ixO^S,mag(2)) = 0.5d0 * dbpole * sin(x(ixO^S,2)) &
           * (drstar/x(ixO^S,1))**3.0d0
      w(ixO^S,mag(3)) = 0.0d0
    endif

    ! Pressure stratification is polytrope
    if (mhd_energy) then
      w(ixO^S,p_) = dasound**2.0 * drhobound &
           * (w(ixO^S,rho_)/drhobound)**mhd_gamma
    endif

    ! If using Dedner+(2002) divergence cleaning
    if (mhd_glm) w(ixO^S,psi_) = 0.0d0

    call mhd_to_conserved(ixI^L,ixO^L,w,x)

    ! Initialise extra vars at 0
    w(ixO^S,nw-nwextra+1:nw) = 0.0d0

  end subroutine initial_conditions

!==============================================================================
! Initialise the vector potential on edges for b_from_vectorpotential().
!==============================================================================
  subroutine initvecpot_usr(ixI^L, ixC^L, xC, A, idir)

    integer, intent(in)  :: ixI^L, ixC^L,idir
    real(8), intent(in)  :: xC(ixI^S,1:ndim)
    real(8), intent(out) :: A(ixI^S)

    if (idir == 3) then
      A(ixC^S) = 0.5d0*dbpole * (drstar/xC(ixC^S,1))**2.0d0 * sin(xC(ixC^S,2))
    else
      A(ixC^S) = 0.d0
    end if

  end subroutine initvecpot_usr

!==============================================================================
! Special user boundary conditions at inner radial boundary:
!   vr, Btheta, Bphi (extrapolated); rho, vtheta, vphi, Br (fixed)
!==============================================================================
  subroutine special_bound(qt, ixI^L, ixB^L, iB, w, x)

    ! Subroutine arguments
    integer, intent(in)    :: ixI^L, ixB^L, iB
    real(8), intent(in)    :: qt, x(ixI^S,1:ndim)
    real(8), intent(inout) :: w(ixI^S,1:nw)

    ! Local variable
    integer :: ir, ixBs^L

    select case (iB)
    case(1)

      call mhd_to_primitive(ixI^L,ixI^L,w,x)

      w(ixB^S,rho_) = drhobound

      ! vr (2nd order accurate constant slope extrapolation in 1st ghost cell)
      do ir = ixBmax1,ixBmin1,-1
        if (ir == ixBmax1) then
          w(ir^%1ixB^S,mom(1)) = 1.0d0/3.0d0 &
               * (-w(ir+2^%1ixB^S,mom(1)) + 4.0d0*w(ir+1^%1ixB^S,mom(1)))
        else
          w(ir^%1ixB^S,mom(1)) = w(ir+1^%1ixB^S,mom(1))
        endif
      enddo

      ! Prohibit radial velocity ghosts to be supersonic, and avoid overloading
      w(ixB^S,mom(1)) = min(w(ixB^S,mom(1)), dasound/5.0d0)
      w(ixB^S,mom(1)) = max(w(ixB^S,mom(1)), -dasound/5.0d0)

      w(ixB^S,mom(3)) = dvrot * sin(x(ixB^S,2)) * w(ixB^S,rho_)

      !******************
      ! Tanaka splitting
      !******************
      if (B0field) then
        do ir = ixBmax1,ixBmin1,-1
          ! r*r*(B0r + delta Br) = constant
          w(ir^%1ixB^S,mag(1)) = ( dbpole * cos(x(ixBmax1+1^%1ixB^S,2))    &
               +                   w(ixBmax1+1^%1ixB^S,mag(1)) )           &
               * (drstar/x(ir^%1ixB^S,1))**2.0d0 - block%B0(ir^%1ixB^S,1,0)

          ! delta Btheta
          w(ir^%1ixB^S,mag(2)) = 1.0d0/3.0d0 &
               * (-w(ir+2^%1ixB^S,mag(2)) + 4.0d0*w(ir+1^%1ixB^S,mag(2)))

          ! delta Bphi
          w(ir^%1ixB^S,mag(3)) = 1.0d0/3.0d0 &
               * (-w(ir+2^%1ixB^S,mag(3)) + 4.0d0*w(ir+1^%1ixB^S,mag(3)))
        enddo
      !***********************
      ! Constrained Transport
      !***********************
      elseif (stagger_grid) then
        ! Br
        ! Index shift for staggered grid in radial direction
        ixBs^L=ixB^L-kr(1,^D);

        block%ws(ixBs^S,1) = dbpole * cos(x(ixB^S,2))      &
             * (drstar - 0.5d0 * block%dx(ixB^S,1))**2.0d0 &
             / (x(ixB^S,1) - 0.5d0 * block%dx(ixB^S,1))**2.0d0

        ! Btheta (Flo: CT contact method only works with etastar distinction)
        ixBsmax^D=ixBmax^D;
        ixBsmin^D=ixBmin^D-kr(^D,1);

        if (etastar >= 1.0d0) then
          do ir = ixBsmax1,ixBsmin1,-1
            block%ws(ir^%1ixBs^S,2) = 1.0d0/3.0d0 &
                 * (- block%ws(ir+2^%1ixBs^S,2)   &
                 +  4.d0 * block%ws(ir+1^%1ixBs^S,2))
          enddo
        else
          block%ws(ixBs^S,2) = 0.0d0
        endif

        call mhd_face_to_center(ixB^L,block)

        ! Bphi
        do ir = ixBmax1,ixBmin1,-1
          w(ir^%1ixB^S,mag(3)) = 1.0d0/3.0d0 * (- w(ir+2^%1ixB^S,mag(3)) &
               + 4.0d0 * w(ir+1^%1ixB^S,mag(3)))
        enddo
      !*********
      ! Regular
      !*********
      else
        do ir = ixBmax1,ixBmin1,-1
          ! r*r*Br = constant
          w(ir^%1ixB^S,mag(1)) = dbpole * cos(x(ixBmax1+1^%1ixB^S,2)) &
               * (drstar/x(ir^%1ixB^S,1))**2.0d0

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

      ! pthermal
      if (mhd_energy) then
        w(ixB^S,p_) = dasound**2.0 * drhobound &
             * (w(ixB^S,rho_)/drhobound)**mhd_gamma
      endif

      if (mhd_glm) w(ixB^S,psi_) = 0.0d0

      call mhd_to_conserved(ixI^L,ixI^L,w,x)

    ! case(2)

    !   do ir = ixBmin1,ixBmax1
    !     ! r*r*rho = constant
    !     w(ir^%1ixB^S,rho_) = w(ixBmin1-1^%1ixB^S,rho_) &
    !          * (x(ixBmin1-1^%1ixB^S,1) / x(ir^%1ixB^S,1))**2.0d0

    !     ! r*r*rho*vr = constant
    !     w(ir^%1ixB^S,mom(1)) = w(ixBmin1-1^%1ixB^S,mom(1)) &
    !          * (x(ixBmin1-1^%1ixB^S,1)/x(ir^%1ixB^S,1))**2.0d0

    !     ! rho*vtheta = constant
    !     w(ir^%1ixB^S,mom(2)) = w(ixBmin1-1^%1ixB^S,mom(2))

    !     ! r*vphi = constant
    !     w(ir^%1ixB^S,mom(3)) = &
    !          (w(ixBmin1-1^%1ixB^S,mom(3)) / w(ixBmin1-1^%1ixB^S,rho_)) &
    !          * (x(ixBmin1-1^%1ixB^S,1) / x(ir^%1ixB^S,1)) * w(ir^%1ixB^S,rho_)

    !     !******************
    !     ! Tanaka splitting
    !     !******************
    !     if (B0field) then
    !       ! r*r*(B0r + delta Br) = constant
    !       w(ir^%1ixB^S,mag(1)) = ( dbpole * cos(x(ixBmin1-1^%1ixB^S,2)) &
    !            * (drstar / x(ixBmin1-1^%1ixB^S,1))**3.0d0               &
    !            + w(ixBmin1-1^%1ixB^S,mag(1)) )                          &
    !            * (x(ixBmin1-1^%1ixB^S,1) / x(ir^%1ixB^S,1))**2.0d0      &
    !            - block%B0(ir^%1ixB^S,1,0)

    !       ! (B0theta + delta Btheta) = constant
    !       w(ir^%1ixB^S,mag(2)) = (0.5d0*dbpole * sin(x(ixBmin1-1^%1ixB^S,2)) &
    !            * (drstar/x(ixBmin1-1^%1ixB^S,1))**3.0d0                      &
    !            + w(ixBmin1-1^%1ixB^S,mag(2)))                                &
    !            - block%B0(ir^%1ixB^S,2,0)

    !       ! r*(B0phi + delta Bphi) = constant
    !       w(ir^%1ixB^S,mag(3)) = w(ixBmin1-1^%1ixB^S,mag(3)) &
    !            * x(ixBmin1-1^%1ixB^S,1) / x(ir^%1ixB^S,1)    &
    !            - block%B0(ir^%1ixB^S,3,0)
    !     !***********************
    !     ! Constrained Transport
    !     !***********************
    !     elseif (stagger_grid) then
    !       ! A*Br = constant
    !       block%ws(ir^%1ixB^S,1) = block%ws(ixBmin1-1^%1ixB^S,1) &
    !            * block%surfaceC(ixBmin1-1^%1ixB^S,1)             &
    !            / block%surfaceC(ir^%1ixB^S,1)

    !       ! Btheta = constant
    !       block%ws(ir^%1ixB^S,2) = block%ws(ixBmin1-1^%1ixB^S,2)

    !       call mhd_face_to_center(ixB^L,block)

    !       ! r*Bphi = constant
    !       w(ir^%1ixB^S,mag(3)) = w(ixBmin1-1^%1ixB^S,mag(3)) &
    !            * x(ixBmin1-1^%1ixB^S,1) / x(ir^%1ixB^S,1)
    !     !*********
    !     ! Regular
    !     !*********
    !     else
    !       ! r*r*Br = constant
    !       w(ir^%1ixB^S,mag(1)) = w(ixBmin1-1^%1ixB^S,mag(1)) &
    !            * (x(ixBmin1-1^%1ixB^S,1) / x(ir^%1ixB^S,1))**2.0d0

    !       ! Btheta = constant
    !       w(ir^%1ixB^S,mag(2)) = w(ixBmin1-1^%1ixB^S,mag(2))

    !       ! r*Bphi = constant
    !       w(ir^%1ixB^S,mag(3)) = w(ixBmin1-1^%1ixB^S,mag(3)) &
    !            * x(ixBmin1-1^%1ixB^S,1) / x(ir^%1ixB^S,1)
    !     endif
    !   enddo

    !   if (mhd_glm) w(ixB^S,psi_) = 0.0d0

    case default
      call mpistop("BC not specified")
    end select

  end subroutine special_bound

!==============================================================================
! Compute stellar gravity
!==============================================================================
  subroutine stellar_gravity(ixI^L, ixO^L, wCT, x, gravity_field)

    ! Subroutine arguments
    integer, intent(in)  :: ixI^L, ixO^L
    real(8), intent(in)  :: x(ixI^S,1:ndim), wCT(ixI^S,1:nw)
    real(8), intent(out) :: gravity_field(ixI^S,ndim)

    gravity_field(ixO^S,:) = 0.0d0

    ! Only in radial direction
    gravity_field(ixO^S,1) = -dGgrav * dmstar/x(ixO^S,1)**2.0d0

  end subroutine stellar_gravity

!==============================================================================
! Routine computes the time-averaged statistical quantity <X> via:
!   <X>_i = <X>_i-1 + dt * X_i
! where <X> is the average of variable X and i','i-1' are the current and
! previous timestep. NOTE: every iteration (un)normalisation has to be done
!==============================================================================
  subroutine compute_stats(igrid, level, ixI^L, ixO^L, qt, w, x)

    ! Subroutine arguments
    integer, intent(in)    :: igrid, level, ixI^L, ixO^L
    real(8), intent(in)    :: qt, x(ixI^S,1:ndim)
    real(8), intent(inout) :: w(ixI^S,1:nw)

    ! Local variables
    real(8) :: tnormp, tnormc

    ! Note: qt is just a placeholder for the 'global_time' variable
    if (qt < dtimestat) RETURN

    call mhd_to_primitive(ixI^L,ixO^L,w,x)

    ! Current ^(n+1) and previous ^(n) timestep normalisation weigths
    tnormc = qt + dt - dtimestat
    tnormp = qt - dtimestat

    ! Average density
    w(ixO^S,my_rhoav) = w(ixO^S,my_rhoav)*tnormp
    w(ixO^S,my_rhoav) = w(ixO^S,my_rhoav) + dt * w(ixO^S,rho_)
    w(ixO^S,my_rhoav) = w(ixO^S,my_rhoav)/tnormc

    ! Average mass density squared
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

    ! Average radial momentum density (correlation mass density-velocity)
    w(ixO^S,my_rhovrav) = w(ixO^S,my_rhovrav)*tnormp
    w(ixO^S,my_rhovrav) = w(ixO^S,my_rhovrav) &
         + dt * w(ixO^S,rho_)*w(ixO^S,mom(1))
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

  end subroutine compute_stats

!==============================================================================
! Computes and stores additional variables of interest in convert stage
!==============================================================================
  subroutine set_extravar_output(ixI^L, ixO^L, w, x, normconv)

    ! Subroutine arguments
    integer, intent(in) :: ixI^L, ixO^L
    real(8), intent(in) :: x(ixI^S,1:ndim)
    real(8)             :: w(ixI^S,nw+nwauxio)
    real(8)             :: normconv(0:nw+nwauxio)

    ! Local variable
    real(8) :: divbboy(ixI^S)

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
! Additional auxiliary io variables: Alfven velocity, divergence B
!==============================================================================
  subroutine set_extravarnames_output(varnames)

    character(len=*) :: varnames
    varnames = 'vAlf divB'

  end subroutine set_extravarnames_output

!==============================================================================
! Add a steady (time-independent) potential dipole background field with polar
! magnetic field strength set by dimensionless 'dbpole' variable
!==============================================================================
  subroutine make_dipoleboy(ixI^L, ixO^L, x, wB0)

    ! Subroutine arguments
    integer, intent(in)    :: ixI^L, ixO^L
    real(8), intent(in)    :: x(ixI^S,1:ndim)
    real(8), intent(inout) :: wB0(ixI^S,1:ndir)

    wB0(ixI^S,1) = dbpole * (drstar/x(ixI^S,1))**3.0d0 * cos(x(ixI^S,2))
    wB0(ixI^S,2) = 0.5d0*dbpole * (drstar/x(ixI^S,1))**3.0d0 * sin(x(ixI^S,2))
    wB0(ixI^S,3) = 0.0d0

  end subroutine make_dipoleboy

end module mod_usr
