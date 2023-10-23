!==============================================================================
! Modified version of AMRVAC mod_oneblock module specifically tailored to my
! problem to read in a simulated 1-D CAK wind model and interpolate it on the
! 2-D or 3-D simulation grid.
! NOTE: if input grid < simulation grid the outer wind will be held constant at
!       values of outermost density and radial velocity of input wind.
!
! Works in 1-D, 2-D, 3-D (M)HD radiation-driven wind models.
!
! To make AMRVAC aware of user-written module include in 'local.make' file:
!   mod_amrvac: mod_initcak.o
!   amrvac: mod_initcak.o
!
! Coded up by Flo for his KU Leuven PhD thesis 2018-2022
!==============================================================================
module mod_initcak

  ! Make some global variables and AMRVAC stuff available to this module
  use mod_amrvac

  implicit none

  ! Global variables to be set from reading input file
  integer              :: nrcells
  real(8), allocatable :: woneblock(:,:), xoneblock(:)

contains

!==============================================================================
! Read in a relaxed 1D CAK profile stored in a .blk file that is produced
! with the CAKwind_1d code -> assumed input format: rho, vr, gcak, fdfac
!==============================================================================
  subroutine read_initial_oned_cakwind(filename)

    ! Subroutine argument
    character(len=*), intent(in) :: filename

    ! Local variables
    integer :: ir, unit=69
    real(8) :: tmp, tmp1
    logical :: alive

    !*************************
    ! Master does the reading
    !*************************
    if (mype == 0) then
      inquire(file=filename, exist=alive)

      if (alive) then
        open(unit, file=filename, status='unknown')
      else
        call mpistop('Input file you want to use cannot be found!')
      endif

      print*,'Wind initialisation with input file: ', filename

      ! The header information:
      read(unit,*) ! skip variables
      read(unit,*) nrcells
      read(unit,*) ! skip time information

      if (nrcells < domain_nx1) then
        print*,'Input grid contains less radial cells than simulation grid'
        print*,'NOTE: constant extrapolation rho, vr done in outer wind'
      endif

      ! Allocate and read the grid and variables
      allocate(xoneblock(nrcells))
      allocate(woneblock(nrcells,1:2))

      do ir = 1,nrcells
        read(unit,*) xoneblock(ir), woneblock(ir,:), tmp, tmp1
      enddo

      close(unit)
    endif

    call MPI_BARRIER(icomm,ierrmpi)

    !****************************
    ! Broadcast what mype=0 read
    !****************************
    if (npe > 1) then
      call MPI_BCAST(nrcells,1,MPI_INTEGER,0,icomm,ierrmpi)

      if (mype /= 0) then
        allocate(xoneblock(nrcells))
        allocate(woneblock(nrcells,1:2))
      endif

      call MPI_BCAST(xoneblock,nrcells,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
      call MPI_BCAST(woneblock,nrcells*2,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
    endif

  end subroutine read_initial_oned_cakwind

!==============================================================================
! Linear interpolation of 1D CAK profile from .blk file onto simulation grid
!==============================================================================
  subroutine interp_oned_cakwind_on_grid(local_r,interp_rho,interp_vr)

    ! Subroutine arguments
    real(8), intent(in)  :: local_r
    real(8), intent(out) :: interp_rho, interp_vr

    ! Local arguments
    integer :: ic1, ic11, ic21
    real(8) :: xd1, rp

    rp = local_r

    ! Find correct block index to get correct value
    ic1 = minloc(abs(rp - xoneblock(:)), dim=1, mask=.true.)

    if (xoneblock(ic1) <= rp) then
      ic11 = ic1
    else
      ic11 = ic1 - 1
    endif
    ic21 = ic11 + 1

    ! Constant extrapolation when outside the range of input radial grid
    if (ic11 <= 1) then
      ic11 = 1
      ic21 = ic11 + 1
      rp   = xoneblock(ic1)
    endif

    if (ic21 >= nrcells) then
      ic21 = ic21
      ic11 = ic21 - 1
      rp   = xoneblock(ic11)
    endif

    ! Interpolate
    xd1        = (rp - xoneblock(ic11)) / (xoneblock(ic21) - xoneblock(ic11))
    interp_rho = woneblock(ic11,1) * (1.0d0 - xd1) + woneblock(ic21,1) * xd1
    interp_vr  = woneblock(ic11,2) * (1.0d0 - xd1) + woneblock(ic21,2) * xd1

  end subroutine interp_oned_cakwind_on_grid

end module mod_initcak
