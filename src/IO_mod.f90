module IO_mod
  use :: types_mod, only: dp
  use netcdf
  implicit none

  contains

   subroutine r8mat_write(output_filename, table, t, x)
    implicit none

    character (len=*), intent(in) :: output_filename
    real (kind=dp), dimension(:,:), intent(in) :: table
    real (kind=dp), dimension(:), intent(in) :: t, x

    integer :: j
    integer :: output_unit_id
    character (len=30) :: string

    integer :: m, n, nt, nx

    integer :: ierr, ncid, varid_table, varid_t, varid_x
    integer :: xtable_dimid, ytable_dimid, t_dimid, x_dimid

    m=size( table(:, :), 1 )
    n=size( table(:, :), 2 )

    nt = size(t)
    nx = size(x)

!    output_unit_id = 10
!    open (unit=output_unit_id, file=output_filename, status='replace')
!    write (string, '(a1,i8,a1,i8,a1,i8,a1)') '(', m, 'g', 24, '.', 16, ')'
!    do j = 1, n
!      write (output_unit_id, string) table(1:m, j)
!    end do
!    close (unit=output_unit_id)

! entering define mode
    ierr = NF90_CREATE(output_filename, NF90_CLOBBER, ncid )

    ierr = NF90_PUT_ATT( ncid, NF90_GLOBAL, "purpose", "Fortran Workshop")
    ierr = NF90_PUT_ATT( ncid, NF90_GLOBAL, "name", "Nicolas Sisourat")

    ierr = NF90_DEF_DIM( ncid,"xtable", m, xtable_dimid)
    ierr = NF90_DEF_DIM( ncid,"ytable", n, ytable_dimid)
    ierr = NF90_DEF_DIM( ncid,"tdim", nt, t_dimid)
    ierr = NF90_DEF_DIM( ncid,"xdim", nx, x_dimid)

    ierr = NF90_DEF_VAR( ncid,"table", NF90_DOUBLE, [ xtable_dimid, ytable_dimid ], varid_table)
    ierr = NF90_PUT_ATT( ncid, varid_table, "units", "celsius,")

    ierr = NF90_DEF_VAR( ncid,"tarr", NF90_DOUBLE, [ t_dimid ], varid_t)
    ierr = NF90_PUT_ATT( ncid, varid_t, "units", "seconds")
    ierr = NF90_DEF_VAR( ncid,"xarr", NF90_DOUBLE, [ x_dimid ], varid_x)
    ierr = NF90_PUT_ATT( ncid, varid_x, "units", "metres")
    ierr = NF90_ENDDEF( ncid )
! end define mode and enter data mode

    ierr = NF90_PUT_VAR( ncid, varid_table, table ) ! write data
    ierr = NF90_PUT_VAR( ncid, varid_t, t ) ! write data
    ierr = NF90_PUT_VAR( ncid, varid_x, x ) ! write data
    ierr = NF90_CLOSE( ncid )

  end subroutine

  subroutine r8vec_linspace(a_first, a_last, a)

    implicit none

    real (kind=dp), intent(in) :: a_first
    real (kind=dp), intent(in) :: a_last
    real (kind=dp), dimension(:), intent(out) :: a

    integer :: n
    integer :: i

    n = size(a)

    do i = 1, n
      a(i) = (real(n-i,kind=dp)*a_first+real(i-1,kind=dp)*a_last)/ &
        real(n-1, kind=dp)
    end do

  end subroutine

  subroutine r8vec_write(output_filename, x)

    implicit none

    character (len=*), intent(in) :: output_filename
    real (kind=dp), dimension(:), intent(in) :: x

    integer :: n
    integer :: m
    integer :: j
    integer :: output_unit_id

    n = size(x)

    output_unit_id = 11
    open (unit=output_unit_id, file=output_filename, status='replace')

    do j = 1, n
      write (output_unit_id, '(2x,g24.16)') x(j)
    end do

    close (unit=output_unit_id)
  end subroutine

end module IO_mod
