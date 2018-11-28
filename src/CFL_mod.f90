module CFL_mod
!> this module calculates the CFL number
  use :: types_mod, only: dp
  implicit none

  contains

!> calculates the CFL number
!> \( \text{CFL} = \kappa\frac{\Delta t}{\Delta x^2} \)
  subroutine fd1d_heat_explicit_cfl(k, t_num, t_min, t_max, x_num, x_min, &
    x_max, cfl)

    implicit none

    real (kind=dp), intent(in) :: k
    integer, intent(in) :: t_num
    real (kind=dp), intent(in) :: t_min
    real (kind=dp), intent(in) :: t_max
    integer, intent(in) :: x_num
    real (kind=dp), intent(in) :: x_min
    real (kind=dp), intent(in) :: x_max
!> calculated CFL number
    real (kind=dp), intent(out) :: cfl
!> the heat constant \( \kappa \)

    real (kind=dp) :: dx
    real (kind=dp) :: dt

    dx = (x_max-x_min)/real(x_num-1, kind=dp)
    dt = (t_max-t_min)/real(t_num-1, kind=dp)

    cfl = k*dt/dx/dx

    write (*, '(a)') ' '
    write (*, '(a,g14.6)') '  CFL stability criterion value = ', cfl

  end subroutine

end module CFL_mod
