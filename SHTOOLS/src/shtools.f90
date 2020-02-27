subroutine C_SHExpandDH ( &
     grid, n, cilm, lmax, norm, sampling, csphase, &
     lmax_calc, exitstatus)
  use shtools
  implicit none
  integer, parameter :: dp = selected_real_kind(p=15)
  integer :: n
  integer :: norm, sampling, csphase, lmax_calc
  real(dp) :: grid(n, sampling*n)
  real(dp) :: cilm(2, lmax_calc+1, lmax_calc+1)
  integer :: lmax
  integer :: exitstatus

  call SHExpandDH( &
       grid, n, cilm, lmax, norm, sampling, csphase, &
       lmax_calc, exitstatus)
end subroutine C_SHExpandDH

subroutine C_MakeGridDH( &
     griddh, n, cilm, lmax, norm, sampling, csphase, &
     lmax_calc, exitstatus)
  use shtools
  implicit none
  integer, parameter :: dp = selected_real_kind(p=15)
  integer :: lmax
  integer :: norm, sampling, csphase, lmax_calc
  real(dp) :: cilm(2, lmax_calc+1, lmax_calc+1)
  real(dp) :: griddh(2*(lmax+1), sampling*(2*(lmax+1)))
  integer :: n
  integer :: exitstatus

  call MakeGridDH( &
       griddh, n, cilm, lmax, norm, sampling, csphase, &
       lmax_calc, exitstatus)
end subroutine C_MakeGridDH
