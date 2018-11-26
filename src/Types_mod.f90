module Types_mod
use, intrinsic :: iso_fortran_env
implicit none
! everything is private unless otherwise stated
private
public :: SP, DP, SI, DI

integer, parameter :: SP = REAL32
integer, parameter :: DP = REAL64
integer, parameter :: SI = INT32
integer, parameter :: DI = INT64

! put other kind parameters here from Listing 1
 contains
   
end module Types_mod
