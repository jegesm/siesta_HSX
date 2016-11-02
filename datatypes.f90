module datatypes

use parameters, only: dp, sp

implicit none

type,public :: hsx_t
  integer :: nspecies
  integer :: no_s
  integer, DIMENSION(:), allocatable :: no 
  integer, DIMENSION(:), allocatable :: iaorb 
  integer, allocatable :: iphorb(:) 
  integer, allocatable :: nquant(:,:)
  integer, allocatable :: lquant(:,:) 
  integer, allocatable :: zeta(:,:) 
  integer :: nspin,nh
  integer :: no_u
  integer :: na_u

  integer, allocatable  :: numh(:)  
  integer, allocatable  :: listhptr(:)  
  integer, allocatable  :: listh(:)   
  integer, allocatable  :: indxuo(:) 
  real(dp), allocatable :: hamilt(:,:) 
  real(dp), allocatable :: Sover(:) 
  real(dp), allocatable :: xij(:,:) 
  integer, allocatable  :: isa(:) 
  real(dp), allocatable :: zval(:) 
  real(dp)          :: qtot, temp           !  fossils
!  character(len=5)  :: label(99)

end type

type,public :: arrayvm
  real(dp), allocatable :: xij(:,:) 
  
end type

contains

!subroutine init_hsx_t(dertype, m, n)
!    type(hsx_t), INTENT(inout) :: dertype
!    INTEGER(4), INTENT(in) :: m, n
!    allocate(dertype%nquant(m,n))
!    allocate(dertype%psi(m,n))
!end subroutine init_hsx_t

!subroutine destroy_hsx_t(dertype)
!    type(hsx_t), INTENT(inout) :: dertype
!    if (allocated(dertype%chi)) deallocate(dertype%chi)
!    if (allocated(dertype%psi)) deallocate(dertype%psi)
!end subroutine destroy_hsx_t


end module datatypes
