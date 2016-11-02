module parameters

implicit none
private
public :: dp, sp

INTEGER, PARAMETER :: dp=kind(1d0) ! double precision, np.float64 equivalent
INTEGER, PARAMETER :: sp=kind(1e0) ! single precision, np.float32 equivalent

contains

end module parameters

