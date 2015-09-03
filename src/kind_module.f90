!*****************************************************************************************
!> author: Jacob Williams
!
!  Kind definitions for the other modules.

	module kind_module
	
	use iso_fortran_env
	
	implicit none
	
	private
	
	! default real precision:
	
	!integer,parameter,public :: wp = real32   !! single precision
	integer,parameter,public :: wp = real64    !! double precision [default]
	!integer,parameter,public :: wp = real128  !! quad precision
	
	end module kind_module
!*****************************************************************************************