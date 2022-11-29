    program diracfinder
    use constants, only : twopi
    use parameters 
    implicit NONE
    call init_mpi(ierr)
    call init_param
    call construct_hamiltonian
    call export_data
    call mpi_finalise(ierr)
    end program diracfinder