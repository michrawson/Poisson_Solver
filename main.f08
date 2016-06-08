
program main
    implicit none
    integer, parameter                       :: n = 128
    real ( kind = 8 ), dimension(n,n,n)      :: rho
    real ( kind = 8 ), dimension(n,n,n)      :: V

    rho=0
    V=0

    call solver( n, rho, V )

    print *,"V",V(1,1,1)

end program main
