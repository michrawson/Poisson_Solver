
program main
    implicit none
    integer, parameter                       :: n = 56
    real ( kind = 8 ), dimension(n,n,n) :: density,potential_true,potential
    real(kind=8) :: hx,hy,hz

    hx=10.0/real(n,kind=8)
    hy=10.0/real(n,kind=8)
    hz=10.0/real(n,kind=8)

    call test_functions(n,n,n,1.0,hx,hy,hz,density,potential_true)

    potential(:,:,:)=density(:,:,:)

    call solver( n, potential )

    print *,"max_diff",maxval(abs(potential(:,:,:)-potential_true(:,:,:)))

end program main
