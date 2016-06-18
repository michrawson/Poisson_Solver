subroutine test_function(n01,n02,n03,a_gauss,hx,hy,hz,density,potential)
  implicit none
  integer, intent(in) :: n01,n02,n03
  real(kind=8), intent(in) :: a_gauss,hx,hy,hz
  real(kind=8), dimension(n01,n02,n03), intent(out) :: density,potential

  !local variables
  integer :: i1,i2,i3
  real(kind=8) :: x1,x2,x3,pi,a2,derf,factor,r,r2

     !grid for the free BC case
     !hgrid=max(hx,hy,hz)

     pi = 4.d0*atan(1.d0)
     a2 = a_gauss**2

     !Normalization
     factor = -4.d0/(a_gauss*a2*sqrt(pi))
     !gaussian function
     do i3=1,n03
        x3 = hz*real(i3-n03/2,kind=8)
        do i2=1,n02
           x2 = hy*real(i2-n02/2,kind=8)
           do i1=1,n01
              x1 = hx*real(i1-n01/2,kind=8)
              r2 = x1*x1+x2*x2+x3*x3
              density(i1,i2,i3) = factor*exp(-r2/a2)
              r = sqrt(r2)
              !Potential from a gaussian
              if (r == 0.d0) then
                 potential(i1,i2,i3) = 2.d0/(sqrt(pi)*a_gauss)
              else
                 potential(i1,i2,i3) = derf(r/a_gauss)/r
              end if
           end do
        end do
     end do
end subroutine test_function


subroutine test_function2(n01,n02,n03,a_gauss,hx,hy,hz,density,potential)
  implicit none
  integer, intent(in) :: n01,n02,n03
  real(kind=8), intent(in) :: a_gauss,hx,hy,hz
  real(kind=8), dimension(n01,n02,n03), intent(out) :: density,potential

  !local variables
  integer :: i1,i2,i3
  real(kind=8) :: x1,x2,x3,pi,a2,derf,factor,r,r2

     !grid for the free BC case
     !hgrid=max(hx,hy,hz)

     pi = 4.d0*atan(1.d0)
     a2 = a_gauss**2

     !Normalization
     factor = 1.0
     !gaussian function
     do i3=1,n03
        x3 = hz*real(i3-n03/2,kind=8)
        do i2=1,n02
           x2 = hy*real(i2-n02/2,kind=8)
           do i1=1,n01
              x1 = hx*real(i1-n01/2,kind=8)
              r2 = x1*x1+x2*x2+x3*x3

              density(i1,i2,i3) = factor*exp(-r2/a2)* &
                        (-2.0/a2*3.0 + 4.0/a2/a2*(x1*x1+x2*x2+x3*x3))

              potential(i1,i2,i3) = factor*exp(-r2/a2)

           end do
        end do
     end do
end subroutine test_function2


program main
    implicit none
    integer, parameter                       :: n = 32
    real ( kind = 8 ), dimension(n,n,n) :: density,potential_true,potential
    real(kind=8) :: hx,hy,hz,rang

    rang = 20.0d0
    hx=rang/real(n,kind=8)
    hy=rang/real(n,kind=8)
    hz=rang/real(n,kind=8)

    call test_function2(n,n,n,2.0,hx,hy,hz,density,potential_true)

    print *,"density edges",density(n/2,n/2,n),density(n,n/2,n/2),&
                density(n/2,n,n/2),density(n/2,n/2,1),&
                density(n/2,1,n/2),density(1,n/2,n/2)

    call solve( n, hx, density, potential )

    print *,"max_diff",maxval(abs(potential(:,:,:)-potential_true(:,:,:)))

    print *,"diff",(abs(potential(n/2,n/2,n/2)-potential_true(n/2,n/2,n/2)))

end program main
