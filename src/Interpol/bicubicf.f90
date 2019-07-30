subroutine bicubic(f, fx, fy, fxy, x0, y0, derivs, res)

  implicit none

  !!-force zero indexing on the arrays for clarity with python
  real, dimension(0:3), intent(in) :: f, fx, fy, fxy
  real, dimension(0:15) :: xall, alpOld
  real, dimension(0:3,0:3) :: alp
  real, dimension(0:15, 0:15) :: A
  real, intent(in) :: x0
  real, intent(in) :: y0
  real, intent(out) :: res
  character(len=3), intent(in) :: derivs
  integer :: i, j

  xall(0:3) = f
  xall(4:7) = fx
  xall(8:11) = fy
  xall(12:15) = fxy
  A(0,:) = (/  1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
  A(1,:) = (/  0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
  A(2,:) = (/ -3, 3, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
  A(3,:) = (/  2,-2, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
  A(4,:) = (/  0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0/)
  A(5,:) = (/  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0/)
  A(6,:) = (/  0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0,-2,-1, 0, 0/)
  A(7,:) = (/  0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 1, 1, 0, 0/)
  A(8,:) = (/ -3, 0, 3, 0, 0, 0, 0, 0,-2, 0,-1, 0, 0, 0, 0, 0/)
  A(9,:) = (/  0, 0, 0, 0,-3, 0, 3, 0, 0, 0, 0, 0,-2, 0,-1, 0/)
  A(10,:) = (/ 9,-9,-9, 9, 6, 3,-6,-3, 6,-6, 3,-3, 4, 2, 2, 1/)
  A(11,:) = (/-6, 6, 6,-6,-3,-3, 3, 3,-4, 4,-2, 2,-2,-2,-1,-1/)
  A(12,:) = (/ 2, 0,-2, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0/)
  A(13,:) = (/ 0, 0, 0, 0, 2, 0,-2, 0, 0, 0, 0, 0, 1, 0, 1, 0/)
  A(14,:) = (/-6, 6, 6,-6,-4,-2, 4, 2,-3, 3,-3, 3,-2,-1,-2,-1/)
  A(15,:) = (/ 4,-4,-4, 4, 2, 2,-2,-2, 2,-2, 2,-2, 1, 1, 1, 1/)

  alpOld = matmul(A, xall)
  alp = reshape(alpOld, (/4, 4/))
  
  res = 0.0
  select case (derivs)
    case ('v')
      do i=0,3
          do j=0,3
              res = res + alp(i, j) * (x0**i) * (y0**j)
          end do
      end do

    case ('vr')
      do i=1,3
        do j=0,3
          res = res + alp(i,j)*(i*x0**(i-1))*(y0**j)
        end do
      end do

    case ('vz')
      do i=0,3
        do j=1,3
           res = res + alp(i,j)*(x0**i)*(j*y0**(j-1))
        end do
      end do
           
    case ('vrz')
      do i=1,3
        do j=1,3
          res = res + alp(i, j)*i*x0**(i-1)*j*y0**(j-1)
        end do
      end do
      
    case default
      print*, "error in Fortran"  
      stop

      
  end select
!  print *, "fortran res = ", res  
end subroutine bicubic


program test_bi3

 implicit none
 
 real :: f(4), fx(4), fy(4), fxy(4)
 real :: x0, y0, res
 character(len=3) :: derivs
 !!-constant field
 f=1.0
 fx=0.0
 fy=0.0
 fxy=0.0


 x0=0.5
 y0=0.5   
 
 derivs = 'v'
 
 call bicubic(f,fx,fy,fxy,x0,y0,derivs,res)
 print *, x0, res

end program
 