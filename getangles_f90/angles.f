      program main
         implicit none

      integer :: i, j, k, l, counter, counter2, counter3
      integer, parameter :: nstart = 1000
      integer, parameter :: nsteps = 20000
      integer, parameter :: natoms = 128
      integer, parameter :: natpermol = 4
      integer, parameter :: nmols = natoms/natpermol
      real, parameter :: bohr2ang = 0.52917721
      real, parameter :: gamma_par = 0.73612
      real, parameter :: charge = 1.1128
      real, parameter :: ebohr2D = 1.0/0.393430307
      real, parameter :: rad2deg = 3.14159265358979/180.0

      real :: cell(3,3)
      character*100 :: input_file
      character*4 :: element
      real :: x, y, z
      real :: r(3, natpermol)
      real :: rdiff(3), rdiff2(3)
      real :: rbar = 0
      real :: rbarold = 0
      real :: dipole = 0
      real :: dipoleold = 0
      real :: angle = 0
      real :: angleold = 0

      cell = 0.0
      do i = 1, 3
         cell(i,i) = 19.105
      end do

      !write(input_file, '(3A)') "test.pos_00.xyz"
      call getarg(1, input_file)

      open(unit=11, file=trim(input_file), action="READ")
      !open(unit=12, file="bonds.output", action="WRITE")
      open(unit=12, file="bonds.output", access="APPEND")
      open(unit=13, file="angles.output", access="APPEND")

      do i = 1, nstart
         read(11,*)
         read(11,*)
         do j = 1, natoms
            read(11,*)
         end do
      end do

      counter = 0
      counter2 = 0
      counter3 = 0
      do i = nstart, nsteps
         read(11,*)
         read(11,*)
         do j = 1, nmols
            do k = 1, natpermol
               read(11,'(A8, E13.5, E13.5, E13.5)') element, x, y, z
               r(1,k) = x
               r(2,k) = y
               r(3,k) = z
            end do

            do k = 2, natpermol
               !calculating bond length
               rbarold = rbar
               rbar = rbar + nearest_distance(cell, r(:,1), r(:,k))
               counter = counter + 1
               if (nearest_distance(cell, r(:,1), r(:,k)) > 4) then
                  write(*,*) i, j, k
                  write(*,*) nearest_distance(cell, r(:,1), r(:,k))
                  stop
               end if
               if (nearest_distance(cell, r(:,1), r(:,k)) < 1.0) then
                  write(*,*) i, j, k
                  write(*,*) nearest_distance(cell, r(:,1), r(:,k))
                  stop
               end if
               if (rbar-rbarold < 1.0) then
                  write(12,*) rbarold/counter*bohr2ang
                  rbar = 0
                  counter = 0
               endif
            end do
!            !calculating dipole moment
!            dipoleold = dipole
!            rdiff = 0
!            do k = 2, natpermol
!               rdiff  = rdiff + r(:,k)
!            end do
!            rdiff = r(:,1)*gamma_par - gamma_par*rdiff/(natpermol-1.0)
!            dipole = dipole + sqrt(dot_product(rdiff, rdiff))
!            counter2 = counter2 + 1
!            if (dipole-dipoleold < 0.01) then
!               write(12,*) "dipole=", dipoleold/counter2*charge*ebohr2D
!               dipole = 0
!               counter2 = 0
!            endif
            !calculating angle
            do k = 2,natpermol
               do l = k+1, natpermol
                  angleold = angle
                  rdiff = nearest_image(cell, r(:,1), r(:,k))
                  rdiff2 = nearest_image(cell, r(:,1), r(:,l))
                  angle = angle + acos(dot_product(rdiff,rdiff2)/
     1                    (sqrt(dot_product(rdiff,rdiff))*
     2                     sqrt(dot_product(rdiff2,rdiff2))))
                  counter3 = counter3 + 1
                  if (angle-angleold < 1.0) then
                     write(13,*) angleold/(counter3*rad2deg)
                     angle = 0
                     counter3 = 0
                  endif
               end do
           end do
         end do
         write(*,*) "Done step ", i-nstart, " of ", nsteps-nstart
         write(*,*) "rbar = ", rbar/counter*bohr2ang, 
     1                  " counter ", counter
         write(*,*) "angle = ", angle/counter3*rad2deg 
      end do

      write(*,*) "r(n-h) = ", rbar/counter*bohr2ang
      write(*,*) "theta(h-n-h) = ", angle/(counter3*rad2deg)
      write(12,*) rbar/counter*bohr2ang
!      write(*,*) "dipole = ", dipole/counter2*charge*ebohr2D
      write(13,*) angle/(counter3*rad2deg)

      contains
         real function nearest_distance(cell, q1, q2)
            real, dimension(3,3), intent(in) :: cell
            real, dimension(3), intent(in) :: q1
            real, dimension(3), intent(in) :: q2
            real, dimension(3) :: s

            integer i

            s = matmul(inv(cell), q1 - q2)
            do i = 1, 3
               s(i) = s(i) - nint(s(i))
            end do
            s = matmul(cell, s)
            nearest_distance = sqrt(dot_product(s,s))
         end function

         function nearest_image(cell, q1, q2)
            real, dimension(3,3), intent(in) :: cell
            real, dimension(3), intent(in) :: q1
            real, dimension(3), intent(in) :: q2
            real, dimension(3) :: nearest_image
            real, dimension(3) :: s

            integer i

            s = matmul(inv(cell), q1 - q2)
            do i = 1, 3
               s(i) = s(i) - nint(s(i))
            end do
            s = matmul(cell, s)
            nearest_image = s
         end function

         function inv(cell)
            real, dimension(3,3), intent(in) :: cell
            real, dimension(3,3) :: inv

            integer i

            do i = 1, 3
               inv(i,i) = 1.0/cell(i,i)
            end do
            inv(1,2) = -inv(1,1)*cell(1,2)*inv(2,2)
            inv(2,3) = -inv(2,2)*cell(2,3)*inv(3,3)
            inv(1,3) = -inv(1,2)*cell(1,2)*inv(1,1) - 
     1                              inv(1,1)*cell(1,3)*inv(3,3)
            inv(2,1) = 0.0
            inv(3,1) = 0.0
            inv(3,2) = 0.0
         end function

      end program
