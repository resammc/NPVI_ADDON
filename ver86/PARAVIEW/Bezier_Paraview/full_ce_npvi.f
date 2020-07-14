!$Id:$
      subroutine full_ce_npvi(c_e1T, c_e2T, c_e3T, p, q, r, nll, 
     &                        type, c_e)

!      * * F E A P * * A Finite Element Analysis Program

!....  Copyright (c) 1984-2018: Regents of the University of California
!                               All rights reserved

!-----[--.----+----.----+----.-----------------------------------------]

!      Written for FEAP 8.5 by Resam Makvandi (resam.makvandi@ovgu.de)
!                              Chair of Computational Mechanics
!                              Institute of Mechanics
!                              Otto von Guericke University Magdeburg
!                              Germany

!-----[--.----+----.----+----.-----------------------------------------]
!     Modification log                                Date (dd/mm/year)
!       Original version                                    30/05/2020
!-----[--.----+----.----+----.-----------------------------------------]
!      Purpose: Calculates the full 2D extraction operator
!-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      integer          :: i, j
      integer          :: p, q, r
      integer          :: nll, type
      integer          :: nn, mm, pp, qq
      real    (kind=8) :: c_e1T(p+1,p+1), c_e2T(q+1,q+1)
      real    (kind=8) :: c_e3T(r+1,r+1)
      real    (kind=8) :: c_e1(p+1,p+1), c_e2(q+1,q+1)
      real    (kind=8) :: c_e3(r+1,r+1)
      real    (kind=8) :: c_e(nll,nll)
      real    (kind=8) :: c_e_temp((p+1)*(q+1),(p+1)*(q+1))

      c_e = 0.0d0
      if (type.eq.1) then  ! univariate case

          do i = 1, p+1
            do j = 1, p+1
              c_e(i,j) = C_e1T(j,i)
            end do
          end do

      else if (type.eq.2) then  ! bivariate case

          ! calcualte the transpose
          do i = 1, p+1
            do j = 1, p+1
              c_e1(i,j) = C_e1T(j,i)
            end do
          end do
          do i = 1, q+1
            do j = 1, q+1
              c_e2(i,j) = C_e2T(j,i)
            end do
          end do
          
          ! Kronecker-product (c_e = c_e2 * c_e1)
          do i = 1,q+1
            do j = 1,q+1
             nn = (i-1)*(p+1) + 1
             mm = nn+ (p+1) - 1
             pp = (j-1)*(p+1) + 1
             qq = pp+ (p+1) - 1
             c_e(nn:mm,pp:qq) = c_e2(i,j)*c_e1
            enddo
          enddo

      else if (type.eq.3) then ! trivariate case

          !calculate the transpose
          do i = 1, p+1
            do j = 1, p+1
              c_e1(i,j) = C_e1T(j,i)
            end do
          end do
          do i = 1, q+1
            do j = 1, q+1
              c_e2(i,j) = C_e2T(j,i)
            end do
          end do
          do i = 1, r+1
            do j = 1, r+1
              c_e3(i,j) = C_e3T(j,i)
            end do
          end do

          ! Kronecker-product (c_e_temp = c_e2 * c_e1)
          do i = 1,q+1
            do j = 1,q+1
             nn = (i-1)*(p+1) + 1
             mm = nn+ (p+1) - 1
             pp = (j-1)*(p+1) + 1
             qq = pp+ (p+1) - 1
             c_e_temp(nn:mm,pp:qq) = c_e2(i,j)*c_e1
            enddo
          enddo

          ! Kronecker-product (c_e = c_e_temp * c_e3)
          do i = 1,(p+1)*(q+1)
            do j = 1,(p+1)*(q+1)
              nn = (i-1)*(r+1) + 1
              mm = nn + (r+1) - 1
              pp = (j-1)*(r+1) + 1
              qq = pp + (r+1) - 1
              c_e(nn:mm,pp:qq) = c_e_temp(i,j)*c_e3
            enddo
          enddo

      end if

      end
