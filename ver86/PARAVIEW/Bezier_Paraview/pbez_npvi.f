!$Id:$
      subroutine pbez_npvi(array, array_local, c_e, p, q, r, nll,
     &                     size, type)

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
!      Purpose: Projects the array on the Bezier mesh
!-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      integer          :: i, j, k
      integer          :: p, q, r, nll
      integer          :: type, size
      real    (kind=8) :: array(size, *), array_local(size,*)
      real    (kind=8) :: c_e(nll,nll)

      do i = 1,nll
        do j = 1,size
            array(j,i) = 0.0d0
        end do ! j
      end do ! i
      do i = 1,nll
        do j = 1,nll
          do k = 1,size
            array(k,j) = array(k,j) + c_e(i,j)*array_local(k,i)
          end do
        end do
      end do

      ! reorder
      do i = 1,size
        call reorder_npvi(array(i,1:nll), nll, p, q, r, type)
      end do

      end