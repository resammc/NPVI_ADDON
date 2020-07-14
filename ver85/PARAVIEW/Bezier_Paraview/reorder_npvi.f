!$Id:$
      subroutine reorder_npvi(array, size, p, q, r, type)
   
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
!      Purpose: Re-orders the Bezier control points for Paraview
!-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      integer          :: i, j, k, pc
      integer          :: p, q, r
      integer          :: type, size
      real    (kind=8) :: array(size), temp_array(size), conn(size)
        
      if (type.eq.1) then  ! univariate case

        ! end cps
        conn(1) = 1
        conn(2) = p+1
        
        ! inner cps
        pc = 2
        do i = 1, p-1
          pc = pc + 1
          conn(pc) = i + 1
        end do
          
      else if (type.eq.2) then  ! bivariate case

        ! vertices
        conn(1) = 1
        conn(2) = p+1
        conn(3) = (p+1)*(q+1)
        conn(4) = q*(p+1) + 1

        ! edges
        pc = 4
        ! edge 1
        do i = 1, p-1
          pc = pc + 1
          conn(pc) = i + 1
        end do
        ! edge 2
        do i = 1, q-1
          pc = pc + 1
          conn(pc) = (p+1)*(i+1)
        end do
        ! edge 3
        do i = 1, p-1
          pc = pc + 1
          conn(pc) = (p+1)*q + i + 1
        end do
        ! edge 4
        do i = 1, q-1
          pc = pc + 1
          conn(pc) = (p+1)*i + 1
        end do
        ! face
        do j = 1, q-1
          do i = 1, p-1
            pc = pc + 1
            conn(pc) = (p+1)*j + i + 1
          end do
        end do


      else if (type.eq.3) then ! trivariate case

        ! vertices
        conn(1) = 1
        conn(2) = p*(r+1) + 1
        conn(3) = q*(r+1)*(p+1) + p*(r+1) + 1
        conn(4) = q*(r+1)*(p+1) + 1
        conn(5) = r + 1
        conn(6) = (p+1)*(r+1)
        conn(7) = (p+1)*(q+1)*(r+1)
        conn(8) = (p+1)*q*(r+1) + r + 1

        ! edges
        pc = 8
        ! edge 1
        do i = 1, p-1
          pc = pc + 1
          conn(pc) = i*(r+1) + 1
        end do
        ! edge 2
        do i = 1, q-1
          pc = pc + 1
          conn(pc) = i*(p+1)*(r+1) + p*(r+1) + 1
        end do
        ! edge 3
        do i = 1, p-1
          pc = pc + 1
          conn(pc) = q*(r+1)*(p+1) + i*(r+1) + 1
        end do
        ! edge 4
        do i = 1, q-1
          pc = pc + 1
          conn(pc) = i*(p+1)*(r+1) + 1
        end do
        ! edge 5
        do i = 1, p-1
          pc = pc + 1
          conn(pc) = (i + 1)*(r + 1)
        end do
        ! edge 6
        do i = 1, q-1
          pc = pc + 1
          conn(pc) = (p + 1)*(r + 1)*(i + 1)
        end do
        ! edge 7
        do i = 1, p-1
          pc = pc + 1
          conn(pc) = (p+1)*(r+1)*q + (i + 1)*(r + 1)
        end do
        ! edge 8
        do i = 1, q-1
          pc = pc + 1
          conn(pc) = (p+1)*(r+1)*i + (r+1)
        end do
        ! edge 9
        do i = 1, r-1
          pc = pc + 1
          conn(pc) = i + 1
        end do
        ! edge 10
        do i = 1, r-1
          pc = pc + 1
          conn(pc) = p*(r+1) + i + 1
        end do
        ! edge 11
        do i = 1, r-1
          pc = pc + 1
          conn(pc) = (p+1)*(r+1)*q + (r+1)*p + i + 1
        end do
        ! edge 12
        do i = 1, r-1
          pc = pc + 1
          conn(pc) = (p+1)*(r+1)*q + i + 1
        end do
        ! face 1 (s-t-1)
        do j = 1, r-1
          do i = 1, q-1
            pc = pc + 1
            conn(pc) = (p+1)*(r+1)*i + j + 1
          end do ! i
        end do ! j
        ! face 2 (s-t-2)
        do j = 1, r-1
          do i = 1, q-1
            pc = pc + 1
            conn(pc) = (p+1)*(r+1)*i + p*(r+1) + j + 1
          end do ! i
        end do ! j
        ! face 3 (r-t-1)
        do j = 1, r-1
          do i = 1, p-1
            pc = pc + 1
            conn(pc) = i*(r+1) + j + 1
          end do ! i
        end do ! j
        ! face 4 (r-t-2)
        do j = 1, r-1
          do i = 1, p-1
            pc = pc + 1
            conn(pc) = (p+1)*(r+1)*q + i*(r+1) + j + 1
          end do ! i
        end do ! j
        ! face 5 (r-s-1)
        do j = 1, q-1
          do i = 1, p-1
            pc = pc + 1
            conn(pc) = j*(p+1)*(r+1) + i*(r+1) + 1
          end do ! i
        end do ! j
        ! face 6 (r-s-2)
        do j = 1, q-1
          do i = 1, p-1
            pc = pc + 1
            conn(pc) = j*(p+1)*(r+1) + (i+1)*(r+1)
          end do ! i
        end do ! j
        ! inside
        do k = 1, r-1
          do j = 1, q-1
            do i = 1, p-1
              pc = pc + 1
              conn(pc) = j*(p+1)*(r+1) + i*(r+1) + k + 1 
            end do ! i
          end do ! j
        end do !k

      end if

      do i = 1, size
        temp_array(i) = array(conn(i))
      end do

      array = temp_array
      
      end
