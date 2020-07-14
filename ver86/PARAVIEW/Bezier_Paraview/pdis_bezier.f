!$Id:$
      subroutine pdis_bezier(ie, ix, lknot, u, v, a, ixbez, u_bez, 
     &                       v_bez, a_bez, velofl, accefl, nenb)

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
!      Purpose: Calculates displacements, velocities, and accelerations
!               on the Bezier mesh
!-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      include 'pointer.h' ! np, up
      include 'comblk.h'  ! hr, mr
      include 'cdata.h'   ! numel, nen
      include 'eldata.h'  ! nel, n_el, ma, eltyp, elty2, elty3
      include 'sdata.h'   ! nen, ndf, nen1
      include 'cdat1.h'   ! nie

      integer          :: i, j
      integer          :: ne, nll, nenb, type
      integer          :: nd_bez, ne_bez
      integer          :: ie(nie,*), ix(nen1,*), lknot(0:4,*)
      integer          :: p, q, r, pp, qq, rr
      integer          :: k1, k2, k3, is, js, ks
      integer          :: ixbez(nenb,*)
      integer (kind=8) :: point1, point2, point3
      real    (kind=8) :: u(ndf,*), v(ndf,*), a(ndf,*)
      real    (kind=8) :: u_bez(ndf,*), v_bez(ndf,*), a_bez(ndf,*)
      real    (kind=8) :: ul(ndf,nen), vl(ndf,nen), al(ndf,nen)
      real    (kind=8) :: ubl(ndf,nen), vbl(ndf,nen), abl(ndf,nen)
      real    (kind=8) :: c_e(nen,nen)
      logical          :: velofl, accefl

      do ne = 1,numel
        n_el = ne ! element number to use
        ma = ix(nen1,n_el)
        eltyp = ix(nen+7,ne)
        if (eltyp.gt.0) then
          elty2 = ix(nen+8,ne)
          elty3 = ix(nen+9,ne)

          type = ie(1,ma)

          if (type.eq.1) then

            k1 = mod(eltyp,500)
      
            p  = lknot(2,k1)
            q  = 0
            r  = 0
      
            pp = p + 1

            is = eltyp/500

            point1 = np(289) + mr(np(273)+k1-1) + pp*pp*is
            point2 = 0
            point3 = 0

          elseif (type.eq.2) then

            k1 = mod(eltyp,500)
            k2 = mod(elty2,500)
        
            p  = lknot(2,k1)
            q  = lknot(2,k2)
            r  = 0
        
            pp = p + 1
            qq = q + 1

            is = eltyp/500
            js = elty2/500

            point1 = np(289) + mr(np(273)+k1-1) + pp*pp*is
            point2 = np(289) + mr(np(273)+k2-1) + qq*qq*js
            point3 = 0

          elseif (type.eq.3) then

            k1 = mod(eltyp,500)
            k2 = mod(elty2,500)
            k3 = mod(elty3,500)
          
            p  = lknot(2,k1)
            q  = lknot(2,k2)
            r  = lknot(2,k3)
          
            pp = p + 1
            qq = q + 1
            rr = r + 1
  
            is = eltyp/500
            js = elty2/500
            ks = elty3/500
  
            point1 = np(289) + mr(np(273)+k1-1) + pp*pp*is
            point2 = np(289) + mr(np(273)+k2-1) + qq*qq*js
            point3 = np(289) + mr(np(273)+k3-1) + rr*rr*ks
              
          endif ! IE

          ! calculate the full extraction operator for the current element
          nll = (p+1)*(q+1)*(r+1)

          c_e = 0.0d0
          call full_ce_npvi(hr(point1), hr(point2), hr(point3), p,
     &                      q, r, nll, type, c_e)

          ! store the displacements, ... on the original mesh
          do i = 1,ndf
            do j = 1,nen
              ul(i,j) = u(i,ix(j,n_el))
              if (velofl) vl(i,j) = v(i,ix(j,n_el))
              if (accefl) al(i,j) = a(i,ix(j,n_el))
            end do
          end do
          
          ubl = 0.0d0
          call pbez_npvi(ubl, ul, c_e, p, q, r, nll, ndf, type)

          if (velofl) then
            vbl = 0.0d0
            call pbez_npvi(vbl, vl, c_e, p, q, r, nll, ndf, type)
          end if
      
          if (accefl) then
            abl = 0.0d0
            call pbez_npvi(abl, al, c_e, p, q, r, nll, ndf, type)
          end if

          ! store the displacements, ... on the Bezier mesh
          do i = 1,ndf
            do j = 1,ixbez(nenb,ne)
              u_bez(i,ixbez(j,n_el)) = ubl(i,j)
              if (velofl) v_bez(i,ixbez(j,n_el)) = vbl(i,j)
              if (accefl) a_bez(i,ixbez(j,n_el)) = abl(i,j)
            end do
          end do

        end if ! ELTYP
      end do ! NE

      end