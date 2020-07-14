!$Id:$
      subroutine pmsh_bezier(ie, ix, x, xl, wt, wtl, lknot,
     &                       x_bez, w_bez, ixbez, nmat, nenb, numpbez,
     &                       mergfl)

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
!      Purpose: Creates the Bezier mesh
!-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      include 'cdata.h'   ! numel, nen
      include 'eldata.h'  ! nel, n_el, ma, eltyp, elty2, elty3
      include 'sdata.h'   ! nen, ndf, nen1
      include 'comblk.h'  ! hr, mr
      include 'pointer.h' ! np, up
      include 'cdat1.h'   ! nie

      integer          :: i, j
      integer          :: ne, nll, nenb, numpbez, type
      integer          :: nd_bez, ne_bez
      integer          :: ie(nie,*), ix(nen1,*), lknot(0:4,*), nmat(*)
      integer          :: p, q, r, pp, qq, rr
      integer          :: k1, k2, k3, is, js, ks
      integer          :: ixbez(nenb,*)
      integer (kind=8) :: point1, point2, point3
      real    (kind=8) :: x(ndm,*), wt(*)
      real    (kind=8) :: xl(ndm,*), wtl(*)
      real    (kind=8) :: c_e(nen,nen)
      real    (kind=8) :: x_bez(ndm,*), w_bez(*)
      real    (kind=8) :: xbl(ndm,nen), wbl(nen)
      logical          :: mergfl


      do ne = 1,numel
        n_el = ne ! element number to use
        ma = ix(nen1,n_el)
        eltyp = ix(nen+7,ne)
        if (eltyp.gt.0) then
          elty2 = ix(nen+8,ne)
          elty3 = ix(nen+9,ne)

          ! Set up local values
          nel = 0
          do i = 1,nen
            if(ix(i,ne).gt.0) then
              do j = 1,ndm
                xl(j,i) = x(j,ix(i,ne))
              end do ! j
              wtl(i)  = wt(ix(i,ne))
              nel     = i
            else
              do j = 1,ndm
                xl(j,i) = 0.0d0
              end do ! j
              wtl(i)  = 1.0d0
            endif
          end do ! i

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

          else

            write(*,*) "Unknown element type. ELTYP = ", eltyp

          endif ! IE
        end if ! ELTYP

        nll = (p+1)*(q+1)*(r+1)

        ! calculate the full extraction operator for the current element
        c_e = 0.0d0
        call full_ce_npvi(hr(point1), hr(point2), hr(point3), p,
     &                    q, r, nll, type, c_e)

        xbl = 0.0d0
        wbl = 0.0d0
        call pbez_npvi(xbl, xl, c_e, p, q, r, nll, ndm, type)
        call pbez_npvi(wbl, wtl, c_e, p, q, r, nll, 1, type)

        ! rearrange the ix array and remove duplicate nodes
        nd_bez = numel*nen
        ne_bez = ne

        call psetbez_npvi(xbl, wbl, x_bez, w_bez, ixbez, nen, 
     &                    numpbez, nd_bez, nmat, ne_bez, ma, 
     &                    nenb, mergfl)

        ixbez(nenb,ne)  = nel
        ixbez(nen+2,ne) = p
        ixbez(nen+3,ne) = q
        ixbez(nen+4,ne) = r

      end do ! NE

      end
