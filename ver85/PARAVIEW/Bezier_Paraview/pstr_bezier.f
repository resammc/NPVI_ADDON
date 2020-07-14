!$Id:$
      subroutine pstr_bezier
     &           (ie, ix, lknot, ixbez, st, ps, s_bez, p_bez,ip,
     &            nenb, numpbez)

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
!      Purpose: Project stresses, etc. on Bezier control points
!-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'cdata.h'
      include   'cdat1.h'
      include   'cnurb.h'
      include   'eldata.h'
      include   'eldatp.h'
      include   'iofile.h'
      include   'pdata3.h'
      include   'qudshp.h'
      include   'sdata.h'
      include   'strnum.h'
      include   'pointer.h'
      include   'comblk.h'

      !include   'cdata.h'
      !include   'cnurb.h'
      !include   'iofile.h'
      !include   'sdata.h'
      include   'setups.h'

      include   'npvi.h'

      integer          :: nenb, numpbez, type
      integer          :: ix(nen1,*), ixbez(nenb,*)
      integer          :: ie(nie,*), node
      integer          :: is, js, ks, nd, ne, nm, nn, ns, nod
      integer          :: i, ii, j, i_bez, d_bez
      integer          :: nll
      integer          :: p, q, r, pp, qq, rr
      integer          :: k1, k2, k3
      integer          :: lknot(0:4,*), nd_bez, ne_bez
      integer          :: ip(0:numnp)
      integer (kind=8) :: point1, point2, point3
      real    (kind=8) :: s_bez(numpbez,*),p_bez(numpbez,*)
      real    (kind=8) :: st(numnm,*), ps(numnm,8)
      real    (kind=8) :: n_stress(nen)
      real    (kind=8) :: b_stress(nen), b_temp(nen)
      real    (kind=8) :: c_e(nen,nen)

      save

!     Loop over NURBS elements

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

          endif ! eltyp

          ! calculate the full extraction operator for the current element
          nll = (p+1)*(q+1)*(r+1)

          c_e = 0.0d0
          call full_ce_npvi(hr(point1), hr(point2), hr(point3), p,
     &                      q, r, nll, type, c_e)

          ! stresses
          do i = 1, npstr-1
            n_stress = 0.0d0
            do j = 1,ixbez(nenb,ne)
              node = ix(j,ne)
              n_stress(j) = st(ip(node),i)
            end do

            ! calculate stresses on the Bezier control points
            b_stress = 0.0d0
            call pbez_npvi(b_stress,n_stress,
     &                     c_e, p, q, r, nll, 1, type)

            do j = 1,ixbez(nenb,ne)
              node = ixbez(j,ne)
              s_bez(node,i) = b_stress(j)
            end do
          end do

          ! principal stresses
          do i = 1, 8
            do j = 1,ixbez(nenb,ne)
              node = ix(j,n_el)
              n_stress(j) = ps(node,i)
            end do

            ! calculate stresses on the Bezier control points
            call pbez_npvi(b_stress,n_stress,
     &                     c_e, p, q, r, nll, 1, type)

            do j = 1,ixbez(nenb,ne)
              node = ixbez(j,ne)
              p_bez(node,i) = b_stress(j)
            end do
          end do

        end if
      end do


      end 

