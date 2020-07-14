!$Id:$
      subroutine pdis_quad_nusd_new
     &           (ix, x,xl, wt,wtl, u, u_quad, ix_quad)

!      * * F E A P * * A Finite Element Analysis Program

!....  Copyright (c) 1984-2017: Regents of the University of California
!                               All rights reserved

!-----[--.----+----.----+----.-----------------------------------------]
!     Modification log                                Date (dd/mm/year)
!       Original version                                    01/11/2006
!       1. Add line plot capability                         17/01/2011
!       2. Add loop over material numbers                   31/03/2012
!       3. Add d_min to separate material nodes             28/09/2012
!       4. Remove loop over material sets and tests on ma   18/12/2014
!       5. Correct index on u(i ... to u(j for 1-d plots    25/03/2016
!       6. Set dist_min in pmsh_lin; add to qudshp.h        23/04/2016
!       7. Recode to use data stored in ix_lin(19,*) array  17/07/2016
!       8. Change 'ne' to 'nn' for set of is in 1-d plots   01/12/2016`
!       9. Change 'n' to 'n_el'                             09/05/2017
!      10. Added calls for npvi and nusd functionality      03/03/2018
!-----[--.----+----.----+----.-----------------------------------------]
!      Purpose: Project solution parameters of linear elements for plots

!      Inputs:
!        ix(nen1,*)   - Mesh global node connections
!        x(ndm,*)     - Nodal coordinates
!        xl(ndm,*)    - Local element coordinate
!        wt(*)        - Nodal weights
!        wtl(*)       - Local weights
!        u(ndf,*)     - Nodal solution
!        ix_quad(31,*)- Quadratic element connection

!      Outputs:
!        u_quad(ndf,*) - Quadrilateral solution
!-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none

      include   'cdata.h'
      include   'cdat1.h'
      include   'cnurb.h'
      include   'eldata.h'
      include   'iofile.h'
      include   'qudshp.h'
      include   'sdata.h'

      include   'npvi.h'

      integer  :: ix(nen1,*), ix_quad(31,*)
      real    (kind=8) :: x(ndm,*), xl(ndm,*), wt(*), wtl(*)
      real    (kind=8) :: u(ndf,*), u_quad(ndf,*)

!     Local variables

      integer  :: is,js,ks, nd, ne, nm, nn, ns, i, j, d_quad
      integer  :: iqx3d(20), iqy3d(20), iqz3d(20), iqx2d(8), iqy2d(8)
      integer  :: iqx1d(3)

      real    (kind=8) :: du
      real    (kind=8) :: gp(3,20)

      save

      data  iqx3d  / 0, 2, 2, 0, 0, 2, 2, 0, 1, 2, 1, 0, 1, 2, 1, 0, 0,
     &               2, 2, 0 /
      data  iqy3d  / 0, 0, 2, 2, 0, 0, 2, 2, 0, 1, 2, 1, 0, 1, 2, 1, 0,
     &               0, 2, 2 /
      data  iqz3d  / 0, 0, 0, 0, 2, 2, 2, 2, 0, 0, 0, 0, 2, 2, 2, 2, 1,
     &               1, 1, 1 /
     
      data  iqx2d   / 0, 2, 2, 0, 1, 2, 1, 0 /
      data  iqy2d   / 0, 0, 2, 2, 0, 1, 2, 1 /
      
      data  iqx1d   / 0, 2, 1 /

!     Get Gaussian points and weights

      lint = 1

      do i = 1,ndm
        du  = 1.0d0/dble(npl_int_nusd(i)) ! half the step length
        gp(i,1)  = -1.000d0
        do j = 1,(2*npl_int_nusd(i)-1)    ! double the number of steps
          gp(i,j+1) = gp(i,j) + du
        end do ! j
        if(tsplfl .or. hsplfl) then
          gp(i,2*npl_int_nusd(i)+1) = 1.d0
        else
          gp(i,2*npl_int_nusd(i)+1) = 0.999999999999d0 ! Keep in elmt
        endif
      end do !i

!     Loop over quadratic element

      do nn = 1,ne_quad

        ne = ix_quad(22,nn)
        nd = ix_quad(23,nn)

        n_el  = ne   ! Element number for use in interpolation
        eltyp = ix(nen+7,ne)
        if(eltyp.gt.0) then

          elty2 = ix(nen+8,ne)
          elty3 = ix(nen+9,ne)

!         Set up local values

          nel = 0
          do i = 1,nen
            if(ix(i,ne).gt.0) then
              xl(1:ndm,i) = x(1:ndm,ix(i,ne))
              wtl(i)      = wt(ix(i,ne))
              nel         = i
            else
              xl(1:ndm,i) = 0.0d0
              wtl(i)      = 1.0d0
            endif
          end do ! i

!         Do 1-d plots

          if(nd.eq.1) then

            is = ix_quad(24,nn)

!           Construct 3-node 'quadratic elements'

            do ns = 1,3

              sg1(1,1) = gp(1,-1+2*is+iqx1d(ns))
              sg1(2,1) = 1.0d0

              call interp1d(1, xl, ndm,nel, .true.)

!             Get local node number

              d_quad = ix_quad(ns,nn)

!             Project dependent variables

              u_quad(1:ndf,d_quad) = 0.0d0
              do nm = 1,nel
                if(ix(nm,ne).gt.0) then
                  u_quad(1:ndf,d_quad) = u_quad(1:ndf,d_quad)
     &                               + shp1(2,nm,1)*u(1:ndf,ix(nm,ne))
                endif
              end do ! nm

            end do ! ns

!         Do 2-d plots

          elseif(nd.eq.2) then

!           For element loop over number of subdivisions 'npl_int'

            is = ix_quad(24,nn)
            js = ix_quad(25,nn)

!           Construct 8-node 'quadratic elements'

            do ns = 1,8

              sg2(1,1) = gp(1,-1+2*is+iqx2d(ns))
              sg2(2,1) = gp(2,-1+2*js+iqy2d(ns))
              sg2(3,1) = 1.0d0

              call interp2d(1, xl, ix(1,ne), ndm,nel, .true.)

!             Set node projection number

              d_quad = ix_quad(ns,nn)

!             Project dependent variables

              u_quad(1:ndf,d_quad) = 0.0d0
              do nm = 1,nel
                if(ix(nm,ne).gt.0) then
                  u_quad(1:ndf,d_quad) = u_quad(1:ndf,d_quad)
     &                               + shp2(3,nm,1)*u(1:ndf,ix(nm,ne))
                endif
              end do ! nm

            end do ! ns

!         Do 3-d plots

          elseif(nd.eq.3) then

!           For element loop over number of subdivisions 'npl_int'

            is = ix_quad(24,nn)
            js = ix_quad(25,nn)
            ks = ix_quad(26,nn)

!           Construct 20-node 'quadratic elements'

            do ns = 1,20

              sg3(1,1) = gp(1,-1+2*is+iqx3d(ns))
              sg3(2,1) = gp(2,-1+2*js+iqy3d(ns))
              sg3(3,1) = gp(3,-1+2*ks+iqz3d(ns))
              sg3(4,1) = 1.0d0

              call interp3d(1, xl, ndm,nel)

              d_quad = ix_quad(ns,nn)

!             Project dependent variables

              u_quad(1:ndf,d_quad) = 0.0d0
              do nm = 1,nel
                if(ix(nm,ne).gt.0) then
                  u_quad(1:ndf,d_quad) = u_quad(1:ndf,d_quad)
     &                               + shp3(4,nm,1)*u(1:ndf,ix(nm,ne))
                endif
              end do ! nn

            end do ! ns

          endif
        endif ! eltyp

      end do ! ne

      end
