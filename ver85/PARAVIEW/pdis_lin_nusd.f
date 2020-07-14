!$Id:$
      subroutine pdis_lin_nusd(ie,ix, x,xl, wt,wtl, u, u_lin, x_lin)

!      for F E A P -- A Finite Element Analysis Program
!
!      coded by:
!                B.Sc. Henning Venghaus
!                Institute of Mechanics
!                Otto von Guericke University
!                Spring 2017

!-----[--.----+----.----+----.-----------------------------------------]
!     Modification log                                Date (dd/mm/year)
!       Original version                                    26/03/2017
!       1. Update to ver8.5                                 03/07/2017
!          Change 'dist' to 'dst' (in qudshp.h); n to n_el
!-----[--.----+----.----+----.-----------------------------------------]
!      Purpose: Project T-spline/NURBS to linear elements for plots with
!               the opportunity of anisotropic element subdivision

!      Inputs:
!        ie(nie,*)    - Element type properties
!        ix(nen1,*)   - Mesh global node connections
!        x(ndm,*)     - Nodal coordinates
!        xl(ndm,*)    - Local element coordinate
!        wt(*)        - Nodal weights
!        wtl(*)       - Local weights
!        u(ndf,*)     - Nodal solution
!        x_lin(ndf,*) - Quadrilateral coordinates

!      Outputs:
!        u_lin(ndf,*) - Quadrilateral solution
!-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none

      include   'cdata.h'
      include   'cdat1.h'
      include   'cnurb.h'
      include   'eldata.h'
      include   'iofile.h'
      include   'pbody.h'
      include   'qudshp.h'
      include   'sdata.h'

      include   'npvi.h'

      integer    ie(nie,*),ix(nen1,*)
      real*8     x(ndm,*), xl(ndm,*), wt(*), wtl(*), x_lin(ndm,*), xx(3)
      real*8     u(ndf,*), u_lin(ndf,*)

      integer    ne, is,js,ks, nn, ns, i,j, na
      integer    ixl(3,8)

      logical    noskip
      integer    e_lin, d_lin, d_min
      real*8     du, dst, dst_min
      real*8     gp(3,9)

      save


      data       ixl / 0,0,0,   1,0,0,  1,1,0,   0,1,0,
     &                 0,0,1,   1,0,1,  1,1,1,   0,1,1 /


!       Set dst_min for mesh searches

        dst_min = 1.0d-03*max(hsize(1),hsize(2)*1.d-3)**2
        if(dst_min.eq.0.0d0) then
          dst_min = 1.d-12
        endif

!     Get Gaussian points and weights

      lint = 1

      do i = 1,ndm
        du  = 2.0d0/dble(npl_int_nusd(i))
        gp(i,1)  = -1.000d0
        do j = 1,npl_int_nusd(i)-1
          gp(i,j+1) = gp(i,j) + du
        end do ! j
        if(tsplfl .or. hsplfl) then
          gp(i,npl_int_nusd(i)+1) = 1.d0
        else
          gp(i,npl_int_nusd(i)+1) = 0.999999999999d0 ! Keep inside elmt
        endif
      end do !i

      e_lin = 0
      d_lin = 0
      d_min = 1
      do na = 1,nummat
        do ne = 1,numel
          n_el  = ne   ! Element number for use in interpolation
          ma = ix(nen1,ne)
          if((ma.eq.na .and. maplt.eq.0) .or.
     &       (ma.eq.na .and. ma.eq.maplt)) then
            eltyp = ix(nen+7,ne)
            if(eltyp.gt.0) then

              elty2 = ix(nen+8,ne)
              elty3 = ix(nen+9,ne)

!             Set up local values

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

!             Do 1-d plots

              if(ie(1,ma).eq.1) then

                do is = 1,npl_int_nusd(1)

                  e_lin = e_lin + 1

!                 Construct 2-node 'linear elements'

                  do ns = 1,2

                    sg1(1,1) = gp(1,is+ixl(1,ns))
                    sg1(2,1) = 1.0d0

                    call interp1d(1, xl, ndm,nel, .true.)

!                   Compute position

                    do i = 1,ndm
                      xx(i) = 0.0d0
                      do nn = 1,nel
                        xx(i) = xx(i) + shp1(2,nn,1)*xl(i,nn)
                      end do ! nn
                    end do ! i

!                   Check for repeated nodes

                    noskip = .true.
                    do i = d_lin,d_min,-1
                      dst = (x_lin(1,i) - xx(1))**2
                      if(ndm.ge.2) then
                        dst = dst + (x_lin(2,i) - xx(2))**2
                      endif
                      if(ndm.eq.3)then
                         dst = dst + (x_lin(3,i) - xx(3))**2
                      endif

                      if(dst.le.dst_min) then
                        noskip = .false.
                        exit
                      endif
                    end do ! i

!                   New node

                    if(noskip) then
                      d_lin = d_lin + 1

!                     Project dependent variables

                      do j = 1,ndf
                        u_lin(j,d_lin) = 0.0d0
                      end do ! j
                      do nn = 1,nel
                        if(ix(nn,ne).gt.0) then
                          do j = 1,ndf
                            u_lin(j,d_lin) = u_lin(j,d_lin)
     &                                 + shp1(2,nn,1)*u(i,ix(nn,ne))
                          end do ! j
                        endif
                      end do ! nn
                    endif

                  end do ! ns
                end do ! is

!             Do 2-d plots

              elseif(ie(1,ma).eq.2) then

!               For element loop over number of subdivisions 'npl_int'

                do js = 1,npl_int_nusd(2)
                  do is = 1,npl_int_nusd(1)

                    e_lin = e_lin + 1

!                   Construct 4-node 'linear elements'

                    do ns = 1,4

                      sg2(1,1) = gp(1,is+ixl(1,ns))
                      sg2(2,1) = gp(2,js+ixl(2,ns))
                      sg2(3,1) = 1.0d0

                      call interp2d(1, xl, ix(1,ne), ndm,nel, .true.)

!                     Compute position

                      do i = 1,ndm
                        xx(i) = 0.0d0
                        do nn = 1,nel
                          xx(i) = xx(i) + shp2(3,nn,1)*xl(i,nn)
                        end do ! nn
                      end do ! i

!                     Check for repeated nodes

                      noskip = .true.
                      do i = d_lin,d_min,-1
                        dst = (x_lin(1,i) - xx(1))**2
     &                       + (x_lin(2,i) - xx(2))**2
                        if(ndm.eq.3) then
                          dst = dst + (x_lin(3,i) - xx(3))**2
                        endif

                        if(dst.le.dst_min) then
                          noskip = .false.
                          exit
                        endif
                      end do ! i

!                     New node

                      if(noskip) then
                        d_lin = d_lin + 1
			write(*,*) "d_lin", d_lin
!                       Project dependent variables

                        do j = 1,ndf
                          u_lin(j,d_lin) = 0.0d0
                        end do ! j
                        do nn = 1,nel
                          if(ix(nn,ne).gt.0) then
                            do j = 1,ndf
                              u_lin(j,d_lin) = u_lin(j,d_lin)
     &                                   + shp2(3,nn,1)*u(j,ix(nn,ne))
                            end do ! j
                          endif
                        end do ! nn
                      endif

                    end do ! ns
                  end do ! is
                end do ! js

!             Do 3-d plots

              elseif(ie(1,ma).eq.3) then

!               For element loop over number of subdivisions 'npl_int'

                do ks = 1,npl_int_nusd(3)
                  do js = 1,npl_int_nusd(2)
                    do is = 1,npl_int_nusd(1)

                      e_lin = e_lin + 1

!                     Construct 8-node 'linear elements'

                      do ns = 1,8

                        sg3(1,1) = gp(1,is+ixl(1,ns))
                        sg3(2,1) = gp(2,js+ixl(2,ns))
                        sg3(3,1) = gp(3,ks+ixl(3,ns))
                        sg3(4,1) = 1.0d0

                        call interp3d(1, xl, ndm,nel)

!                       Compute position

                        do i = 1,ndm
                          xx(i) = 0.0d0
                          do nn = 1,nel
                            xx(i) = xx(i)+shp3(4,nn,1)*xl(i,nn)
                          end do ! nn
                        end do ! i

!                       Check for repeated nodes

                        noskip = .true.
                        do i = d_lin,d_min,-1
                          dst = (x_lin(1,i) - xx(1))**2
     &                         + (x_lin(2,i) - xx(2))**2
     &                         + (x_lin(3,i) - xx(3))**2

                          if(dst.le.dst_min) then
                            noskip = .false.
                            exit
                          endif
                        end do ! i

!                       New node

                        if(noskip) then
                          d_lin = d_lin + 1

!                         Project dependent variables

                          do j = 1,ndf
                            u_lin(j,d_lin) = 0.0d0
                          end do ! j
                          do nn = 1,nel
                            if(ix(nn,ne).gt.0) then
                              do j = 1,ndf
                                u_lin(j,d_lin) = u_lin(j,d_lin)
     &                                     + shp3(4,nn,1)*u(j,ix(nn,ne))
                              end do ! j
                            endif
                          end do ! nn

                        endif

                      end do ! ns

                    end do ! is
                  end do ! js
                end do ! ks

              endif
            endif ! eltyp
          endif ! ma tests

        end do ! ne
        d_min = d_lin + 1
      end do ! na

      end subroutine pdis_lin_nusd
