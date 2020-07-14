      subroutine pdis_quad_nusd(ie,ix, x,xl,wt,wtl,u,u_quad,x_quad)

c      for F E A P -- A Finite Element Analysis Program
c
c      coded by:
c                B.Sc. Henning Venghaus
c                Institute of Mechanics
c                Otto von Guericke University
c                Spring 2017

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    02/03/2017
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Project T-spline/NURBS to linear elements for plots
c               works like pdis_lin.f but for quad. elements

c      Inputs:
c        ie(nie,*)    - Element type properties
c        ix(nen1,*)   - Mesh global node connections
c        x(ndm,*)     - Nodal coordinates
c        xl(ndm,*)    - Local element coordinate
c        wt(*)        - Nodal weights
c        wtl(*)       - Local weights
c        u(ndf,*)     - Nodal solution
c        x_quad(ndf,*) - quadratic coordinates

c      Outputs:
c        u_quad(ndf,*) - quadratic solution
c-----[--+---------+---------+---------+---------+---------+---------+-]
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
      real*8     x(ndm,*), xl(ndm,*), wt(*), wtl(*), x_quad(ndm,*),xx(3)
      real*8     u(ndf,*), u_quad(ndf,*)

      integer    ne, is,js,ks, nn, ns, i,j, na
      integer    ixq3d(3,20),ixq2d(2,8),ixq1d(1,3)

      logical    noskip
      integer    e_quad, d_quad, d_min
      real*8     du, dist, dist_min
      real*8     gp(3,20)

      save

      data       dist_min / 1.0d-12 /

      data ixq3d /0,0,0,
     &            2,0,0,
     &            2,2,0,
     &            0,2,0,
     &            0,0,2,
     &            2,0,2,
     &            2,2,2,
     &            0,2,2,
     &            1,0,0,
     &            2,1,0,
     &            1,2,0,
     &            0,1,0,
     &            1,0,2,
     &            2,1,2,
     &            1,2,2,
     &            0,1,2,
     &            0,0,1,
     &            2,0,1,
     &            2,2,1,
     &            0,2,1/
 
      data ixq2d /0,0,
     &            2,0,
     &            2,2,
     &            0,2,
     &            1,0,
     &            2,1,
     &            1,2,
     &            0,1/

      data ixq1d /0,
     &            2,
     &            1/

c     Get Gaussian points and weights

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
          gp(i,2*npl_int_nusd(i)+1) = 0.999999999999d0 ! Keep point inside element
        endif
      end do !i
      
      e_quad = 0
      d_quad = 0
      d_min = 1
      
      do na = 1,nummat
        do ne = 1,numel
          n  = ne   ! Element number for use in interpolation
          ma = ix(nen1,ne)
          if((ma.eq.na .and. maplt.eq.0) .or.
     &       (ma.eq.na .and. ma.eq.maplt)) then
            eltyp = ix(nen+7,ne)
            if(eltyp.gt.0) then

              elty2 = ix(nen+8,ne)
              elty3 = ix(nen+9,ne)

c             Set up local values

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

c             Do 1-d plots

              if(ie(1,ma).eq.1) then

                do is = 1,npl_int_nusd(1)

                  e_quad = e_quad + 1

c                 Construct 3-node 'quadratic elements'

                  do ns = 1,3

                    sg1(1,1) = gp(1,-1+2*is+ixq1d(1,ns))
                    sg1(2,1) = 1.0d0

                    call interp1d(1, xl, ndm,nel, .true.)

c                   Compute position

                    do i = 1,ndm
                      xx(i) = 0.0d0
                      do nn = 1,nel
                        xx(i) = xx(i) + shp1(2,nn,1)*xl(i,nn)
                      end do ! nn
                    end do ! i

c                   Check for repeated nodes

                    noskip = .true.
                    do i = d_quad,d_min,-1
                      dist = (x_quad(1,i) - xx(1))**2
                      if(ndm.ge.2) then
                        dist = dist + (x_quad(2,i) - xx(2))**2
                      endif
                      if(ndm.eq.3)then
                         dist = dist + (x_quad(3,i) - xx(3))**2
                      endif

                      if(dist.le.dist_min) then
                        noskip = .false.
                        exit
                      endif
                    end do ! i

c                   New node

                    if(noskip) then
                      d_quad = d_quad + 1

c                     Project dependent variables

                      do j = 1,ndf
                        u_quad(j,d_quad) = 0.0d0
                      end do ! j
                      do nn = 1,nel
                        if(ix(nn,ne).gt.0) then
                          do j = 1,ndf
                            u_quad(j,d_quad) = u_quad(j,d_quad)
     &                                 + shp1(2,nn,1)*u(i,ix(nn,ne))
                          end do ! j
                        endif
                      end do ! nn
                    endif

                  end do ! ns
                end do ! is

c             Do 2-d plots

              elseif(ie(1,ma).eq.2) then

c               For element loop over number of subdivisions 'npl_int'

                do js = 1,npl_int_nusd(2)
                  do is = 1,npl_int_nusd(1)

                    e_quad = e_quad + 1

c                   Construct 8-node 'quadratic elements'

                    do ns = 1,8

                      sg2(1,1) = gp(1,-1+2*is+ixq2d(1,ns))
                      sg2(2,1) = gp(2,-1+2*js+ixq2d(2,ns))
                      sg2(3,1) = 1.0d0
                      
                      call interp2d(1, xl, ix(1,ne), ndm,nel, .true.)

c                     Compute position

                      do i = 1,ndm
                        xx(i) = 0.0d0
                        do nn = 1,nel
                          xx(i) = xx(i) + shp2(3,nn,1)*xl(i,nn)
                        end do ! nn
                      end do ! i

c                     Check for repeated nodes

                      noskip = .true.
                      do i = d_quad,d_min,-1
                        dist = (x_quad(1,i) - xx(1))**2
     &                       + (x_quad(2,i) - xx(2))**2
                        if(ndm.eq.3) then
                          dist = dist + (x_quad(3,i) - xx(3))**2
                        endif

                        if(dist.le.dist_min) then
                          noskip = .false.
                          exit
                        endif
                      end do ! i

c                     New node

                      if(noskip) then
                        d_quad = d_quad + 1

c                       Project dependent variables

                        do j = 1,ndf
                          u_quad(j,d_quad) = 0.0d0
                        end do ! j
                        do nn = 1,nel
                          if(ix(nn,ne).gt.0) then
                            do j = 1,ndf
                              u_quad(j,d_quad) = u_quad(j,d_quad)
     &                                   + shp2(3,nn,1)*u(j,ix(nn,ne))
                            end do ! j
                          endif
                        end do ! nn
                      endif

                    end do ! ns
                  end do ! is
                end do ! js

c             Do 3-d plots

              elseif(ie(1,ma).eq.3) then
c           CANT USE THIS DUE TO CK_LIN8 !! 
c               For element loop over number of subdivisions 'npl_int'

                do ks = 1,npl_int_nusd(3)
                  do js = 1,npl_int_nusd(2)
                    do is = 1,npl_int_nusd(1)

                      e_quad = e_quad + 1

c                     Construct 20-node 'quadratic elements'

                      do ns = 1,20

                        sg3(1,1) = gp(1,-1+2*is+ixq3d(1,ns))
                        sg3(2,1) = gp(2,-1+2*js+ixq3d(2,ns))
                        sg3(3,1) = gp(3,-1+2*ks+ixq3d(3,ns))
                        sg3(4,1) = 1.0d0
!                         sg3(1,1) = gp(is+ixl(1,ns))
!                         sg3(2,1) = gp(js+ixl(2,ns))
!                         sg3(3,1) = gp(ks+ixl(3,ns))
!                         sg3(4,1) = 1.0d0

                        call interp3d(1, xl, ndm,nel)

c                       Compute position

                        do i = 1,ndm
                          xx(i) = 0.0d0
                          do nn = 1,nel
                            xx(i) = xx(i)+shp3(4,nn,1)*xl(i,nn)
                          end do ! nn
                        end do ! i

c                       Check for repeated nodes

                        noskip = .true.
                        do i = d_quad,d_min,-1
                          dist = (x_quad(1,i) - xx(1))**2
     &                         + (x_quad(2,i) - xx(2))**2
     &                         + (x_quad(3,i) - xx(3))**2

                          if(dist.le.dist_min) then
                            noskip = .false.
                            exit
                          endif
                        end do ! i

c                       New node

                        if(noskip) then
                          d_quad = d_quad + 1

c                         Project dependent variables

                          do j = 1,ndf
                            u_quad(j,d_quad) = 0.0d0
                          end do ! j
                          do nn = 1,nel
                            if(ix(nn,ne).gt.0) then
                              do j = 1,ndf
                                u_quad(j,d_quad) = u_quad(j,d_quad)
     &                                     + shp3(4,nn,1)*u(j,ix(nn,ne))
                              end do ! j
                            endif
                          end do ! nn

                        endif

                      end do ! ns

                    end do ! is
                  end do ! js
                end do ! ks

              endif ! ie(1,ma) test
            endif ! eltyp
          endif ! ma tests

        end do ! ne
        d_min = d_quad + 1
      end do ! na

      end