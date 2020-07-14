!$Id:$
      subroutine pmsh_quad_lgrn_nusd(ie,ix,x,xl,wt,wtl,ix_quad,x_quad,
     &                          nd_l,ne_l)

!      for F E A P -- A Finite Element Analysis Program
!
!      coded by:
!                B.Sc. Henning Venghaus
!                Institute of Mechanics
!                Otto von Guericke University
!                Spring 2017

!-----[--.----+----.----+----.-----------------------------------------]
!     Modification log                                Date (dd/mm/year)
!       Original version                                    02/03/2017
!       1. Update to ver8.5                                 03/07/2017
!          Change xl8(3,8) to xl8(3,27); dist to dst;
!          n to n_el
!       2. Fixed a bug in nodal coordiantes 21-27           23/06/2020
!       3. Commented out the call to ck_lin8 for quad elems 23/06/2020
!-----[--.----+----.----+----.-----------------------------------------]
!      Purpose: Prepare mesh of quadratic elements for 3-d plots
!               works like pmsh_lin but for quad. elements

!      Inputs:
!        ie
!        ix
!        x
!        xl
!        wt
!        wtl
!        ix_quad
!        x_quad
!        nd_l
!        ne_l

!      Outputs:
!-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none

      include   'cdata.h'
      include   'iofile.h'
      include   'cdat1.h'
      include   'cnurb.h'
      include   'eldata.h'
      include   'pbody.h'
      include   'qudshp.h'
      include   'sdata.h'
      include   'pointer.h'

      include   'npvi.h'

      integer    ie(nie,*),ix(nen1,*), ix_quad(38,*)
      real*8     x(ndm,*), xl(ndm,*), x_quad(ndm,*), wt(*), wtl(*)

      logical    noskip
      integer    nd_l,ne_l, nd_min, ne, is,js,ks, nn, ns, i,j, na
      integer    ixq3d(3,27),ixq2d(2,9),ixq1d(1,3)
      real*8     dst, dst_min, du
      real*8     gp(3,20), xx(3), xl8(3,27)

      save

      data       dst_min / 1.0d-12 /

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
     &            0,2,1,
     &            1,1,0,
     &            1,1,2,
     &            0,1,1,
     &            2,1,1,
     &            1,0,1,
     &            1,2,1,
     &            1,1,1/

      data ixq2d /0,0,
     &            2,0,
     &            2,2,
     &            0,2,
     &            1,0,
     &            2,1,
     &            1,2,
     &            0,1,
     &            1,1/

      data ixq1d /0,
     &            2,
     &            1/

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

!       write(*,*) gp

      nd_quad = 0
      nd_min = 1
      ne_quad = 0
      do na = 1,nummat
        do ne = 1,numel
          n_el  = ne   ! Element number for use
          ma = ix(nen1,ne)
          if((ma.eq.na .and. maplt.eq.0 ) .or.
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

!             For element loop over number of subdivisions 'npl_int'

!             Line mesh

              if(ie(1,ma).eq.1) then

                do is = 1,npl_int_nusd(1)

                  ne_quad = ne_quad + 1

                  if(ne_quad.gt.ne_l) then
                    write(*,*) ' NE_QUAD TOO LARGE',n_el,ne_quad,ne_l
                    call plstop(.true.)
                  endif

                  ix_quad(38,ne_quad) = ma

!                 Construct 3-node 'quadratic elements'

                  do ns = 1,3

                    sg1(1,1) = gp(1,-1+2*is+ixq1d(1,ns))
                    sg1(2,1) = 1.0d0

                    call interp1d(1, xl, ndm,nel, .false.)

                    do j = 1,ndm
                      xx(j) = 0.0d0
                      do nn = 1,nel
                        xx(j) = xx(j) + shp1(2,nn,1)*xl(j,nn)
                      end do ! nn
                    end do ! j

!                   Check for repeated nodes

                    noskip = .true.
                    do i = nd_quad,nd_min,-1
                      dst = (x_quad(1,i) - xx(1))**2
                      if(ndm.ge.2) dst = dst + (x_quad(2,i)-xx(2))**2
                      if(ndm.eq.3) dst = dst + (x_quad(3,i)-xx(3))**2

                      if(dst.le.dst_min) then
                        ix_quad(ns,ne_quad) = i
                        noskip = .false.
                        exit
                      endif
                    end do ! i

!                   New node

                    if(noskip) then
                      nd_quad            = nd_quad + 1
                      if(nd_quad.gt.nd_l) then
                        write(*,*) 'ND_QUAD TOO LARGE',
     &                              n_el,ns,nd_quad,nd_l
                        call plstop(.true.)
                      endif
                      ix_quad(ns,ne_quad) = nd_quad
                      ix_quad(34,ne_quad) = (-1)
                      do j = 1,ndm
                        x_quad(j,nd_quad)   = xx(j)
                      end do ! j
                    endif

                  end do ! ns
                end do ! is

!             Surface mesh

              elseif(ie(1,ma).eq.2) then
                do js = 1,npl_int_nusd(2)
                  do is = 1,npl_int_nusd(1)

                    ne_quad = ne_quad + 1

                    if(ne_quad.gt.ne_l) then
                      write(*,*) 'NE_QUAD TOO LARGE',n_el,ne_quad,ne_l
                      call plstop(.true.)
                    endif

                    ix_quad(38,ne_quad) = ma

!                   Construct 8-node 'quadratic elements'

                    do ns = 1,9

                      sg2(1,1) = gp(1,-1+2*is+ixq2d(1,ns))
                      sg2(2,1) = gp(2,-1+2*js+ixq2d(2,ns))
                      sg2(3,1) = 1.0d0

                      call interp2d(1, xl, ix(1,ne), ndm,nel, .false.)

                      do j = 1,ndm
                        xx(j) = 0.0d0
                        do nn = 1,nel
                          xx(j) = xx(j) + shp2(3,nn,1)*xl(j,nn)
                        end do ! nn
                      end do ! j


!                     Check for repeated nodes

                      noskip = .true.
                      do i = nd_quad,nd_min,-1
                        dst = (x_quad(1,i) - xx(1))**2
     &                       + (x_quad(2,i) - xx(2))**2
                        if(ndm.eq.3) then
                          dst = dst + (x_quad(3,i)-xx(3))**2
                        endif

                        if(dst.le.dst_min) then
                          ix_quad(ns,ne_quad) = i
                          noskip = .false.
                          exit
                        endif
                      end do ! i

!                     New node
!                       write(*,*) xx, noskip
                      if(noskip) then
                        nd_quad            = nd_quad + 1
                        if(nd_quad.gt.nd_l) then
                          write(*,*) ' ND_QUAD TOO LARGE',
     &                                 n_el,ns,nd_quad,nd_l
                          call plstop(.true.)
                        endif
                        ix_quad(ns,ne_quad) = nd_quad
                        ix_quad(34,ne_quad) = (-3)
                        do j = 1,ndm
                          x_quad(j,nd_quad)   = xx(j)
                        end do ! j
                      endif

                    end do ! ns
                  end do ! is
                end do ! js

!             Solid mesh

              elseif(ie(1,ma).eq.3) then

                do ks = 1,npl_int_nusd(3)
                  do js = 1,npl_int_nusd(2)
                    do is = 1,npl_int_nusd(1)

                      ne_quad = ne_quad + 1

                      if(ne_quad.gt.ne_l) then
                        write(*,*) ' NE_QUAD TOO LARGE',
     &                              n_el,ne_quad,ne_l
                        call plstop(.true.)
                      endif
                      ix_quad(38,ne_quad) = ma
                      ix_quad(1:37,ne_quad) = 0

!                     Construct 27-node 'quadratic elements'

                      do ns = 1,27

                        sg3(1,1) = gp(1,-1+2*is+ixq3d(1,ns))
                        sg3(2,1) = gp(2,-1+2*js+ixq3d(2,ns))
                        sg3(3,1) = gp(3,-1+2*ks+ixq3d(3,ns))
                        sg3(4,1) = 1.0d0

                        call interp3d(1, xl, ndm,nel)

                        do j = 1,ndm
                          xx(j) = 0.0d0
                          do nn = 1,nel
                            xx(j) = xx(j) + shp3(4,nn,1)*xl(j,nn)
                          end do ! nn
                        end do ! j

!                       Check for repeated nodes

                        noskip = .true.
                        do i = nd_quad,nd_min,-1
                          dst = (x_quad(1,i) - xx(1))**2
     &                         + (x_quad(2,i) - xx(2))**2
     &                         + (x_quad(3,i) - xx(3))**2

                          if(dst.le.dst_min) then
                            ix_quad(ns,ne_quad) = i
                            noskip = .false.
                            exit
                          endif
                        end do ! i

!                       New node

                        if(noskip) then
                          nd_quad            = nd_quad + 1
                          if(nd_quad.gt.nd_l) then
                            write(*,*) ' ND_QUAD TOO LARGE',
     &                                   n_el,ns,nd_quad,nd_l
                            call plstop(.true.)
                          endif
                          ix_quad(ns,ne_quad) = nd_quad
                          ix_quad(34,ne_quad) = (-5)
                          do j = 1,ndm
                            x_quad(j,nd_quad) = xx(j)
                          end do ! j
                        endif

                      end do ! ns

!                     Check for positive jacobian

                      do i = 1,27
                        if(ix_quad(i,ne_quad).gt.0) then
                          do j = 1,ndm
                            xl8(j,i) = x_quad(j,ix_quad(i,ne_quad))
                          end do ! j
                        else
                          do j = 1,ndm
                            xl8(j,i) = 0.0d0
                          end do ! j
                        endif
                      end do ! i
                      ! call ck_lin8(ix_quad(1,ne_quad),xl8)

                    end do ! is
                  end do ! js
                end do ! ks

              endif ! ie(1,ma) test
            else  ! eltyp < 0
!             write(*,*) ' ELTYP = ',eltyp
            endif ! eltyp
          endif ! ma tests

        end do ! ne
        nd_min = nd_quad + 1
      end do ! na

      end subroutine pmsh_quad_lgrn_nusd
