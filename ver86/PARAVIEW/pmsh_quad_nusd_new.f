!$Id:$
      subroutine pmsh_quad_nusd_new(ie,ix, x,xl, wt,wtl,
     &                         ix_quad, x_quad,nd_l,ne_l)

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
!          Change "dist" to "dst" (in qudshp); n to  n_el
!       2. Update to ver8.5 (new ix_quad array)             02/03/2018
!       3. Commented out the call to ck_lin8 for quad elems 23/06/2020
!-----[--.----+----.----+----.-----------------------------------------]
!      Purpose: Prepare mesh of quadratic elements for 3d plots with the
!               opportunity of anisotropic element subdivision

!      Inputs:

!      Outputs:
!-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none

      include   'cdata.h'
      include   'cdat1.h'
      include   'cnurb.h'
      include   'eldata.h'
      include   'iofile.h'
      include   'qudshp.h'
      include   'sdata.h'
      include   'setups.h'

      include   'npvi.h'

      integer    ie(nie,*),ix(nen1,*), ix_quad(31,*)
      real*8     x(ndm,*), xl(ndm,*), x_quad(ndm,*), wt(*), wtl(*)

      logical    noskip
      integer    nd_l,ne_l, nd_min, ne, is,js,ks, nn, ns, i,j, na
      integer    iqx3d(20), iqy3d(20), iqz3d(20), iqx2d(8), iqy2d(8)
      integer    iqx1d(3)
      real*8     dst, du, dst_min
      real*8     gp(3,20), xx(3), xl8(3,20)

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


!       Set dst_min for mesh searches

        dst_min = 1.0d-03*max(hsize(1),hsize(2)*1.d-3)**2
        if(dst_min.eq.0.0d0) then
          dst_min = 1.d-12
        endif

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
      
      nd_quad = 0
      nd_min = 1
      ne_quad = 0
!      do na = 1,nummat
        do ne = 1,numel
          n_el  = ne   ! Element number for use
          ma = ix(nen1,ne)
          eltyp = ix(nen+7,ne)
          if(eltyp.gt.0) then

            elty2 = ix(nen+8,ne)
            elty3 = ix(nen+9,ne)

!           Set up local values

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

!           For element loop over number of subdivisions 'npl_int'

!           Line mesh
            if(ie(1,ma).eq.1) then

              do is = 1,npl_int_nusd(1)

                ne_quad = ne_quad + 1

                if(ne_quad.gt.ne_l) then
                  if(rank.eq.0) then
                    write(*,*) ' NE_QUAD TOO LARGE',n_el,ne_quad,ne_l
                  endif
                  write(iow,*) ' NE_QUAD TOO LARGE',n_el,ne_quad,ne_l
                  write(ilg,*) ' NE_QUAD TOO LARGE',n_el,ne_quad,ne_l
                  call plstop(.true.)
                endif

                ix_quad(1:30,ne_quad) = 0

!             Save element, dimension, vertex point, material no.

                ix_quad(22,ne_quad) = ne
                ix_quad(23,ne_quad) = 1
                ix_quad(24,ne_quad) = is
                ix_quad(31,ne_quad) = ma
                ix_quad(27,ne_quad) = -1
                
!               Construct 3-node 'quadratic elements'

                do ns = 1,3

                  sg1(1,1) = gp(1,-1+2*is+iqx1d(ns))
                  sg1(2,1) = 1.0d0

                  call interp1d(1, xl, ndm,nel, .false.)

                  do j = 1,ndm
                    xx(j) = 0.0d0
                    do nn = 1,nel
                      xx(j) = xx(j) + shp1(2,nn,1)*xl(j,nn)
                    end do ! nn
                  end do ! j

!             Check for repeated nodes

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

!             New node

              if(noskip) then
                   nd_quad          = nd_quad + 1
               if(nd_quad.gt.nd_l) then
                 if(rank.eq.0) then
                   write(*,*) ' ND_QUAD TOO LARGE',n_el,ns,nd_quad,nd_l
                 endif
                 write(iow,*) ' ND_QUAD TOO LARGE',n_el,ns,nd_quad,nd_l
                 write(ilg,*) ' ND_QUAD TOO LARGE',n_el,ns,nd_quad,nd_l
                 call plstop(.true.)
               endif
               ix_quad(ns,ne_quad) = nd_quad
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
                    if(rank.eq.0) then
                      write(*,*) ' NE_QUAD TOO LARGE',n_el,ne_quad,ne_l
                    endif
                    write(iow,*) ' NE_QUAD TOO LARGE',n_el,ne_quad,ne_l
                    write(ilg,*) ' NE_QUAD TOO LARGE',n_el,ne_quad,ne_l
                    call plstop(.true.)
                  endif

                  ix_quad(1:30,ne_quad) = 0

!                 Save element, dimension, vertex point, material no.

                  ix_quad(22,ne_quad) = ne
                  ix_quad(23,ne_quad) = 2
                  ix_quad(24,ne_quad) = is
                  ix_quad(25,ne_quad) = js
                  ix_quad(31,ne_quad) = ma
                  ix_quad(27,ne_quad) = -3
                  
!                 Construct 8-node 'quadratic elements'

                    do ns = 1,8

                      sg2(1,1) = gp(1,-1+2*is+iqx2d(ns))
                      sg2(2,1) = gp(2,-1+2*js+iqy2d(ns))
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
     &                      + (x_quad(2,i) - xx(2))**2
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

                      if(noskip) then
                        nd_quad            = nd_quad + 1
                      if(nd_quad.gt.nd_l) then
                        if(rank.eq.0) then
                          write(*,*) ' ND_QUAD TOO LARGE',
     &                                 n_el,ns,nd_quad,nd_l
                        endif
                        write(iow,*) ' ND_QUAD TOO LARGE',
     &                                 n_el,ns,nd_quad,nd_l
                        write(ilg,*) ' ND_QUAD TOO LARGE',
     &                                 n_el,ns,nd_quad,nd_l
                        call plstop(.true.)
                      endif

                        ix_quad(ns,ne_quad) = nd_quad
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
                    if(rank.eq.0) then
                      write(*,*) ' NE_QUAD TOO LARGE',n_el,ne_quad,ne_l
                    endif
                    write(iow,*) ' NE_QUAD TOO LARGE',n_el,ne_quad,ne_l
                    write(ilg,*) ' NE_QUAD TOO LARGE',n_el,ne_quad,ne_l
                    call plstop(.true.)
                  endif

                  ix_quad(1:30,ne_quad) = 0

!                 Save element, dimension, vertex point, material no.

                    ix_quad(22,ne_quad) = ne
                    ix_quad(23,ne_quad) = 3
                    ix_quad(24,ne_quad) = is
                    ix_quad(25,ne_quad) = js
                    ix_quad(26,ne_quad) = ks
                    ix_quad(31,ne_quad) = ma
                    ix_quad(27,ne_quad) = -5

!                 Construct 20-node 'quadratic elements'

                    do ns = 1,20

                      sg3(1,1) = gp(1,-1+2*is+iqx3d(ns))
                      sg3(2,1) = gp(2,-1+2*js+iqy3d(ns))
                      sg3(3,1) = gp(3,-1+2*ks+iqz3d(ns))
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
     &                        + (x_quad(2,i) - xx(2))**2
     &                        + (x_quad(3,i) - xx(3))**2

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
                            if(rank.eq.0) then
                              write(*,*) ' ND_QUAD TOO LARGE',
     &                                     n_el,ns,nd_quad,nd_l
                            endif
                            write(iow,*) ' ND_QUAD TOO LARGE',
     &                                     n_el,ns,nd_quad,nd_l
                            write(ilg,*) ' ND_QUAD TOO LARGE',
     &                                     n_el,ns,nd_quad,nd_l
                            call plstop(.true.)
                          endif
                          ix_quad(ns,ne_quad) = nd_quad
                          do j = 1,ndm
                            x_quad(j,nd_quad) = xx(j)
                          end do ! j
                        endif

                      end do ! ns

!                     Check for positive jacobian

                      do i = 1,20
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

        end do ! ne
        nd_min = nd_quad + 1
!      end do ! na

      end subroutine pmsh_quad_nusd_new
