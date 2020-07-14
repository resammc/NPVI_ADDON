      subroutine pmsh_lin_nusd(ie,ix, x,xl, wt,wtl,
     &                         ix_lin, x_lin,nd_l,ne_l)

c      for F E A P -- A Finite Element Analysis Program
c
c      coded by:
c                B.Sc. Henning Venghaus
c                Institute of Mechanics
c                Otto von Guericke University
c                Spring 2017

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    26/03/2017
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Prepare mesh of linear elements for 3-d plots with the
c               opportunity of anisotropic element subdivision

c      Inputs:

c      Outputs:
c-----[--+---------+---------+---------+---------+---------+---------+-]
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

      integer    ie(nie,*),ix(nen1,*), ix_lin(19,*)
      real*8     x(ndm,*), xl(ndm,*), x_lin(ndm,*), wt(*), wtl(*)

      logical    noskip
      integer    nd_l,ne_l, nd_min, ne, is,js,ks, nn, ns, i,j, na
      integer    ixl(3,8)
      real*8     dist, dist_min, du
      real*8     gp(3,9), xx(3), xl8(3,8)

      save


      data       ixl / 0,0,0,   1,0,0,  1,1,0,   0,1,0,
     &                 0,0,1,   1,0,1,  1,1,1,   0,1,1 /

!       Set dist_min for mesh searches
	 
	    dist_min = 1.0d-03*max(hsize(1),hsize(2)*1.d-3)**2
        if(dist_min.eq.0.0d0) then
          dist_min = 1.d-12
        endif
		
		
c     Get Gaussian points and weights

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
          gp(i,npl_int_nusd(i)+1) = 0.999999999999d0 ! Keep point inside element
        endif
      end do !i
      
      

      nd_lin = 0
      nd_min = 1
      ne_lin = 0
      do na = 1,nummat
        do ne = 1,numel
          n  = ne   ! Element number for use
          ma = ix(nen1,ne)
          if((ma.eq.na .and. maplt.eq.0 ) .or.
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

c             For element loop over number of subdivisions 'npl_int'

c             Line mesh

              if(ie(1,ma).eq.1) then

                do is = 1,npl_int_nusd(1)

                  ne_lin = ne_lin + 1

                  if(ne_lin.gt.ne_l) then
                    write(*,*) ' NE_LIN TOO LARGE',n,ne_lin,ne_l
                    call plstop()
                  endif

                  ix_lin(19,ne_lin) = ma

c                 Construct 2-node 'linear elements'

                  do ns = 1,2

                    sg1(1,1) = gp(1,is+ixl(1,ns))
                    sg1(2,1) = 1.0d0

                    call interp1d(1, xl, ndm,nel, .false.)

                    do j = 1,ndm
                      xx(j) = 0.0d0
                      do nn = 1,nel
                        xx(j) = xx(j) + shp1(2,nn,1)*xl(j,nn)
                      end do ! nn
                    end do ! j

c                   Check for repeated nodes

                    noskip = .true.
                    do i = nd_lin,nd_min,-1
                      dist = (x_lin(1,i) - xx(1))**2
                      if(ndm.ge.2) dist = dist + (x_lin(2,i)-xx(2))**2
                      if(ndm.eq.3) dist = dist + (x_lin(3,i)-xx(3))**2

                      if(dist.le.dist_min) then
                        ix_lin(ns,ne_lin) = i
                        noskip = .false.
                        exit
                      endif
                    end do ! i

c                   New node

                    if(noskip) then
                      nd_lin            = nd_lin + 1
                      if(nd_lin.gt.nd_l) then
                        write(*,*) ' ND_LIN TOO LARGE',n,ns,nd_lin,nd_l
                        call plstop()
                      endif
                      ix_lin(ns,ne_lin) = nd_lin
					  ix_lin(15,ne_lin) = (-1)
                      do j = 1,ndm
                        x_lin(j,nd_lin)   = xx(j)
                      end do ! j
                    endif

                  end do ! ns
                end do ! is

c             Surface mesh

              elseif(ie(1,ma).eq.2) then
                do js = 1,npl_int_nusd(2)
                  do is = 1,npl_int_nusd(1)

                    ne_lin = ne_lin + 1

                    if(ne_lin.gt.ne_l) then
                      write(*,*) ' NE_LIN TOO LARGE',n,ne_lin,ne_l
                      call plstop()
                    endif

                    ix_lin(19,ne_lin) = ma

c                   Construct 4-node 'linear elements'

                    do ns = 1,4

                      sg2(1,1) = gp(1,is+ixl(1,ns))
                      sg2(2,1) = gp(2,js+ixl(2,ns))
                      sg2(3,1) = 1.0d0

                      call interp2d(1, xl, ix(1,ne), ndm,nel, .false.)

                      do j = 1,ndm
                        xx(j) = 0.0d0
                        do nn = 1,nel
                          xx(j) = xx(j) + shp2(3,nn,1)*xl(j,nn)
                        end do ! nn
                      end do ! j

c                     Check for repeated nodes

                      noskip = .true.
                      do i = nd_lin,nd_min,-1
                        dist = (x_lin(1,i) - xx(1))**2
     &                       + (x_lin(2,i) - xx(2))**2
                        if(ndm.eq.3) then
                          dist = dist + (x_lin(3,i)-xx(3))**2
                        endif

                        if(dist.le.dist_min) then
                          ix_lin(ns,ne_lin) = i
                          noskip = .false.
                          exit
                        endif
                      end do ! i

c                     New node

                      if(noskip) then
                        nd_lin            = nd_lin + 1
                        if(nd_lin.gt.nd_l) then
                          write(*,*) ' ND_LIN TOO LARGE',
     &                                 n,ns,nd_lin,nd_l
                          call plstop()
                        endif
                        ix_lin(ns,ne_lin) = nd_lin
						ix_lin(15,ne_lin) = (-3)
                        do j = 1,ndm
                          x_lin(j,nd_lin)   = xx(j)
                        end do ! j
                      endif

                    end do ! ns
                  end do ! is
                end do ! js

c             Solid mesh

              elseif(ie(1,ma).eq.3) then

                do ks = 1,npl_int_nusd(3)
                  do js = 1,npl_int_nusd(2)
                    do is = 1,npl_int_nusd(1)

                      ne_lin = ne_lin + 1

                      if(ne_lin.gt.ne_l) then
                        write(*,*) ' NE_LIN TOO LARGE',n,ne_lin,ne_l
                        call plstop()
                      endif
                      ix_lin(19,ne_lin) = ma
                      ix_lin(1:18,ne_lin) = 0

c                     Construct 8-node 'linear elements'

                      do ns = 1,8

                        sg3(1,1) = gp(1,is+ixl(1,ns))
                        sg3(2,1) = gp(2,js+ixl(2,ns))
                        sg3(3,1) = gp(3,ks+ixl(3,ns))
                        sg3(4,1) = 1.0d0

                        call interp3d(1, xl, ndm,nel)

                        do j = 1,ndm
                          xx(j) = 0.0d0
                          do nn = 1,nel
                            xx(j) = xx(j) + shp3(4,nn,1)*xl(j,nn)
                          end do ! nn
                        end do ! j

c                       Check for repeated nodes

                        noskip = .true.
                        do i = nd_lin,nd_min,-1
                          dist = (x_lin(1,i) - xx(1))**2
     &                         + (x_lin(2,i) - xx(2))**2
     &                         + (x_lin(3,i) - xx(3))**2

                          if(dist.le.dist_min) then
                            ix_lin(ns,ne_lin) = i
                            noskip = .false.
                            exit
                          endif
                        end do ! i

c                       New node

                        if(noskip) then
                          nd_lin            = nd_lin + 1
                          if(nd_lin.gt.nd_l) then
                            write(*,*) ' ND_LIN TOO LARGE',
     &                                   n,ns,nd_lin,nd_l
                            call plstop()
                          endif
                          ix_lin(ns,ne_lin) = nd_lin
						  ix_lin(15,ne_lin) = (-5)
                          do j = 1,ndm
                            x_lin(j,nd_lin) = xx(j)
                          end do ! j
                        endif

                      end do ! ns

c                     Check for positive jacobian

                      do i = 1,8
                        if(ix_lin(i,ne_lin).gt.0) then
                          do j = 1,ndm
                            xl8(j,i) = x_lin(j,ix_lin(i,ne_lin))
                          end do ! j
                        else
                          do j = 1,ndm
                            xl8(j,i) = 0.0d0
                          end do ! j
                        endif
                      end do ! i
                      call ck_lin8(ix_lin(1,ne_lin),xl8)

                    end do ! is
                  end do ! js
                end do ! ks

              endif ! ie(1,ma) test
            else  ! eltyp < 0
c             write(*,*) ' ELTYP = ',eltyp
            endif ! eltyp
          endif ! ma tests

        end do ! ne
        nd_min = nd_lin + 1
      end do ! na

      end
