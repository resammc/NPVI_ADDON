      subroutine pstr_lin_nusd(ie,ix,x,xl,wt,wtl,ix_lin,st,s_lin,
     &                         p_lin,u,u_lin,x_lin)

c      for F E A P -- A Finite Element Analysis Program
c
c      coded by:
c                B.Sc. Henning Venghaus
c                Institute of Mechanics
c                Otto von Guericke University
c                Spring 2017

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    27/03/2017
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Projects element variables for T-spline/NURBS plots

c      Inputs:
c        ie(nie,*)    - Element type properties
c        ix(nen1,*)   - Mesh connection list
c        x(ndm,*)     - Mesh control points
c        xl(ndm,*)    - Element control points
c        wt(*)        - Mesh weights
c        wtl(*)       - Element weights
c        ix_lin(19,*) - Linear element connectionlist
c        st(numnp,*)  - Mesh projected values
c        u(ndf,*)     - Mesh solution values
c        x_lin(ndm,*) - Linear mesh coordinates

c      Outputs:
c        s_lin(nd_lin,*) - Plot projected values
c        p_lin(nd_lin,*) - Plot projected values
c        u_lin(ndf,*)    - Plot projected values
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none

      include   'cdata.h'
      include   'cdat1.h'
      include   'cnurb.h'
      include   'eldata.h'
      include   'eldatp.h'
      include   'iofile.h'
      include   'pbody.h'
      include   'pdata3.h'
      include   'qudshp.h'
      include   'sdata.h'
      include   'pointer.h'
      include   'comblk.h'
      
      include   'npvi.h'

      integer    ie(nie,*), ix(nen1,*), ix_lin(19,*)
      real*8     x(ndm,*), xl(ndm,*), wt(*), wtl(*), st(numnp,*), xx(3)
      real*8     s_lin(nd_lin,*),p_lin(nd_lin,*),u_lin(ndf,*),u(ndf,*)
      real*8     x_lin(ndm,*)

      logical    noskip
      integer    nd, ne, is,js,ks, nn, ns, i,j, na
      integer    ixl(3,8)
      real*8     gp(3,9), du, dist, dist_min

      integer    e_lin, d_lin, d_min

      save


      data       ixl / 0,0,0,   1,0,0,  1,1,0,   0,1,0,
     &                 0,0,1,   1,0,1,  1,1,1,   0,1,1 /

!       Set dist_min for mesh searches

	    dist_min = 1.0d-03*max(hsize(1),hsize(2)*1.d-3)**2
        if(dist_min.eq.0.0d0) then
          dist_min = 1.d-12
        endif
		
c     Get Gaussian points and weights

      lint   = 1
      
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

      d_lin = 0
      d_min = 1
      e_lin = 0
      do na = 1,nummat
        do ne = 1,numel
          n  = ne   ! Element number for use
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
                  xl(1,i) = x(1,ix(i,ne))
                  xl(2,i) = x(2,ix(i,ne))
                  if(ndm.eq.3) xl(3,i) = x(3,ix(i,ne))
                  wtl(i)  = wt(ix(i,ne))
                  nel     = i
                else
                  xl(1,i) = 0.0d0
                  xl(2,i) = 0.0d0
                  if(ndm.eq.3) xl(3,i) = 0.0d0
                  wtl(i)  = 1.0d0
                endif
              end do ! i

c             Do 1-d plots

              if(ie(1,ma).eq.1) then

c               For element loop over number of subdivisions 'npl_int'

                do is = 1,npl_int_nusd(1)

                  e_lin = e_lin + 1

c                 Construct 2-node 'linear elements'

                  do ns = 1,2

                    nd = ix_lin(ns,e_lin)

                    sg1(1,1) = gp(1,is+ixl(1,ns))
                    sg1(2,1) = 1.0d0

                    call interp1d(1, xl, ndm,nel, .false.)

c                   Compute position

                    do i = 1,ndm
                      xx(i) = 0.0d0
                      do nn = 1,nel
                        xx(i) = xx(i) + shp1(2,nn,1)*xl(i,nn)
                      end do ! nn
                    end do ! i

c                   Check for repeated nodes

                    noskip = .true.
                    do i = d_lin,d_min,-1
                      dist = (x_lin(1,i) - xx(1))**2
                      if(ndm.ge.2) dist = dist + (x_lin(2,i)-xx(2))**2
                      if(ndm.eq.3) dist = dist + (x_lin(3,i)-xx(3))**2

                      if(dist.le.dist_min) then
                        noskip = .false.
                        exit
                      endif
                    end do ! i

c                   New node

                    if(noskip) then
                      d_lin = d_lin + 1

c                     Project dependent variables

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

c                     Project stress variables

                      s_lin(d_lin,1) = 1.0d0
                      do i = 2,npstr
                        s_lin(d_lin,i) = 0.0d0
                      end do ! i
                      do nn = 1,nel
                        if(ix(nn,ne).gt.0) then
                          do i = 2,npstr
                            s_lin(d_lin,i) = s_lin(d_lin,i)
     &                              + shp1(2,nn,1)*st(ix(nn,ne),i)
                          end do ! i
                        endif ! ix > 0
                      end do ! nn
                      if(histpltfl) then
                        call hp_lin(ix(1,ne),shp1,hr(np(305)),2,nel,
     &                              numnp,hr(np(306)),d_lin,nd_lin)
                      endif

                    endif ! noskip

                  end do ! ns
                end do ! is

c             Do 2-d plots

              elseif(ie(1,ma).eq.2) then

c               For element loop over number of subdivisions 'npl_int'

                do js = 1,npl_int_nusd(2)
                  do is = 1,npl_int_nusd(1)

                    e_lin = e_lin + 1

c                   Construct 4-node 'linear elements'

                    do ns = 1,4

                      sg2(1,1) = gp(1,is+ixl(1,ns))
                      sg2(2,1) = gp(2,js+ixl(2,ns))
                      sg2(3,1) = 1.0d0

                      call interp2d(1, xl, ix(1,ne), ndm,nel, .false.)

                      if(tsplfl .or. hsplfl) then

c                       Project dependent variables

                        nd = ix_lin(ns,e_lin)
                        do i = 1,ndf
                          u_lin(i,nd) = 0.0d0
                        end do ! i
                        do nn = 1,nel
                          do i = 1,ndf
                            u_lin(i,nd) = u_lin(i,nd)
     &                              + shp2(3,nn,1)*u(i,ix(nn,ne))
                          end do ! nn
                        end do ! i

c                       Project stress variables

                        call tsp_lin2d(n,nd,sg2,   mr(np(295)),st(1,1),
     &                             mr(np(297)),hr(np(289)),mr(np(290)),
     &                             mr(np(291)),mr(np(292)),hr(np(264)),
     &                             s_lin(1,1) ,nel)

c                     NURBS block

                      else

c                       Compute position

                        do i = 1,ndm
                          xx(i) = 0.0d0
                          do nn = 1,nel
                            xx(i) = xx(i) + shp2(3,nn,1)*xl(i,nn)
                          end do ! nn
                        end do ! i

c                       Check for repeated nodes

                        noskip = .true.
                        do i = d_lin,d_min,-1
                          dist = (x_lin(1,i) - xx(1))**2
     &                         + (x_lin(2,i) - xx(2))**2
                          if(ndm.eq.3) then
                            dist = dist + (x_lin(3,i) - xx(3))**2
                          endif

                          if(dist.le.dist_min) then
                            noskip = .false.
                            exit
                          endif
                        end do ! i

c                       New node

                        if(noskip) then
                          d_lin = d_lin + 1

c                         Project dependent variables

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

                          s_lin(d_lin,1) = 1.0d0
                          do i = 2,npstr
                            s_lin(d_lin,i) = 0.0d0
                          end do ! i
                          do nn = 1,nel
                            if(ix(nn,ne).gt.0) then
                              do i = 2,npstr
                                s_lin(d_lin,i) = s_lin(d_lin,i)
     &                                  + shp2(3,nn,1)*st(ix(nn,ne),i)
                              end do ! i
                            endif ! ix > 0
                          end do ! nn
                          if(histpltfl) then
                            call hp_lin(ix(1,ne),shp2,hr(np(305)),3,nel,
     &                                  numnp,hr(np(306)),d_lin,nd_lin)
                          endif

                        endif ! noskip

                      endif ! tsplfl .or. hsplfl

                    end do ! ns

                  end do ! is
                end do ! js

c             Do 3-d plots

              elseif(ie(1,ma).eq.3) then

c               For element loop over number of subdivisions 'npl_int'

                do ks = 1,npl_int_nusd(3)
                  do js = 1,npl_int_nusd(2)
                    do is = 1,npl_int_nusd(1)

                      e_lin = e_lin + 1

c                     Construct 8-node 'linear elements'

                      do ns = 1,8

                        sg3(1,1) = gp(1,is+ixl(1,ns))
                        sg3(2,1) = gp(2,js+ixl(2,ns))
                        sg3(3,1) = gp(3,ks+ixl(3,ns))
                        sg3(4,1) = 1.0d0

                        call interp3d(1, xl, ndm,nel)

c                       Project stress variables

                        if(tsplfl .or. hsplfl) then

c                         Project dependent variables

                          nd = ix_lin(ns,ne_lin)
                          do i = 1,ndf
                            u_lin(i,nd) = 0.0d0
                          end do ! i
                          do nn = 1,nel
                            do i = 1,ndf
                              u_lin(i,nd) = u_lin(i,nd)
     &                                + shp3(4,nn,1)*u(i,ix(nn,ne))
                            end do ! nn
                          end do ! i

                          call tsp_lin3d(n,nd,sg3,  mr(np(295)),st(1,1),
     &                              mr(np(297)),hr(np(289)),mr(np(290)),
     &                              mr(np(291)),mr(np(292)),hr(np(264)),
     &                              s_lin(1,1),nel)

c                       NURBS block

                        else

c                         Compute position

                          do i = 1,ndm
                            xx(i) = 0.0d0
                            do nn = 1,nel
                              xx(i) = xx(i)+shp3(4,nn,1)*xl(i,nn)
                            end do ! nn
                          end do ! i

c                         Check for repeated nodes

                          noskip = .true.
                          do i = d_lin,d_min,-1
                            dist = (x_lin(1,i) - xx(1))**2
     &                           + (x_lin(2,i) - xx(2))**2
     &                           + (x_lin(3,i) - xx(3))**2

                            if(dist.le.dist_min) then
                              noskip = .false.
                              exit
                            endif
                          end do ! i

c                         New node

                          if(noskip) then
                            d_lin = d_lin + 1

c                           Project dependent variables

                            do j = 1,ndf
                              u_lin(j,d_lin) = 0.0d0
                            end do ! j
                            do nn = 1,nel
                              if(ix(nn,ne).gt.0) then
                                do j = 1,ndf
                                  u_lin(j,d_lin) = u_lin(j,d_lin)
     &                                    + shp3(4,nn,1)*u(j,ix(nn,ne))
                                end do ! j
                              endif
                            end do ! nn

c                           Project stress variables

                            s_lin(d_lin,1) = 1.0d0
                            do i = 2,npstr
                              s_lin(d_lin,i) = 0.0d0
                            end do ! i
                            do nn = 1,nel
                              if(ix(nn,ne).gt.0) then
                                do i = 2,npstr
                                  s_lin(d_lin,i) = s_lin(d_lin,i)
     &                                    + shp3(4,nn,1)*st(ix(nn,ne),i)
                                end do ! i
                              endif ! ix > 0
                            end do ! nn
                            if(histpltfl) then
                              call hp_lin(ix(1,ne),shp3,hr(np(305)),4,
     &                                    nel,numnp,hr(np(306)),
     &                                    d_lin,nd_lin)
                            endif

                          endif

                        endif ! tsplfl .or. hsplfl

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

      call pltstr(s_lin,p_lin(1,2),s_lin(1,2),nd_lin,ie(1,ma),.false.)

      end
      
!  usually there would be other routines used in pstr_lin but they are
!  not effected by the anisotropic subdivision change so they dont 
!  appear here - they are used from the old file /plot/pst_lin