      subroutine pstr_quad_nusd(ie,ix, x,xl, wt,wtl, ix_quad, st, 
     &                          s_quad,p_quad, u, u_quad, x_quad)

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
c      Purpose: Projects element variables for T-spline/NURBS plots
c               routine works like pstr_lin.f but for quad. elements

c      Inputs:
c        ie(nie,*)    - Element type properties
c        ix(nen1,*)   - Mesh connection list
c        x(ndm,*)     - Mesh control points
c        xl(ndm,*)    - Element control points
c        wt(*)        - Mesh weights
c        wtl(*)       - Element weights
c        ix_quad(31,*)- Linear element connectionlist
c        st(numnp,*)  - Mesh projected values
c        u(ndf,*)     - Mesh solution values
c        x_quad(ndm,*)- Linear mesh coordinates

c      Outputs:
c        s_quad(nd_quad,*) - Plot projected values
c        p_quad(nd_quad,*) - Plot projected values
c        u_quad(ndf,*)     - Plot projected values
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

      integer    ie(nie,*), ix(nen1,*), ix_quad(31,*)
      real*8     x(ndm,*), xl(ndm,*), wt(*), wtl(*), st(numnp,*), xx(3)
      real*8     s_quad(nd_quad,*),p_quad(nd_quad,*),u_quad(ndf,*)
      real*8     u(ndf,*)
      real*8     x_quad(ndm,*)

      logical    noskip
      integer    nd, ne, is,js,ks, nn, ns, i,j, na
      integer    ixq3d(3,20),ixq2d(2,8),ixq1d(1,3)
      real*8     gp(3,20), du, dist, dist_min

      integer    e_quad, d_quad, d_min

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

      d_quad = 0
      d_min = 1
      e_quad = 0
      
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

                  e_quad = e_quad + 1

c                 Construct 3-node 'quadratic elements'

                  do ns = 1,3

                    nd = ix_quad(ns,e_quad)

                    sg1(1,1) = gp(1,-1+2*is+ixq1d(1,ns))
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
                    do i = d_quad,d_min,-1
                      dist = (x_quad(1,i) - xx(1))**2
                      if(ndm.ge.2) dist = dist + (x_quad(2,i)-xx(2))**2
                      if(ndm.eq.3) dist = dist + (x_quad(3,i)-xx(3))**2

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

c                     Project stress variables

                      s_quad(d_quad,1) = 1.0d0
                      do i = 2,npstr
                        s_quad(d_quad,i) = 0.0d0
                      end do ! i
                      do nn = 1,nel
                        if(ix(nn,ne).gt.0) then
                          do i = 2,npstr
                            s_quad(d_quad,i) = s_quad(d_quad,i)
     &                              + shp1(2,nn,1)*st(ix(nn,ne),i)
                          end do ! i
                        endif ! ix > 0
                      end do ! nn
                      if(histpltfl) then
                        call hp_lin(ix(1,ne),shp1,hr(np(305)),2,nel,
     &                              numnp,hr(np(306)),d_quad,nd_quad)
                      endif

                    endif ! noskip

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

                      call interp2d(1, xl, ix(1,ne), ndm,nel, .false.)

                      if(tsplfl .or. hsplfl) then

c                       Project dependent variables

                        nd = ix_quad(ns,e_quad)
                        do i = 1,ndf
                          u_quad(i,nd) = 0.0d0
                        end do ! i
                        do nn = 1,nel
                          do i = 1,ndf
                            u_quad(i,nd) = u_quad(i,nd)
     &                              + shp2(3,nn,1)*u(i,ix(nn,ne))
                          end do ! nn
                        end do ! i

c                       Project stress variables

                        call tsp_quad2d(n,nd,sg2, mr(np(295)),st(1,1),
     &                             mr(np(297)),hr(np(289)),mr(np(290)),
     &                             mr(np(291)),mr(np(292)),hr(np(264)),
     &                             s_quad(1,1) ,nel)

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
                        do i = d_quad,d_min,-1
                          dist = (x_quad(1,i) - xx(1))**2
     &                         + (x_quad(2,i) - xx(2))**2
                          if(ndm.eq.3) then
                            dist = dist + (x_quad(3,i) - xx(3))**2
                          endif

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
     &                                   + shp2(3,nn,1)*u(j,ix(nn,ne))
                              end do ! j
                            endif
                          end do ! nn

                          s_quad(d_quad,1) = 1.0d0
                          do i = 2,npstr
                            s_quad(d_quad,i) = 0.0d0
                          end do ! i
                          do nn = 1,nel
                            if(ix(nn,ne).gt.0) then
                              do i = 2,npstr
                                s_quad(d_quad,i) = s_quad(d_quad,i)
     &                                  + shp2(3,nn,1)*st(ix(nn,ne),i)
                              end do ! i
                            endif ! ix > 0
                          end do ! nn
                          if(histpltfl) then
                            call hp_lin(ix(1,ne),shp2,hr(np(305)),3,nel,
     &                               numnp,hr(np(306)),d_quad,nd_quad)
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

                      e_quad = e_quad + 1

c                     Construct 20-node 'quadratic elements'

                      do ns = 1,20

                        sg3(1,1) = gp(1,-1+2*is+ixq3d(1,ns))
                        sg3(2,1) = gp(2,-1+2*js+ixq3d(2,ns))
                        sg3(3,1) = gp(3,-1+2*ks+ixq3d(3,ns))
                        sg3(4,1) = 1.0d0

                        call interp3d(1, xl, ndm,nel)

c                       Project stress variables

                        if(tsplfl .or. hsplfl) then

c                         Project dependent variables

                          nd = ix_quad(ns,ne_quad)
                          do i = 1,ndf
                            u_quad(i,nd) = 0.0d0
                          end do ! i
                          do nn = 1,nel
                            do i = 1,ndf
                              u_quad(i,nd) = u_quad(i,nd)
     &                                + shp3(4,nn,1)*u(i,ix(nn,ne))
                            end do ! nn
                          end do ! i

                          call tsp_quad3d(n,nd,sg3,mr(np(295)),st(1,1),
     &                              mr(np(297)),hr(np(289)),mr(np(290)),
     &                              mr(np(291)),mr(np(292)),hr(np(264)),
     &                              s_quad(1,1),nel)

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
                          do i = d_quad,d_min,-1
                            dist = (x_quad(1,i) - xx(1))**2
     &                           + (x_quad(2,i) - xx(2))**2
     &                           + (x_quad(3,i) - xx(3))**2

                            if(dist.le.dist_min) then
                              noskip = .false.
                              exit
                            endif
                          end do ! i

c                         New node

                          if(noskip) then
                            d_quad = d_quad + 1

c                           Project dependent variables

                            do j = 1,ndf
                              u_quad(j,d_quad) = 0.0d0
                            end do ! j
                            do nn = 1,nel
                              if(ix(nn,ne).gt.0) then
                                do j = 1,ndf
                                  u_quad(j,d_quad) = u_quad(j,d_quad)
     &                                    + shp3(4,nn,1)*u(j,ix(nn,ne))
                                end do ! j
                              endif
                            end do ! nn

c                           Project stress variables

                            s_quad(d_quad,1) = 1.0d0
                            do i = 2,npstr
                              s_quad(d_quad,i) = 0.0d0
                            end do ! i
                            do nn = 1,nel
                              if(ix(nn,ne).gt.0) then
                                do i = 2,npstr
                                  s_quad(d_quad,i) = s_quad(d_quad,i)
     &                                    + shp3(4,nn,1)*st(ix(nn,ne),i)
                                end do ! i
                              endif ! ix > 0
                            end do ! nn
                            if(histpltfl) then
                              call hp_lin(ix(1,ne),shp3,hr(np(305)),4,
     &                                    nel,numnp,hr(np(306)),
     &                                    d_quad,nd_quad)
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
        d_min = d_quad + 1
      end do ! na

      call pltstr(s_quad,p_quad(1,2),s_quad(1,2),nd_quad,ie(1,ma),
     &            .false.)

      end

      subroutine tsp_quad2d(n,nd,sg,ixbez,st,el_e,c_e,rc_e,ic_e,ip,wt,
     &                     s_quad,nel)

      implicit   none

      include   'cnurb.h'
      include   'pdata3.h'
      
      include   'npvi.h'

      integer    n,nd,nel
      integer    ixbez(nnpl+1,*),el_e(2,*),rc_e(*),ic_e(*),ip(3,*)
      real*8     sg(2),c_e(*),wt(*),st(numpln,*), s_quad(nd_quad,*)

      integer    pp,qq,nll, isb,jsb, i,j, i1
      real*8     shpl(3,64), wb(64)

      pp  = ip(1,n)
      qq  = ip(2,n)
      nll = (pp+1)*(qq+1)

c     Multiply by extraction operator to compute Bezier weights

      isb = el_e(1,n) - 1
      jsb = el_e(2,n) - 1

      do j = 1,nll
        wb(j)   = 0.0d0
      end do ! j

c     Compute Bezier weights

      if(cetype.eq.1) then

        do i = 1,nel
          i1 = isb + i
          do j = ic_e(i1),ic_e(i1+1)-1
            wb(rc_e(jsb+j)) = wb(rc_e(jsb+j)) + c_e(jsb+j)*wt(i)
          end do ! j
        end do ! i

      elseif(cetype.eq.2) then

        i1 = 0
        do i = 1,nel
          do j = 1,nll
            wb(j) = wb(j) + wt(i)*c_e(jsb+j+i1)
          end do ! j
          i1 = i1 + nll
        end do ! i

      endif

      call shp2d_bez(sg,wb, pp,qq, nll, shpl)

      s_quad(nd,1) = 1.0d0
      do i = 2,npstr
        s_quad(nd,i) = 0.0d0
        do j = 1,nll
          s_quad(nd,i) = s_quad(nd,i) + shpl(3,j)*st(ixbez(j,n),i)
        end do ! j
      end do ! i

      end

      subroutine tsp_quad3d(n,nd,sg,ixbez,st,el_e,c_e,rc_e,ic_e,
     &                     ip,wt,s_quad,nel)

      implicit   none

      include   'cnurb.h'
      include   'pdata3.h'
      
      include   'npvi.h'

      integer    n,nd,nel,ixbez(nnpl+1,*),rc_e(*),ic_e(*),ip(3,*)
      integer    el_e(2,*)
      real*8     sg(3),c_e(*),wt(*),st(numpln,*)
      real*8     s_quad(nd_quad,*)

      integer    pp,qq,rr,nll, i,j
      integer    isu,jsu,liu, isb,jsb,lib,nuni,nbiv
      real*8     shpl(4,100), wb(100)

      pp  = ip(1,n)
      qq  = ip(2,n)
      rr  = ip(3,n)
      nll = (pp+1)*(qq+1)*(rr+1)

      nuni = rr + 1
      nbiv = nel/nuni

c     Project weights to Bezier ones

      do j = 1,nll
        wb(j)   = 0.0d0
      end do ! j

c     Multiply by extraction operator

      isu = el_e(1,n)
      jsu = el_e(2,n)

      if(cetype.eq.1) then
        liu = nuni

        isb = isu + liu + 1
        jsb = jsu + ic_e(isu+nuni) - 1
        lib = el_e(1,n+1) - isb - 1

        call sparse_wb(ic_e(isu),ic_e(isb),rc_e(jsu),rc_e(jsb),
     &                 c_e(jsu),c_e(jsb), nuni,nbiv, lib, wt,wb)

      elseif(cetype.eq.2) then
        isb = nll/nuni
        jsb = jsu + nuni*nuni

        call full3d_wb(c_e(jsu),c_e(jsb),nuni,nbiv,isb,wt,wb)

      elseif(cetype.eq.3) then

        call full3d_w(c_e(jsu),nel,nll,wt,wb)

      endif

      call shp3d_bez(sg,wb, pp,qq,rr, nll, shpl)

      s_quad(nd,1) = 1.0d0
      do i = 2,npstr
        s_quad(nd,i) = 0.0d0
        do j = 1,nll
          s_quad(nd,i) = s_quad(nd,i) + shpl(4,j)*st(ixbez(j,n),i)
        end do ! j
      end do ! i

      end