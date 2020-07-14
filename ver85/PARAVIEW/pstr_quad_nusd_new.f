!$Id:$
      subroutine pstr_quad_nusd_new
     &           (ix, x,xl, wt,wtl, ix_quad, st, s_quad, p_quad,
     &                    u, u_quad, ip,imat, jp,jmat)

!      * * F E A P * * A Finite Element Analysis Program

!....  Copyright (c) 1984-2017: Regents of the University of California
!                               All rights reserved

!-----[--.----+----.----+----.-----------------------------------------]
!     Modification log                                Date (dd/mm/year)
!       Original version                                    01/11/2006
!       1. Add line plot capability                         17/01/2011
!       2. Dimension ip(sdim,*) to ip(3,*)                  08/17/2011
!       3. Correct projection procedure to be same as for   05/03/2012
!          coordinate choice.  Add x_lin to argument list.
!       4. Add loop over material numbers                   31/03/2012
!       5. Add d_min to separate nodes by materials         28/09/2012
!       6. Remove loop over material sets and tests on ma   18/12/2014
!       7. change ie(1,ma) to ndm for call to pltstr        23/01/2015
!       8. Set dist_min in pmsh_lin; add to qudshp.h        23/04/2016
!       9. Correct index on u(i ... to u(j for 1-d plots    29/04/2016
!      10. Modify to allow discontinuous stress projections 14/07/2016
!      11. Modify 3-d stress projection for discontinuous   25/02/2017
!      13. updated for npvi and nusd functionality          02/03/2018
!-----[--.----+----.----+----.-----------------------------------------]
!      Purpose: Projects element variables for T-spline/NURBS plots

!      Inputs:
!        ix(nen1,*)     - Mesh connection list
!        x(ndm,*)       - Mesh control points
!        xl(ndm,*)      - Element control points
!        wt(*)          - Mesh weights
!        wtl(*)         - Element weights
!        ix_quad(31,*)  - Quadratic element connectionlist
!        st(numnm,*)    - Mesh projected values
!        u(ndf,*)       - Mesh solution values
!        x_quad(ndm,*)  - Quadratic mesh coordinates
!        ip(0:numnp)    - Pointer for materials
!        imat(*)        - Materials for each node
!        jp(0:numni)    - Quadratic mesh pointer for materials
!        jmat(*)        - Quadratic mesh materials for each node

!      Outputs:
!        s_quad(numni,*) - Plot projected values
!        p_quad(numni,*) - Plot projected values
!        u_quad(ndf,*)   - Plot projected values
!-----[--+---------+---------+---------+---------+---------+---------+-]
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

      include   'npvi.h'

      integer  :: ix(nen1,*), ix_quad(31,*)
      integer  :: ip(0:numnp), imat(*),  jp(0:numni), jmat(*)
      real    (kind=8) :: x(ndm,*), xl(ndm,*), wt(*), wtl(*)
      real    (kind=8) :: s_quad(numni,*),p_quad(numni,*), st(numnm,*)
      real    (kind=8) :: u_quad(ndf,*),u(ndf,*)

      integer  :: is,js,ks, nd,ne,nm,nn,ns, nod
      integer  :: i, ii, j, i_quad, d_quad
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

      lint   = 1

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

!     Loop over quadratic elements

      do nn = 1,ne_quad

        ne = ix_quad(22,nn)
        nd = ix_quad(23,nn)
        ma = ix_quad(31,nn)

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

            is = ix_quad(24,ne)

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

!             Project stress variables

              do i_quad = jp(d_quad-1)+1,jp(d_quad)
                if(jmat(i_quad).eq.ma) then
                  s_quad(i_quad,1) = 1.0d0
                  do i = 2,npstr
                    s_quad(i_quad,i) = 0.0d0
                  end do ! i
                  do nm = 1,nel
                    nod = ix(nm,ne)
                    if(nod.gt.0) then

                      do ii = ip(nod-1)+1,ip(nod)
                        if(ma.eq.imat(ii)) then
                          do i = 2,npstr
                            s_quad(i_quad,i) = s_quad(i_quad,i)
     &                                     + shp1(2,nm,1)*st(ii,i)
                          end do ! i
                          endif ! ma check
                      end do ! ii
                    endif ! nod > 0
                  end do ! nm

                  if(histpltfl) then
                    call hp_lin(ix(1,ne),shp1,hr(np(305)),2,nel,
     &                          hr(np(306)),d_quad,ip,imat, jp,jmat)
                  endif
                endif ! jmat(i_quad)
              end do ! i_quad

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

!             Project stress variables

              if(tsplfl .or. hsplfl) then

                call tsp_quad2d(n_el,nd,sg2,   mr(np(295)),st(1,1),
     &                     mr(np(297)),hr(np(289)),mr(np(290)),
     &                     mr(np(291)),mr(np(292)),hr(np(264)),
     &                     s_quad(1,1) ,nel)

!             NURBS block

              else

                do i_quad = jp(d_quad-1)+1,jp(d_quad)
                  if(jmat(i_quad).eq.ma) then
                    s_quad(i_quad,1) = 1.0d0
                    do i = 2,npstr
                      s_quad(i_quad,i) = 0.0d0
                    end do ! i
                    do nm = 1,nel
                      nod = ix(nm,ne)
                      if(nod.gt.0) then
                        do ii = ip(nod-1)+1,ip(nod)
                          if(ma.eq.imat(ii)) then
                            do i = 2,npstr
                              s_quad(i_quad,i) = s_quad(i_quad,i)
     &                                       + shp2(3,nm,1)*st(ii,i)
                            end do ! i
                          endif ! ma check
                        end do ! ii
                      endif ! nod > 0
                    end do ! nm

                    if(histpltfl) then
                      call hp_lin(ix(1,ne),shp2,hr(np(305)),3,nel,
     &                            hr(np(306)),d_quad,
     &                            ip,imat, jp,jmat)
                    endif
                  endif ! jmat(i_quad)
                end do ! i_quad

              endif ! tsplfl .or. hsplfl

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

              d_quad = ix_quad(ns,nn)

!             Project dependent variables

              u_quad(1:ndf,d_quad) = 0.0d0
              do nm = 1,nel
                if(ix(nm,ne).gt.0) then
                  u_quad(1:ndf,d_quad) = u_quad(1:ndf,d_quad)
     &                               + shp3(4,nm,1)*u(1:ndf,ix(nm,ne))
                endif
              end do ! nm

!             Project stress variables

              if(tsplfl .or. hsplfl) then

                call tsp_quad3d(n_el,nd,sg3,  mr(np(295)),st(1,1),
     &                         mr(np(297)),hr(np(289)),mr(np(290)),
     &                         mr(np(291)),mr(np(292)),hr(np(264)),
     &                         s_quad(1,1),nel)

!             NURBS block

              else

                do i_quad = jp(d_quad-1)+1,jp(d_quad)
                  if(jmat(i_quad).eq.ma) then
                    s_quad(i_quad,1) = 1.0d0
                    do i = 2,npstr
                      s_quad(i_quad,i) = 0.0d0
                    end do ! i
                    do nm = 1,nel
                      nod = ix(nm,ne)
                      if(nod.gt.0) then
                        do ii = ip(nod-1)+1,ip(nod)
!               write(*,*) 'NOD;',nod,ma,imat(ii)
                          if(ma.eq.imat(ii)) then
                            do i = 2,npstr
                              s_quad(i_quad,i) = s_quad(i_quad,i)
     &                                       + shp3(4,nm,1)*st(ii,i)
                            end do ! i
                          endif ! ma check
                        end do ! ii
                      endif ! nod > 0
                    end do ! nm

                    if(histpltfl) then
                      call hp_lin(ix(1,ne),shp3,hr(np(305)),4,nel,
     &                            hr(np(306)),d_quad,
     &                            ip,imat, jp,jmat)
                    endif
                  endif ! jmat(i_quad)
                end do ! i_quad

!               s_lin(d_lin,1) = 1.0d0
!               do i = 2,npstr
!                 s_lin(d_lin,i) = 0.0d0
!               end do ! i
!               do nm = 1,nel
!                 if(ix(nm,ne).gt.0) then
!                   do i = 2,npstr
!                     s_lin(d_lin,i) = s_lin(d_lin,i)
!     &                              + shp3(4,nm,1)*st(ix(nm,ne),i)
!                   end do ! i
!                 endif ! ix > 0
!               end do ! nm
!               if(histpltfl) then
!                 call hp_lin(ix(1,ne),shp3,hr(np(305)),4,nel,
!     &                                hr(np(306)),d_lin,
!     &                                ip,imat, jp,jmat)
!               endif

              endif ! tsplfl .or. hsplfl

            end do ! ns

          endif

        endif ! eltyp

      end do ! ne

      call pltstr(s_quad,p_quad(1,2),s_quad(1,2),jp(nd_quad),ndm,
     &           .false.)

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

      save

      pp  = ip(1,n)
      qq  = ip(2,n)
      nll = (pp+1)*(qq+1)

!     Multiply by extraction operator to compute Bezier weights

      isb = el_e(1,n) - 1
      jsb = el_e(2,n) - 1

      do j = 1,nll
        wb(j)   = 0.0d0
      end do ! j

!     Compute Bezier weights

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
      real*8     sg(3),c_e(*),wt(*),st(numpln,*), s_quad(nd_quad,*)

      integer    pp,qq,rr,nll, i,j
      integer    isu,jsu,liu, isb,jsb,lib,nuni,nbiv
      real*8     shpl(4,100), wb(100)

      save

      pp  = ip(1,n)
      qq  = ip(2,n)
      rr  = ip(3,n)
      nll = (pp+1)*(qq+1)*(rr+1)

      nuni = rr + 1
      nbiv = nel/nuni

!     Project weights to Bezier ones

      do j = 1,nll
        wb(j)   = 0.0d0
      end do ! j

!     Multiply by extraction operator

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

