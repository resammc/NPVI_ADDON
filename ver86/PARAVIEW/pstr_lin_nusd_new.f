!$Id:$
      subroutine pstr_lin_nusd_new
     &           (ix, x,xl, wt,wtl, ix_lin, st, s_lin, p_lin,
     &                    u, u_lin, ip,imat, jp,jmat)

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
!        ix_lin(19,*)   - Linear element connectionlist
!        st(numnm,*)    - Mesh projected values
!        u(ndf,*)       - Mesh solution values
!        x_lin(ndm,*)   - Linear mesh coordinates
!        ip(0:numnp)    - Pointer for materials
!        imat(*)        - Materials for each node
!        jp(0:numni)    - Linear mesh pointer for materials
!        jmat(*)        - Linear mesh materials for each node

!      Outputs:
!        s_lin(numni,*) - Plot projected values
!        p_lin(numni,*) - Plot projected values
!        u_lin(ndf,*)   - Plot projected values
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

      integer  :: ix(nen1,*), ix_lin(19,*)
      integer  :: ip(0:numnp), imat(*),  jp(0:numni), jmat(*)
      real    (kind=8) :: x(ndm,*), xl(ndm,*), wt(*), wtl(*)
      real    (kind=8) :: s_lin(numni,*),p_lin(numni,*), st(numnm,*)
      real    (kind=8) :: u_lin(ndf,*),u(ndf,*)

      integer  :: is,js,ks, nd,ne,nm,nn,ns, nod
      integer  :: i, ii, j, i_lin, d_lin
      integer  :: iqx(8), iqy(8), iqz(8)

      real    (kind=8) :: du
      real    (kind=8) :: gp(3,9)

      save

      data    iqx  / 0, 1, 1, 0, 0, 1, 1, 0 /
      data    iqy  / 0, 0, 1, 1, 0, 0, 1, 1 /
      data    iqz  / 0, 0, 0, 0, 1, 1, 1, 1 /

!     Get Gaussian points and weights

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
          gp(i,npl_int_nusd(i)+1) = 0.999999999999d0 ! Keep in elmt
        endif
      end do !i


!     Loop over linear element

      do nn = 1,ne_lin

        ne = ix_lin(10,nn)
        nd = ix_lin(11,nn)
        ma = ix_lin(19,nn)

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

            is = ix_lin(12,nn)

!           Construct 2-node 'linear elements'

            do ns = 1,2

              sg1(1,1) = gp(1,is+iqx(ns))
              sg1(2,1) = 1.0d0

              call interp1d(1, xl, ndm,nel, .true.)

!             Get local node number

              d_lin = ix_lin(ns,nn)

!             Project dependent variables

              u_lin(1:ndf,d_lin) = 0.0d0
              do nm = 1,nel
                if(ix(nm,ne).gt.0) then
                  u_lin(1:ndf,d_lin) = u_lin(1:ndf,d_lin)
     &                               + shp1(2,nm,1)*u(1:ndf,ix(nm,ne))
                endif
              end do ! nm

!             Project stress variables

              do i_lin = jp(d_lin-1)+1,jp(d_lin)
                if(jmat(i_lin).eq.ma) then
                  s_lin(i_lin,1) = 1.0d0
                  do i = 2,npstr
                    s_lin(i_lin,i) = 0.0d0
                  end do ! i
                  do nm = 1,nel
                    nod = ix(nm,ne)
                    if(nod.gt.0) then

                      do ii = ip(nod-1)+1,ip(nod)
                        if(ma.eq.imat(ii)) then
                          do i = 2,npstr
                            s_lin(i_lin,i) = s_lin(i_lin,i)
     &                                     + shp1(2,nm,1)*st(ii,i)
                          end do ! i
                          endif ! ma check
                      end do ! ii
                    endif ! nod > 0
                  end do ! nm

                  if(histpltfl) then
                    call hp_lin(ix(1,ne),shp1,hr(np(305)),2,nel,
     &                          hr(np(306)),d_lin,ip,imat, jp,jmat)
                  endif
                endif ! jmat(i_lin)
              end do ! i_lin

            end do ! ns

!         Do 2-d plots

          elseif(nd.eq.2) then

!           For element loop over number of subdivisions 'npl_int'

            is = ix_lin(12,nn)
            js = ix_lin(13,nn)

!           Construct 4-node 'linear elements'

            do ns = 1,4

              sg2(1,1) = gp(1,is+iqx(ns))
              sg2(2,1) = gp(2,js+iqy(ns))
              sg2(3,1) = 1.0d0

              call interp2d(1, xl, ix(1,ne), ndm,nel, .true.)

!             Set node projection number

              d_lin = ix_lin(ns,nn)

!             Project dependent variables

              u_lin(1:ndf,d_lin) = 0.0d0
              do nm = 1,nel
                if(ix(nm,ne).gt.0) then
                  u_lin(1:ndf,d_lin) = u_lin(1:ndf,d_lin)
     &                               + shp2(3,nm,1)*u(1:ndf,ix(nm,ne))
                endif
              end do ! nm

!             Project stress variables

!              if(tsplfl .or. hsplfl) then

!                call tsp_lin2d(n_el,nd,sg2,   mr(np(295)),st(1,1),
!     &                     mr(np(297)),hr(np(289)),mr(np(290)),
!     &                     mr(np(291)),mr(np(292)),hr(np(264)),
!     &                     s_lin(1,1) ,nel)

!             NURBS block

!              else

                do i_lin = jp(d_lin-1)+1,jp(d_lin)
                  if(jmat(i_lin).eq.ma) then
                    s_lin(i_lin,1) = 1.0d0
                    do i = 2,npstr
                      s_lin(i_lin,i) = 0.0d0
                    end do ! i
                    do nm = 1,nel
                      nod = ix(nm,ne)
                      if(nod.gt.0) then
                        do ii = ip(nod-1)+1,ip(nod)
                          if(ma.eq.imat(ii)) then
                            do i = 2,npstr
                              s_lin(i_lin,i) = s_lin(i_lin,i)
     &                                       + shp2(3,nm,1)*st(ii,i)
                            end do ! i
                          endif ! ma check
                        end do ! ii
                      endif ! nod > 0
                    end do ! nm

                    if(histpltfl) then
                      call hp_lin(ix(1,ne),shp2,hr(np(305)),3,nel,
     &                            hr(np(306)),d_lin,
     &                            ip,imat, jp,jmat)
                    endif
                  endif ! jmat(i_lin)
                end do ! i_lin

!              endif ! tsplfl .or. hsplfl

            end do ! ns

!         Do 3-d plots

          elseif(nd.eq.3) then

!           For element loop over number of subdivisions 'npl_int'

            is = ix_lin(12,nn)
            js = ix_lin(13,nn)
            ks = ix_lin(14,nn)

!           Construct 8-node 'linear elements'

            do ns = 1,8

              sg3(1,1) = gp(1,is+iqx(ns))
              sg3(2,1) = gp(2,js+iqy(ns))
              sg3(3,1) = gp(3,ks+iqz(ns))
              sg3(4,1) = 1.0d0

              call interp3d(1, xl, ndm,nel)

              d_lin = ix_lin(ns,nn)

!             Project dependent variables

              d_lin = ix_lin(ns,nn)

!             Project dependent variables

              u_lin(1:ndf,d_lin) = 0.0d0
              do nm = 1,nel
                if(ix(nm,ne).gt.0) then
                  u_lin(1:ndf,d_lin) = u_lin(1:ndf,d_lin)
     &                               + shp3(4,nm,1)*u(1:ndf,ix(nm,ne))
                endif
              end do ! nm

!             Project stress variables

!              if(tsplfl .or. hsplfl) then

!                call tsp_lin3d(n_el,nd,sg3,  mr(np(295)),st(1,1),
!     &                         mr(np(297)),hr(np(289)),mr(np(290)),
!     &                         mr(np(291)),mr(np(292)),hr(np(264)),
!     &                         s_lin(1,1),nel)

!             NURBS block

!              else

                do i_lin = jp(d_lin-1)+1,jp(d_lin)
                  if(jmat(i_lin).eq.ma) then
                    s_lin(i_lin,1) = 1.0d0
                    do i = 2,npstr
                      s_lin(i_lin,i) = 0.0d0
                    end do ! i
                    do nm = 1,nel
                      nod = ix(nm,ne)
                      if(nod.gt.0) then
                        do ii = ip(nod-1)+1,ip(nod)
!               write(*,*) 'NOD;',nod,ma,imat(ii)
                          if(ma.eq.imat(ii)) then
                            do i = 2,npstr
                              s_lin(i_lin,i) = s_lin(i_lin,i)
     &                                       + shp3(4,nm,1)*st(ii,i)
                            end do ! i
                          endif ! ma check
                        end do ! ii
                      endif ! nod > 0
                    end do ! nm

                    if(histpltfl) then
                      call hp_lin(ix(1,ne),shp3,hr(np(305)),4,nel,
     &                            hr(np(306)),d_lin,
     &                            ip,imat, jp,jmat)
                    endif
                  endif ! jmat(i_lin)
                end do ! i_lin

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

!              endif ! tsplfl .or. hsplfl

            end do ! ns

          endif

        endif ! eltyp

      end do ! ne

      call pltstr(s_lin,p_lin(1,2),s_lin(1,2),jp(nd_lin),ndm,.false.)

      end

