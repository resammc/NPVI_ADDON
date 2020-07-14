!$Id:$
      subroutine psetbez_npvi(xb, wb, xbez, wbez, ixbez, nll, numpbez,
     &                        nd_bez, nmat, ne_bez, ma, nenb, mergfl)

!      * * F E A P * * A Finite Element Analysis Program

!....  Copyright (c) 1984-2018: Regents of the University of California
!                               All rights reserved

!-----[--.----+----.----+----.-----------------------------------------]

!          Adapted for NPVI by Resam Makvandi (resam.makvandi@ovgu.de)
!                              Chair of Computational Mechanics
!                              Institute of Mechanics
!                              Otto von Guericke University Magdeburg
!                              Germany

!-----[--.----+----.----+----.-----------------------------------------]
!     Modification log                                Date (dd/mm/year)
!       Original version                                    01/11/2006
!       Customized for NPVI (originally psetbez.f)          30/05/2020
!-----[--.----+----.----+----.-----------------------------------------]
!      Purpose: Compute Bezier mesh connections from T-spline inputs

!      Inputs:
!         xb(ndm,*)  -
!         wb(*)      -
!         nd_bez     -
!         ne_bez     = Number of element

!      Working
!         nmat(*)    - Mark material number for nodes

!      Outputs:
!         xbez(ndm,*)-
!         wbez(*)    -
!         ixbez(*,*) - Bezier connection data
!-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'cdata.h'
      include   'cnurb.h'
      include   'iofile.h'
      include   'sdata.h'
      include   'setups.h'

      integer    ixbez(nenb,*), nmat(*), nll, nenb, numpbez
      real*8     xb(ndm,*), wb(*), xbez(ndm,*), wbez(*)

      logical    pflag, mergfl
      integer    ma, nd_bez, nn, i,ne_bez
      real*8     ptol, xb3

      save

      data       ptol / 1.d-6 /

      if(ndm.eq.2) then
        xb3   = 0.0d0
      endif

      do i = 1,nll
        pflag = .true.
        if (mergfl) then
          do nn = 1,numpbez
           ! if(nmat(nn).eq.ma) then
              if(ndm.eq.3) then
                xb3 = abs(xb(3,i) - xbez(3,nn))
              endif
              if(max(abs(xb(1,i)-xbez(1,nn)),
     &               abs(xb(2,i)-xbez(2,nn)),xb3) .lt. ptol) then
                pflag = .false.
                exit
              endif
           ! endif
          end do ! nn
        end if ! mergfl

        if(pflag) then
            numpbez = numpbez + 1

          if(numpbez.gt.nd_bez) then
            if(rank.eq.0) then
              write(*,*) ' Number of computed Bezier points too large'
              write(*,*) ' NUMPBEZ = ',numpbez,' ND_BEZ =',nd_bez,
     &                   ' ELEMENT = ',ne_bez
            endif
            write(iow,*) ' Number of computed Bezier points too large'
            write(iow,*) ' NUMPBEZ = ',numpbez,' ND_BEZ =',nd_bez,
     &                   ' ELEMENT = ',ne_bez
            write(ilg,*) ' Number of computed Bezier points too large'
            write(ilg,*) ' NUMPBEZ = ',numpbez,' ND_BEZ =',nd_bez,
     &                   ' ELEMENT = ',ne_bez
            call plstop(.true.)
          endif
          do nn = 1,ndm
            xbez(nn,numpbez) = xb(nn,i)
          end do ! nn
          wbez(numpbez) = wb(i)
          ixbez(i,ne_bez)   = numpbez
          nmat(numpbez)      = ma
        else
          ixbez(i,ne_bez)   = nn
        endif
      end do ! i
      ixbez(nll+1,ne_bez)  = ma

      end
