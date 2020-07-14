!$Id:$
      subroutine pjstrs_lin(trifl)

!      * * F E A P * * A Finite Element Analysis Program

!....  Copyright (c) 1984-2017: Regents of the University of California
!                               All rights reserved

!-----[--.----+----.----+----.-----------------------------------------]
!     Modification log                                Date (dd/mm/year)
!       Original version                                    01/11/2006
!       1. Set 'pltmfl' to true for reactions               22/03/2009
!       2. Set history projection array to zero             11/01/2012
!       3. ADD PTEM3 to store active material numbers       10/07/2016
!          Change numnp to numnm for projection quantities.
!       4. Set allocation of projection of stress variables 20/04/2017
!          Moved from pplotf.f.  Add 'cnurb.h'
!       5. Added calls for npvi and nusd functionality      02/03/2018
!-----[--.----+----.----+----.-----------------------------------------]
!      Purpose: Project nodal stresses

!      Inputs:
!         trifl      - Flag, generate element size for tri2d if true

!      Outputs:
!         none       - Output stored in pointers to arrays
!-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'cnurb.h'
      include  'comblk.h'
      include  'elcapt.h'
      include  'eldatp.h'
      include  'hdatam.h'
      include  'qudshp.h'
      include  'pdata3.h'
      include  'pointer.h'
      include  'prstrs.h'
      include  'sdata.h'
      include  'strnum.h'
      include  'fdata.h'

      include  'p_int.h'

      logical   trifl, setvar, palloc
      integer    nd_pre, ne_pre
      integer   i,ii

      save

!     Allocate storage for element variable projections
	if(plfl) then
        if(tsplfl .or. hsplfl) then
          setvar = palloc( 57,'NDER',numpln*8    ,2)
          setvar = palloc( 58,'NDNP',numpln*npstr,2)
          call pzero (hr(np(57)),numpln*8)
          call pzero (hr(np(58)),numpln*npstr)
          if(histpltfl) then
            setvar = palloc(305,'HDNP ',numpln*hplmax,2)
          endif
        else
          setvar = palloc( 57,'NDER',numnm*8    ,2)
          setvar = palloc( 58,'NDNP',numnm*npstr,2)
          if(histpltfl) then
            setvar = palloc(305,'HDNP ',numnm*hplmax,2)
          endif
        endif
        setvar = palloc( 60,'NDNS',max(nen*npstr,nst*nst),2)
        setvar = palloc(207,'NSCR',numel      ,2)
        nper   = np(57)
        npnp   = np(58)
        plfl   =.false.
        if(histpltfl) then
          setvar = palloc(304,'HSELM',nen*hplmax  ,2)
        endif
	endif !plfl

      ner   = nper
      nph   = npnp

!     Clear caption array

      do i = 1,50
        ecapt(i) = '  '
      end do ! i

!     Stress projections

      istv = npstr - 1

!     Save element material status
      call pset_projs_nusd(mr(np(32)),mr(np(33)),mr(np(128)))

c     setvar = palloc(323,'PTEM3',numel,1)
c     ii = 0
c     do i = nen1-1,nen1*numel-1,nen1
c       mr(np(323)+ii) =  mr(np(33)+i)
c       write(*,*) 'IX:MA(1):33 =',ii,mr(np(33)+i),nen1,i+1
c       ii = ii + 1
c     end do ! i

c     ii = 0
c     do i = nen1-1,nen1*numel-1,nen1
c       if(mr(np(128)+ii).lt.0) then
c         mr(np(33)+i) = -abs(mr(np(33)+i))
c       else
c         mr(np(33)+i) =  abs(mr(np(33)+i))
c       endif
c       write(*,*) 'IX:MA(2):33:128 =',ii,mr(np(33)+i),mr(np(128)+ii)
c       ii = ii + 1
c     end do ! i

      call pzero(hr(npnp), npstr*numnm)
      call pzero(hr(nper),     8*numnm)
      if(histpltfl) then
        call pzero(hr(np(305)),numnm*hplmax)
      endif
      if(.not.trifl) call pzero(hr(np(207)),numel)

!     Compute element stress values for projections

      fp(1)  = np(36)
      np(36) = np(60)
      pltmfl = .true.
      call formfe(np(40),np(26),np(26),np(26),
     &           .false.,.false.,.false.,.false.,8,1,numel,1)
      pltmfl = .false.
      np(36) = fp(1)

!     Compute nodal stress values by average for projections
      if(tsplfl .or. hsplfl) then
        call pltstr(hr(npnp),hr(nper+numnm),hr(npnp+numnm),
     &              numnm,ndm,.true.)
      else
        call pltstr(hr(npnp),hr(nper+numnm),hr(npnp+numnm),
     &              numnm,ndm,.true.)
      endif

!     Reset element material status

      ii = 0
      do i = nen1-1,nen1*numel-1,nen1
        mr(np(33)+i) = mr(np(323)+ii)
        ii           = ii + 1
      end do ! i
      setvar = palloc(323,'PTEM3',0,1)

!     If NURBS or T-spline set up plot region & convert to quads

      if(nurbfl) then
!        if(np(281).eq.0) then
          setvar = palloc(281, 'S_LIN', numni*npstr, 2) ! s_lin
          setvar = palloc(282, 'P_LIN', numni*8    , 2) ! p_lin
          if(histpltfl) then
            setvar = palloc(306, 'H_LIN', numni*hplmax, 2) ! h_lin
          endif
!        endif
        call pstr_lin_nusd_new(mr(np(33)), hr(np(43)), hr(np(44)),
     &                hr(np(263)),hr(np(264)),mr(np(276)),
     &                hr(np(58)), hr(np(281)),hr(np(282)),
     &                hr(np(40)), hr(np(278)),mr(np(339)),
     &                mr(np(340)),mr(np(341)),mr(np(342)))
!       call iprint(mr(np(342)),1,numni,1,'IMAT_S_LIN')
!       call mprint(hr(np(281)),numni,npstr,numni,'S_LIN')
      endif

      end subroutine pjstrs_lin

      subroutine pset_projs_nusd(ie,ix,ia)

      implicit   none

      include   'cdata.h'   ! numel
      include   'cdat1.h'   ! nie
      include   'cnurb.h'   ! npl_int
      include   'sdata.h'   ! nen1
      include   'pointer.h' ! np(*)
      include   'comblk.h'  ! mr(*)

      include   'npvi.h'

      logical          :: setvar, palloc
      integer  :: ie(nie,*)
      integer  :: ix(nen1,*)
      integer  :: ia(*)

      integer  :: i, ii, ma

      save

      setvar = palloc(323,'PTEM3',numel,1)
      do i = 1,numel
        mr(np(323)+i-1) =  ix(nen1,i)
!       write(*,*) 'IX:MA(1):33 =',i,ix(nen1,i)
      end do ! i

      ii = 1
      do i = 1,numel
        ma = ix(nen1,i)
        if(ie(1,ma).le.0) then
          ix(nen1,i) = -abs(ix(nen1,i))
        else
          if(ia(ii).lt.0) then
            ix(nen1,i) = -abs(ix(nen1,i))
          endif
!         write(*,*) 'IX:MA(2):33:128 =',ii,mr(np(33)+i),mr(np(128)+ii)
          ii = ii + product(npl_int_nusd)!npl_int**ie(1,ma)
        endif
      end do ! i

      end subroutine pset_projs_nusd
