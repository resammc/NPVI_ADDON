!$Id:$
      subroutine pstr_bezier_prepare(nenb,numpbez)

!      * * F E A P * * A Finite Element Analysis Program

!....  Copyright (c) 1984-2018: Regents of the University of California
!                               All rights reserved

!-----[--.----+----.----+----.-----------------------------------------]

!      Written for FEAP 8.5 by Resam Makvandi (resam.makvandi@ovgu.de)
!                              Chair of Computational Mechanics
!                              Institute of Mechanics
!                              Otto von Guericke University Magdeburg
!                              Germany

!-----[--.----+----.----+----.-----------------------------------------]
!     Modification log                                Date (dd/mm/year)
!       Original version                                    30/05/2020
!-----[--.----+----.----+----.-----------------------------------------]
!      Purpose: Project stresses
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
      include  'eldata.h'

      include  'p_int.h'

      logical  :: setvar, npvi_alloc, palloc, ualloc
      integer  :: nd_pre, ne_pre
      integer  :: i, ii
      integer  :: nenb, numpbez, numni_bez

      save

!     Allocate storage for element variable projections
      setvar = ualloc(15,'UTEM5',2*numpbez,1)

      call proj_matc(mr(up(13)),nen,nenb,numel,mr(up(15)),
     &               numpbez,numni_bez)

      ! setvar = npvi_alloc(109,'NTEM9',numel*nen+1,1)  ! 341 : MIGPT
      ! setvar = npvi_alloc(110,'NTE10',numni  ,1)      ! 342 : MIGVA  
      setvar = palloc(341,'MIGPT',numpbez+1,1)  ! 341 : MIGPT
      setvar = palloc(342,'MIGVA',numni_bez  ,1)      ! 342 : MIGVA  

      call proj_mats(mr(up(13)),nen,nenb,numel, mr(up(15)),
     &               mr(np(341)),numpbez, mr(np(342)))

      setvar = ualloc(15,'UTEM5',0,1)
      
      setvar = ualloc(15,'UTEM5',2*numnp,1)
      call proj_matc(mr(np(33)),nen,nen1,numel,mr(up(15)),numnp,numnm)

      ! setvar = npvi_alloc(111,'NTE11',numnp+1,1) ! 339 : MAPTR
      ! setvar = npvi_alloc(112,'NTE12',numnp  ,1) ! 340 : MAVAL
      setvar = palloc(339,'MAPTR',numnp+1,1) ! 339 : MAPTR
      setvar = palloc(340,'MAVAL',numnm  ,1) ! 340 : MAVAL

      call proj_mats(mr(np(33)),nen,nen1,numel,mr(up(15)),
     &               mr(np(339)),numnp,mr(np(340)))
      !mr(up(112)) = 1

      setvar = ualloc(15,'UTEM5',0,1)
      
	    if(plfl) then
        !if(tsplfl .or. hsplfl) then
        !  setvar = palloc( 57,'NDER',numpln*8    ,2)
        !  setvar = palloc( 58,'NDNP',numpln*npstr,2)
        !  call pzero (hr(np(57)),numpln*8)
        !  call pzero (hr(np(58)),numpln*npstr)
        !  if(histpltfl) then
        !    setvar = palloc(305,'HDNP ',numpln*hplmax,2)
        !  endif
        !else
          setvar = palloc(57,'NDER',numnm*8    ,2) ! 57 : NDER
          setvar = palloc(58,'NDNP',numnm*npstr,2) ! 58 : NDNP
          !write(*,*) "npstr",npstr
          if(histpltfl) then
            ! setvar = npvi_alloc(115,'NTE15 ',numnm*hplmax,2) ! 305 : HDNP
            setvar = palloc(305,'HDNP ',numnm*hplmax,2) ! 305 : HDNP
          endif
        !endif
        setvar = palloc(60,'NDNS',max(nen*npstr,nst*nst),2) ! 60 : NDNS
        setvar = palloc(207,'NSCR',numel      ,2) ! 207 : NSCR
        nper   = np(57)
        npnp   = np(58)
        !plfl   =.false.
        if(histpltfl) then
          ! setvar = npvi_alloc(118,'NTE18',nen*hplmax  ,2) ! 304 : HSELM
          setvar = palloc(304,'HSELM',nen*hplmax  ,2) ! 304 : HSELM
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

      call pzero(hr(npnp), npstr*numnm)
      call pzero(hr(nper),     8*numnm)

      if(histpltfl) then
        call pzero(hr(np(305)),numnm*hplmax)
      endif
      call pzero(hr(np(207)),numel)

!     Compute element stress values for projections
      !setvar = ualloc(-num_nps,'     ', 1, 2) ! A_BEZ
      fp(5)  = np(36)
      np(36) = np(60)
      pltmfl = .true.
      call formfe(np(40),np(26),np(26),np(26),
     &           .false.,.false.,.false.,.false.,8,1,numel,1)
      pltmfl = .false.
      np(36) = fp(5)
      !setvar = npvi_alloc(-num_nps,'    ', 1, 2) ! A_BEZ
!     Compute nodal stress values by average for projections
!      if(tsplfl .or. hsplfl) then
!        call pltstr(hr(npnp),hr(nper+numnm),hr(npnp+numnm),
!     &              numnm,ndm,.true.)
!      else
      call pltstr(hr(npnp),hr(nper+numnm),hr(npnp+numnm),
     &              numnm,ndm,.true.)
!      endif

!     Reset element material status

      ii = 0
      do i = nen1-1,nen1*numel-1,nen1
        mr(np(33)+i) = mr(np(323)+ii)
        ii           = ii + 1
      end do ! i
      setvar = palloc(323,'PTEM3',0,1)

!     If NURBS or T-spline set up plot region & convert to quads

      if(nurbfl) then
       ! if(np(119).eq.0) then
          setvar = ualloc(18, 'UTEM8', numpbez*npstr, 2) ! s_bez
          setvar = ualloc(19, 'UTEM9', numpbez*8    , 2) ! p_bez
       !   if (histpltfl) then
       !     setvar = npvi_alloc(121, 'NTE21', numpbez*hplmax, 2) ! h_bez
       !   endif
       !endif
      endif

      end