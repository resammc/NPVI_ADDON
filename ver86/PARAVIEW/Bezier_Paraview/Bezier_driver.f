!$Id:$
      subroutine bezier_driver(timefl, strefl, velofl, accefl, histfl,
     &                         mergfl, helpfl, dbglvl)

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
!      Purpose: Post-process IGA data using Bezier cells
!-----[--.----+----.----+----.-----------------------------------------]

      implicit  none
      
      include 'sdata.h'   ! ndm, ndf
      include 'cdata.h'   ! numel, nen
      include 'pointer.h' ! num_nps
      include 'comblk.h'  ! hr, mr
      include 'fdata.h'   ! fl
      include 'strnum.h'  ! numnm
      include 'eldatp.h'  ! hplmax

      integer    :: dbglvl
      integer    :: numpbez, nenb
      logical    :: timefl, strefl, velofl, accefl, histfl, helpfl
      logical    :: mergfl
      logical    :: setvar, npvi_alloc, ualloc
      logical    :: fl11_backup

      if (dbglvl.ge.1) write(*,4400) 'Started creating the Bezier mesh'

      ! initialize numpbez (stores the total number of Bezier cps)
      numpbez = 0

      nenb = nen + 5 ! nen+1: ma, nen+2: p, nen+3: q, nen+4: r, nen+5 : nel
      
      ! allocate storage for Bezier nodes and elements
      setvar = ualloc(11,'UTEM1',numel*nen*ndm , 2)    ! X_BEZ
      setvar = ualloc(12,'UTEM2',numel*nen     , 2)    ! W_BEZ
      setvar = ualloc(13,'UTEM3',numel*(nen+5) , 1)    ! IXBEZ
      setvar = ualloc(15,'UTEM5',numel*nen     , 1)    ! Temporary

      ! activate projection
      fl11_backup = fl(11)
      fl(11) = .true.

      ! create the Bezier mesh
      call pmsh_bezier( mr(np(32))  ,   ! IE    (FEAP)
     &                  mr(np(33))  ,   ! IX    (FEAP)
     &                  hr(np(43))  ,   ! X     (FEAP)
     &                  hr(np(44))  ,   ! XL    (FEAP)
     &                  hr(np(263)) ,   ! WT    (FEAP)
     &                  hr(np(264)) ,   ! WTL   (FEAP)
     &                  mr(np(308)) ,   ! LKNOT (FEAP)
     &                  hr(up(11))  ,   ! X_BEZ (USER)
     &                  hr(up(12))  ,   ! W_BEZ (USER)
     &                  mr(up(13))  ,   ! IXBEZ (USER)
     &                  mr(up(15))  ,   ! NMAT  (USER)
     &                  nenb        ,   ! like nen1 for bezier
     &                  numpbez     ,   ! number of Bezier cps
     &                  mergfl      )   ! merging flag

      ! deallocate the temporary array
      setvar = ualloc(15,'UTEM5',    0         , 1)    ! Temporary

      ! Set actual storage for Bezier mesh data
      setvar = ualloc(11,'UTEM1',numpbez*ndm   , 2)    ! X_BEZ
      setvar = ualloc(12,'UTEM2',numpbez       , 2)    ! W_BEZ
      setvar = ualloc(13,'UTEM3',numel*nenb    , 1)    ! IXBEZ

      ! calculate displacements, velocities and accelerations
      ! on the Bezier mesh

      if (dbglvl.ge.1) write(*,4400) 'Started calculating DOFs'

      ! allocate storage for displacements, velocities and accelerations
      ! on the Bezier cps
        setvar = ualloc(14,'UTEM4',numpbez*ndf , 2)    ! U_BEZ
      if (velofl) then
        setvar = ualloc(16,'UTEM6',numpbez*ndf , 2)    ! V_BEZ
      end if
      if (accefl) then
        setvar = ualloc(17,'UTEM7',numpbez*ndf , 2)    ! A_BEZ
      end if

      call pdis_bezier( mr(np(32))      ,   ! IE    (FEAP)
     &                  mr(np(33))      ,   ! IX    (FEAP)
     &                  mr(np(308))     ,   ! LKNOT (FEAP)
     &                  hr(np(40))      ,   ! U     (FEAP)
     &                  hr(np(42))      ,   ! V     (FEAP)
     &                  hr(np(42))+nneq ,   ! A     (FEAP)
     &                  mr(up(13))      ,   ! IXBEZ (USER)
     &                  hr(up(14))      ,   ! U_BEZ (USER)
     &                  hr(up(16))      ,   ! V_BEZ (USER)
     &                  hr(up(17))      ,   ! A_BEZ (USER)
     &                  velofl          ,
     &                  accefl          ,
     &                  nenb            )

      ! calculate stresses, strains, ... on the Bezier mesh
      if (strefl) then
        if (dbglvl.ge.1) write(*,4400) 'Start projection of STRE'
        if (dbglvl.ge.1) write(*,4400) 'Start projection of PSTR'

        ! prepare program and user arrays
        call pstr_bezier_prepare(nenb,numpbez)

        ! project the stresses on the Bezier mesh
        call pstr_bezier( mr(np(32))       ,   ! IE    (FEAP)
     &                    mr(np(33))       ,   ! IX    (FEAP)
     &                    mr(np(308))      ,   ! LKNOT (FEAP)
     &                    mr(up(13))       ,   ! IXBEZ (USER)
     &                    hr(np(58)+numnm) ,   ! ST    (USER)
     &                    hr(np(57))       ,   ! PS    (USER)
     &                    hr(up(18))       ,   ! S_BEZ (USER)
     &                    hr(up(19))       ,   ! P_BEZ (USER)
     &                    mr(np(339))      ,
     &                    nenb             ,
     &                    numpbez          )

        if (dbglvl.ge.1) write(*,4401) 'Finished projection of STRE'
        if (dbglvl.ge.1) write(*,4401) 'Finished projection of PSTR'

      end if ! strefl

      ! Allocate array to store element offset pointers
      setvar = ualloc(15,'UTEM5',numel+1,1)

              call uparaview_bez(
     &              hr(up(11)),      !x(ndm,numnp)     - Node coordinat
     &              hr(up(12)),      !w(numnp)   - weights
     &              mr(up(13)),      !ix(nen1,numel)   - Element connec
     &              hr(up(14)),      !u(ndf,numnp)     - Node displacem
     &              hr(up(16)),      !ud(ndf,numel)  - Node velocity
     &              hr(up(17)),      !ud(ndf,numel)  - Node acceleration
     &              hr(up(18)),      !sig(numnp,*)     - Node stre
     &              hr(up(19)),      !psig(numnp,*)    - Node prin
     &              hr(np(305)),     !his(numnp,*)     - History variab
     &              mr(up(15)),      !el_ptr(numel+1)  - Pointer to ele
     &              mr(np(341)),     !ip(0:numnp)      - Pointer to sig
     &              mr(np(32)),      !ie(nie,*)
     &              ndm,             !ndm              - Mesh dimension
     &              ndf,             !ndf              - DOF's/node
     &              nenb,            !nen1             - Dimension of i
     &              nen,             !nen              - Max nodes/elem
     &              hplmax,          !hplmax           - Number of hist
     &              numpbez,         !nd_export        - Number of node
     &              numel,           !ne_export        - Number of elem
     &              timefl,          !                 - Time Flag
     &              strefl,          !                 - Stress Flag
     &              velofl,          !                 - Velocity Flag
     &              accefl,          !                 - Acceleration F
     &              histfl,          !                 - History Flag
     &              dbglvl)          !                 - debug level

      ! restore fl(11) flag
      fl(11) = fl11_backup

      ! deallocate all the storage arrays
      setvar = ualloc(11,'UTEM1',    0         , 2)
      setvar = ualloc(12,'UTEM2',    0         , 2)
      setvar = ualloc(13,'UTEM3',    0         , 1)
      setvar = ualloc(14,'UTEM4',    0         , 2)
      if (velofl) then
        setvar = ualloc(16,'UTEM6',  0         , 2)
      end if
      if (accefl) then
        setvar = ualloc(17,'UTEM7',  0         , 2)
      end if

      if (strefl) then
        setvar = ualloc(18,'UTEM8',  0         , 2)
        if (histpltfl) then
          setvar = ualloc(20,'UTEM0',0         , 2)
        end if
      end if
      setvar = ualloc(15,'UTEM5',    0         , 1)

4400  format('[....] ',a)
4401  format('[ OK ] ',a)
4402  format('[info] ',a,L2)
4403  format('[info] ',a,i7)
4404  format('[warn] ',a)

      end
      
