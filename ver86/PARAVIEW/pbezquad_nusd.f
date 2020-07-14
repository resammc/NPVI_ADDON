!$Id:$
      subroutine pbezquad_nusd()

!      for F E A P -- A Finite Element Analysis Program
!
!      coded by:
!                B.Sc. Henning Venghaus
!                Institute of Mechanics
!                Otto von Guericke University
!                Spring 2017

!-----[--.----+----.----+----.-----------------------------------------]
!     Modification log                                Date (dd/mm/year)
!       Original version                                    02/03/2017
!       1. Update to ver8.5                                 03/07/2017
!       2. added calls to new routines                      02/03/2018
!-----[--.----+----.----+----.-----------------------------------------]
!      Purpose: Construct Bezier mesh and Quadratic meshes for plots
!               works like pbez_lin.f but for quad. elements

!      Inputs:

!      Outputs:
!         none   - Users are responsible for generating outputs
!                  through common blocks, etc.  See programmer
!                  manual for example.
!-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'cdata.h'
      include   'cnurb.h'
      include   'qudshp.h'
      include   'sdata.h'

      include   'pointer.h'
      include   'comblk.h'

      include   'npvi.h'

      logical    setvar,palloc
      integer    nd_pre, ne_pre, nd_order

!     Set quadratic mesh for plot use

!     nd_order
!     linear = 2 | quadratic = 3 | cubic = 4 | quatic = 5
      nd_order = 3

      ne_quad = numel*product(npl_int_nusd)
      ne_pre = ne_quad
      nd_pre = (nd_order**ndm)*ne_quad

!     S E R E N D I P I T Y   E L E M E N T (element for vtu file)
      if (serefl) then
        setvar = palloc(276, 'I_LIN', 31*ne_quad , 1)    ! IX_LIN
        setvar = palloc(277, 'X_LIN', ndm*nd_pre, 2)     ! X_LIN
        call pzeroi(mr(np(276)),31*ne_quad)              ! Zero IX_LIN
        call pzero (hr(np(277)),ndm*nd_pre)              ! Zero X_LIN
        call pmsh_quad_nusd_new(mr(np(32)) , mr(np(33)) , hr(np(43)),
     &                 hr(np(44)) , hr(np(263)), hr(np(264)),
     &                 mr(np(276)), hr(np(277)),nd_pre,ne_pre )
        if(ne_quad.gt.ne_pre .or. nd_quad.gt.nd_pre) then
          write(*,*) ' Element Storage =',ne_quad,ne_pre
          write(*,*) ' Coord.  Storage =',nd_quad,nd_pre
        endif
        setvar = palloc(276, 'I_LIN',  31*ne_quad, 1)     ! IX_LIN
        setvar = palloc(277, 'X_LIN', ndm*nd_quad, 2)     ! X_LIN

      else
!       L A G R A N G N E   E L E M E N T (for FE output)
        setvar = palloc(276, 'I_LIN', 38*ne_quad , 1)    ! IX_LIN
        setvar = palloc(277, 'X_LIN', ndm*nd_pre, 2)     ! X_LIN
        call pzeroi(mr(np(276)),38*ne_quad)              ! Zero IX_LIN
        call pzero (hr(np(277)),ndm*nd_pre)              ! Zero X_LIN
        call pmsh_quad_lgrn_nusd(mr(np(32)),mr(np(33)),hr(np(43)),
     &                 hr(np(44)),hr(np(263)),hr(np(264)),
     &                 mr(np(276)),hr(np(277)),nd_pre,ne_pre )
        if(ne_quad.gt.ne_pre .or. nd_quad.gt.nd_pre) then
          write(*,*) ' Element Storage =',ne_quad,ne_pre
          write(*,*) ' Coord.  Storage =',nd_quad,nd_pre
        endif
        setvar = palloc(276, 'I_LIN',  38*ne_quad, 1)    ! IX_LIN
        setvar = palloc(277, 'X_LIN', ndm*nd_quad, 2)    ! X_LIN
      end if !serefl

      end subroutine pbezquad_nusd
