      subroutine pbezquad_nusd()

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
c      Purpose: Construct Bezier mesh and Quadratic meshes for plots
c               works like pbez_lin.f but for quad. elements

c      Inputs:

c      Outputs:
c         none   - Users are responsible for generating outputs
c                  through common blocks, etc.  See programmer
c                  manual for example.
c-----[--.----+----.----+----.-----------------------------------------]
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

c     Set quadratic mesh for plot use
      
c     nd_order
c     linear = 2 | quadratic = 3 | cubic = 4 | quatic = 5      
      nd_order = 3
           
      ne_quad = numel*product(npl_int_nusd)
      ne_pre = ne_quad
      nd_pre = (nd_order**ndm)*ne_quad

      if (serefl) then
!       S E R E N D I P I T Y   E L E M E N T (common element for vtu file)
        setvar = palloc(276, 'I_LIN', 31*ne_quad , 1)     ! IX_LIN
        setvar = palloc(277, 'X_LIN', ndm*nd_pre, 2)     ! X_LIN
        call pzeroi(mr(np(276)),31*ne_quad)               ! Zero IX_LIN
        call pzero (hr(np(277)),ndm*nd_pre)              ! Zero X_LIN
        call pmsh_quad_nusd(mr(np(32)) , mr(np(33)) , hr(np(43)),
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
        setvar = palloc(276, 'I_LIN', 38*ne_quad , 1)     ! IX_LIN
        setvar = palloc(277, 'X_LIN', ndm*nd_pre, 2)     ! X_LIN
        call pzeroi(mr(np(276)),38*ne_quad)               ! Zero IX_LIN
        call pzero (hr(np(277)),ndm*nd_pre)              ! Zero X_LIN
        call pmsh_quad_lgrn_nusd(mr(np(32)),mr(np(33)),hr(np(43)),
     &                 hr(np(44)),hr(np(263)),hr(np(264)),
     &                 mr(np(276)),hr(np(277)),nd_pre,ne_pre )
        if(ne_quad.gt.ne_pre .or. nd_quad.gt.nd_pre) then
          write(*,*) ' Element Storage =',ne_quad,ne_pre
          write(*,*) ' Coord.  Storage =',nd_quad,nd_pre
        endif
        setvar = palloc(276, 'I_LIN',  38*ne_quad, 1)     ! IX_LIN
        setvar = palloc(277, 'X_LIN', ndm*nd_quad, 2)     ! X_LIN
      end if !serefl
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      end