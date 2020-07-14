      subroutine pbezlin_nusd()

c      for F E A P -- A Finite Element Analysis Program
c
c      coded by:
c                B.Sc. Henning Venghaus
c                Institute of Mechanics
c                Otto von Guericke University
c                Spring 2017

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    26/03/2017
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Construct Bezier mesh and Linear meshes with the
c               opportunity of anisotropic nurb element subdivision

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
      integer    nd_pre, ne_pre

c     Set linear mesh for plot use

      ne_lin = numel*product(npl_int_nusd)
      ne_pre = ne_lin
      nd_pre = max(2*ndm*ne_lin,4*(ndm-1)*ne_lin)
      setvar = palloc(276, 'I_LIN',19*ne_lin  , 1)     ! IX_LIN
      setvar = palloc(277, 'X_LIN', ndm*nd_pre, 2)     ! X_LIN
      call pzeroi(mr(np(276)),19*ne_lin)               ! Zero IX_LIN
      call pzero (hr(np(277)),ndm*nd_pre)              ! Zero X_LIN
      call pmsh_lin_nusd(mr(np(32)) , mr(np(33)) , hr(np(43)),
     &                   hr(np(44)) , hr(np(263)), hr(np(264)),
     &                   mr(np(276)), hr(np(277)),nd_pre,ne_pre )
      if(ne_lin.gt.ne_pre .or. nd_lin.gt.nd_pre) then
        write(*,*) ' Element Storage =',ne_lin,ne_pre
        write(*,*) ' Coord.  Storage =',nd_lin,nd_pre
      endif
      setvar = palloc(276, 'I_LIN',  19*ne_lin, 1)     ! IX_LIN
      setvar = palloc(277, 'X_LIN', ndm*nd_lin, 2)     ! X_LIN

      end
