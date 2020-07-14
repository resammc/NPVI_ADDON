!$Id:$
      subroutine FEoutput(x,ix,ndm,nen1,nen,numnp,numel,quadfl)

!      for F E A P -- A Finite Element Analysis Program
!
!      modified by:
!                B.Sc. Henning Venghaus
!                Institute of Mechanics
!                Otto von Guericke University
!                Spring 2017

!-----[--.----+----.----+----.-----------------------------------------]
!     Modification log                                Date (dd/mm/year)
!       Original version                                    02/03/2017
!       1. Modify for ver8.5; delete timefl argument        03/07/2017
!       2. Fixed a bug in models with multiple 
!          material cards                                   23/06/2020
!       3. Fixed a bug in writing out the conn. array       23/06/2020
!-----[--.----+----.----+----.-----------------------------------------]
!      Purpose: write FE input file from nurbfeap input file

!      Inputs:
!        x(ndm,*)     - Mesh control points
!        ix(nen1,*)   - Mesh connection list
!        ndm          - number of dimensions
!        nen1         - size of ix array
!        nen          - number of nodes per element
!        numnp        - number of nodes
!        numel        - number of elements
!        quadfl       - true if elements are quadratic

!      Outputs:
!        to file
!-----[--+---------+---------+---------+---------+---------+---------+-]

      implicit none

      include 'iofile.h'
      include 'npvi.h'

      character  name*8
      integer    ndm,nen,nen1,numnp,numel,numma, i,n
      integer    ix(nen1,numel)
      real*8     x(ndm,numnp)
      logical    quadfl


      if (quadfl .and. serefl) then
        nen = (2+ndm)*2**(ndm-1)
      else if (quadfl .and. (.not.serefl)) then
        nen = 3**ndm
      else
        nen = (2**ndm)
      end if !(quadfl)
!     Determine maximum number of materials

      numma = 0
      do n = 1,numel
        numma = max(ix(nen1,n),numma)
      end do ! n

      write(*,*) ndm,nen,nen1,numnp,numel,numma

!     Set output filename
      name = 'FEoutput'

      if(ior.lt.0) then
        write(*,3000) ' Output FE-element file: ',name
      endif
      write(iow,3000) ' Output FE-element file: ',name

!     Open file and write control information

!       open(unit=3,file = name,status = 'unknown')
      open(unit=3,file = name,access='sequential')
      write(3,2000) 'feap * * File=',name,numnp,numel,numma,ndm,ndm,nen

!     Output a simple material
      do i = 1, numma
       write(3,2006) '    '
       write(3,2007) 'MATE,',i
       write(3,2006) '  SOLI'
       write(3,2006) '    PLAN STRA'
       write(3,2006) '    ELAS ISOT 1.e5 0.25'
      end do

!     Output coordinates

      write(3,2001) 'COORdinates all'
      do n = 1,numnp
        write(3,2002) n,0,(x(i,n),i=1,ndm)
      end do ! n

!     Output elements

      write(3,2001) 'ELEMents all'
      do n = 1,numel
        if (ndm.eq.3 .and. quadfl) then
        write(3,2003) n,0,ix(nen1,n),(ix(i,n),i=1,13)
        write(3,2003) (ix(i,n),i=14,nen)
      else
        write(3,2003) n,0,ix(nen1,n),(ix(i,n),i=1,nen)
      end if !(ndm.eq.3 .and. quadfl)
      end do ! n

!     Write end mesh

      write(3,2001) 'END MESH'
      write(3,2001) 'INTER'
      write(3,2001) 'STOP'


      close(unit=3,status='KEEP')
!     Formats

2000  format(a,a/6i8)
2001  format(/a)
2002  format(2i8,1p,3e16.7)
2003  format(16i8)
2006  format(a)
2007  format(a,2i4)
3000  format(a,a/1x)

      end subroutine FEoutput
