      subroutine FEoutput(x,ix,ndm,nen1,nen,numnp,numel,quadfl,timefl)
      
c      for F E A P -- A Finite Element Analysis Program
c
c      modified by:
c                B.Sc. Henning Venghaus
c                Institute of Mechanics
c                Otto von Guericke University
c                Spring 2017

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    02/03/2017
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: write FE input file from nurbfeap input file

c      Inputs:
c        x(ndm,*)     - Mesh control points
c        ix(nen1,*)   - Mesh connection list
c        ndm          - number of dimensions
c        nen1         - size of ix array
c        nen          - number of nodes per element
c        numnp        - number of nodes
c        numel        - number of elements
c        quadfl       - true if elements are quadratic
c        timefl       - true if calculation is transient

c      Outputs:
c        to file
c-----[--+---------+---------+---------+---------+---------+---------+-]
      
      implicit none
      
      include 'iofile.h'
      include 'npvi.h'
      
      character  name*8
      integer    ndm,nen,nen1,numnp,numel,numma, i,n
      integer    ix(nen1,numel)
      real*8     x(ndm,numnp)
      logical    quadfl,timefl
      
      
      if (quadfl .and. serefl) then
        nen = (2+ndm)*2**(ndm-1)
      else if (quadfl .and. (.not.serefl)) then
        nen = 3**ndm
      else
        nen = (2**ndm)
      end if !(quadfl)
c     Determine maximum number of materials

      numma = 0
      do n = 1,numel
        numma = max(ix(nen1,n),numma)
      end do ! n

      write(*,*) ndm,nen,nen1,numnp,numel,numma
      
c     Set output filename
      name = 'FEoutput'

      if(ior.lt.0) then
        write(*,3000) ' Output FE-element file: ',name
      endif
      write(iow,3000) ' Output FE-element file: ',name

c     Open file and write control information

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
      
c     Output coordinates

      write(3,2001) 'COORdinates all'
      do n = 1,numnp
        write(3,2002) n,0,(x(i,n),i=1,ndm)
      end do ! n

c     Output elements

      write(3,2001) 'ELEMents all'
      do n = 1,numel
        if (ndm.eq.3 .and. quadfl) then
        write(3,2003) n,0,ix(nen1,n),(ix(i,n),i=1,13)
        write(3,2003) (ix(i,n),i=14,nen)
      else
        write(3,2003) n,0,ix(nen1,n),(ix(i,n),i=1,nen)
      end if !(ndm.eq.3 .and. quadfl)
      end do ! n

c     Write end mesh

      write(3,2001) 'END MESH'
      write(3,2001) 'INTER'
      write(3,2001) 'STOP'

      
      close(unit=3,status='KEEP')
c     Formats

2000  format(a,a/6i8)
2001  format(/a)
2002  format(2i8,1p,3e16.7)
2003  format(16i8)
2005  format(11i8)
2006  format(a)
2007  format(a,2i4)
3000  format(a,a/1x)

      end
