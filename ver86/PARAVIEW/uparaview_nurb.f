!$Id:$
      subroutine uparaview_nurb(x,ix,u,ud,sig,psig,his, el_ptr,ip,
     &                          ndm,ndf,nen1,nen,hplmax, numnp,numel,
     &                          timefl, strefl, velofl, accefl, histfl,
     &                          dbglvl,quadfl)

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
!       1. Modified for ver8.5                              03/07/2017
!          use numnm, istv (strnum.h); add ip(0:numnp)
!-----[--.----+----.----+----.-----------------------------------------]
!      Purpose:  Interface to PARAVIEW

!      Inputs:
!         x(ndm,numnp)     - Node coordinates
!         ix(nen1,numel)   - Element connections
!         u(ndf,numnp)     - Node displacements
!         ud(ndf,numel,2)  - Node velocity and acceleration
!         sig(numnm,*)     - Node stresses
!         psig(numnm,*)    - Node principal stresses
!         his(numnm,*)     - History variables
!         ip(0:numnp)      - Pointer for material nodal array
!         ndm              - Mesh dimension
!         ndf              - DOF's/node
!         nen1             - Dimension of ix array
!         nen              - Max nodes/element
!         istv             - Number stresses/node
!         nplmax           - Number of history variables

!      Working array
!         el_ptr(numel+1) - Pointer to elements

!      Outputs:
!         To file with appender *.vtu
!-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'comfil.h'
      include  'counts.h'
      include  'iofile.h'
      include  'pointer.h'
      include  'strnum.h'
      include  'tdata.h'

      character parafile*128,parext*4
      integer   ndm,ndf,nen1,nen,hplmax, numnp,numel

      integer   ix(nen1,*), el_ptr(0:numel)
      real*8    x(ndm,numnp),u(ndf,numnp),ud(ndf,numnp,2)
      real*8    sig(numnm,*),psig(numnm,*),his(numnm,*)
      integer   i,ii,node, plu, dbglvl, etyp

      integer   ip(0:numnp)
      integer   jj,ns
      real*8    sigl(50)

      logical   timefl, strefl, velofl, accefl, histfl, quadfl

      integer   linestart

      save

      data      plu / 99 /

      if (dbglvl.ge.1) write(*,4401) 'VTU file routine called'
      if (dbglvl.ge.2) write(*,4403) 'imported nen1:' ,nen1
      if (dbglvl.ge.2) write(*,4403) 'imported nen: ' ,nen

!      call pltgmv(2)    ! Parameters Set for IgA

!     Set extender name
      parext = 'vtu '

!     Assign name
      parafile = '    '
      if(.not.timefl) then
        parafile(1:13) = 'feap_paraview'
!     Assign name from 'fplt': Allows for multiple file outputs

      else
        i = index(fplt,' ')
        parafile(1:128) = ' '
        parafile(1:i-1)=fplt(1:i-1)
        parafile(i:i+4) = '00000'
        if (nstep.le.9) then
          write(parafile(i+4:i+4),'(i1)') nstep
        elseif (nstep.le.99) then
          write(parafile(i+3:i+4),'(i2)') nstep
        elseif (nstep.le.999) then
          write(parafile(i+2:i+4),'(i3)') nstep
        elseif (nstep.le.9999) then
          write(parafile(i+1:i+4),'(i4)') nstep
        elseif (nstep.le.99999) then
          write(parafile(i:i+4),'(i5)') nstep
        endif
      end if
      call addext(parafile,parext,128,4)

      if (dbglvl.ge.1) write(*,4400) 'Start writing file'
      open(unit=plu,file=parafile,access='sequential')
      if (dbglvl.ge.0) write(*,4404)'Saving PARAVIEW data to "'//
     & trim(parafile)//'"'

!     Write out top header
      write(plu,1000)

      linestart=14
      if (strefl) then
        linestart = linestart + 2
      endif
      if (velofl) then
        linestart = linestart + 1
      endif
      if (accefl) then
        linestart = linestart + 1
      endif

      write(plu,1002) 'step:',nstep
      write(plu,1003) 'time:',ttim
      write(plu,1003) 'dt:  ',dt
      write(plu,1001) 'nodes',linestart
      linestart = linestart + numnp
      linestart = linestart + 4
      write(plu,1001) 'conne',linestart
      linestart = linestart + numel
      linestart = linestart + 2
      write(plu,1001) 'offse',linestart
      linestart = linestart + numel
      linestart = linestart + 2
      write(plu,1001) 'types',linestart
      linestart = linestart + numel
      linestart = linestart + 4
      write(plu,1001) 'displ',linestart
      linestart = linestart + numnp
      linestart = linestart + 2
      write(plu,1001) 'stres',linestart
      linestart = linestart + numnp
      linestart = linestart + 2
      write(plu,1001) 'pstre',linestart

      write(plu,1020) numnp,numel               ! Start Mesh/Piece Sect.
      write(plu,1010) '<Points>'                ! Start Point/Node data

      write(plu,1030) 'Float64','nodes',3
      do i = 1,numnp
        write(plu,2010) (x(ii,i),ii = 1,ndm) ,(0.0d0,ii = ndm+1,3)
        write(plu,*)
      end do ! i
      write(plu,1010) '</DataArray>'            ! Close Node data

      write(plu,1010) '</Points>'               ! Close Points section

      write(plu,1010) '<Cells>'                 ! Start Cell Section
      write(plu,1030) 'Int32','connectivity',1  ! Start Elements


!     Offsets memory allocation

      el_ptr(0) = 0
      do i = 1,numel
        if (i.le.3) then
        end if !i.le.3
        el_ptr(i) = el_ptr(i-1)
        do ii = 1,nen
          node = ix(ii,i)

          if (node .gt. 0) then
            write(plu,2020) node-1
            el_ptr(i) = el_ptr(i) + 1
          endif
        end do ! ii
        write(plu,*)
      end do ! i

      write(plu,1010) '</DataArray>'            ! Close Elements

!     Output offsets for elements

      write(plu,1030) 'Int32','offsets',1       ! Start Offsets
      do i = 1,numel
        write(plu,2021) el_ptr(i)
      end do ! i
      write(plu,1010) '</DataArray>'            ! Close Offsets




!     Output element connectivity type

      write(plu,1030) 'UInt8','types',1         ! Start Element types
      do i = 1,numel
        if (quadfl) then
          etyp = ix(27,i)
        else
          etyp = ix(15,i)
        end if !(quadfl)
        if(etyp.eq.-1) then                          ! LINE
          if (    el_ptr(i)-el_ptr(i-1).eq.2) then   ! 2 node
            write(plu,2021) 3
          elseif (el_ptr(i)-el_ptr(i-1).eq.3) then   ! 3 node
            write(plu,2021) 21
          endif
        elseif(etyp.eq.-2) then                      ! TRIANGLE
          if (    el_ptr(i)-el_ptr(i-1).eq.3) then   ! 3 node
            write(plu,2021) 5
          elseif (el_ptr(i)-el_ptr(i-1).eq.6) then   ! 6 node
            write(plu,2021) 22
          endif
        elseif(etyp.eq.-3) then                      ! QUADRILATERAL
          if (    el_ptr(i)-el_ptr(i-1).eq.4) then   ! 4 node
            write(plu,2021) 9
          elseif (el_ptr(i)-el_ptr(i-1).eq.8) then   ! 8 node
            write(plu,2021) 23
          elseif (el_ptr(i)-el_ptr(i-1).eq.9) then   ! 9 node
            write(plu,2021) 23
          endif
        elseif(etyp.eq.-4) then                      ! TETRAHEDRON
          if (    el_ptr(i)-el_ptr(i-1).eq.4) then   ! 4 node
            write(plu,2021) 10
          elseif (el_ptr(i)-el_ptr(i-1).eq.10) then  ! 10 node
            write(plu,2021) 24
          endif
        elseif(etyp.eq.-5) then                      ! HEXAHEDRON
          if (    el_ptr(i)-el_ptr(i-1).eq.8) then   ! 8 node
            write(plu,2021) 12
          elseif (el_ptr(i)-el_ptr(i-1).eq.20) then  ! 20 node
            write(plu,2021) 25
          elseif (el_ptr(i)-el_ptr(i-1).eq.27) then  ! 27 node
            write(plu,2021) 25
          endif
        elseif(etyp.eq.-6) then                      ! WEDGE
          if (    el_ptr(i)-el_ptr(i-1).eq.6) then   ! 6 node
            write(plu,2021) 13
          endif
        elseif(etyp.eq.-7) then                      ! PYRAMID
          if (    el_ptr(i)-el_ptr(i-1).eq.5) then   ! 5 node
            write(plu,2021) 14
          endif
        else
        if (el_ptr(i)-el_ptr(i-1)    .eq.2) then   ! 2 node line
          write(plu,2021) 3
        elseif (el_ptr(i)-el_ptr(i-1).eq.3) then   ! 3 node triangle
          write(plu,2021) 5
        elseif (el_ptr(i)-el_ptr(i-1).eq.4.and.ndm.eq.3) then
          write(plu,2021) 10                       ! 4 node tet
        elseif (el_ptr(i)-el_ptr(i-1).eq.4) then   ! 4 node quad
          write(plu,2021) 9
        elseif (el_ptr(i)-el_ptr(i-1).eq.5.and.ndm.eq.3) then
          write(plu,2021) 14                       ! 5 node pyramid
        elseif (el_ptr(i)-el_ptr(i-1).eq.6.and.ndm.eq.3) then
          write(plu,2021) 13                       ! 6 node wedge
        elseif (el_ptr(i)-el_ptr(i-1).eq.6) then   ! 6 node triangle
          write(plu,2021) 22
        elseif (el_ptr(i)-el_ptr(i-1).eq.8.and.ndm.eq.3) then
          write(plu,2021) 12                       ! 8 node brick
        elseif (el_ptr(i)-el_ptr(i-1).eq.8) then   ! 8 node quad
          write(plu,2021) 23
        elseif (el_ptr(i)-el_ptr(i-1).eq.9) then   ! 9 node quad
          write(plu,2021) 23
        elseif (el_ptr(i)-el_ptr(i-1).eq.10) then  ! 10 node tet
          write(plu,2021) 24
        elseif (el_ptr(i)-el_ptr(i-1).eq.20) then  ! 20 node brick
          write(plu,2021) 25
        elseif (el_ptr(i)-el_ptr(i-1).eq.27) then  ! 27 node brick
          write(plu,2021) 25
        endif
        endif ! etyp
      end do ! i

      write(plu,1010) '</DataArray>'            ! Close Element types
      write(plu,1010) '</Cells>'                ! Close Cell Section


      write(plu,1010) '<PointData>'             ! Start Point Data


!     Output displacements
! ======================================================================
      write(plu,1030) 'Float64','Displacements',ndf  ! Start Displ.
      do i = 1,numnp
        do ii = 1,ndf
        write(plu,2000) u(ii,i)
        end do ! ii
        write(plu,*)
      end do ! i
      write(plu,1010) '</DataArray>'                 ! Close Displ.

!     Output Velocity
! ======================================================================
      if (velofl) then
        write(plu,1030) 'Float64','Velocity',ndf     ! Start Velocity
        do i = 1,numnp
          do ii = 1,ndf
          write(plu,2000) ud(ii,i,1)
          end do ! ii
          write(plu,*)
        end do ! i
        write(plu,1010) '</DataArray>'               ! Close Velocity
      end if ! velofl

!     Output Acceleration
! ======================================================================
      if (accefl) then
        write(plu,1030) 'Float64','Acceleration',ndf ! Start Accel.
        do i = 1,numnp
          do ii = 1,ndf
          write(plu,2000) ud(ii,i,2)
          end do ! ii
          write(plu,*)
        end do ! i
        write(plu,1010) '</DataArray>'               ! Close Accel.
      end if ! accefl


!     Output stresses
! ======================================================================
      if(strefl) then
        write(plu,1030) 'Float64','Stress/Strain',istv ! Start Stresses
!        write(*,*) numnm,numni,numnp
!        read(*,*)
        do i = 1,numnp
          sigl = 0.0d0
          do ii =1,istv
            ns = 0
            do jj = ip(i-1)+1,ip(i)
              if(sig(jj,ii).ne.0.0d0) then
                ns = ns + 1
                sigl(ii) = sigl(ii) + sig(jj,ii)
              endif
            end do ! jj
            if(ns.gt.1) sigl(ii) = sigl(ii)/dble(ns)  
            write(plu,2000) sigl(ii)
!            write(*,*) "node number:",i,ip(i+1),ii
!	    write(*,*) sigl(ii)
!	    read(*,*)
          end do ! ii
          write(plu,*)
        end do ! i
        write(plu,1010) '</DataArray>'               ! Close Stresses

        write(plu,1030) 'Float64','Principal Stress',7 ! Start Prin Str.
        do i = 1,numnp
          sigl(1:7) = 0.0d0
          do ii =1,7
            ns = 0
            do jj = ip(i-1)+1,ip(i)
              if(psig(jj,ii).ne.0.0d0) then
                ns       = ns + 1
                sigl(ii) = sigl(ii) + psig(jj,ii)
              endif
            end do ! jj
            if(ns.gt.1) sigl(ii) = sigl(ii)/dble(ns)
            write(plu,2000) sigl(ii)
          end do ! ii
          write(plu,*)
        end do ! i
        write(plu,1010) '</DataArray>'               ! Close Prin Str.
      else
        write(*,*) ' No stresses output to PARAVIEW file'
      end if !

!     Output history variables
! ======================================================================
      if(histfl) then
        write(plu,1030) 'Float64','History',hplmax   ! Start History
        do i = 1,numnp
          do ii =1,hplmax
            ns = 0
            do jj = ip(i-1)+1,ip(i)
              if(his(jj,ii).ne.0.0d0) then
                ns       = ns + 1
                sigl(ii) = sigl(ii) + his(jj,ii)
              endif
            end do ! jj
            if(ns.gt.1) sigl(ii) = sigl(ii)/dble(ns)
            write(plu,2000) sigl(ii)
          end do ! ii
          write(plu,*)
        end do ! i
        write(plu,1010) '</DataArray>'               ! Close History
      else
        if (dbglvl.ge.2) write(*,4404) 'No history variables'//
     &  ' output to PARAVIEW file'
      endif


      write(plu,1010) '</PointData>'        ! Close Point Data Section

      write(plu,1010) '</Piece>'            ! Close Mesh/Piece

!     Close the XML file

      write(plu,1010) '</UnstructuredGrid> </VTKFile>'
      close(plu, status = 'keep')
      if (dbglvl.ge.1) write(*,4401) 'Finished writing file'
      if (dbglvl.ge.1) write(*,*)

!     Formats

1000  format('<?xml version="1.0"?>',/
     &       '<VTKFile type="UnstructuredGrid" version="0.1">',/
     &       '<UnstructuredGrid>')
1001  format('<!--',a,' Line:',i10,'-->')
1002  format('<!--',a,i10,'-->')
1003  format('<!--',a,1e14.5,'-->')

1010  format(a)

1020  format('<Piece NumberOfPoints="',i10,
     &       '" NumberOfCells="',i10,'">')

1030  format('<DataArray type="',a,'" Name="',a,
     &       '" NumberOfComponents="',i2,'" format="ascii">')

2000  format(1p,1e14.5E3,$)
2010  format(3(1p,1e14.5),$)
2020  format(i8,$)
2021  format(i8)

4400  format('[....] ',a)
4401  format('[ OK ] ',a)
4403  format('[info] ',a,i4)
4404  format('[info] ',a)

      end subroutine uparaview_nurb
