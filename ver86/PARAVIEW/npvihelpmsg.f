!$Id:$
      subroutine npvihelpmsg()

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
!       2. Added Bezier elements and merge flags            23/06/2020
!-----[--.----+----.----+----.-----------------------------------------]
!      Purpose: Display help message for npvi routine

!      Inputs:

!      Outputs: On stdout
!-----[--+---------+---------+---------+---------+---------+---------+-]

      implicit none

      write(*,*) 'Usage: npvi,specifier,n1,n2,n3'
      write(*,*) ' specifier: can contain certin substings'
      write(*,*) '  "help" -> display this help'
      write(*,*) '  "t"    -> time flag, used for transient problems'
      write(*,*) '  "s"    -> stress flag, stress output if active'
      write(*,*) '  "v"    -> velocity flag, velocity output if active'
      write(*,*) '  "a"    -> acceleration flag, acc. output if active'
      write(*,*) '  "m"    -> merge flag for repeated nodes/control '//
     &           'points'
      write(*,*) '  "nusd" -> input parameters for non-uniform '//
     &           'subdivision'
      write(*,*) '  "feop" -> generate Standard Feap input file'
      write(*,*) ' n1: FE-element type'
      write(*,*) '  1 -> linear FE element will be used (default)'
      write(*,*) '  2 -> quadratic Serendipity FE element will be'//
     &           ' used'
      write(*,*) '  3 -> quadratic Lagranian FE element will be'//
     &           ' used (only Standard Feap output)'
      write(*,*) '  4 -> Bezier elements will be used'
      write(*,*) ' n2: subdivision factor'
      write(*,*) '  0 -> value for n2 if left blank, indicates'//
     &           ' *non-uniform* subdivision'
      write(*,*) '  1 -> one NURBS element = one FE element (default)'
      write(*,*) '  2+-> one NURBS element = multiple FE elements'
      write(*,*) ' n3: output level'
      write(*,*) '  0 -> least output, only .vtu file name (default)'
      write(*,*) '  1 -> status information'
      write(*,*) '  2 -> additional information (e.g.'//
     &           'mesh dimension or number of elements)'
      write(*,*) 'Report Bugs to:'
      write(*,*) ' "resam.makvandi@ovgu.de" or'
      write(*,*) ' "henning.venghaus@st.ovgu.de"'

      end subroutine npvihelpmsg
