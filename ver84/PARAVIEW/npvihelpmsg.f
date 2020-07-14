      subroutine npvihelpmsg()
      
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
c      Purpose: display help message for npvi routine

c      Inputs:

c      Outputs: on stdout
c-----[--+---------+---------+---------+---------+---------+---------+-]      
            
      implicit none
      write(*,*) 'Usage: npvi,specifier,n1,n2,n3'
      write(*,*) ' specifier: can contain certin substings'
      write(*,*) '  "help" -> display this help  *OR*'
      write(*,*) '  "t"    -> time flag, used for transient problems'
      write(*,*) '  "s"    -> stress flag, stress output if active'
      write(*,*) '  "v"    -> velocity flag, velocity output if active'
      write(*,*) '  "a"    -> acceleration flag, acc. output if active'
      write(*,*) '  "nusd" -> input parameters for non-uniform'//
     & 'subdivision'
      write(*,*) '  "feop" -> generate Standard Feap input file'
      write(*,*) ' n1: FE-element type'
      write(*,*) '  1 -> linear FE element will be used (default)'
      write(*,*) '  2 -> quadratic Serendipity FE element will be'//
     & ' used'
      write(*,*) '  3 -> quadratic Lagranian FE element will be'// 
     & ' used (only Standard Feap output)'
      write(*,*) ' n2: subdivision factor'
      write(*,*) '  0 -> value for n2 if left blank, indicates'//
     & ' *non-uniform* subdivision'
      write(*,*) '  1 -> one NURBS element = one FE element (default)'
      write(*,*) '  2+-> one NURBS element = multiple FE elements'
      write(*,*) ' n3: output level'
      write(*,*) '  0 -> least output, only .vtu file name (default)'
      write(*,*) '  1 -> status information'
      write(*,*) '  2 -> additional information (e.g.'//
     &'mesh dimension or number of elements)'
      write(*,*) 'Report Bugs to:'
      write(*,*) ' "henning.venghaus@st.ovgu.de" or'
      write(*,*) ' "resam.makvandi@ovgu.de"'
      end subroutine