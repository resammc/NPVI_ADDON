!$Id:$
      subroutine npvi_python()

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
!       Original version                                    14/07/2020
!-----[--.----+----.----+----.-----------------------------------------]
!      Purpose: Calculates the full 2D extraction operator
!-----[--.----+----.----+----.-----------------------------------------]
      implicit   none
      
      character(:), allocatable :: filename_prefix, compressor_type, compression_level
      
      filename_prefix   = 'Comp_'
      ! Compressor types:
      ! 1: vtkZLibDataCompressor, 2: vtkLZ4DataCompressor , 3: vtkLZMADataCompressor
      compressor_type   = '3'
      compression_level = '9'

      open(9001, file = 'npvi_vtkcomp.py') 
      
      write(9001,'(a)') 'from paraview.simple import *'
      write(9001,'(a)') 'import os'
      write(9001,'(a)') 'import fnmatch'
      write(9001,'(a)') ' '
      write(9001,'(a)') '# get the current directory path'
      write(9001,'(a)') 'cwd = os.getcwd()'
      write(9001,'(a)') ' '
      write(9001,'(a)') '# find all the vtu files in the current directory'
      write(9001,'(a)') 'files = fnmatch.filter(os.listdir(cwd), "*.vtu")'
      write(9001,'(a)') 'files.sort()'
      write(9001,'(a)') '# print the names of the files (to be converted) on the screen'
      write(9001,'(a)') 'print("Found:")'
      write(9001,'(a)') 'print(files)'
      write(9001,'(a)') ' '
      write(9001,'(a)') 'for fileno, entry in enumerate(files):'
      write(9001,'(a)') '  print("==========")'
      write(9001,'(a)') '  print("Processing File", fileno + 1, "of", len(files))'
      write(9001,'(a)') '  # check the file to see if it is already compressed or not'
      write(9001,'(a)') '  print("Checking: ", entry)'
      write(9001,'(a)') '  fline = open(entry).readline().rstrip()'
      write(9001, 1000) '  if ("compressor" in fline or os.path.isfile("',filename_prefix,'"+entry)):'
      write(9001,'(a)') '    print(entry, "is already compressed! skipping to the next file...")'
      write(9001,'(a)') '    # if already compressed, skip the file'
      write(9001,'(a)') '    continue'
      write(9001,'(a)') '  print("Converting: ", entry)'
      write(9001,'(a)') '  vtu_file = XMLUnstructuredGridReader(FileName=[entry])'
      write(9001,'(a)') '  writer = XMLUnstructuredGridWriter \'
      write(9001, 1001) '  (FileName=["',filename_prefix,'"+entry],CompressorType=',compressor_type,',CompressionLevel=',compression_level,') '
      write(9001,'(a)') '  writer.DataMode="Binary" # appended, ascii'
      write(9001,'(a)') '  writer.UpdatePipeline()'
      write(9001,'(a)') '  print("Done.")'
      write(9001,'(a)') ' '
      write(9001,'(a)') 'print("==========")'
      write(9001,'(a)') 'print("***FINISHED***")'
      write(9001,'(a)') ''
      
1000  format('',3a)
1001  format('',7a)
      

      end
