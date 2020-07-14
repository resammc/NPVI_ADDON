c-----[--.----+----.----+----.-----------------------------------------]

c      Date             : 25 October 2016
c      Last modified:   : 07 April 2017
c      Release          : 3.0

c-----[--.----+----.----+----.-----------------------------------------]
c                  NURBS PARAVIEW EXPORT (NPVI)
c                    
C              BY:            HENNING VENGHAUS
C              EMAIL ADDRESS: henning.venghaus@st.ovgu.de
C-----[--.----+----.----+----.-----------------------------------------]
C     quadratic paraview export by B.Sc. Henning Venghaus
C     Plot     values

      integer        nd_quad, ne_quad
      common /npvi/  nd_quad, ne_quad 
C-----[--.----+----.----+----.-----------------------------------------]
C     non-uniform subdivision of NURBS elements
      integer        npl_int_nusd
      common /npvi/  npl_int_nusd(3)
      
      logical        serefl, nusdfl
      common /npvi/  serefl, nusdfl

