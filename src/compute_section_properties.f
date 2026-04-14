c=======================================================================
c compute_section_properties.f
c=======================================================================
c MIT License
c Copyright (c) 2024 Bruno Zilli & DeepSeek
c
c Permission is hereby granted, free of charge, to any person obtaining a copy
c of this software and associated documentation files (the "Software"), to deal
c in the Software without restriction, including without limitation the rights
c to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
c copies of the Software, and to permit persons to whom the Software is
c furnished to do so, subject to the following conditions:
c
c The above copyright notice and this permission notice shall be included in all
c copies or substantial portions of the Software.
c
c THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
c IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
c FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
c AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
c LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
c OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
c SOFTWARE.
c=======================================================================

      subroutine compute_section_properties(nnode_sec, ntri_sec,
     &                                      conn_sec, x_sec, y_sec,
     &                                      A, x_c, y_c,
     &                                      Ixx, Iyy, Ixy, J_polar)

      implicit none

c     INPUT
      integer, intent(in) :: nnode_sec, ntri_sec
      integer, intent(in) :: conn_sec(3, ntri_sec)
      double precision, intent(in) :: x_sec(nnode_sec)
      double precision, intent(in) :: y_sec(nnode_sec)

c     OUTPUT
      double precision, intent(out) :: A, x_c, y_c
      double precision, intent(out) :: Ixx, Iyy, Ixy
      double precision, intent(out) :: J_polar

c     LOCALS
      integer :: i, i1, i2, i3
      double precision :: x1, x2, x3, y1, y2, y3
      double precision :: area_tri
      double precision :: S_x, S_y
      double precision :: Ixx0, Iyy0, Ixy0
      double precision :: det

c     ==========================================================
c     CHECK FOR 1D MESH (no triangles)
c     ==========================================================
      if (ntri_sec .le. 0) then
         write(*,*) 'INFO: 1D mesh detected - skipping area calculation'
         A = 0.0d0
         x_c = 0.0d0
         y_c = 0.0d0
         Ixx = 0.0d0
         Iyy = 0.0d0
         Ixy = 0.0d0
         J_polar = 0.0d0
         return
      endif

c     Initialise output variables
      A   = 0.0d0
      S_x = 0.0d0
      S_y = 0.0d0
      Ixx = 0.0d0
      Iyy = 0.0d0
      Ixy = 0.0d0

c     Loop over all triangular elements
      do i = 1, ntri_sec

        i1 = conn_sec(1,i)
        i2 = conn_sec(2,i)
        i3 = conn_sec(3,i)

        x1 = x_sec(i1)
        y1 = y_sec(i1)
        x2 = x_sec(i2)
        y2 = y_sec(i2)
        x3 = x_sec(i3)
        y3 = y_sec(i3)

c       Oriented area (2 x area) - sign indicates orientation
        det = (x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1)
        area_tri = 0.5d0 * det

c       Accumulate total area and first moments
        A = A + area_tri
        S_x = S_x + area_tri * (x1 + x2 + x3) / 3.0d0
        S_y = S_y + area_tri * (y1 + y2 + y3) / 3.0d0

c       Second moments about origin (0,0) using exact triangle formulae
c       Ixx = ∫ y² dA  (moment about X axis)
        Ixx0 = area_tri *
     &        (y1*y1 + y2*y2 + y3*y3 +
     &         y1*y2 + y2*y3 + y3*y1) / 6.0d0

c       Iyy = ∫ x² dA  (moment about Y axis)
        Iyy0 = area_tri *
     &        (x1*x1 + x2*x2 + x3*x3 +
     &         x1*x2 + x2*x3 + x3*x1) / 6.0d0

c       Ixy = ∫ x*y dA  (product of inertia)
        Ixy0 = area_tri *
     &        (2.0d0*(x1*y1 + x2*y2 + x3*y3) +
     &         x1*y2 + x2*y1 +
     &         x2*y3 + x3*y2 +
     &         x3*y1 + x1*y3) / 12.0d0

        Ixx = Ixx + Ixx0
        Iyy = Iyy + Iyy0
        Ixy = Ixy + Ixy0

      end do

c     Check for valid area
      if (dabs(A) .lt. 1.0d-16) then
         write(*,*) 'ERROR: zero or invalid section area'
         return
      endif

c     Section centroid
      x_c = S_x / A
      y_c = S_y / A

c     Parallel axis theorem (Steiner) - shift to centroid
      Ixx = Ixx - A * y_c*y_c
      Iyy = Iyy - A * x_c*x_c
      Ixy = Ixy - A * x_c*y_c

c     Polar moment of inertia (not the torsion constant!)
      J_polar = Ixx + Iyy

      return
      end subroutine compute_section_properties
