c     ------------------------------------------------------------------
c     shear_center.f - Shear center calculation for beam sections
c     ------------------------------------------------------------------
c     Copyright (c) 2024 Bruno Zilli & DeepSeek
c     
c     Permission is hereby granted, free of charge, to any person
c     obtaining a copy of this software and associated documentation
c     files (the "Software"), to deal in the Software without
c     restriction, including without limitation the rights to use,
c     copy, modify, merge, publish, distribute, sublicense, and/or sell
c     copies of the Software, and to permit persons to whom the
c     Software is furnished to do so, subject to the following
c     conditions:
c     
c     The above copyright notice and this permission notice shall be
c     included in all copies or substantial portions of the Software.
c     
c     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
c     EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
c     OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
c     NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
c     HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
c     WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
c     FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
c     OTHER DEALINGS IN THE SOFTWARE.
c     ------------------------------------------------------------------
c     Authors: Bruno Zilli & DeepSeek
c     ------------------------------------------------------------------
c     
c     Final stable version (validated on HEB200):
c     - 1 fixed node to remove singularity
c     - RHS normalization (Neumann compatibility)
c     - DGELS solver (least squares, tolerant to near-singular)
c     - Integration of gradient of phi (shear flow)
c     
c     ------------------------------------------------------------------
      
      subroutine compute_shear_center(nn, ne, nodes, elements,
     +                                 Iy, Iz, y_s, z_s)
      
      implicit none
      
c     ---------------- INPUT ----------------
      integer, intent(in) :: nn, ne
      double precision, intent(in) :: nodes(2, nn)
      integer, intent(in) :: elements(3, ne)
      double precision, intent(in) :: Iy, Iz
      
c     ---------------- OUTPUT ---------------
      double precision, intent(out) :: y_s, z_s
      
c     ---------------- LOCAL ----------------
      integer :: i, j, e, info, lwork
      integer, allocatable :: ipiv(:)
      double precision, allocatable :: stiff(:,:), rhs(:,:), phi(:,:)
      double precision, allocatable :: work(:)
      double precision :: area, yc, zc
      double precision :: b(3), c(3)
      double precision :: ke(3,3)
      integer :: n1, n2, n3
      double precision :: fmean1, fmean2, total_area
      double precision :: dphidy, dphidz
      double precision :: y1, y2, y3, z1, z2, z3
      double precision :: eps
      
      write(*,*) '  DEBUG: compute_shear_center started'
      write(*,*) '  DEBUG: nn=', nn, ' ne=', ne
      write(*,*) '  DEBUG: Iy=', Iy, ' Iz=', Iz
      
c     allocate system
      allocate(stiff(nn, nn))
      allocate(rhs(nn, 2))
      allocate(phi(nn, 2))
      allocate(ipiv(nn))
      
c     init
      do i = 1, nn
         do j = 1, nn
            stiff(i,j) = 0.0d0
         end do
         rhs(i,1) = 0.0d0
         rhs(i,2) = 0.0d0
      end do
      
c     ==============================
c     ASSEMBLY
c     ==============================
      do e = 1, ne
         
         n1 = elements(1, e)
         n2 = elements(2, e)
         n3 = elements(3, e)
         
         y1 = nodes(1, n1)
         z1 = nodes(2, n1)
         y2 = nodes(1, n2)
         z2 = nodes(2, n2)
         y3 = nodes(1, n3)
         z3 = nodes(2, n3)
         
         yc = (y1 + y2 + y3) / 3.0d0
         zc = (z1 + z2 + z3) / 3.0d0
         
         area = 0.5d0 * ((y2 - y1)*(z3 - z1) - (y3 - y1)*(z2 - z1))
         if (area .lt. 0.0d0) area = -area
         
c        shape function gradients
         b(1) = (z2 - z3) / (2.0d0 * area)
         b(2) = (z3 - z1) / (2.0d0 * area)
         b(3) = (z1 - z2) / (2.0d0 * area)
         
         c(1) = (y3 - y2) / (2.0d0 * area)
         c(2) = (y1 - y3) / (2.0d0 * area)
         c(3) = (y2 - y1) / (2.0d0 * area)
         
c        element stiffness
         do i = 1, 3
            do j = 1, 3
               ke(i, j) = area * (b(i)*b(j) + c(i)*c(j))
            end do
         end do
         
c        assemble stiff
         do i = 1, 3
            do j = 1, 3
               stiff(elements(i,e), elements(j,e)) =
     +         stiff(elements(i,e), elements(j,e)) + ke(i, j)
            end do
         end do
         
c        RHS: Case 1 (f = -2z) and Case 2 (f = 2y)
         do i = 1, 3
            rhs(elements(i,e), 1) = rhs(elements(i,e), 1)
     +                              - 2.0d0 * area * zc / 3.0d0
            
            rhs(elements(i,e), 2) = rhs(elements(i,e), 2)
     +                              + 2.0d0 * area * yc / 3.0d0
         end do
         
      end do
      
      write(*,*) '  DEBUG: assembly completed'
      
c     ==============================
c     RHS NORMALIZATION (Neumann compatibility)
c     ==============================
      fmean1 = 0.0d0
      fmean2 = 0.0d0
      total_area = 0.0d0
      
      do e = 1, ne
         n1 = elements(1, e)
         n2 = elements(2, e)
         n3 = elements(3, e)
         
         y1 = nodes(1, n1)
         z1 = nodes(2, n1)
         y2 = nodes(1, n2)
         z2 = nodes(2, n2)
         y3 = nodes(1, n3)
         z3 = nodes(2, n3)
         
         area = 0.5d0 * ((y2 - y1)*(z3 - z1) - (y3 - y1)*(z2 - z1))
         if (area .lt. 0.0d0) area = -area
         
         total_area = total_area + area
         
         fmean1 = fmean1 + area * 
     +           (rhs(n1,1) + rhs(n2,1) + rhs(n3,1)) / 3.0d0
         
         fmean2 = fmean2 + area * 
     +           (rhs(n1,2) + rhs(n2,2) + rhs(n3,2)) / 3.0d0
      end do
      
      if (total_area .gt. 0.0d0) then
         fmean1 = fmean1 / total_area
         fmean2 = fmean2 / total_area
      end if
      
      do i = 1, nn
         rhs(i,1) = rhs(i,1) - fmean1
         rhs(i,2) = rhs(i,2) - fmean2
      end do
      
      write(*,*) '  DEBUG: RHS normalized'
      
c     ==============================
c     FIX SINGULARITY (1 node)
c     ==============================
      i = 1
      do j = 1, nn
         stiff(i, j) = 0.0d0
         stiff(j, i) = 0.0d0
      end do
      stiff(i, i) = 1.0d0
      rhs(i, 1) = 0.0d0
      rhs(i, 2) = 0.0d0
      
      eps = 1.0d-12
      do i = 1, nn
         stiff(i, i) = stiff(i, i) + eps
      end do
      
      write(*,*) '  DEBUG: singularity fixed (1 node + epsilon)'
      
c     ==============================
c     SOLVE with DGELS (least squares)
c     ==============================
      do i = 1, nn
         phi(i, 1) = rhs(i, 1)
         phi(i, 2) = rhs(i, 2)
      end do
      
      lwork = -1
      allocate(work(1))
      call DGELS('N', nn, nn, 2, stiff, nn, phi, nn, work, lwork, info)
      lwork = int(work(1))
      deallocate(work)
      allocate(work(lwork))
      
      call DGELS('N', nn, nn, 2, stiff, nn, phi, nn, work, lwork, info)
      
      if (info .ne. 0) then
         write(*,*) 'ERROR: DGELS failed with info =', info
         y_s = 0.0d0
         z_s = 0.0d0
         deallocate(stiff, rhs, phi, ipiv, work)
         return
      end if
      
      write(*,*) '  DEBUG: DGELS completed successfully'
      
c     ==============================
c     GRADIENT INTEGRATION (shear flow)
c     ==============================
      y_s = 0.0d0
      z_s = 0.0d0
      
      do e = 1, ne
         n1 = elements(1, e)
         n2 = elements(2, e)
         n3 = elements(3, e)
         
         y1 = nodes(1, n1)
         z1 = nodes(2, n1)
         y2 = nodes(1, n2)
         z2 = nodes(2, n2)
         y3 = nodes(1, n3)
         z3 = nodes(2, n3)
         
         area = 0.5d0 * ((y2 - y1)*(z3 - z1) - (y3 - y1)*(z2 - z1))
         if (area .lt. 0.0d0) area = -area
         
         yc = (y1 + y2 + y3) / 3.0d0
         zc = (z1 + z2 + z3) / 3.0d0
         
c        shape function gradients
         b(1) = (z2 - z3) / (2.0d0 * area)
         b(2) = (z3 - z1) / (2.0d0 * area)
         b(3) = (z1 - z2) / (2.0d0 * area)
         
         c(1) = (y3 - y2) / (2.0d0 * area)
         c(2) = (y1 - y3) / (2.0d0 * area)
         c(3) = (y2 - y1) / (2.0d0 * area)
         
c        Gradient of φ₁ (case 1)
         dphidz = phi(n1,1)*c(1) + phi(n2,1)*c(2) + phi(n3,1)*c(3)
         y_s = y_s + area * dphidz * yc
         
c        Gradient of φ₂ (case 2)
         dphidy = phi(n1,2)*b(1) + phi(n2,2)*b(2) + phi(n3,2)*b(3)
         z_s = z_s + area * dphidy * zc
         
      end do
      
c     Final normalization
      if (Iz .gt. 0.0d0) y_s = y_s / Iz
      if (Iy .gt. 0.0d0) z_s = z_s / Iy
      
      write(*,*) 'DEBUG: y_s =', y_s
      write(*,*) 'DEBUG: z_s =', z_s
      
c     cleanup
      deallocate(stiff, rhs, phi, ipiv, work)
      
      return
      end
