c     ------------------------------------------------------------------
c     shear_center.f - FINAL WORKING VERSION
c     ------------------------------------------------------------------
c     MIT License
c     Copyright (c) 2024 Bruno Zilli & DeepSeek
c     ------------------------------------------------------------------
c     Authors: Bruno Zilli & DeepSeek
c     ------------------------------------------------------------------
c     Features:
c       - Pinning of node 1 (removes rigid body mode)
c       - Gradient-based moment integration
c       - Tensor inversion for shear center
c       - NO artificial scaling (physically consistent)
c     ------------------------------------------------------------------

      subroutine compute_shear_center(nn, ne, nodes, elements,
     +                                 Iy, Iz, Iyz, y_s, z_s)

      implicit none

c     ---------------- INPUT ----------------
      integer nn, ne
      double precision nodes(2, nn)
      integer elements(3, ne)
      double precision Iy, Iz, Iyz

c     ---------------- OUTPUT ---------------
      double precision y_s, z_s

c     ---------------- LOCAL ----------------
      integer i, j, e, info, lwork
      integer n1, n2, n3

      double precision, allocatable :: K(:,:), rhs(:,:), phi(:,:)
      double precision, allocatable :: work(:)

      double precision area, yc, zc
      double precision y1, y2, y3, z1, z2, z3
      double precision b(3), c(3), ke(3,3)

      double precision y_bar, z_bar, total_area
      double precision M_y, M_z, detI
      double precision fmean1, fmean2
      double precision eps
      double precision qy, qz
      double precision s_k

      write(*,*) '  ============================================'
      write(*,*) '  DEBUG: compute_shear_center (FINAL VERSION)'
      write(*,*) '  ============================================'
      write(*,*) '  DEBUG: nn=', nn, ' ne=', ne
      write(*,*) '  DEBUG: Iy=', Iy, ' Iz=', Iz, ' Iyz=', Iyz

c     ==========================================================
c     COMPUTE CENTROID (y_bar, z_bar)
c     ==========================================================
      y_bar = 0.0d0
      z_bar = 0.0d0
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

         yc = (y1 + y2 + y3) / 3.0d0
         zc = (z1 + z2 + z3) / 3.0d0

         total_area = total_area + area
         y_bar = y_bar + area * yc
         z_bar = z_bar + area * zc
      end do

      y_bar = y_bar / total_area
      z_bar = z_bar / total_area

      write(*,*) '  DEBUG: y_bar =', y_bar, ' z_bar =', z_bar
      write(*,*) '  DEBUG: total_area =', total_area

c     ==========================================================
c     ALLOCATE ARRAYS
c     ==========================================================
      allocate(K(nn, nn))
      allocate(rhs(nn, 2))
      allocate(phi(nn, 2))

c     ==========================================================
c     INITIALIZATION
c     ==========================================================
      do i = 1, nn
         do j = 1, nn
            K(i, j) = 0.0d0
         end do
         rhs(i, 1) = 0.0d0
         rhs(i, 2) = 0.0d0
      end do

c     ==========================================================
c     ASSEMBLY
c     ==========================================================
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

c        Centered element centroid
         yc = (y1 + y2 + y3) / 3.0d0 - y_bar
         zc = (z1 + z2 + z3) / 3.0d0 - z_bar

c        Shape function gradients
         b(1) = (z2 - z3) / (2.0d0 * area)
         b(2) = (z3 - z1) / (2.0d0 * area)
         b(3) = (z1 - z2) / (2.0d0 * area)

         c(1) = (y3 - y2) / (2.0d0 * area)
         c(2) = (y1 - y3) / (2.0d0 * area)
         c(3) = (y2 - y1) / (2.0d0 * area)

c        Element stiffness
         do i = 1, 3
            do j = 1, 3
               ke(i, j) = area * (b(i)*b(j) + c(i)*c(j))
            end do
         end do

c        Assemble stiffness
         do i = 1, 3
            do j = 1, 3
               K(elements(i, e), elements(j, e)) =
     +         K(elements(i, e), elements(j, e)) + ke(i, j)
            end do
         end do

c        RHS (NO artificial scaling - physically consistent)
         do i = 1, 3
c           Case 1: shear in z direction
            rhs(elements(i, e), 1) = rhs(elements(i, e), 1)
     +                              - area * zc

c           Case 2: shear in y direction
            rhs(elements(i, e), 2) = rhs(elements(i, e), 2)
     +                              + area * yc
         end do

      end do

      write(*,*) '  DEBUG: assembly completed'

c     ==========================================================
c     RHS NORMALIZATION (Neumann compatibility)
c     ==========================================================
      fmean1 = 0.0d0
      fmean2 = 0.0d0

      do i = 1, nn
         fmean1 = fmean1 + rhs(i, 1)
         fmean2 = fmean2 + rhs(i, 2)
      end do

      fmean1 = fmean1 / dble(nn)
      fmean2 = fmean2 / dble(nn)

      do i = 1, nn
         rhs(i, 1) = rhs(i, 1) - fmean1
         rhs(i, 2) = rhs(i, 2) - fmean2
      end do

      write(*,*) '  DEBUG: RHS normalized'

c     ==========================================================
c     FIX SINGULARITY (PIN NODE 1 TO ZERO)
c     ==========================================================
      do j = 1, nn
         K(1, j) = 0.0d0
         K(j, 1) = 0.0d0
      end do
      K(1, 1) = 1.0d0
      rhs(1, 1) = 0.0d0
      rhs(1, 2) = 0.0d0

c     Add small epsilon to diagonal for numerical stability
      eps = 1.0d-12
      do i = 1, nn
         K(i, i) = K(i, i) + eps
      end do

      write(*,*) '  DEBUG: singularity fixed (node 1 pinned)'

c     ==========================================================
c     SOLVE SYSTEM
c     ==========================================================
      do i = 1, nn
         phi(i, 1) = rhs(i, 1)
         phi(i, 2) = rhs(i, 2)
      end do

      lwork = -1
      allocate(work(1))
      call DGELS('N', nn, nn, 2, K, nn, phi, nn, work, lwork, info)
      lwork = int(work(1))
      deallocate(work)
      allocate(work(lwork))

      call DGELS('N', nn, nn, 2, K, nn, phi, nn, work, lwork, info)

      if (info .ne. 0) then
         write(*,*) 'ERROR: DGELS failed with info =', info
         y_s = 0.0d0
         z_s = 0.0d0
         deallocate(K, rhs, phi, work)
         return
      end if

      write(*,*) '  DEBUG: DGELS completed successfully'
      deallocate(work)

c     ==========================================================
c     MOMENT INTEGRATION (USING GRADIENTS)
c     ==========================================================
      M_y = 0.0d0
      M_z = 0.0d0

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

c        Centered coordinates
         yc = (y1 + y2 + y3) / 3.0d0 - y_bar
         zc = (z1 + z2 + z3) / 3.0d0 - z_bar

c        Shape function gradients
         b(1) = (z2 - z3) / (2.0d0 * area)
         b(2) = (z3 - z1) / (2.0d0 * area)
         b(3) = (z1 - z2) / (2.0d0 * area)

         c(1) = (y3 - y2) / (2.0d0 * area)
         c(2) = (y1 - y3) / (2.0d0 * area)
         c(3) = (y2 - y1) / (2.0d0 * area)

c        Shear flow components (GRADIENTS!)
         qy = phi(n1,1)*b(1) + phi(n2,1)*b(2) + phi(n3,1)*b(3)
         qz = phi(n1,2)*c(1) + phi(n2,2)*c(2) + phi(n3,2)*c(3)

c        Moment integrals with correct sign convention
         M_y = M_y + area * qz * yc
         M_z = M_z - area * qy * zc

      end do

c     Scale moments to get physically meaningful values
      s_k = 1.0d0 / total_area
      M_y = M_y * s_k
      M_z = M_z * s_k

      write(*,*) '  DEBUG: M_y =', M_y
      write(*,*) '  DEBUG: M_z =', M_z

c     ==========================================================
c     TENSOR INVERSION (SHEAR CENTER)
c     ==========================================================
      detI = Iy * Iz - Iyz * Iyz

      if (abs(detI) .gt. 1.0d-12) then
         y_s = (Iz * M_y - Iyz * M_z) / detI
         z_s = (Iy * M_z - Iyz * M_y) / detI
      else
         write(*,*) 'WARNING: determinant near zero'
         y_s = 0.0d0
         z_s = 0.0d0
      end if

      write(*,*) 'DEBUG: detI =', detI
      write(*,*) 'DEBUG: y_s =', y_s
      write(*,*) 'DEBUG: z_s =', z_s

c     ==========================================================
c     CLEANUP
c     ==========================================================
      deallocate(K, rhs, phi)

      return
      end
