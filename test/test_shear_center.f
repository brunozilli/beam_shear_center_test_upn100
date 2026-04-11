c     ------------------------------------------------------------------
c     test_shear_center.f - Test program for shear center calculation
c     ------------------------------------------------------------------
c     Copyright (c) 2024 Bruno Zilli & DeepSeek
c     ------------------------------------------------------------------
c     Authors: Bruno Zilli & DeepSeek
c     ------------------------------------------------------------------
      
      program test_shear_center
      
      implicit none
      
c     .. Parameters ..
      integer, parameter :: MAX_NODES = 10000
      integer, parameter :: MAX_ELEMENTS = 20000
      
c     .. Local variables ..
      integer :: nn, ne, i
      double precision, allocatable :: y(:), z(:)
      integer, allocatable :: conn(:,:)
      double precision, allocatable :: nodes(:,:)
      integer, allocatable :: elements(:,:)
      double precision :: y_s, z_s
      double precision :: y_c, z_c, area, Iy, Iz, Iyz, J_polar
      character*80 :: filename
      
c     Temporary arrays for reading
      double precision :: ytmp(MAX_NODES), ztmp(MAX_NODES)
      integer :: conntmp(3, MAX_ELEMENTS)
      
c     Test: L100x100x10
      write(*,*) '========================================'
      write(*,*) 'Test: L100x100x10 45° rotated and CG centered'
      write(*,*) '========================================'
      
      filename = 'meshes/Mesh_test_L.unv'
      call read_section_mesh_unv(filename, MAX_NODES, MAX_ELEMENTS,
     +                           nn, ne, conntmp, ytmp, ztmp)
      
      if (nn .eq. 0 .or. ne .eq. 0) then
         write(*,*) 'ERROR: Failed to read mesh'
         stop
      end if
      
      write(*,*) 'Mesh read: nodes =', nn, ' elements =', ne
      
c     Allocate arrays
      allocate(y(nn))
      allocate(z(nn))
      allocate(conn(3, ne))
      allocate(nodes(2, nn))
      allocate(elements(3, ne))
      
c     Copy data
      do i = 1, nn
         y(i) = ytmp(i)
         z(i) = ztmp(i)
         nodes(1, i) = ytmp(i)
         nodes(2, i) = ztmp(i)
      end do
      
      do i = 1, ne
         conn(1, i) = conntmp(1, i)
         conn(2, i) = conntmp(2, i)
         conn(3, i) = conntmp(3, i)
         elements(1, i) = conntmp(1, i)
         elements(2, i) = conntmp(2, i)
         elements(3, i) = conntmp(3, i)
      end do
      
c     Compute section properties
      call compute_section_properties(nn, ne, conn, y, z,
     +                                area, y_c, z_c, Iy, Iz, Iyz, 
     +                                J_polar)
      
      write(*,*) 'Section Properties:'
      write(*,*) '  Area:', area
      write(*,*) '  Centroid (y_c, z_c):', y_c, z_c
      write(*,*) '  Iy:', Iy, ' Iz:', Iz, ' Iyz:', Iyz
      
c     Compute shear center
      call compute_shear_center(nn, ne, nodes, elements, Iy, Iz, Iyz,
     +                          y_s, z_s)
      
      write(*,*) 'Shear Center (y_s, z_s):', y_s, z_s
      
c     Verification for L100x100x10 (mesh centered at centroid)
      write(*,*) ''
      write(*,*) 'Expected values (from theory):'
      write(*,*) '  y_s ≈ 6.8 mm (from centroid)'
      write(*,*) '  z_s ≈ 6.8 mm (from centroid)'
      write(*,*) '  (absolute from corner: ~35.7 mm)'
      
c     Check symmetry (for equal legs)
      if (abs(y_s - z_s) .lt. 1.0d0) then
         write(*,*) 'PASS: y_s ≈ z_s (symmetry)'
      else
         write(*,*) 'WARNING: y_s != z_s, diff =', y_s - z_s
      end if
      
c     Check that shear center is not at centroid (for L-section)
      if (abs(y_s) .gt. 1.0d0 .and. abs(z_s) .gt. 1.0d0) then
         write(*,*) 'PASS: Shear center is NOT at centroid (expected)'
      else
         write(*,*) 'WARNING: Shear center too close to centroid'
      end if
      
      deallocate(y, z, conn, nodes, elements)
      
      stop
      end
