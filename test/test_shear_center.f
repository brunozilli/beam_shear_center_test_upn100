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
      
c     Test: HEB200 (Doubly Symmetric)
      write(*,*) '========================================'
      write(*,*) 'Test: HEB200 (Doubly Symmetric)'
      write(*,*) '========================================'
      
      filename = 'meshes/HEB200_mm.unv'
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
      call compute_shear_center(nn, ne, nodes, elements, Iy, Iz,
     +                          y_s, z_s)
      
      write(*,*) 'Shear Center (y_s, z_s):', y_s, z_s
      
c     Verification - should be at centroid (0,0) for HEB200
      if (abs(y_s) .lt. 1.0d-1) then
         write(*,*) 'PASS: y_s ≈ 0'
      else
         write(*,*) 'WARNING: y_s != 0, diff =', y_s
      end if
      
      if (abs(z_s) .lt. 1.0d-1) then
         write(*,*) 'PASS: z_s ≈ 0'
      else
         write(*,*) 'WARNING: z_s != 0, diff =', z_s
      end if
      
      deallocate(y, z, conn, nodes, elements)
      
      stop
      end
