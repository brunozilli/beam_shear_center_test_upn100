c     ------------------------------------------------------------------
c     test_shear_center.f - Test shear centre calculation
c     ------------------------------------------------------------------
c     MIT License
c     Copyright (c) 2024 Bruno Zilli & DeepSeek
c     ------------------------------------------------------------------

      program test_shear_center
      
      implicit none
      
c     Parameters
      integer, parameter :: max_nodes = 50000
      integer, parameter :: max_triangles = 100000
      
c     Variables
      character(len=256) :: filename
      integer :: n_nodes, n_triangles
      double precision :: y_arr(max_nodes), z_arr(max_nodes)
      integer :: conn(3, max_triangles)
      
c     Additional variables for section properties
      double precision :: A, y_c, z_c, I_y, I_z, I_yz, J_polar
      double precision :: y_s, z_s
      
c     Converted arrays for shear_center (nodes format)
      double precision, allocatable :: nodes(:,:)
      integer, allocatable :: elements(:,:)
      integer :: i
      
c     Get filename from command line
      if (iargc() .lt. 1) then
         write(*,*) 'Usage: test_shear_center <meshfile.unv>'
         write(*,*) 'Example:'
         write(*,*) '  ./test/test_shear_center meshes/HEB200_mm.unv'
         stop
      end if
      
      call getarg(1, filename)
      
      write(*,*) '========================================='
      write(*,*) '  SHEAR CENTRE CALCULATION'
      write(*,*) '========================================='
      write(*,*) 'File: ', trim(filename)
      write(*,*)
      
c     Read mesh using working subroutine
      call read_section_mesh_unv(trim(filename), max_nodes, 
     +                           max_triangles, n_nodes, n_triangles,
     +                           conn, y_arr, z_arr)
      
      write(*,*) 'Nodes: ', n_nodes
      write(*,*) 'Triangles: ', n_triangles
      write(*,*)
      
c     Compute section properties (using the correct interface)
      call compute_section_properties(n_nodes, n_triangles, conn,
     +                                y_arr, z_arr,
     +                                A, y_c, z_c, I_y, I_z, I_yz, 
     +                                J_polar)
      
      write(*,*) 'Section Properties:'
      write(*,*) '  Area:      ', A
      write(*,*) '  Centroid:  (', y_c, ', ', z_c, ')'
      write(*,*) '  I_y:       ', I_y
      write(*,*) '  I_z:       ', I_z
      write(*,*) '  I_yz:      ', I_yz
      write(*,*) '  J_polar:   ', J_polar
      write(*,*)
      
c     Convert to format expected by shear_center
c     shear_center expects nodes(2, nn) and elements(3, ne)
      allocate(nodes(2, n_nodes))
      allocate(elements(3, n_triangles))
      
      do i = 1, n_nodes
         nodes(1, i) = y_arr(i)
         nodes(2, i) = z_arr(i)
      end do
      
      do i = 1, n_triangles
         elements(1, i) = conn(1, i)
         elements(2, i) = conn(2, i)
         elements(3, i) = conn(3, i)
      end do
      
c     Compute shear centre
      call compute_shear_center(n_nodes, n_triangles, nodes, elements,
     +                          I_y, I_z, I_yz, y_s, z_s)
      
      write(*,*) 'Shear Centre:'
      write(*,*) '  y_s: ', y_s
      write(*,*) '  z_s: ', z_s
      write(*,*)
      
c     For doubly symmetric sections (like HEB200), shear centre should
c     coincide with centroid
      if (abs(y_s - y_c) .lt. 1.0d-6 .and. 
     +    abs(z_s - z_c) .lt. 1.0d-6) then
         write(*,*) 'PASS: Shear centre matches centroid'
      else
         write(*,*) 'NOTE: Section may be asymmetric'
      end if
      
c     Clean up
      deallocate(nodes, elements)
      
      stop
      end
