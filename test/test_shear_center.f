c     ------------------------------------------------------------------
c     test_shear_center.f - Test shear centre calculation
c     ------------------------------------------------------------------
c     MIT License
c     Copyright (c) 2024 Bruno Zilli & DeepSeek
c     ------------------------------------------------------------------
c     C section with opening towards +X, CG centered in plane X-Y,
c     X symmetry axis.
c     ------------------------------------------------------------------

      program test_shear_center
      
      implicit none
      
c     Parameters
      integer, parameter :: max_nodes = 50000
      integer, parameter :: max_triangles = 100000
      
c     Variables
      character(len=256) :: filename
      integer :: n_nodes, n_triangles
      double precision :: x_arr(max_nodes), y_arr(max_nodes)
      integer :: conn(3, max_triangles)
      
c     Additional variables for section properties
      double precision :: A, x_c, y_c, Ixx, Iyy, Ixy, J_polar
      double precision :: x_sc, y_sc
      
c     Converted arrays for shear_center (nodes format)
      double precision, allocatable :: nodes(:,:)
      integer, allocatable :: elements(:,:)
      integer :: i
      
c     Get filename from command line
      if (iargc() .lt. 1) then
         write(*,*) 'Usage: test_shear_center <meshfile.unv>'
         write(*,*) 'Example:'
         write(*,*) '  ./test/test_shear_center meshes/upn100_1D.unv'
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
     +                           conn, x_arr, y_arr)
      
      write(*,*) 'Nodes: ', n_nodes
      write(*,*) 'Triangles: ', n_triangles
      write(*,*)
      
c     Compute section properties (using X and Y axes)
      call compute_section_properties(n_nodes, n_triangles, conn,
     +                                x_arr, y_arr,
     +                                A, x_c, y_c, Ixx, Iyy, Ixy, 
     +                                J_polar)
      
      write(*,*) 'Section Properties:'
      write(*,*) '  Area:      ', A
      write(*,*) '  Centroid:  (', x_c, ', ', y_c, ')'
      write(*,*) '  Ixx:       ', Ixx
      write(*,*) '  Iyy:       ', Iyy
      write(*,*) '  Ixy:       ', Ixy
      write(*,*) '  J_polar:   ', J_polar
      write(*,*)
      
c     Convert to format expected by shear_center
c     shear_center expects nodes(2, nn) and elements(3, ne)
      allocate(nodes(2, n_nodes))
      allocate(elements(3, n_triangles))
      
      do i = 1, n_nodes
         nodes(1, i) = x_arr(i)
         nodes(2, i) = y_arr(i)
      end do
      
      do i = 1, n_triangles
         elements(1, i) = conn(1, i)
         elements(2, i) = conn(2, i)
         elements(3, i) = conn(3, i)
      end do
      
c     Compute shear centre
      call compute_shear_center(n_nodes, n_triangles, nodes, elements,
     +                          Ixx, Iyy, Ixy, x_sc, y_sc)
      
      write(*,*) 'Shear Centre:'
      write(*,*) '  X_sc: ', x_sc
      write(*,*) '  Y_sc: ', y_sc
      write(*,*)
      
c     For doubly symmetric sections (like HEB), shear centre should
c     coincide with centroid
      if (abs(x_sc - x_c) .lt. 1.0d-6 .and. 
     +    abs(y_sc - y_c) .lt. 1.0d-6) then
         write(*,*) 'PASS: Shear centre matches centroid'
      else
         write(*,*) 'NOTE: Section may be asymmetric'
      end if
      
c     Clean up
      deallocate(nodes, elements)
      
      stop
      end
