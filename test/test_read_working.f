c     test_read_working.f - Test read_section_mesh_unv
c     MIT License - Copyright (c) 2024 Bruno Zilli & DeepSeek

      program test_read_working
      
      implicit none
      
c     Parameters
      integer, parameter :: max_nodes = 50000
      integer, parameter :: max_triangles = 100000
      
c     Variables
      character(len=256) :: filename
      integer :: n_nodes, n_triangles
      double precision :: y(max_nodes), z(max_nodes)
      integer :: conn(3, max_triangles)
      integer :: i
      
c     Get filename from command line
      if (iargc() .lt. 1) then
         write(*,*) 'Usage: test_read_working <meshfile.unv>'
         stop
      end if
      
      call getarg(1, filename)
      
      write(*,*) '========================================='
      write(*,*) '  TEST MESH READER'
      write(*,*) '========================================='
      write(*,*) 'File: ', trim(filename)
      write(*,*)
      
c     Read mesh using working subroutine
      call read_section_mesh_unv(filename, max_nodes, max_triangles,
     +                           n_nodes, n_triangles, conn, y, z)
      
      write(*,*)
      write(*,*) '-----------------------------------------'
      write(*,*) 'Mesh summary:'
      write(*,*) '  Nodes:    ', n_nodes
      write(*,*) '  Triangles:', n_triangles
      write(*,*) '-----------------------------------------'
      
      if (n_nodes .gt. 0) then
         write(*,*) 'First 5 nodes:'
         do i = 1, min(5, n_nodes)
            write(*,*) '  Node', i, ': y=', y(i), ' z=', z(i)
         end do
      end if
      
      if (n_triangles .gt. 0) then
         write(*,*) 'First 5 triangles:'
         do i = 1, min(5, n_triangles)
            write(*,*) '  Tri', i, ': ', conn(1,i), conn(2,i),
     +                 conn(3,i)
         end do
      end if
      
      stop
      end
