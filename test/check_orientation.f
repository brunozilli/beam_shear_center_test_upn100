      program check_orientation
      implicit none
      integer, parameter :: max_nodes = 50000, max_tri = 100000
      integer :: n_nodes, n_tri, i, n1, n2, n3
      double precision :: y(max_nodes), z(max_nodes)
      integer :: conn(3, max_tri)
      double precision :: y_min, y_max, z_min, z_max
      double precision :: y_sum, z_sum
      
      call read_section_mesh_unv('meshes/Mesh_test_L.unv', 
     +     max_nodes, max_tri, n_nodes, n_tri, conn, y, z)
      
      y_min = y(1)
      y_max = y(1)
      z_min = z(1)
      z_max = z(1)
      y_sum = 0.0d0
      z_sum = 0.0d0
      
      do i = 1, n_nodes
         y_min = min(y_min, y(i))
         y_max = max(y_max, y(i))
         z_min = min(z_min, z(i))
         z_max = max(z_max, z(i))
         y_sum = y_sum + y(i)
         z_sum = z_sum + z(i)
      end do
      
      write(*,*) '========================================='
      write(*,*) '  MESH ORIENTATION ANALYSIS'
      write(*,*) '========================================='
      write(*,*) 'Nodes:', n_nodes
      write(*,*) 'y range: ', y_min, ' to ', y_max
      write(*,*) 'z range: ', z_min, ' to ', z_max
      write(*,*) 'Mean y:  ', y_sum / n_nodes
      write(*,*) 'Mean z:  ', z_sum / n_nodes
      write(*,*)
      write(*,*) 'Interpretation:'
      if (y_max > 0 .and. y_min < 0) then
         write(*,*) '  L extends in BOTH y directions'
      else if (y_max > 0) then
         write(*,*) '  L extends in +y direction only'
      else
         write(*,*) '  L extends in -y direction only'
      end if
      
      if (z_max > 0 .and. z_min < 0) then
         write(*,*) '  L extends in BOTH z directions'
      else if (z_max > 0) then
         write(*,*) '  L extends in +z direction only'
      else
         write(*,*) '  L extends in -z direction only'
      end if
      
      write(*,*)
      write(*,*) 'Theoretical shear centre for L-section:'
      write(*,*) '  Located on the bisector of the internal angle'
      write(*,*) '  Outside the profile, in the quadrant where'
      write(*,*) '  both extensions go'
      end
